"""Deterministic polymer–adsorbate pose-set construction.

Pose models are saved as ordinary ``mdl_*`` children so the existing xTB
runner and NGL endpoints can consume them without a parallel file format.
The pose-set manifest only groups those immutable child models and records how
their rigid-body placements were generated.
"""

from __future__ import annotations

import json
import math
import platform
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors

from .builder import RUNS_DIR, _safe_model_dir, _write_outputs, get_model, model_file


POSE_SETS_DIR = RUNS_DIR / "pose_sets"
POSE_SETS_DIR.mkdir(parents=True, exist_ok=True)


def _safe_pose_set_path(pose_set_id: str) -> Path:
    if not pose_set_id.startswith("pset_") or not pose_set_id.replace("_", "").isalnum():
        raise FileNotFoundError(pose_set_id)
    path = (POSE_SETS_DIR / f"{pose_set_id}.json").resolve()
    if path.parent != POSE_SETS_DIR.resolve():
        raise FileNotFoundError(pose_set_id)
    return path


def _load_sdf(model_id: str) -> Chem.Mol:
    path = model_file(model_id, "structure.sdf")
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    mol = next((item for item in supplier if item is not None), None)
    if mol is None or not mol.GetNumConformers():
        raise ValueError(f"Saved model has no usable SDF coordinates: {model_id}")
    return Chem.Mol(mol)


def _unit(vector: np.ndarray, fallback: np.ndarray | None = None) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-10:
        if fallback is None:
            raise ValueError("Cannot normalize a zero-length geometry vector")
        return _unit(fallback)
    return vector / norm


def _rotation_axis(axis: np.ndarray, angle_deg: float) -> np.ndarray:
    axis = _unit(axis)
    angle = math.radians(angle_deg)
    cross = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    return np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * (cross @ cross)


def _rotation_from_to(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    source = _unit(source)
    target = _unit(target)
    cross = np.cross(source, target)
    dot = float(np.clip(np.dot(source, target), -1.0, 1.0))
    if np.linalg.norm(cross) < 1e-10:
        if dot > 0:
            return np.eye(3)
        helper = np.array([1.0, 0.0, 0.0]) if abs(source[0]) < 0.8 else np.array([0.0, 1.0, 0.0])
        return _rotation_axis(np.cross(source, helper), 180.0)
    axis = _unit(cross)
    return _rotation_axis(axis, math.degrees(math.acos(dot)))


def _surface_frame(mol: Chem.Mol) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    coords = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    heavy = np.array([atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1], dtype=int)
    if len(heavy) < 2:
        heavy = np.arange(mol.GetNumAtoms(), dtype=int)
    heavy_coords = coords[heavy]
    center = heavy_coords.mean(axis=0)
    covariance = (heavy_coords - center).T @ (heavy_coords - center)
    values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(values)[::-1]
    axis_x = _unit(vectors[:, order[0]], np.array([1.0, 0.0, 0.0]))
    normal = _unit(vectors[:, order[-1]], np.array([0.0, 0.0, 1.0]))
    if float(np.dot(normal, np.array([0.0, 0.0, 1.0]))) < 0:
        normal *= -1.0
    axis_y = _unit(np.cross(normal, axis_x), np.array([0.0, 1.0, 0.0]))
    normal = _unit(np.cross(axis_x, axis_y), normal)
    return center, axis_x, axis_y, normal


def _ligand_features(mol: Chem.Mol) -> dict[str, np.ndarray]:
    coords = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    heavy = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    center = coords[heavy].mean(axis=0)
    aromatic = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    ring = aromatic if len(aromatic) >= 3 else heavy[: min(len(heavy), 6)]
    ring_coords = coords[ring]
    ring_center = ring_coords.mean(axis=0)
    _, _, vh = np.linalg.svd(ring_coords - ring_center, full_matrices=False)
    ring_normal = _unit(vh[-1], np.array([0.0, 0.0, 1.0]))
    ring_x = _unit(ring_coords[0] - ring_center, np.array([1.0, 0.0, 0.0]))
    ring_x = _unit(ring_x - ring_normal * np.dot(ring_x, ring_normal), np.array([1.0, 0.0, 0.0]))
    ring_y = _unit(np.cross(ring_normal, ring_x), np.array([0.0, 1.0, 0.0]))
    nitrogen = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    oxygens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    n_vector = coords[nitrogen].mean(axis=0) - ring_center if nitrogen else ring_x
    o_vector = coords[oxygens].mean(axis=0) - ring_center if oxygens else -ring_x
    return {
        "coords": coords,
        "center": center,
        "ring_center": ring_center,
        "ring_frame": np.column_stack([ring_x, ring_y, ring_normal]),
        "n_vector": n_vector,
        "o_vector": o_vector,
    }


def _pose_templates(count: int) -> list[dict[str, Any]]:
    base = [
        {"name": "parallel", "label": "Ring parallel · 0°", "kind": "ring", "tilt": 0.0, "azimuth": 0.0},
        {"name": "parallel_90", "label": "Ring parallel · 90°", "kind": "ring", "tilt": 0.0, "azimuth": 90.0},
        {"name": "tilted_plus", "label": "Ring tilted · +45°", "kind": "ring", "tilt": 45.0, "azimuth": 0.0},
        {"name": "tilted_minus", "label": "Ring tilted · −45°", "kind": "ring", "tilt": -45.0, "azimuth": 0.0},
        {"name": "edge", "label": "Ring edge-on · 0°", "kind": "ring", "tilt": 90.0, "azimuth": 0.0},
        {"name": "edge_90", "label": "Ring edge-on · 90°", "kind": "ring", "tilt": 90.0, "azimuth": 90.0},
        {"name": "amine_facing", "label": "Amine-facing", "kind": "amine", "tilt": None, "azimuth": 0.0},
        {"name": "catechol_facing", "label": "Catechol-facing", "kind": "catechol", "tilt": None, "azimuth": 0.0},
    ]
    result = []
    for index in range(count):
        template = dict(base[index % len(base)])
        cycle = index // len(base)
        if cycle:
            template["azimuth"] = float(template.get("azimuth", 0.0)) + cycle * 37.0
            template["name"] = f"{template['name']}_{cycle + 1}"
            template["label"] = f"{template['label']} · az {template['azimuth']:.0f}°"
        result.append(template)
    return result


def _fibonacci_directions(count: int = 96) -> np.ndarray:
    """Return deterministic, approximately equal-area unit-sphere points."""

    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    directions = []
    for index in range(count):
        z = 1.0 - 2.0 * (index + 0.5) / count
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        angle = golden_angle * index
        directions.append([radius * math.cos(angle), radius * math.sin(angle), z])
    return np.asarray(directions, dtype=float)


def _surface_region(
    point: np.ndarray,
    center: np.ndarray,
    axis_x: np.ndarray,
    axis_y: np.ndarray,
    normal: np.ndarray,
    heavy_coords: np.ndarray,
) -> str:
    relative = point - center
    extents = np.array([
        max(float(np.max(np.abs((heavy_coords - center) @ axis))), 1e-6)
        for axis in (axis_x, axis_y, normal)
    ])
    scaled = np.array([
        float(np.dot(relative, axis_x)) / extents[0],
        float(np.dot(relative, axis_y)) / extents[1],
        float(np.dot(relative, normal)) / extents[2],
    ])
    if abs(scaled[0]) >= 0.72:
        return "terminal +" if scaled[0] > 0 else "terminal −"
    if abs(scaled[2]) >= 0.48:
        return "upper face" if scaled[2] > 0 else "lower face"
    return "lateral +" if scaled[1] > 0 else "lateral −"


def _sample_surface_sites(
    mol: Chem.Mol,
    site_count: int,
    probe_radius_A: float = 1.4,
    seed: int = 20260827,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Sample solvent-accessible sites, then spread them by farthest points.

    Candidate points are generated on expanded van-der-Waals spheres and
    rejected when buried inside another expanded atom.  Greedy farthest-point
    selection prevents the old global-centre bias and gives deterministic
    coverage of faces, sides, and chain ends for a finite oligomer.
    """

    if not 1 <= site_count <= 25:
        raise ValueError("Surface site count must be 1-25")
    if not 0.0 <= probe_radius_A <= 3.0:
        raise ValueError("SASA probe radius must be 0.0-3.0 Å")
    coords = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    heavy_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_indices:
        raise ValueError("Surface model has no heavy atoms")
    periodic_table = Chem.GetPeriodicTable()
    heavy_coords = coords[heavy_indices]
    expanded_radii = np.array([
        float(periodic_table.GetRvdw(mol.GetAtomWithIdx(index).GetAtomicNum())) + probe_radius_A
        for index in heavy_indices
    ])
    contact_radii = expanded_radii - probe_radius_A
    directions = _fibonacci_directions(96)
    candidates: list[dict[str, Any]] = []
    for local_index, atom_index in enumerate(heavy_indices):
        atom_coord = heavy_coords[local_index]
        for direction in directions:
            accessible_point = atom_coord + direction * expanded_radii[local_index]
            distances = np.linalg.norm(heavy_coords - accessible_point, axis=1)
            mask = np.arange(len(heavy_indices)) != local_index
            if np.any(distances[mask] < expanded_radii[mask] - 0.05):
                continue
            contact_point = atom_coord + direction * contact_radii[local_index]
            candidates.append({
                "point": contact_point,
                "accessible_point": accessible_point,
                "normal": _unit(direction),
                "anchor_atom_index": atom_index,
                "anchor_element": mol.GetAtomWithIdx(atom_index).GetSymbol(),
            })
    if len(candidates) < site_count:
        raise ValueError(f"Only {len(candidates)} solvent-accessible candidates were found")

    accessible = np.asarray([item["accessible_point"] for item in candidates], dtype=float)
    rng = np.random.default_rng(seed)
    first = int(rng.integers(0, len(candidates)))
    selected_indices = [first]
    minimum_distances = np.linalg.norm(accessible - accessible[first], axis=1)
    while len(selected_indices) < site_count:
        minimum_distances[selected_indices] = -1.0
        next_index = int(np.argmax(minimum_distances))
        selected_indices.append(next_index)
        minimum_distances = np.minimum(
            minimum_distances,
            np.linalg.norm(accessible - accessible[next_index], axis=1),
        )

    center, axis_x, axis_y, global_normal = _surface_frame(mol)
    sites = []
    for number, candidate_index in enumerate(selected_indices, start=1):
        item = candidates[candidate_index]
        site = {
            "site_id": f"S{number:02d}",
            "region": _surface_region(item["point"], center, axis_x, axis_y, global_normal, heavy_coords),
            "point_A": [round(float(value), 6) for value in item["point"]],
            "accessible_point_A": [round(float(value), 6) for value in item["accessible_point"]],
            "normal": [round(float(value), 8) for value in item["normal"]],
            "anchor_atom_index": int(item["anchor_atom_index"]),
            "anchor_element": item["anchor_element"],
        }
        sites.append(site)

    selected_points = np.asarray([sites[index]["accessible_point_A"] for index in range(len(sites))], dtype=float)
    pairwise = np.linalg.norm(selected_points[:, None, :] - selected_points[None, :, :], axis=2)
    pairwise[pairwise < 1e-9] = np.inf
    min_spacing = None if len(sites) == 1 else float(np.min(pairwise))
    summary = {
        "algorithm": "atomic SASA candidates + greedy farthest-point sampling",
        "probe_radius_A": probe_radius_A,
        "candidate_count": len(candidates),
        "site_count": len(sites),
        "minimum_site_spacing_A": None if min_spacing is None else round(min_spacing, 4),
    }
    return sites, summary


def _minimum_heavy_clearance(surface: Chem.Mol, ligand: Chem.Mol, ligand_coords: np.ndarray) -> float:
    surface_coords = np.asarray(surface.GetConformer().GetPositions(), dtype=float)
    surface_heavy = [atom.GetIdx() for atom in surface.GetAtoms() if atom.GetAtomicNum() > 1]
    ligand_heavy = [atom.GetIdx() for atom in ligand.GetAtoms() if atom.GetAtomicNum() > 1]
    deltas = surface_coords[surface_heavy, None, :] - ligand_coords[ligand_heavy][None, :, :]
    return float(np.linalg.norm(deltas, axis=2).min())


def _oriented_ligand_coords(
    ligand: Chem.Mol,
    surface: Chem.Mol,
    template: dict[str, Any],
    distance: float,
    site: dict[str, Any],
) -> np.ndarray:
    surface_coords = np.asarray(surface.GetConformer().GetPositions(), dtype=float)
    center = np.asarray(site["point_A"], dtype=float)
    normal = _unit(np.asarray(site["normal"], dtype=float), np.array([0.0, 0.0, 1.0]))
    helper = np.array([1.0, 0.0, 0.0]) if abs(normal[0]) < 0.82 else np.array([0.0, 1.0, 0.0])
    sx = _unit(helper - normal * np.dot(helper, normal), np.array([1.0, 0.0, 0.0]))
    sy = _unit(np.cross(normal, sx), np.array([0.0, 1.0, 0.0]))
    features = _ligand_features(ligand)
    ligand_centered = features["coords"] - features["center"]

    if template["kind"] == "ring":
        tilt_matrix = _rotation_axis(sx, float(template["tilt"]))
        target_normal = tilt_matrix @ normal
        target_x = sx
        target_y = _unit(np.cross(target_normal, target_x), sy)
        target_x = _unit(np.cross(target_y, target_normal), sx)
        target_frame = np.column_stack([target_x, target_y, target_normal])
        rotated = ligand_centered @ features["ring_frame"] @ target_frame.T
    else:
        source = features["n_vector"] if template["kind"] == "amine" else features["o_vector"]
        rotation = _rotation_from_to(source, -normal)
        rotated = ligand_centered @ rotation.T

    azimuth = float(template.get("azimuth", 0.0))
    if abs(azimuth) > 1e-9:
        rotated = rotated @ _rotation_axis(normal, azimuth).T

    ligand_heavy = [atom.GetIdx() for atom in ligand.GetAtoms() if atom.GetAtomicNum() > 1]
    ligand_bottom = float(np.min(rotated[ligand_heavy] @ normal))
    translation = center + normal * (distance - ligand_bottom)
    placed = rotated + translation

    # The plane offset above is orientation-stable, but an irregular finite
    # oligomer can still approach the ligand from the side.  Move only along
    # the surface normal until the true heavy-atom clearance equals the user
    # requested value.  This prevents hidden steric clashes while preserving
    # the selected tilt and azimuth.
    def minimum_clearance(extra: float) -> float:
        shifted = placed[ligand_heavy] + normal * extra
        surface_heavy = [atom.GetIdx() for atom in surface.GetAtoms() if atom.GetAtomicNum() > 1]
        deltas = surface_coords[surface_heavy, None, :] - shifted[None, :, :]
        return float(np.linalg.norm(deltas, axis=2).min())

    clearance = minimum_clearance(0.0)
    if abs(clearance - distance) > 1e-4:
        if clearance < distance:
            low, high = 0.0, max(1.0, distance)
            while minimum_clearance(high) < distance:
                high *= 2.0
                if high > 100.0:
                    raise ValueError("Could not place adsorbate without a steric clash")
        else:
            high, low = 0.0, -max(1.0, distance / 2.0)
            while minimum_clearance(low) > distance:
                low *= 2.0
                if low < -100.0:
                    # A highly concave site can be non-monotonic along its
                    # normal.  Keep the safe, non-clashing placement rather
                    # than forcing the ligand through the surface.
                    low = high
                    break
        if low != high:
            for _ in range(40):
                middle = (low + high) / 2.0
                if minimum_clearance(middle) < distance:
                    low = middle
                else:
                    high = middle
            placed = placed + normal * high
    return placed


def _set_pdb_info(mol: Chem.Mol, surface_atoms: int, pose_number: int = 1) -> None:
    for index, atom in enumerate(mol.GetAtoms()):
        is_surface = index < surface_atoms
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{(index + 1) % 1000:03d}"[-4:])
        info.SetResidueName("POL" if is_surface else "DAH")
        info.SetResidueNumber(1 if is_surface else pose_number)
        info.SetChainId("A" if is_surface else "B")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)


def _combined_pose(surface: Chem.Mol, ligand: Chem.Mol, ligand_coords: np.ndarray, pose_number: int) -> Chem.Mol:
    ligand_copy = Chem.Mol(ligand)
    conf = Chem.Conformer(ligand_copy.GetNumAtoms())
    for index, xyz in enumerate(ligand_coords):
        conf.SetAtomPosition(index, tuple(float(value) for value in xyz))
    ligand_copy.RemoveAllConformers()
    ligand_copy.AddConformer(conf, assignId=True)
    combined = Chem.Mol(Chem.CombineMols(surface, ligand_copy))
    _set_pdb_info(combined, surface.GetNumAtoms(), pose_number)
    return combined


def _write_child_manifest(
    combined: Chem.Mol,
    surface_manifest: dict[str, Any],
    adsorbate_manifest: dict[str, Any],
    pose_set_id: str,
    pose_number: int,
    pose_meta: dict[str, Any],
    seed: int,
    parent_pose_model_id: str | None = None,
) -> dict[str, Any]:
    model_id = f"mdl_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    directory = _safe_model_dir(model_id)
    directory.mkdir(parents=True, exist_ok=False)
    pose_id = f"P{pose_number:03d}"
    title = f"{surface_manifest['title']} + {adsorbate_manifest['title']} · {pose_id}"
    files = _write_outputs(combined, directory, title)
    formal_charge = int(Chem.GetFormalCharge(combined))
    manifest = {
        "id": model_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "title": title,
        "engine": "SASA farthest-point pose builder" if parent_pose_model_id is None else "CREST GFN-FF refinement",
        "optimization": "Surface-site rigid-body placement" if parent_pose_model_id is None else "GFN-FF conformer sampling",
        "model_kind": "complex_pose",
        "smiles": f"{surface_manifest.get('smiles', '')}.{adsorbate_manifest.get('smiles', '')}",
        "formula": rdMolDescriptors.CalcMolFormula(combined),
        "formal_charge": formal_charge,
        "atom_count": combined.GetNumAtoms(),
        "heavy_atom_count": combined.GetNumHeavyAtoms(),
        "molecular_weight": round(float(Descriptors.MolWt(combined)), 4),
        "conformers_requested": 1,
        "conformers_generated": 1,
        "conformer_results": [],
        "seed": seed,
        "files": files,
        "source": {
            "family": "complex",
            "pose_set_id": pose_set_id,
            "pose_id": pose_id,
            "pose_number": pose_number,
            "surface_model_id": surface_manifest["id"],
            "adsorbate_model_id": adsorbate_manifest["id"],
            "surface_atom_count": surface_manifest["atom_count"],
            "adsorbate_atom_count": adsorbate_manifest["atom_count"],
            "parent_pose_model_id": parent_pose_model_id,
            "pose": pose_meta,
        },
        "warnings": [
            "Finite oligomer complex pose; this is not a periodic nanoplastic surface.",
            "Starting-pose label is a sampling hypothesis, not a final binding-mode assignment.",
        ],
        "validation": {
            "rdkit_sanitized": True,
            "has_3d_coordinates": True,
            "formal_charge_checked": True,
            "periodic": False,
            "pose_atom_mapping_preserved": True,
        },
        "versions": {"rdkit": rdBase.rdkitVersion, "python": platform.python_version()},
        "wall_time_seconds": 0.0,
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def create_pose_set(
    surface_model_id: str,
    adsorbate_model_id: str,
    count: int | None = None,
    distance_A: float = 3.2,
    seed: int = 20260827,
    strategy: str = "surface_scan",
    site_count: int | None = None,
    orientations_per_site: int | None = None,
    probe_radius_A: float = 1.4,
    orientation_names: list[str] | None = None,
) -> dict[str, Any]:
    if strategy not in {"surface_scan", "systematic", "crest", "crest_assisted", "manual_seed"}:
        raise ValueError("Unknown pose strategy")
    # ``count`` is retained for older CLI/API clients.  New clients define a
    # surface grid explicitly as sites × orientations.
    if site_count is None and orientations_per_site is None and count is not None:
        legacy_count = int(count)
        if not 1 <= legacy_count <= 100:
            raise ValueError("Pose count must be 1-100")
        site_count = min(5, legacy_count)
        orientations_per_site = max(1, math.ceil(legacy_count / site_count))
        total_limit = legacy_count
    else:
        site_count = 5 if site_count is None else int(site_count)
        orientations_per_site = 4 if orientations_per_site is None else int(orientations_per_site)
        total_limit = site_count * orientations_per_site
    if not 1 <= site_count <= 25:
        raise ValueError("Surface site count must be 1-25")
    if not 1 <= orientations_per_site <= 8:
        raise ValueError("Orientations per site must be 1-8")
    if not 1 <= total_limit <= 100:
        raise ValueError("Pose set must contain 1-100 poses")
    if not 2.5 <= distance_A <= 8.0:
        raise ValueError("Initial separation must be 2.5-8.0 Å")
    surface_manifest = get_model(surface_model_id)
    adsorbate_manifest = get_model(adsorbate_model_id)
    if (surface_manifest.get("source") or {}).get("family") != "polymer":
        raise ValueError("Surface model must be a saved polymer model")
    if (adsorbate_manifest.get("source") or {}).get("compound") != "dopamine":
        raise ValueError("Adsorbate model must be a saved dopamine model")

    surface = _load_sdf(surface_model_id)
    ligand = _load_sdf(adsorbate_model_id)
    pose_set_id = f"pset_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    sites, sampling_summary = _sample_surface_sites(surface, site_count, probe_radius_A, seed)
    templates = _pose_templates(orientations_per_site)
    if orientation_names is not None:
        available = {t["name"]: t for t in _pose_templates(8)}
        if len(orientation_names) != orientations_per_site or any(n not in available for n in orientation_names):
            raise ValueError("orientation_names must match orientations_per_site and known template names")
        templates = [dict(available[n]) for n in orientation_names]
    poses = []
    pose_number = 0
    for site in sites:
        for template in templates:
            if pose_number >= total_limit:
                break
            pose_number += 1
            coords = _oriented_ligand_coords(ligand, surface, template, distance_A, site)
            combined = _combined_pose(surface, ligand, coords, pose_number)
            clearance = _minimum_heavy_clearance(surface, ligand, coords)
            pose_meta = {
                **template,
                "site_id": site["site_id"],
                "surface_region": site["region"],
                "surface_point_A": site["point_A"],
                "surface_normal": site["normal"],
                "anchor_atom_index": site["anchor_atom_index"],
                "anchor_element": site["anchor_element"],
                "initial_distance_A": distance_A,
                "actual_clearance_A": round(clearance, 5),
                "probe_radius_A": probe_radius_A,
                "strategy": "surface_scan",
            }
            child = _write_child_manifest(
                combined, surface_manifest, adsorbate_manifest, pose_set_id,
                pose_number, pose_meta, seed + pose_number,
            )
            seed_method = "manual_seed" if strategy == "manual_seed" else "surface_scan"
            label = f"{site['site_id']} · {site['region']} · {template['label']}"
            poses.append({
                "pose_id": child["source"]["pose_id"],
                "model_id": child["id"],
                "label": label,
                "method": seed_method,
                "parent_model_id": None,
                "parameters": pose_meta,
            })

    manifest = {
        "id": pose_set_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "title": f"{surface_manifest['title']} + {adsorbate_manifest['title']}",
        "surface_model_id": surface_model_id,
        "adsorbate_model_id": adsorbate_model_id,
        "strategy": "surface_scan",
        "initial_distance_A": distance_A,
        "probe_radius_A": probe_radius_A,
        "site_count": len(sites),
        "orientations_per_site": orientations_per_site,
        "surface_sampling": {**sampling_summary, "sites": sites},
        "seed": seed,
        "formal_charge": int(surface_manifest.get("formal_charge", 0)) + int(adsorbate_manifest.get("formal_charge", 0)),
        "surface_atom_count": int(surface_manifest["atom_count"]),
        "adsorbate_atom_count": int(adsorbate_manifest["atom_count"]),
        "poses": poses,
        "warnings": [
            "SASA/farthest-point poses are reproducible surface hypotheses, not a thermodynamic ensemble.",
            "CREST refinement is optional and local; it does not replace whole-surface coverage.",
        ],
    }
    _safe_pose_set_path(pose_set_id).write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def get_pose_set(pose_set_id: str) -> dict[str, Any]:
    path = _safe_pose_set_path(pose_set_id)
    if not path.is_file():
        raise FileNotFoundError(pose_set_id)
    return json.loads(path.read_text(encoding="utf-8"))


def list_pose_sets(limit: int = 50) -> list[dict[str, Any]]:
    paths = sorted(POSE_SETS_DIR.glob("pset_*.json"), key=lambda item: item.stat().st_mtime, reverse=True)
    records = []
    for path in paths[: max(1, min(limit, 100))]:
        try:
            records.append(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError):
            continue
    return records


def add_coordinate_variants(pose_set_id: str, parent_model_id: str, coordinates: list[np.ndarray], method: str = "crest_gfnff") -> list[dict[str, Any]]:
    pose_set = get_pose_set(pose_set_id)
    parent = get_model(parent_model_id)
    if (parent.get("source") or {}).get("pose_set_id") != pose_set_id:
        raise ValueError("Parent pose does not belong to the pose set")
    surface_manifest = get_model(pose_set["surface_model_id"])
    adsorbate_manifest = get_model(pose_set["adsorbate_model_id"])
    template_mol = _load_sdf(parent_model_id)
    template_coords = np.asarray(template_mol.GetConformer().GetPositions(), dtype=float)
    surface_atoms = int(pose_set["surface_atom_count"])
    parent_record = next(
        (record for record in pose_set["poses"] if record.get("model_id") == parent_model_id),
        {},
    )
    parent_parameters = dict(parent_record.get("parameters") or {})
    added = []
    for coords in coordinates:
        if len(coords) != template_mol.GetNumAtoms():
            continue
        # CREST is free to translate/rotate the complete constrained complex.
        # Superpose its fixed polymer atoms onto the parent pose before saving,
        # so ensemble overlays compare ligand motion in one common frame.
        moving_surface = np.asarray(coords[:surface_atoms], dtype=float)
        target_surface = template_coords[:surface_atoms]
        moving_center = moving_surface.mean(axis=0)
        target_center = target_surface.mean(axis=0)
        covariance = (moving_surface - moving_center).T @ (target_surface - target_center)
        left, _, right = np.linalg.svd(covariance)
        rotation = left @ right
        if np.linalg.det(rotation) < 0:
            left[:, -1] *= -1
            rotation = left @ right
        coords = (np.asarray(coords, dtype=float) - moving_center) @ rotation + target_center
        mol = Chem.Mol(template_mol)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for atom_index, xyz in enumerate(coords):
            conf.SetAtomPosition(atom_index, tuple(float(value) for value in xyz))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        number = len(pose_set["poses"]) + 1
        parent_pose = (parent.get("source") or {}).get("pose_id", "pose")
        # Preserve the sampled site, region, surface normal and starting
        # orientation.  A CREST child is a local refinement of that hypothesis,
        # not a newly discovered whole-surface site.
        meta = {
            **parent_parameters,
            "name": f"{method}_{parent_pose}",
            "label": f"CREST refinement of {parent_pose}",
            "kind": "crest",
            "strategy": method,
            "parent_pose_id": parent_pose,
        }
        child = _write_child_manifest(mol, surface_manifest, adsorbate_manifest, pose_set_id, number, meta, int(pose_set["seed"]) + number, parent_model_id)
        record = {"pose_id": child["source"]["pose_id"], "model_id": child["id"], "label": meta["label"], "method": method, "parent_model_id": parent_model_id, "parameters": meta}
        pose_set["poses"].append(record)
        added.append(record)
    _safe_pose_set_path(pose_set_id).write_text(json.dumps(pose_set, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return added


def parse_xyz_ensemble(path: Path, limit: int) -> list[np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    structures: list[np.ndarray] = []
    cursor = 0
    while cursor < len(lines) and len(structures) < limit:
        try:
            count = int(lines[cursor].strip())
        except (ValueError, IndexError):
            break
        block = []
        for line in lines[cursor + 2: cursor + 2 + count]:
            fields = line.split()
            if len(fields) < 4:
                break
            block.append([float(fields[1]), float(fields[2]), float(fields[3])])
        if len(block) == count:
            structures.append(np.asarray(block, dtype=float))
        cursor += count + 2
    return structures


def overlay_pdb(
    pose_set_id: str,
    offset: int = 0,
    limit: int = 8,
    pose_ids: list[str] | None = None,
) -> str:
    pose_set = get_pose_set(pose_set_id)
    if pose_ids:
        requested = set(pose_ids)
        selected = [
            record for record in pose_set["poses"]
            if record.get("pose_id") in requested or record.get("model_id") in requested
        ][:100]
    else:
        selected = pose_set["poses"][offset: offset + min(max(1, limit), 100)]
    if not selected:
        raise FileNotFoundError("No poses in requested page")
    # Keep a pose's chain/color identity stable even when the user changes
    # the selected subset (P003 should not become P002's colour on the next
    # click).  The chain is therefore based on its position in the full set,
    # not its position in the current filtered overlay.
    pose_order = {record.get("model_id"): index for index, record in enumerate(pose_set["poses"])}
    surface = _load_sdf(pose_set["surface_model_id"])
    adsorbate_atoms = int(pose_set["adsorbate_atom_count"])
    combined = Chem.Mol(surface)
    chains = "BCDEFGHIJKLMNOPQRSTUVWXYZbcdefghijklmnopqrstuvwxyz0123456789"
    for index, record in enumerate(selected, start=1):
        pose_mol = _load_sdf(record["model_id"])
        coords = np.asarray(pose_mol.GetConformer().GetPositions(), dtype=float)[-adsorbate_atoms:]
        ligand = _load_sdf(pose_set["adsorbate_model_id"])
        conf = Chem.Conformer(ligand.GetNumAtoms())
        for atom_index, xyz in enumerate(coords):
            conf.SetAtomPosition(atom_index, tuple(float(value) for value in xyz))
        ligand.RemoveAllConformers()
        ligand.AddConformer(conf, assignId=True)
        combined = Chem.Mol(Chem.CombineMols(combined, ligand))
    surface_count = surface.GetNumAtoms()
    for atom_index, atom in enumerate(combined.GetAtoms()):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{(atom_index + 1) % 1000:03d}"[-4:])
        if atom_index < surface_count:
            info.SetResidueName("POL"); info.SetResidueNumber(1); info.SetChainId("A")
        else:
            selected_index = (atom_index - surface_count) // adsorbate_atoms
            pose_index = pose_order.get(selected[selected_index].get("model_id"), selected_index)
            info.SetResidueName("DAH"); info.SetResidueNumber(offset + pose_index + 1); info.SetChainId(chains[pose_index % len(chains)])
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    return Chem.MolToPDBBlock(combined)
