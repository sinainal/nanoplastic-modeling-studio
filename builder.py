from __future__ import annotations

import json
import math
import os
import platform
import re
import time
import uuid
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from .catalog import DOPAMINE_STATES, POLYMERS, polymer_smiles, resolve_dopamine_state
from .external_engines import build_with_polyply, build_with_psp


_LOCAL_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("MODELING_DATA_DIR", str(_LOCAL_ROOT if (_LOCAL_ROOT / "runs").is_dir() else Path.home() / ".local" / "share" / "nanoplastic-studio"))).expanduser().resolve()
RUNS_DIR = ROOT / "runs"
RUNS_DIR.mkdir(parents=True, exist_ok=True)


def _safe_model_dir(model_id: str) -> Path:
    if not model_id.startswith("mdl_") or not model_id.replace("_", "").isalnum():
        raise FileNotFoundError(model_id)
    target = (RUNS_DIR / model_id).resolve()
    if target.parent != RUNS_DIR.resolve():
        raise FileNotFoundError(model_id)
    return target


def _best_conformer(mol: Chem.Mol, count: int, seed: int) -> tuple[int, list[dict[str, Any]], str]:
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useRandomCoords = True
    params.pruneRmsThresh = 0.3
    params.numThreads = 0
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=count, params=params))
    if not ids:
        raise RuntimeError("RDKit could not generate a 3D conformer")

    force_field = "MMFF94"
    try:
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=1500)
    except Exception:
        force_field = "UFF"
        results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=1500)

    records: list[dict[str, Any]] = []
    for index, conf_id in enumerate(ids):
        status, energy = results[index]
        records.append({"id": int(conf_id), "converged": int(status) == 0, "energy": float(energy)})
    converged = [record for record in records if record["converged"] and math.isfinite(record["energy"])]
    if not converged:
        raise RuntimeError("No converged force-field conformer; cannot select a valid starting geometry")
    pool = converged
    best = min(pool, key=lambda record: record["energy"])
    return int(best["id"]), records, force_field


def _canonicalize(coords: np.ndarray) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    covariance = centered.T @ centered
    _, vectors = np.linalg.eigh(covariance)
    basis = vectors[:, ::-1]
    if np.linalg.det(basis) < 0:
        basis[:, 2] *= -1
    return centered @ basis


def _surface_patch(mol: Chem.Mol, conf_id: int) -> Chem.Mol:
    coords = np.array(mol.GetConformer(conf_id).GetPositions(), dtype=float)
    coords = _canonicalize(coords)
    extent = np.ptp(coords, axis=0)
    spacing_x = max(float(extent[0]) + 3.0, 6.0)
    spacing_y = max(float(extent[1]) + 3.0, 6.0)

    combined = Chem.Mol(mol)
    for _ in range(3):
        combined = Chem.CombineMols(combined, mol)
    combined = Chem.Mol(combined)
    conformer = Chem.Conformer(combined.GetNumAtoms())
    atom_count = mol.GetNumAtoms()
    shifts = [
        np.array([-spacing_x / 2, -spacing_y / 2, 0.0]),
        np.array([spacing_x / 2, -spacing_y / 2, 0.0]),
        np.array([-spacing_x / 2, spacing_y / 2, 0.0]),
        np.array([spacing_x / 2, spacing_y / 2, 0.0]),
    ]
    for copy_index, shift in enumerate(shifts):
        for atom_index, xyz in enumerate(coords + shift):
            conformer.SetAtomPosition(copy_index * atom_count + atom_index, tuple(float(v) for v in xyz))
    combined.RemoveAllConformers()
    combined.AddConformer(conformer, assignId=True)
    return combined


def _write_xyz(mol: Chem.Mol, path: Path, comment: str) -> None:
    conformer = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), comment]
    for atom_index, atom in enumerate(mol.GetAtoms()):
        point = conformer.GetAtomPosition(atom_index)
        lines.append(f"{atom.GetSymbol():2s} {point.x: .8f} {point.y: .8f} {point.z: .8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_outputs(mol: Chem.Mol, directory: Path, title: str) -> list[str]:
    pdb_path = directory / "structure.pdb"
    sdf_path = directory / "structure.sdf"
    xyz_path = directory / "structure.xyz"
    mol.SetProp("_Name", title)
    Chem.MolToPDBFile(mol, str(pdb_path), confId=0)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol, confId=0)
    writer.close()
    _write_xyz(mol, xyz_path, title)
    return [pdb_path.name, sdf_path.name, xyz_path.name]


def build_model(request: dict[str, Any]) -> dict[str, Any]:
    started = time.perf_counter()
    family = request.get("family", "small_molecule")
    model_kind = request.get("model_kind", "molecule")
    engine = str(request.get("engine", "rdkit_etkdg"))
    conformer_count = max(1, min(int(request.get("conformers", 6)), 32))
    seed = int(request.get("seed", 20260826))
    model_id = f"mdl_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    directory = _safe_model_dir(model_id)
    directory.mkdir(parents=True, exist_ok=False)
    warnings: list[str] = []

    if family == "small_molecule":
        pH = float(request.get("pH", 7.0))
        custom_smiles = str(request.get("smiles") or "").strip()
        chembl_id = str(request.get("chembl_id") or "").strip().upper()
        if custom_smiles or chembl_id:
            if custom_smiles and chembl_id:
                warnings.append("Both custom SMILES and ChEMBL ID were provided; the SMILES input was used.")
            if not custom_smiles:
                if not re.fullmatch(r"CHEMBL\d+", chembl_id):
                    raise ValueError("ChEMBL ID must look like CHEMBL followed by digits")
                url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
                try:
                    with urllib.request.urlopen(url, timeout=12) as response:
                        record = json.loads(response.read().decode("utf-8"))
                    custom_smiles = str(((record.get("molecule_structures") or {}).get("canonical_smiles") or "")).strip()
                except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
                    # PubChem's PUG endpoint exposes the same cross-reference
                    # and is a reliable fallback when the ChEMBL service is
                    # temporarily slow or unavailable.
                    try:
                        fallback = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chembl_id}/property/CanonicalSMILES/JSON"
                        with urllib.request.urlopen(fallback, timeout=12) as response:
                            record = json.loads(response.read().decode("utf-8"))
                        properties = (record.get("PropertyTable") or {}).get("Properties") or []
                        custom_smiles = str((properties[0] if properties else {}).get("ConnectivitySMILES") or (properties[0] if properties else {}).get("CanonicalSMILES") or "").strip()
                    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as fallback_exc:
                        raise RuntimeError(f"Could not retrieve {chembl_id} from ChEMBL or fallback resolver: {fallback_exc}") from exc
                if not custom_smiles:
                    raise ValueError(f"ChEMBL record {chembl_id} has no canonical SMILES")
            title = str(request.get("name") or chembl_id or "Custom small molecule")
            smiles = custom_smiles
            source_metadata = {"family": family, "compound": "custom", "pH": pH, "name": title, "chembl_id": chembl_id or None, "input": "SMILES" if not chembl_id else "SMILES/ChEMBL"}
        else:
            microstate = resolve_dopamine_state(pH, str(request.get("microstate", "auto")))
            definition = DOPAMINE_STATES[microstate]
            title = definition["label"]
            smiles = definition["smiles"]
            source_metadata = {"family": family, "compound": "dopamine", "pH": pH, "microstate": microstate, "microstate_resolution": "curated fixed-state selection"}
            if request.get("microstate", "auto") == "auto":
                warnings.append("pH selects a representative fixed microstate; it is not constant-pH sampling.")
    elif family == "polymer":
        polymer = str(request.get("polymer", "PE")).upper()
        repeats = int(request.get("repeats", POLYMERS.get(polymer, {}).get("default_repeats", 4)))
        pH = float(request.get("pH", 7.0))
        temperature_K = float(request.get("temperature_K", 310.15))
        smiles = polymer_smiles(polymer, repeats)
        title = f"{POLYMERS[polymer]['label']} · {repeats} repeat pilot"
        source_metadata = {
            "family": family,
            "polymer": polymer,
            "repeats": repeats,
            "repeat_formula": POLYMERS[polymer]["repeat_formula"],
            "tacticity": "unspecified" if polymer in {"PP", "PS"} else "not_set",
            "pH": pH,
            "temperature_K": temperature_K,
            "solvent": "water",
            "condition_label": "physiological pH 7 pilot" if abs(pH - 7.0) < 1e-9 else "fixed-pH pilot",
        }
        warnings.append("Capped oligomer model; it is not a density-matched amorphous or crystalline polymer.")
    else:
        raise ValueError(f"Unknown model family: {family}")

    if engine in {"psp", "polyply"}:
        if family != "polymer":
            raise ValueError(f"{engine} is currently available for polymer models only")
        if model_kind == "surface_patch" and engine == "polyply":
            raise ValueError("Polyply generates a condensed-phase coordinate box; choose Single oligomer for this engine")
        if engine == "psp":
            external = build_with_psp(polymer, repeats, conformer_count, seed, directory)
        else:
            external = build_with_polyply(polymer, repeats, directory, seed)
        best = external["mol"]
        title = f"{title} · {external['engine']}"
        warnings.extend(external["warnings"])
        source_metadata["engine"] = engine
        source_metadata["engine_metadata"] = external["engine_metadata"]
        if model_kind == "surface_patch":
            best = _surface_patch(best, 0)
            title += " · 2×2 surface proxy"
            warnings.append("Finite 2×2 proxy without periodic boundary conditions; intended for visualization and pilot work.")
        elif model_kind != "molecule":
            raise ValueError(f"Unknown model kind: {model_kind}")
        files = _write_outputs(best, directory, title)
        formal_charge = int(Chem.GetFormalCharge(best))
        formula = rdMolDescriptors.CalcMolFormula(best)
        manifest = {
            "id": model_id,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "title": title,
            "engine": external["engine"],
            "optimization": external["optimization"],
            "model_kind": model_kind,
            "smiles": smiles,
            "formula": formula,
            "formal_charge": formal_charge,
            "atom_count": best.GetNumAtoms(),
            "heavy_atom_count": best.GetNumHeavyAtoms(),
            "molecular_weight": round(float(Descriptors.MolWt(best)), 4),
            "conformers_requested": conformer_count,
            "conformers_generated": 1,
            "conformer_results": [],
            "seed": seed,
            "files": files,
            "source": source_metadata,
            "warnings": warnings,
            "validation": {
                "rdkit_sanitized": best.GetNumAtoms() > 0,
                "has_3d_coordinates": best.GetNumConformers() == 1,
                "formal_charge_checked": True,
                "density_validated": False,
                "periodic": False,
            },
            "versions": {
                "rdkit": rdBase.rdkitVersion,
                "python": platform.python_version(),
            },
            "wall_time_seconds": round(time.perf_counter() - started, 4),
        }
        (directory / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return manifest

    if engine != "rdkit_etkdg":
        raise ValueError(f"Unknown modeling engine: {engine}")

    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError("The generated SMILES could not be parsed")
    Chem.SanitizeMol(base)
    mol = Chem.AddHs(base)
    best_conf, conformers, force_field = _best_conformer(mol, conformer_count, seed)
    best = Chem.Mol(mol)
    selected = Chem.Conformer(mol.GetConformer(best_conf))
    best.RemoveAllConformers()
    best.AddConformer(selected, assignId=True)

    if model_kind == "surface_patch":
        if family != "polymer":
            raise ValueError("Surface patch is currently available for polymer models only")
        best = _surface_patch(best, 0)
        title += " · 2×2 surface proxy"
        warnings.append("Finite 2×2 proxy without periodic boundary conditions; intended for visualization and xTB pilot work.")
    elif model_kind != "molecule":
        raise ValueError(f"Unknown model kind: {model_kind}")

    files = _write_outputs(best, directory, title)
    formal_charge = int(Chem.GetFormalCharge(best))
    formula = rdMolDescriptors.CalcMolFormula(best)
    manifest = {
        "id": model_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "title": title,
        "engine": "RDKit ETKDGv3",
        "optimization": force_field,
        "model_kind": model_kind,
        "smiles": smiles,
        "formula": formula,
        "formal_charge": formal_charge,
        "atom_count": best.GetNumAtoms(),
        "heavy_atom_count": best.GetNumHeavyAtoms(),
        "molecular_weight": round(float(Descriptors.MolWt(best)), 4),
        "conformers_requested": conformer_count,
        "conformers_generated": len(conformers),
        "conformer_results": conformers,
        "seed": seed,
        "files": files,
        "source": source_metadata,
        "warnings": warnings,
        "validation": {
            "rdkit_sanitized": True,
            "has_3d_coordinates": best.GetNumConformers() == 1,
            "formal_charge_checked": True,
            "density_validated": False,
            "periodic": False,
        },
        "versions": {
            "rdkit": rdBase.rdkitVersion,
            "python": platform.python_version(),
        },
        "wall_time_seconds": round(time.perf_counter() - started, 4),
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def get_model(model_id: str) -> dict[str, Any]:
    manifest = _safe_model_dir(model_id) / "manifest.json"
    if not manifest.is_file():
        raise FileNotFoundError(model_id)
    return json.loads(manifest.read_text(encoding="utf-8"))


def model_file(model_id: str, filename: str = "structure.pdb") -> Path:
    if filename not in {"structure.pdb", "structure.sdf", "structure.xyz", "manifest.json"}:
        raise FileNotFoundError(filename)
    path = _safe_model_dir(model_id) / filename
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def list_models(limit: int = 30) -> list[dict[str, Any]]:
    manifests = sorted(RUNS_DIR.glob("mdl_*/manifest.json"), key=lambda path: path.stat().st_mtime, reverse=True)
    records = []
    # Simulation needs both the recently generated complex poses and their
    # shared parent dopamine/polymer models.  A 100-record cap could hide the
    # parent adsorbate after a large pilot run, leaving the pose builder with
    # an apparently empty dopamine selector.
    for path in manifests[: max(1, min(limit, 500))]:
        try:
            records.append(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError):
            continue
    return records
