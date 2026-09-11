"""Adapters for the locally vendored polymer builders.

The adapters deliberately keep the web layer small.  PSP and Polyply are
optional engines: their source is present in ``modeling/tools`` and their
Python dependencies are installed in the modeling virtual environment, but a
request still fails clearly if a required engine cannot be imported.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem


APP_ROOT = Path(__file__).resolve().parents[1]
TOOLS_ROOT = APP_ROOT / "tools"
PSP_ROOT = TOOLS_ROOT / "PSP"
POLYPLY_ROOT = TOOLS_ROOT / "polyply"
_ENGINE_LOCK = threading.Lock()


def _with_psp_compatibility() -> tuple[Any, Any]:
    """Import PSP while keeping it compatible with pandas 3.x.

    PSP is an older research codebase.  It uses ``DataFrame.append`` and the
    removed ``delim_whitespace`` argument; both are patched only for the
    duration of a PSP call and restored afterwards.
    """

    import pandas as pd

    old_append = getattr(pd.DataFrame, "append", None)
    old_read_csv = pd.read_csv

    if old_append is None:
        pd.DataFrame.append = lambda self, other, ignore_index=False, **kwargs: pd.concat(  # type: ignore[attr-defined]
            [self, other], ignore_index=ignore_index
        )

    def read_csv_compat(*args: Any, **kwargs: Any) -> Any:
        if kwargs.pop("delim_whitespace", False):
            kwargs.setdefault("sep", r"\s+")
        return old_read_csv(*args, **kwargs)

    pd.read_csv = read_csv_compat  # type: ignore[assignment]
    try:
        import psp.MoleculeBuilder as molecule_builder
    except Exception:
        if old_append is None:
            delattr(pd.DataFrame, "append")
        else:
            pd.DataFrame.append = old_append  # type: ignore[attr-defined]
        pd.read_csv = old_read_csv  # type: ignore[assignment]
        raise
    return pd, (old_append, old_read_csv, molecule_builder)


def _restore_psp_compatibility(pd: Any, state: tuple[Any, Any, Any]) -> None:
    old_append, old_read_csv, _ = state
    if old_append is None:
        if hasattr(pd.DataFrame, "append"):
            delattr(pd.DataFrame, "append")
    else:
        pd.DataFrame.append = old_append  # type: ignore[attr-defined]
    pd.read_csv = old_read_csv  # type: ignore[assignment]


def engine_status() -> dict[str, dict[str, Any]]:
    """Return runtime availability without starting a model job."""

    status: dict[str, dict[str, Any]] = {
        "rdkit_etkdg": {
            "id": "rdkit_etkdg",
            "label": "RDKit ETKDGv3 + MMFF",
            "status": "ready",
            "scope": "small molecules and pilot oligomers",
        },
        "psp": {
            "id": "psp",
            "label": "PSP · all-atom polymer builder",
            "status": "unavailable",
            "scope": "all-atom oligomers and polymer chains",
        },
        "polyply": {
            "id": "polyply",
            "label": "Polyply · condensed-phase builder",
            "status": "unavailable",
            "scope": "coarse-grained / force-field polymer coordinates",
        },
    }
    try:
        import psp  # noqa: F401
        import openbabel  # noqa: F401
    except Exception as exc:
        status["psp"]["reason"] = f"PSP dependencies unavailable: {exc}"
    else:
        status["psp"]["status"] = "ready"
        status["psp"]["reason"] = "Local PSP source and Open Babel available"
    try:
        import polyply  # noqa: F401
    except Exception as exc:
        status["polyply"]["reason"] = f"Polyply dependencies unavailable: {exc}"
    else:
        status["polyply"]["status"] = "ready"
        status["polyply"]["reason"] = "Local Polyply source available"
    return status


def _psp_repeat_smiles(polymer: str) -> str:
    repeats = {
        "PE": "[*]CC[*]",
        "PP": "C(C(C)[*])[*]",
        "PS": "C(C(c1ccccc1)[*])[*]",
    }
    try:
        return repeats[polymer.upper()]
    except KeyError as exc:
        raise ValueError("PSP currently supports PE, PP and PS; PET needs a validated repeat definition") from exc


def _load_pdb(path: Path) -> Chem.Mol:
    molecule = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)
    if molecule is None:
        raise RuntimeError(f"External engine produced an unreadable PDB: {path.name}")
    try:
        Chem.SanitizeMol(molecule)
    except Exception:
        # PDB bond records can be incomplete for coarse-grained outputs.  The
        # coordinates remain useful for visualization, so keep the molecule.
        pass
    return molecule


def build_with_psp(polymer: str, repeats: int, conformers: int, seed: int, workdir: Path) -> dict[str, Any]:
    """Build a finite oligomer with PSP's MoleculeBuilder."""

    if engine_status()["psp"]["status"] != "ready":
        raise RuntimeError(engine_status()["psp"].get("reason", "PSP is not available"))
    import pandas as pd

    data = pd.DataFrame(
        [[polymer, _psp_repeat_smiles(polymer), "C[*]", "C[*]"]],
        columns=["ID", "smiles", "LeftCap", "RightCap"],
    )
    output_dir = workdir / "psp_models"
    output_dir.mkdir(parents=True, exist_ok=True)
    old_cwd = Path.cwd()
    with _ENGINE_LOCK:
        pd_module, psp_state = _with_psp_compatibility()
        try:
            os.chdir(workdir)
            builder = psp_state[2].Builder(
                data,
                OutDir=str(output_dir),
                Length=[int(repeats)],
                NumConf=max(1, min(int(conformers), 4)),
                NCores=1,
                IrrStruc=False,
                OPLS=False,
                GAFF2=False,
                Subscript=True,
            )
            builder.Build()
        finally:
            os.chdir(old_cwd)
            _restore_psp_compatibility(pd_module, psp_state)

    candidates = sorted(output_dir.glob(f"{polymer}_N{int(repeats)}_C*.pdb"))
    if not candidates:
        raise RuntimeError("PSP finished without a PDB oligomer output")
    pdb_path = candidates[0]
    return {
        "mol": _load_pdb(pdb_path),
        "engine": "PSP MoleculeBuilder",
        "optimization": "PSP UFF / constrained local optimization",
        "warnings": [
            "PSP output is an all-atom capped oligomer; it is not an amorphous or density-equilibrated nanoparticle.",
            "PSP polymer tacticity is not specified for this pilot model.",
        ],
        "engine_metadata": {"source_output": str(pdb_path.relative_to(workdir)), "seed": int(seed)},
    }


def _write_polyply_topology(itp: Path, top: Path, molecule_name: str) -> None:
    top.write_text(
        "[ defaults ]\n"
        "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n"
        "1 2 no 1.0 1.0\n\n"
        f'#include "{itp.name}"\n\n'
        "[ system ]\nPolyply modeling pilot\n\n"
        "[ molecules ]\n"
        f"{molecule_name} 1\n",
        encoding="utf-8",
    )


def _gro_to_pdb(gro: Path, pdb: Path) -> None:
    from openbabel import openbabel as ob

    molecule = ob.OBMol()
    conversion = ob.OBConversion()
    conversion.SetInAndOutFormats("gro", "pdb")
    if not conversion.ReadFile(molecule, str(gro)):
        raise RuntimeError("Open Babel could not read Polyply GRO output")
    if not conversion.WriteFile(molecule, str(pdb)):
        raise RuntimeError("Open Babel could not convert Polyply GRO output to PDB")


def _fallback_rdkit_coordinates(polymer: str, repeats: int, seed: int) -> Chem.Mol:
    """Create a visible pilot coordinate set if Polyply lacks base FF files.

    The vendored Polyply source can generate the polymer ``.itp`` but its
    coordinate optimizer needs the full GROMOS/Martini base force-field files.
    Keeping this fallback explicit lets the UI show the requested third
    pipeline without claiming that a force-field-equilibrated box was made.
    """

    repeat_smiles = {
        "PE": "C" * (2 * repeats),
        "PP": "CC(C)" * repeats,
        "PS": "CC(c1ccccc1)" * repeats,
    }
    base = Chem.MolFromSmiles(repeat_smiles[polymer.upper()])
    if base is None:
        raise RuntimeError("Polyply fallback could not parse the repeat structure")
    molecule = Chem.AddHs(base)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(molecule, params=params) < 0:
        raise RuntimeError("Polyply fallback could not generate 3D coordinates")
    try:
        AllChem.MMFFOptimizeMolecule(molecule, maxIters=1500)
    except Exception:
        AllChem.UFFOptimizeMolecule(molecule, maxIters=1500)
    return molecule


def build_with_polyply(polymer: str, repeats: int, workdir: Path, seed: int = 20260826) -> dict[str, Any]:
    """Generate a small condensed-phase coordinate box with Polyply.

    PE, PP and PS have local 2016H66 definitions.  PET is intentionally
    rejected until a validated force-field block is added rather than silently
    substituting a different chemistry.
    """

    if polymer.upper() not in {"PE", "PP", "PS"}:
        raise ValueError("Polyply currently supports PE, PP and PS; PET needs a validated force-field block")
    if engine_status()["polyply"]["status"] != "ready":
        raise RuntimeError(engine_status()["polyply"].get("reason", "Polyply is not available"))

    script = POLYPLY_ROOT / "bin" / "polyply"
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(TOOLS_ROOT), str(POLYPLY_ROOT), env.get("PYTHONPATH", "")])
    itp = workdir / f"{polymer.lower()}_{int(repeats)}.itp"
    top = workdir / "system.top"
    gro = workdir / "structure.gro"
    molecule_name = polymer.lower()
    commands = [
        [sys.executable, str(script), "gen_params", "-lib", "2016H66", "-name", molecule_name, "-seq", f"{polymer}:{int(repeats)}", "-o", str(itp)],
    ]
    box = max(3.0, min(8.0, 2.5 + 0.3 * int(repeats)))
    _write_polyply_topology(itp, top, molecule_name)
    commands.append([
        sys.executable, str(script), "gen_coords", "-p", str(top), "-o", str(gro),
        "-name", molecule_name, "-box", f"{box:.3f}", f"{box:.3f}", f"{box:.3f}",
    ])
    polyply_error = None
    with _ENGINE_LOCK:
        try:
            for command in commands:
                completed = subprocess.run(command, cwd=workdir, env=env, capture_output=True, text=True, timeout=180)
                if completed.returncode != 0:
                    detail = (completed.stderr or completed.stdout).strip().splitlines()
                    raise RuntimeError(f"Polyply failed: {detail[-1] if detail else 'unknown error'}")
            pdb = workdir / "structure_polyply.pdb"
            _gro_to_pdb(gro, pdb)
            molecule = _load_pdb(pdb)
            engine = "Polyply gen_params + gen_coords"
            optimization = "Polyply random-walk placement"
            engine_metadata = {"topology": top.name, "coordinate_output": gro.name, "box_nm": [box, box, box]}
        except (RuntimeError, OSError, ValueError) as exc:
            polyply_error = str(exc)
            molecule = _fallback_rdkit_coordinates(polymer, repeats, seed)
            engine = "Polyply gen_params + RDKit coordinate fallback"
            optimization = "Polyply topology; RDKit pilot coordinates"
            engine_metadata = {
                "topology": top.name,
                "box_nm": [box, box, box],
                "coordinate_builder": "RDKit ETKDGv3 fallback",
                "polyply_coordinate_error": polyply_error,
            }
    return {
        "mol": molecule,
        "engine": engine,
        "optimization": optimization,
        "warnings": [
            "Polyply output is a force-field-based condensed-phase coordinate pilot; run minimization and equilibration before MD.",
            "Polyply uses a coarse-grained/force-field representation and is not an all-atom DFT geometry.",
            *( [f"Polyply gen_coords was unavailable in this local checkout; coordinates came from an explicit RDKit pilot fallback: {polyply_error}"] if polyply_error else [] ),
        ],
        "engine_metadata": engine_metadata,
    }
