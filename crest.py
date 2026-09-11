"""Optional constrained CREST/GFN-FF refinement for saved complex poses."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .builder import _safe_model_dir, get_model
from .pose_sets import add_coordinate_variants, parse_xyz_ensemble
from .xtb import _execution_environment, _find_xtb


_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_BUNDLED_CREST = _PROJECT_ROOT / "modeling" / "tools" / "crest" / "crest"
_CONSTRAINT_COMPAT_CREST = _PROJECT_ROOT / "modeling" / "tools" / "crest" / "versions" / "v2.12" / "crest"


def _find_crest(constraint_compatible: bool = False) -> str | None:
    if constraint_compatible and _CONSTRAINT_COMPAT_CREST.is_file() and os.access(_CONSTRAINT_COMPAT_CREST, os.X_OK):
        return str(_CONSTRAINT_COMPAT_CREST.resolve())
    configured = os.environ.get("CREST_BIN")
    candidates = [configured] if configured else []
    candidates.extend(["crest", "/usr/local/bin/crest", "/usr/bin/crest", str(_BUNDLED_CREST)])
    for candidate in candidates:
        if not candidate:
            continue
        resolved = shutil.which(candidate) if os.path.basename(candidate) == candidate else candidate
        if resolved and Path(resolved).is_file() and os.access(resolved, os.X_OK):
            return str(Path(resolved).resolve())
    return None


def _crest_environment() -> dict[str, str]:
    xtb = _find_xtb()
    env = _execution_environment(xtb) if xtb else os.environ.copy()
    if xtb:
        xtb_path = str(Path(xtb).resolve().parent)
        env["PATH"] = f"{xtb_path}:{env.get('PATH', '')}"
        bundled_share = Path(xtb).resolve().parents[1] / "share" / "xtb"
        if bundled_share.is_dir():
            env["XTBPATH"] = str(bundled_share)
    env.setdefault("OMP_NUM_THREADS", "4")
    return env


def get_crest_status() -> dict[str, Any]:
    crest = _find_crest()
    constrained = _find_crest(constraint_compatible=True)
    xtb = _find_xtb()
    if crest and constrained and xtb:
        return {"status": "ready", "binary": crest, "constraint_binary": constrained, "xtb_binary": xtb, "method": "CREST 2.12 constrained / GFN-FF", "note": "CREST 2.12 is used for fixed-substructure jobs because the official 3.0.2 binary has a known --cinp parser defect."}
    reasons = []
    if not crest:
        reasons.append("CREST executable was not found")
    if not xtb:
        reasons.append("xTB executable was not found")
    return {"status": "unavailable", "binary": crest, "xtb_binary": xtb, "method": "CREST / GFN-FF", "reason": "; ".join(reasons)}


def run_crest_refinement(model_id: str, top_k: int = 3, solvent: str = "water", timeout_seconds: int = 600) -> dict[str, Any]:
    if not 1 <= top_k <= 10:
        raise ValueError("CREST top_k must be 1-10")
    if solvent not in {"water", "none"}:
        raise ValueError("CREST solvent must be water or none")
    model = get_model(model_id)
    source = model.get("source") or {}
    if source.get("family") != "complex":
        raise ValueError("CREST refinement requires a saved complex pose")
    surface_atoms = int(source.get("surface_atom_count", 0))
    if surface_atoms < 1 or surface_atoms >= int(model.get("atom_count", 0)):
        raise ValueError("Complex pose has invalid surface atom mapping")
    # CREST 3.0.2 has a reproducible --cinp allocation bug for whole-
    # substructure constraints.  The official 2.12 binary remains compatible
    # with the same xTB constraint syntax and is used only for this workflow.
    binary = _find_crest(constraint_compatible=True)
    if not binary:
        raise RuntimeError(get_crest_status().get("reason", "CREST unavailable"))

    model_dir = _safe_model_dir(model_id)
    run_id = f"run_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"
    output_dir = model_dir / "crest" / run_id
    output_dir.mkdir(parents=True, exist_ok=True)
    input_xyz = output_dir / "input.xyz"
    shutil.copy2(model_dir / "structure.xyz", input_xyz)
    env = _crest_environment()

    constraint_command = [binary, input_xyz.name, "--constrain", f"1-{surface_atoms}"]
    constraint = subprocess.run(constraint_command, cwd=output_dir, env=env, capture_output=True, text=True, timeout=60, check=False)
    (output_dir / "constraint.out").write_text(constraint.stdout + "\n" + constraint.stderr, encoding="utf-8")
    control = output_dir / ".xcontrol.sample"
    if constraint.returncode != 0 or not control.is_file():
        raise RuntimeError("CREST could not create the fixed-surface constraint file")

    command = [binary, input_xyz.name, "--gfnff", "--quick", "--cinp", control.name, "-chrg", str(int(model.get("formal_charge", 0)))]
    if solvent == "water":
        command.extend(["-alpb", "water"])
    started = time.perf_counter()
    try:
        completed = subprocess.run(command, cwd=output_dir, env=env, capture_output=True, text=True, timeout=timeout_seconds, check=False)
    except subprocess.TimeoutExpired:
        result = {"status": "failed", "model_id": model_id, "pose_set_id": source.get("pose_set_id"), "method": "CREST 2.12 GFN-FF quick", "run_id": run_id, "reason": f"CREST timed out after {timeout_seconds} seconds"}
        (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return result
    elapsed = time.perf_counter() - started
    (output_dir / "crest.out").write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
    ensemble = output_dir / "crest_conformers.xyz"
    if completed.returncode != 0 or not ensemble.is_file():
        detail = (completed.stderr or completed.stdout).strip().splitlines()
        result = {"status": "failed", "model_id": model_id, "pose_set_id": source.get("pose_set_id"), "method": "CREST 2.12 GFN-FF quick", "run_id": run_id, "returncode": completed.returncode, "wall_time_seconds": round(elapsed, 3), "reason": detail[-1] if detail else "CREST failed", "files": ["crest.out", "constraint.out"]}
        (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return result

    coordinates = parse_xyz_ensemble(ensemble, top_k)
    variants = add_coordinate_variants(str(source["pose_set_id"]), model_id, coordinates, method="crest_gfnff")
    result = {
        "status": "completed",
        "model_id": model_id,
        "pose_set_id": source.get("pose_set_id"),
        "method": "CREST 2.12 GFN-FF quick",
        "run_id": run_id,
        "surface_mode": "fixed positional constraints",
        "solvent": solvent,
        "requested_top_k": top_k,
        "imported_variants": variants,
        "wall_time_seconds": round(elapsed, 3),
        "files": ["crest.out", "crest_conformers.xyz", "crest_best.xyz", "crest.energies", "constraint.out"],
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return result
