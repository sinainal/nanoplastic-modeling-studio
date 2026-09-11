"""Small, optional GFN2-xTB runner for saved modeling structures.

The app must remain usable when the external ``xtb`` executable is not
installed, so status detection is explicit and a missing binary is returned as
an actionable result instead of crashing the API process.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import time
import uuid
from pathlib import Path
from typing import Any

from .builder import RUNS_DIR, _safe_model_dir, get_model


_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_BUNDLED_XTB_ROOT = _PROJECT_ROOT / "artice" / "execute_e1" / "tools" / "xtb_root"


def _find_xtb() -> str | None:
    configured = os.environ.get("XTB_BIN")
    candidates = [configured] if configured else []
    # E1 ships a self-contained xTB 6.6.1 installation.  Prefer an explicit
    # XTB_BIN override, then PATH/system installs, and finally this repository
    # copy so the app and the original E1 runner use the same executable.
    candidates.extend([
        "xtb",
        "/usr/local/bin/xtb",
        "/usr/bin/xtb",
        str(_BUNDLED_XTB_ROOT / "usr" / "bin" / "xtb"),
    ])
    for candidate in candidates:
        if not candidate:
            continue
        resolved = shutil.which(candidate) if os.path.basename(candidate) == candidate else candidate
        if resolved and Path(resolved).is_file() and os.access(resolved, os.X_OK):
            return str(Path(resolved).resolve())
    return None


def _execution_environment(binary: str) -> dict[str, str]:
    """Return runtime variables required by a system or bundled xTB binary."""

    env = os.environ.copy()
    binary_path = Path(binary).resolve()
    # The bundled layout is .../xtb_root/usr/bin/xtb.  Setting these variables
    # is harmless for a normal system install and avoids loader/resource errors
    # when the app discovers the repository-local binary automatically.
    bundled_usr = _BUNDLED_XTB_ROOT / "usr"
    if binary_path == (bundled_usr / "bin" / "xtb").resolve():
        lib_dir = bundled_usr / "lib" / "x86_64-linux-gnu"
        share_dir = bundled_usr / "share" / "xtb"
        if lib_dir.is_dir():
            existing = env.get("LD_LIBRARY_PATH", "")
            env["LD_LIBRARY_PATH"] = f"{lib_dir}:{existing}" if existing else str(lib_dir)
        if share_dir.is_dir():
            env["XTBHOME"] = str(share_dir)
    # Keep the defaults conservative; callers can tune them without changing
    # the API/CLI.  The old E1 execution used these same values.
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("MKL_NUM_THREADS", "4")
    return env


def get_xtb_status() -> dict[str, Any]:
    binary = _find_xtb()
    if binary:
        return {"status": "ready", "binary": binary, "method": "GFN2-xTB"}
    return {
        "status": "unavailable",
        "binary": None,
        "method": "GFN2-xTB",
        "reason": "The xtb executable is not installed or XTB_BIN is not configured.",
    }


def _result_from_output(output: str, elapsed: float, command: list[str], solvent: str, optimize: bool) -> dict[str, Any]:
    # xTB prints the final value in a bordered summary line (``| TOTAL
    # ENERGY ... |``); older builds also emit a plain ``TOTAL ENERGY`` line.
    # Parse both forms instead of relying on a particular table decoration.
    energy_matches = re.findall(r"TOTAL ENERGY\s+(-?(?:\d+(?:\.\d*)?|\.\d+))\s+Eh", output, flags=re.IGNORECASE)
    energy = float(energy_matches[-1]) if energy_matches else None
    gap_matches = re.findall(r"HOMO-LUMO gap\s+(-?(?:\d+(?:\.\d*)?|\.\d+))\s+eV", output, flags=re.IGNORECASE)
    gradient_matches = re.findall(r"gradient norm\s+(-?(?:\d+(?:\.\d*)?|\.\d+))\s+Eh", output, flags=re.IGNORECASE)
    return {
        "status": "completed",
        "method": "GFN2-xTB",
        "solvent": solvent,
        "optimize": optimize,
        "simulation_type": "optimization" if optimize else "singlepoint",
        "energy_hartree": energy,
        "homo_lumo_gap_ev": float(gap_matches[-1]) if gap_matches else None,
        "gradient_norm_hartree": float(gradient_matches[-1]) if gradient_matches else None,
        "wall_time_seconds": round(elapsed, 3),
        "command": command,
    }


def calculation_converged(output: str, optimize: bool, output_dir: Path, energy: float | None) -> bool:
    """Require positive optimizer evidence; a failure also contains 'converg'."""
    if energy is None or not math.isfinite(energy):
        return False
    if any((output_dir / name).exists() for name in ("NOT_CONVERGED", ".sccnotconverged")):
        return False
    if re.search(r"failed to converge|SCC.*not converged|abnormal termination", output, re.I):
        return False
    if not optimize:
        return True
    return bool(re.search(r"GEOMETRY OPTIMIZATION CONVERGED", output, re.I)) and (output_dir / "xtbopt.xyz").is_file()


def _write_result(output_dir: Path, result: dict[str, Any]) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def run_xtb(
    model_id: str,
    solvent: str = "water",
    optimize: bool = True,
    timeout_seconds: int = 300,
    simulation_type: str | None = None,
    charge: int | None = None,
    uhf: int = 0,
    optimization_mode: str = "full",
) -> dict[str, Any]:
    """Run an optional xTB optimization or single-point calculation.

    ``optimize`` remains for backwards compatibility with the first UI. New
    callers should use ``simulation_type`` explicitly.
    """

    if solvent not in {"none", "water"}:
        raise ValueError("xTB solvent must be 'none' or 'water'")
    model = get_model(model_id)
    model_dir = _safe_model_dir(model_id)
    xyz = model_dir / "structure.xyz"
    if not xyz.is_file():
        raise FileNotFoundError(xyz)
    if simulation_type is None:
        simulation_type = "optimization" if optimize else "singlepoint"
    if simulation_type not in {"optimization", "singlepoint"}:
        raise ValueError("simulation_type must be 'optimization' or 'singlepoint'")
    if optimization_mode not in {"surface_fixed", "full"}:
        raise ValueError("optimization_mode must be 'surface_fixed' or 'full'")
    optimize = simulation_type == "optimization"
    if charge is None:
        charge = int(model.get("formal_charge", 0))
    if uhf < 0:
        raise ValueError("uhf must be non-negative")
    binary = _find_xtb()
    output_dir = model_dir / "xtb"
    if not binary:
        return _write_result(output_dir, {
            "status": "unavailable",
            "method": "GFN2-xTB",
            "model_id": model_id,
            "model_title": model.get("title"),
            "simulation_type": simulation_type,
            "solvent": solvent,
            "optimize": optimize,
            "charge": charge,
            "uhf": uhf,
            "optimization_mode": optimization_mode,
            "reason": get_xtb_status()["reason"],
            "install_hint": "Install xtb 6.7+ and set XTB_BIN if it is outside PATH.",
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    # Archive before execution, never while publishing the current result.
    previous = [p for p in output_dir.iterdir() if p.is_file()]
    if previous:
        archive = output_dir / "history" / uuid.uuid4().hex
        archive.mkdir(parents=True)
        for path in previous:
            shutil.move(str(path), str(archive / path.name))
    input_xyz = output_dir / "input.xyz"
    shutil.copy2(xyz, input_xyz)
    command = [binary, input_xyz.name]
    if optimize:
        command.append("--opt")
    constrained_atoms = 0
    constraint_file: Path | None = None
    if optimize and optimization_mode == "surface_fixed":
        source = model.get("source") or {}
        if source.get("family") != "complex":
            raise ValueError("Surface-fixed xTB requires a saved complex pose")
        constrained_atoms = int(source.get("surface_atom_count", 0))
        if constrained_atoms < 1 or constrained_atoms >= int(model.get("atom_count", 0)):
            raise ValueError("Complex pose has invalid surface atom mapping")
        constraint_file = output_dir / "xcontrol.inp"
        constraint_file.write_text(
            f"$fix\n  atoms: 1-{constrained_atoms}\n$end\n",
            encoding="utf-8",
        )
        command.extend(["--input", constraint_file.name])
    if solvent == "water":
        command.extend(["--alpb", "water"])
    command.extend(["--chrg", str(int(charge)), "--uhf", str(int(uhf))])
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=output_dir,
            env=_execution_environment(binary),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        def decoded(value):
            return value.decode(errors="replace") if isinstance(value, bytes) else (value or "")
        (output_dir / "xtb.out").write_text(decoded(exc.stdout) + "\n" + decoded(exc.stderr), encoding="utf-8")
        return _write_result(output_dir, {
            "status": "failed",
            "method": "GFN2-xTB",
            "model_id": model_id,
            "model_title": model.get("title"),
            "simulation_type": simulation_type,
            "solvent": solvent,
            "optimize": optimize,
            "charge": charge,
            "uhf": uhf,
            "optimization_mode": optimization_mode,
            "surface_atoms_fixed": constrained_atoms,
            "reason": f"xTB timed out after {timeout_seconds} seconds",
            "wall_time_seconds": time.perf_counter() - started,
            "converged": False,
            "command": command,
            "files": ["xtb.out"],
        })
    elapsed = time.perf_counter() - started
    (output_dir / "xtb.out").write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip().splitlines()
        return _write_result(output_dir, {
            "status": "failed",
            "method": "GFN2-xTB",
            "model_id": model_id,
            "model_title": model.get("title"),
            "simulation_type": simulation_type,
            "solvent": solvent,
            "optimize": optimize,
            "charge": charge,
            "uhf": uhf,
            "optimization_mode": optimization_mode,
            "surface_atoms_fixed": constrained_atoms,
            "wall_time_seconds": round(elapsed, 3),
            "returncode": completed.returncode,
            "reason": detail[-1] if detail else "unknown xTB error",
            "files": ["xtb.out"],
        })
    result = _result_from_output(completed.stdout, elapsed, command, solvent, optimize)
    result.update({
        "model_id": model_id,
        "model_title": model.get("title"),
        "charge": charge,
        "uhf": uhf,
        "optimization_mode": optimization_mode,
        "surface_atoms_fixed": constrained_atoms,
    })
    result["converged"] = calculation_converged(completed.stdout + completed.stderr, optimize, output_dir, result["energy_hartree"])
    if not result["converged"]:
        result.update(status="failed", reason="Convergence or finite-energy validation failed; inspect xtb.out")
    # Do not report an old optimization geometry as an output of a later
    # single-point run.  The file is retained on disk for inspection, but the
    # result manifest lists only artifacts produced by this calculation type.
    result["files"] = ["xtb.out"]
    if optimize and (output_dir / "xtbopt.xyz").is_file():
        result["files"].append("xtbopt.xyz")
    if (output_dir / "charges").is_file():
        result["files"].append("charges")
    if constraint_file is not None:
        result["files"].append(constraint_file.name)
    return _write_result(output_dir, result)


def get_xtb_result(model_id: str) -> dict[str, Any] | None:
    """Return the last persisted xTB result for one model, if present."""

    path = _safe_model_dir(model_id) / "xtb" / "result.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def list_xtb_results(limit: int = 100) -> list[dict[str, Any]]:
    """List persisted xTB results for the Analysis page and CLI."""

    paths = sorted(RUNS_DIR.glob("mdl_*/xtb/result.json"), key=lambda path: path.stat().st_mtime, reverse=True)
    results: list[dict[str, Any]] = []
    for path in paths[: max(1, min(limit, 500))]:
        try:
            results.append(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError):
            continue
    return results
