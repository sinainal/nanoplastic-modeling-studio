from __future__ import annotations

import threading
import uuid
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from . import __version__
from .builder import build_model, get_model, list_models, model_file, RUNS_DIR
from .xtb import get_xtb_result, get_xtb_status, list_xtb_results, run_xtb
from .catalog import public_catalog
from .crest import get_crest_status, run_crest_refinement
from .pose_sets import create_pose_set, get_pose_set, list_pose_sets, overlay_pdb
from .experiments import create_experiment, get_experiment, list_experiments


ROOT = Path(__file__).resolve().parent
STATIC_DIR = ROOT / "static"
WORKSPACE_ROOT = ROOT.parents[1]

app = FastAPI(title="Nanoplastic Modeling Studio", version=__version__)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

_XTB_EXECUTOR = ThreadPoolExecutor(max_workers=1, thread_name_prefix="xtb-job")
_XTB_JOBS: dict[str, dict] = {}
_XTB_JOBS_LOCK = threading.Lock()
_CREST_JOBS: dict[str, dict] = {}
_CREST_JOBS_LOCK = threading.Lock()
_EXPERIMENT_EXECUTOR = ThreadPoolExecutor(max_workers=1, thread_name_prefix="experiment-job")
_EXPERIMENT_JOBS: dict[str, dict] = {}
_EXPERIMENT_JOBS_LOCK = threading.Lock()


class ModelRequest(BaseModel):
    family: Literal["small_molecule", "polymer"] = "small_molecule"
    model_kind: Literal["molecule", "surface_patch"] = "molecule"
    engine: Literal["rdkit_etkdg", "psp", "polyply"] = "rdkit_etkdg"
    pH: float = Field(default=7.0, ge=0.0, le=14.0)
    temperature_K: float = Field(default=310.15, ge=200.0, le=500.0)
    microstate: str = "auto"
    polymer: str = "PE"
    repeats: int = Field(default=6, ge=1, le=24)
    conformers: int = Field(default=6, ge=1, le=32)
    seed: int = 20260826
    # Optional custom small-molecule input.  If supplied, it replaces the
    # curated dopamine catalog entry while retaining the same RDKit pipeline.
    smiles: str | None = Field(default=None, max_length=4000)
    chembl_id: str | None = Field(default=None, max_length=32)
    name: str | None = Field(default=None, max_length=160)


class XTBRequest(BaseModel):
    solvent: Literal["none", "water"] = "water"
    simulation_type: Literal["optimization", "singlepoint"] = "optimization"
    optimize: bool | None = None  # compatibility with the first UI/API
    charge: int | None = Field(default=None, ge=-20, le=20)
    uhf: int = Field(default=0, ge=0, le=20)
    optimization_mode: Literal["surface_fixed", "full"] = "full"
    timeout_seconds: int = Field(default=300, ge=30, le=1800)


class XTBJobRequest(XTBRequest):
    model_ids: list[str] = Field(min_length=1, max_length=100)


class PoseSetRequest(BaseModel):
    surface_model_id: str
    adsorbate_model_id: str
    count: int | None = Field(default=None, ge=1, le=100)
    # Optional values preserve the first API's ``count``-only contract.  The
    # current UI always submits both values and therefore defines an explicit
    # surface grid.
    site_count: int | None = Field(default=None, ge=1, le=25)
    orientations_per_site: int | None = Field(default=None, ge=1, le=8)
    distance_A: float = Field(default=3.2, ge=2.5, le=8.0)
    probe_radius_A: float = Field(default=1.4, ge=0.0, le=3.0)
    seed: int = 20260827
    strategy: Literal["surface_scan", "systematic", "crest_assisted", "manual_seed"] = "surface_scan"


class CrestJobRequest(BaseModel):
    model_ids: list[str] = Field(min_length=1, max_length=20)
    top_k: int = Field(default=3, ge=1, le=10)
    solvent: Literal["none", "water"] = "water"
    timeout_seconds: int = Field(default=600, ge=60, le=3600)


class ExperimentPlanRequest(BaseModel):
    title: str = Field(min_length=1, max_length=120)
    stage: Literal["model_generation", "pose_generation", "simulation"]
    design: dict
    notes: str = Field(default="", max_length=1000)


class ExperimentRunRequest(BaseModel):
    mode: Literal["models_only", "end_to_end", "ready_pose_sets", "composed"]
    model_plan: dict = Field(default_factory=dict)
    pose_plan: dict = Field(default_factory=dict)
    simulation_plan: dict = Field(default_factory=dict)
    # A composed experiment is an ordered collection of independent work
    # items.  The three legacy fields are retained so existing UI/CLI calls
    # remain valid.
    work_items: list[dict] = Field(default_factory=list)


@app.get("/", include_in_schema=False)
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/simulation", include_in_schema=False)
def simulation() -> FileResponse:
    return FileResponse(STATIC_DIR / "simulation.html")


@app.get("/simulation/run", include_in_schema=False)
def simulation_run() -> FileResponse:
    return FileResponse(STATIC_DIR / "simulation-run.html")


@app.get("/experiments", include_in_schema=False)
def experiments_workspace() -> FileResponse:
    return FileResponse(STATIC_DIR / "experiments.html")


@app.get("/analysis", include_in_schema=False)
def analysis() -> FileResponse:
    return FileResponse(STATIC_DIR / "analysis.html")


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok", "service": "modeling-studio", "version": __version__, "xtb": get_xtb_status(), "crest": get_crest_status()}


@app.get("/api/crest/status")
def crest_status() -> dict:
    return get_crest_status()


@app.get("/api/pose-sets")
def pose_sets(limit: int = Query(default=50, ge=1, le=100)) -> list[dict]:
    return list_pose_sets(limit)


@app.get("/api/experiments")
def experiments(limit: int = Query(default=100, ge=1, le=200)) -> list[dict]:
    return list_experiments(limit)


@app.post("/api/experiments", status_code=201)
def make_experiment(request: ExperimentPlanRequest) -> dict:
    try:
        return create_experiment(**request.model_dump())
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/api/experiments/{experiment_id}")
def experiment(experiment_id: str) -> dict:
    try:
        return get_experiment(experiment_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Experiment plan not found") from exc


def _experiment_event(job_id: str, *, stage: str, label: str, status: str, **extra: object) -> None:
    with _EXPERIMENT_JOBS_LOCK:
        job = _EXPERIMENT_JOBS[job_id]
        event = {"index": len(job["events"]) + 1, "stage": stage, "label": label, "status": status, **extra}
        job["events"].append(event)
        if status == "running":
            job["total"] += 1
        if status in {"completed", "failed"}:
            job["completed"] += 1
        if extra.get("model_id"):
            job["latest_model_id"] = extra["model_id"]
        if extra.get("pose_set_id"):
            job["latest_pose_set_id"] = extra["pose_set_id"]


def _as_int(value: object, default: int, *, minimum: int, maximum: int) -> int:
    try:
        parsed = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        parsed = default
    return max(minimum, min(maximum, parsed))


def _work_prefix(index: int, title: str = "") -> str:
    return f"W{index + 1:02d}" + (f" · {title}" if title else "")


def _run_model_matrix(job_id: str, plan: dict, *, prefix: str = "") -> list[dict]:
    polymers = [item for item in plan.get("polymers", []) if item in {"PE", "PP", "PS", "PET"}]
    repeats = [_as_int(item, 6, minimum=1, maximum=24) for item in plan.get("repeats", [])]
    engines = [item for item in plan.get("engines", []) if item in {"rdkit_etkdg", "psp", "polyply"}]
    if not engines:
        engines = [str(plan.get("engine", "rdkit_etkdg"))]
    seeds = _as_int(plan.get("independent_seeds"), 1, minimum=1, maximum=999)
    conformers = _as_int(plan.get("conformers_per_composition"), 3, minimum=1, maximum=32)
    model_kind = "surface_patch" if plan.get("model_kind") == "surface_patch" else "molecule"
    pH = float(plan.get("pH", 7.0))
    temperature_K = float(plan.get("temperature_K", 310.15))
    if not polymers or not repeats:
        raise ValueError("Select at least one polymer and oligomer size for the model matrix")
    models: list[dict] = []
    for polymer in polymers:
        for repeat in repeats:
            for engine in engines:
                for replica in range(seeds):
                    label = f"{prefix} Build {polymer} {repeat}-mer · {engine} · seed {replica + 1}/{seeds}".strip()
                    _experiment_event(job_id, stage="model", label=label, status="running")
                    try:
                        model = build_model({"family": "polymer", "model_kind": model_kind, "engine": engine, "pH": pH, "temperature_K": temperature_K, "microstate": "auto", "polymer": polymer, "repeats": repeat, "conformers": conformers, "seed": 20260826 + replica})
                        models.append(model)
                        _experiment_event(job_id, stage="model", label=label, status="completed", model_id=model["id"])
                    except (ValueError, RuntimeError) as exc:
                        _experiment_event(job_id, stage="model", label=label, status="failed", reason=str(exc))
    return models


def _run_pose_generation(job_id: str, surface_ids: list[str], pose_plan: dict, *, prefix: str = "") -> list[dict]:
    adsorbate_ids = [str(item) for item in pose_plan.get("adsorbate_model_ids", [])]
    if not surface_ids or not adsorbate_ids:
        raise ValueError("Select at least one polymer structure and one adsorbate for pose generation")
    site_count = _as_int(pose_plan.get("site_count"), 5, minimum=1, maximum=25)
    orientations = _as_int(pose_plan.get("orientations_per_site"), 4, minimum=1, maximum=8)
    strategy = str(pose_plan.get("strategy", "surface_scan"))
    if strategy not in {"surface_scan", "systematic", "crest_assisted", "manual_seed"}:
        strategy = "surface_scan"
    pose_sets: list[dict] = []
    for surface_id in surface_ids:
        for adsorbate_id in adsorbate_ids:
            label = f"{prefix} Generate poses · {surface_id[-8:]} + {adsorbate_id[-8:]}".strip()
            _experiment_event(job_id, stage="pose", label=label, status="running")
            try:
                pose_set = create_pose_set(surface_model_id=surface_id, adsorbate_model_id=adsorbate_id, site_count=site_count, orientations_per_site=orientations, distance_A=float(pose_plan.get("distance_A", 3.2)), probe_radius_A=float(pose_plan.get("probe_radius_A", 1.4)), seed=_as_int(pose_plan.get("seed"), 20260827, minimum=1, maximum=2_147_483_647), strategy=strategy)
                pose_sets.append(pose_set)
                first_pose_id = next((pose.get("model_id") for pose in pose_set.get("poses", []) if pose.get("model_id")), None)
                _experiment_event(job_id, stage="pose", label=label, status="completed", pose_set_id=pose_set["id"], model_id=first_pose_id)
            except (FileNotFoundError, ValueError, RuntimeError) as exc:
                _experiment_event(job_id, stage="pose", label=label, status="failed", reason=str(exc))
    return pose_sets


def _load_pose_sets(job_id: str, ids: list[str], *, prefix: str = "") -> list[dict]:
    pose_sets: list[dict] = []
    for pose_set_id in ids:
        try:
            pose_sets.append(get_pose_set(pose_set_id))
        except FileNotFoundError:
            _experiment_event(job_id, stage="pose", label=f"{prefix} Load pose set · {pose_set_id[-8:]}", status="failed", reason="Pose set not found")
    return pose_sets


def _run_composed_experiment(job_id: str, items: list[dict]) -> None:
    for index, item in enumerate(items):
        kind = str(item.get("kind", ""))
        prefix = _work_prefix(index, str(item.get("title", "")).strip())
        try:
            if kind == "model_matrix":
                _run_model_matrix(job_id, item.get("model_plan", {}), prefix=prefix)
            elif kind == "pose_generation":
                _run_pose_generation(job_id, [str(value) for value in item.get("surface_model_ids", [])], item.get("pose_plan", {}), prefix=prefix)
            elif kind == "end_to_end":
                surfaces = _run_model_matrix(job_id, item.get("model_plan", {}), prefix=prefix)
                pose_sets = _run_pose_generation(job_id, [model["id"] for model in surfaces], item.get("pose_plan", {}), prefix=prefix)
                _run_experiment_xtb(job_id, pose_sets, item.get("simulation_plan", {}), prefix=prefix)
            elif kind == "end_to_end_ready":
                pose_sets = _run_pose_generation(job_id, [str(value) for value in item.get("surface_model_ids", [])], item.get("pose_plan", {}), prefix=prefix)
                _run_experiment_xtb(job_id, pose_sets, item.get("simulation_plan", {}), prefix=prefix)
            elif kind == "xtb_ready":
                pose_sets = _load_pose_sets(job_id, [str(value) for value in item.get("pose_set_ids", [])], prefix=prefix)
                if not pose_sets:
                    raise ValueError("Select at least one ready pose set for the xTB work item")
                _run_experiment_xtb(job_id, pose_sets, item.get("simulation_plan", {}), prefix=prefix)
            else:
                raise ValueError(f"Unknown experiment work item: {kind or 'missing kind'}")
        except (ValueError, RuntimeError) as exc:
            _experiment_event(job_id, stage="workflow", label=prefix, status="failed", reason=str(exc))


def _run_experiment_job(job_id: str, request: ExperimentRunRequest) -> None:
    with _EXPERIMENT_JOBS_LOCK:
        _EXPERIMENT_JOBS[job_id]["status"] = "running"
    try:
        if request.work_items:
            _run_composed_experiment(job_id, request.work_items)
        elif request.mode in {"models_only", "end_to_end"}:
            surfaces = _run_model_matrix(job_id, request.model_plan)
            if request.mode == "end_to_end":
                _run_experiment_xtb(job_id, _run_pose_generation(job_id, [model["id"] for model in surfaces], request.pose_plan), request.simulation_plan)
        elif request.mode == "ready_pose_sets":
            _run_experiment_xtb(job_id, _load_pose_sets(job_id, [str(item) for item in request.simulation_plan.get("pose_set_ids", [])]), request.simulation_plan)
    except (ValueError, RuntimeError) as exc:
        _experiment_event(job_id, stage="workflow", label="Workflow validation", status="failed", reason=str(exc))
    finally:
        with _EXPERIMENT_JOBS_LOCK:
            _EXPERIMENT_JOBS[job_id]["status"] = "completed"


def _run_experiment_xtb(job_id: str, pose_sets: list[dict], plan: dict, *, prefix: str = "") -> None:
    representatives = _as_int(plan.get("representatives_per_set"), 8, minimum=1, maximum=100)
    solvent = "water" if plan.get("solvent", "water") == "water" else "none"
    optimization_mode = "surface_fixed" if plan.get("optimization_mode", "surface_fixed") == "surface_fixed" else "full"
    timeout = _as_int(plan.get("timeout_seconds"), 300, minimum=30, maximum=1800)
    for pose_set in pose_sets:
        poses = (pose_set.get("poses") or [])[:representatives]
        for pose in poses:
            model_id = pose.get("model_id")
            if not model_id:
                continue
            label = f"{prefix} GFN2-xTB · {pose_set['id'][-8:]} / {pose.get('pose_id', model_id[-8:])}".strip()
            _experiment_event(job_id, stage="xtb", label=label, status="running")
            try:
                result = run_xtb(model_id, solvent=solvent, optimize=True, timeout_seconds=timeout, simulation_type="optimization", charge=pose_set.get("formal_charge"), uhf=0, optimization_mode=optimization_mode)
                if result.get("status") == "completed":
                    _experiment_event(job_id, stage="xtb", label=label, status="completed", model_id=model_id, pose_set_id=pose_set["id"], optimized=bool("xtbopt.xyz" in result.get("files", [])))
                else:
                    _experiment_event(job_id, stage="xtb", label=label, status="failed", model_id=model_id, pose_set_id=pose_set["id"], reason=result.get("reason", "xTB failed"))
            except (FileNotFoundError, ValueError, RuntimeError) as exc:
                _experiment_event(job_id, stage="xtb", label=label, status="failed", model_id=model_id, pose_set_id=pose_set["id"], reason=str(exc))


@app.post("/api/experiments/run", status_code=202)
def run_experiment(request: ExperimentRunRequest) -> dict:
    job_id = f"experiment_{uuid.uuid4().hex[:12]}"
    record = {"job_id": job_id, "status": "queued", "mode": request.mode, "total": 0, "completed": 0, "events": [], "latest_model_id": None, "latest_pose_set_id": None}
    with _EXPERIMENT_JOBS_LOCK:
        _EXPERIMENT_JOBS[job_id] = record
    _EXPERIMENT_EXECUTOR.submit(_run_experiment_job, job_id, request)
    return dict(record)


@app.get("/api/experiments/jobs/{job_id}")
def experiment_job(job_id: str) -> dict:
    with _EXPERIMENT_JOBS_LOCK:
        record = _EXPERIMENT_JOBS.get(job_id)
        if record is None:
            raise HTTPException(status_code=404, detail="Experiment job not found")
        return dict(record)


@app.post("/api/pose-sets", status_code=201)
def make_pose_set(request: PoseSetRequest) -> dict:
    try:
        return create_pose_set(**request.model_dump())
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Surface or adsorbate model not found") from exc
    except (ValueError, RuntimeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/api/pose-sets/{pose_set_id}")
def pose_set(pose_set_id: str) -> dict:
    try:
        return get_pose_set(pose_set_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Pose set not found") from exc


@app.get("/api/pose-sets/{pose_set_id}/overlay")
def pose_set_overlay(
    pose_set_id: str,
    offset: int = Query(default=0, ge=0),
    limit: int = Query(default=8, ge=1, le=100),
    pose_ids: str | None = Query(default=None, description="Comma-separated pose or model IDs to overlay"),
) -> Response:
    try:
        requested = [item.strip() for item in (pose_ids or "").split(",") if item.strip()]
        content = overlay_pdb(pose_set_id, offset=offset, limit=limit, pose_ids=requested or None)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Pose page not found") from exc
    return Response(content=content, media_type="chemical/x-pdb")


@app.get("/api/xtb/status")
def xtb_status() -> dict:
    return get_xtb_status()


@app.get("/api/xtb/results")
def xtb_results(limit: int = Query(default=100, ge=1, le=500)) -> list[dict]:
    return list_xtb_results(limit)


@app.get("/api/analysis/cohorts")
def analysis_cohorts() -> list[dict]:
    return []


@app.get("/api/analysis/cohorts/{cohort_id}")
def analysis_cohort(cohort_id: str) -> dict:
    raise HTTPException(status_code=404, detail="No analysis cohort imported")


@app.get("/api/models/{model_id}/xtb/result")
def xtb_result(model_id: str) -> dict:
    try:
        get_model(model_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Model not found") from exc
    result = get_xtb_result(model_id)
    if result is None:
        raise HTTPException(status_code=404, detail="No xTB result recorded for this model")
    return result


def _optimized_xyz_as_pdb(model_id: str, xyz_path: Path) -> str:
    """Return xTB coordinates in PDB form while retaining the input topology.

    NGL 2.0 (the bundled browser viewer) does not ship an XYZ parser.  xTB
    writes optimized coordinates as XYZ, whereas the saved model PDB also
    contains the atom records and CONECT topology needed for a faithful view.
    Replacing only the coordinate columns gives the viewer an interoperable
    representation without changing the raw xTB artifact on disk.
    """
    structure_path = model_file(model_id, "structure.pdb")
    xyz_lines = xyz_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not xyz_lines:
        raise HTTPException(status_code=404, detail="Optimized coordinates are empty")
    try:
        atom_count = int(xyz_lines[0].strip())
    except ValueError as exc:
        raise HTTPException(status_code=422, detail="Invalid optimized XYZ atom count") from exc

    coordinates: list[tuple[float, float, float]] = []
    for line in xyz_lines[2:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            coordinates.append((float(fields[1]), float(fields[2]), float(fields[3])))
        except ValueError:
            continue
    if atom_count != len(coordinates):
        raise HTTPException(status_code=422, detail="Optimized XYZ atom count does not match its coordinates")

    pdb_lines = structure_path.read_text(encoding="utf-8", errors="replace").splitlines()
    atom_indices = [index for index, line in enumerate(pdb_lines) if line.startswith(("ATOM  ", "HETATM"))]
    if len(atom_indices) != atom_count:
        raise HTTPException(status_code=422, detail="Optimized coordinates do not match the model topology")
    for index, (x, y, z) in zip(atom_indices, coordinates):
        line = pdb_lines[index]
        # PDB Cartesian coordinates occupy columns 31–54 (1-indexed).
        padded = line.ljust(54)
        pdb_lines[index] = f"{padded[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{padded[54:]}"
    return "\n".join(pdb_lines) + "\n"


@app.get("/api/models/{model_id}/xtb/output")
def xtb_output(model_id: str, file: Literal["xtb.out", "xtbopt.xyz", "xtbopt.pdb", "charges", "input.xyz"] = "xtb.out") -> Response:
    """Serve a whitelisted raw xTB artifact for inspection in the Analysis UI."""

    try:
        get_model(model_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Model not found") from exc
    if file == "xtbopt.pdb":
        xyz_path = (RUNS_DIR / model_id / "xtb" / "xtbopt.xyz").resolve()
        runs_root = RUNS_DIR.resolve()
        if runs_root not in xyz_path.parents or not xyz_path.is_file():
            raise HTTPException(status_code=404, detail="Optimized xTB coordinates not found")
        return Response(content=_optimized_xyz_as_pdb(model_id, xyz_path), media_type="chemical/x-pdb")

    path = (RUNS_DIR / model_id / "xtb" / file).resolve()
    runs_root = RUNS_DIR.resolve()
    if runs_root not in path.parents or not path.is_file():
        raise HTTPException(status_code=404, detail="xTB output file not found")
    media_type = "chemical/x-xyz" if file in {"input.xyz", "xtbopt.xyz"} else "text/plain"
    return FileResponse(path, media_type=media_type, filename=f"{model_id}_{file}")


@app.get("/api/catalog")
def catalog() -> dict:
    return public_catalog()


@app.get("/api/models")
def models(limit: int = Query(default=30, ge=1, le=500)) -> list[dict]:
    return list_models(limit)


@app.post("/api/models", status_code=201)
def create_model(request: ModelRequest) -> dict:
    try:
        return build_model(request.model_dump())
    except (ValueError, RuntimeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/api/models/{model_id}")
def model(model_id: str) -> dict:
    try:
        return get_model(model_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Model not found") from exc


@app.get("/api/models/{model_id}/structure")
def structure(model_id: str, format: Literal["pdb", "sdf", "xyz"] = "pdb") -> FileResponse:
    filename = f"structure.{format}"
    media_types = {"pdb": "chemical/x-pdb", "sdf": "chemical/x-mdl-sdfile", "xyz": "chemical/x-xyz"}
    try:
        path = model_file(model_id, filename)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Structure not found") from exc
    return FileResponse(path, media_type=media_types[format], filename=f"{model_id}.{format}")


@app.post("/api/models/{model_id}/xtb")
def xtb(model_id: str, request: XTBRequest) -> dict:
    try:
        # Validate the model id before creating an output directory.
        get_model(model_id)
        optimize = request.optimize if request.optimize is not None else request.simulation_type == "optimization"
        return run_xtb(model_id, solvent=request.solvent, optimize=optimize, timeout_seconds=request.timeout_seconds,
                       simulation_type=request.simulation_type, charge=request.charge, uhf=request.uhf,
                       optimization_mode=request.optimization_mode)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail="Model not found") from exc
    except (ValueError, RuntimeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


def _run_xtb_job(job_id: str, request: XTBJobRequest) -> None:
    with _XTB_JOBS_LOCK:
        job = _XTB_JOBS[job_id]
        job["status"] = "running"
    results = []
    for index, model_id in enumerate(request.model_ids, start=1):
        try:
            optimize = request.optimize if request.optimize is not None else request.simulation_type == "optimization"
            result = run_xtb(model_id, solvent=request.solvent, optimize=optimize,
                             timeout_seconds=request.timeout_seconds, simulation_type=request.simulation_type,
                             charge=request.charge, uhf=request.uhf,
                             optimization_mode=request.optimization_mode)
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            result = {"status": "failed", "model_id": model_id, "reason": str(exc)}
        results.append(result)
        with _XTB_JOBS_LOCK:
            job["completed"] = index
            job["results"] = results
    with _XTB_JOBS_LOCK:
        job["status"] = "completed"


@app.post("/api/xtb/jobs", status_code=202)
def create_xtb_job(request: XTBJobRequest) -> dict:
    # Validate all targets before creating a background job.
    for model_id in request.model_ids:
        try:
            get_model(model_id)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=f"Model not found: {model_id}") from exc
    job_id = f"xtb_{uuid.uuid4().hex[:12]}"
    record = {
        "job_id": job_id,
        "status": "queued",
        "total": len(request.model_ids),
        "completed": 0,
        "model_ids": request.model_ids,
        "simulation_type": request.simulation_type,
        "optimization_mode": request.optimization_mode,
        "solvent": request.solvent,
        "results": [],
    }
    with _XTB_JOBS_LOCK:
        _XTB_JOBS[job_id] = record
    _XTB_EXECUTOR.submit(_run_xtb_job, job_id, request)
    return record


@app.get("/api/xtb/jobs/{job_id}")
def xtb_job(job_id: str) -> dict:
    with _XTB_JOBS_LOCK:
        record = _XTB_JOBS.get(job_id)
        if record is None:
            raise HTTPException(status_code=404, detail="xTB job not found")
        return dict(record)


def _run_crest_job(job_id: str, request: CrestJobRequest) -> None:
    with _CREST_JOBS_LOCK:
        _CREST_JOBS[job_id]["status"] = "running"
    results = []
    for index, model_id in enumerate(request.model_ids, start=1):
        try:
            result = run_crest_refinement(model_id, top_k=request.top_k, solvent=request.solvent,
                                          timeout_seconds=request.timeout_seconds)
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            result = {"status": "failed", "model_id": model_id, "reason": str(exc)}
        results.append(result)
        with _CREST_JOBS_LOCK:
            job = _CREST_JOBS[job_id]
            job["completed"] = index
            job["results"] = results
    with _CREST_JOBS_LOCK:
        _CREST_JOBS[job_id]["status"] = "completed"


@app.post("/api/crest/jobs", status_code=202)
def create_crest_job(request: CrestJobRequest) -> dict:
    if get_crest_status().get("status") != "ready":
        raise HTTPException(status_code=422, detail=get_crest_status().get("reason", "CREST unavailable"))
    for model_id in request.model_ids:
        try:
            model = get_model(model_id)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=f"Model not found: {model_id}") from exc
        if (model.get("source") or {}).get("family") != "complex":
            raise HTTPException(status_code=422, detail=f"CREST target is not a complex pose: {model_id}")
    job_id = f"crest_{uuid.uuid4().hex[:12]}"
    record = {"job_id": job_id, "status": "queued", "total": len(request.model_ids), "completed": 0,
              "model_ids": request.model_ids, "top_k": request.top_k, "solvent": request.solvent, "results": []}
    with _CREST_JOBS_LOCK:
        _CREST_JOBS[job_id] = record
    # Share the single simulation executor with xTB so CPU-heavy methods do
    # not start concurrently on this workstation.
    _XTB_EXECUTOR.submit(_run_crest_job, job_id, request)
    return record


@app.get("/api/crest/jobs/{job_id}")
def crest_job(job_id: str) -> dict:
    with _CREST_JOBS_LOCK:
        record = _CREST_JOBS.get(job_id)
        if record is None:
            raise HTTPException(status_code=404, detail="CREST job not found")
        return dict(record)
