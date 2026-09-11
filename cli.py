"""Command-line entry points for the modeling pilot and xTB runs.

Examples::

    python -m modeling.app.cli xtb-run --model-id mdl_...
    python -m modeling.app.cli xtb-batch --polymer PE --repeats 9 15
    python -m modeling.app.cli xtb-dopamine
    python -m modeling.app.cli xtb-results
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .builder import list_models, RUNS_DIR
from .crest import get_crest_status, run_crest_refinement
from .pose_sets import create_pose_set, get_pose_set, list_pose_sets
from .xtb import get_xtb_status, list_xtb_results, run_xtb


def _json(value: object) -> None:
    print(json.dumps(value, indent=2, ensure_ascii=False))


def _run_one(args: argparse.Namespace) -> int:
    result = run_xtb(args.model_id, solvent=args.solvent, optimize=args.type == "optimization",
                     simulation_type=args.type, charge=args.charge, uhf=args.uhf,
                     timeout_seconds=args.timeout, optimization_mode=args.optimization_mode)
    _json(result)
    return 0 if result.get("status") in {"completed", "unavailable"} else 1


def _batch(args: argparse.Namespace) -> int:
    repeats = set(args.repeats or [9, 15])
    models = []
    for model in list_models(500):
        source = model.get("source") or {}
        if source.get("polymer") != args.polymer.upper():
            continue
        if int(source.get("repeats", -1)) not in repeats:
            continue
        if args.engine != "all" and model.get("engine") != args.engine:
            continue
        models.append(model)
    results = []
    for model in models:
        results.append(run_xtb(model["id"], solvent=args.solvent, optimize=args.type == "optimization",
                               simulation_type=args.type, charge=args.charge, uhf=args.uhf,
                               timeout_seconds=args.timeout))
    report = {"polymer": args.polymer.upper(), "repeats": sorted(repeats), "simulation_type": args.type,
              "solvent": args.solvent, "model_count": len(models), "results": results}
    output = RUNS_DIR / "xtb_cli_latest.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    _json({"status": "completed", "report": str(output), "model_count": len(models), "results": results})
    return 0 if all(item.get("status") in {"completed", "unavailable"} for item in results) else 1


def _dopamine_batch(args: argparse.Namespace) -> int:
    """Run the saved dopamine microstate models as a reproducible CLI batch."""

    models = []
    for model in list_models(500):
        source = model.get("source") or {}
        if source.get("compound") != "dopamine":
            continue
        if args.pH is not None and abs(float(source.get("pH", -1)) - args.pH) > 1e-9:
            continue
        if args.microstate != "all" and source.get("microstate") != args.microstate:
            continue
        models.append(model)
    results = []
    for model in models:
        results.append(run_xtb(model["id"], solvent=args.solvent, optimize=args.type == "optimization",
                               simulation_type=args.type, charge=args.charge, uhf=args.uhf,
                               timeout_seconds=args.timeout))
    report = {"compound": "dopamine", "pH": args.pH, "microstate": args.microstate,
              "simulation_type": args.type, "solvent": args.solvent,
              "model_count": len(models), "results": results}
    output = RUNS_DIR / "xtb_dopamine_latest.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    _json({"status": "completed", "report": str(output), "model_count": len(models), "results": results})
    return 0 if all(item.get("status") in {"completed", "unavailable"} for item in results) else 1


def _pose_create(args: argparse.Namespace) -> int:
    result = create_pose_set(
        args.surface_model,
        args.adsorbate_model,
        count=args.count,
        site_count=args.sites,
        orientations_per_site=args.orientations,
        distance_A=args.distance,
        probe_radius_A=args.probe_radius,
        seed=args.seed,
        strategy=args.strategy,
    )
    _json(result)
    return 0


def _pose_xtb(args: argparse.Namespace) -> int:
    pose_set = get_pose_set(args.pose_set_id)
    requested = set(args.poses or [])
    poses = [pose for pose in pose_set["poses"] if not requested or pose["pose_id"] in requested]
    results = [run_xtb(pose["model_id"], solvent=args.solvent, optimize=True,
                       simulation_type="optimization", charge=pose_set["formal_charge"],
                       uhf=args.uhf, timeout_seconds=args.timeout,
                       optimization_mode=args.optimization_mode) for pose in poses]
    _json({"status": "completed", "pose_set_id": args.pose_set_id, "pose_count": len(poses), "results": results})
    return 0 if all(item.get("status") == "completed" for item in results) else 1


def _crest_refine(args: argparse.Namespace) -> int:
    result = run_crest_refinement(args.model_id, top_k=args.top_k, solvent=args.solvent,
                                  timeout_seconds=args.timeout)
    _json(result)
    return 0 if result.get("status") == "completed" else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Nanoplastic modeling/xTB pilot CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    status = sub.add_parser("xtb-status", help="show GFN2-xTB executable status")
    status.set_defaults(handler=lambda _: (_json(get_xtb_status()) or 0))

    run = sub.add_parser("xtb-run", help="run xTB for one saved model")
    run.add_argument("--model-id", required=True)
    run.add_argument("--type", choices=["optimization", "singlepoint"], default="optimization")
    run.add_argument("--solvent", choices=["water", "none"], default="water")
    run.add_argument("--charge", type=int, default=None)
    run.add_argument("--uhf", type=int, default=0)
    run.add_argument("--optimization-mode", choices=["surface_fixed", "full"], default="full")
    run.add_argument("--timeout", type=int, default=300)
    run.set_defaults(handler=_run_one)

    batch = sub.add_parser("xtb-batch", help="run xTB for matching polymer models")
    batch.add_argument("--polymer", default="PE")
    batch.add_argument("--repeats", type=int, nargs="+")
    batch.add_argument("--engine", default="all")
    batch.add_argument("--type", choices=["optimization", "singlepoint"], default="optimization")
    batch.add_argument("--solvent", choices=["water", "none"], default="water")
    batch.add_argument("--charge", type=int, default=None)
    batch.add_argument("--uhf", type=int, default=0)
    batch.add_argument("--timeout", type=int, default=300)
    batch.set_defaults(handler=_batch)

    dopamine = sub.add_parser("xtb-dopamine", help="run xTB for saved dopamine microstate models")
    dopamine.add_argument("--pH", type=float, default=None)
    dopamine.add_argument("--microstate", choices=["all", "DAH_plus", "DA_neutral", "DA_zwitterion", "DA_anion"], default="all")
    dopamine.add_argument("--type", choices=["optimization", "singlepoint"], default="optimization")
    dopamine.add_argument("--solvent", choices=["water", "none"], default="water")
    dopamine.add_argument("--charge", type=int, default=None)
    dopamine.add_argument("--uhf", type=int, default=0)
    dopamine.add_argument("--timeout", type=int, default=300)
    dopamine.set_defaults(handler=_dopamine_batch)

    results = sub.add_parser("xtb-results", help="list persisted xTB results")
    results.add_argument("--limit", type=int, default=100)
    results.set_defaults(handler=lambda args: (_json(list_xtb_results(args.limit)) or 0))

    pose_create = sub.add_parser("pose-create", help="create a systematic polymer–dopamine pose set")
    pose_create.add_argument("--surface-model", required=True)
    pose_create.add_argument("--adsorbate-model", required=True)
    pose_create.add_argument("--count", type=int, default=None, help="legacy total-pose limit")
    pose_create.add_argument("--sites", type=int, default=5)
    pose_create.add_argument("--orientations", type=int, default=4)
    pose_create.add_argument("--distance", type=float, default=3.2)
    pose_create.add_argument("--probe-radius", type=float, default=1.4)
    pose_create.add_argument("--seed", type=int, default=20260827)
    pose_create.add_argument("--strategy", choices=["surface_scan", "systematic", "crest_assisted", "manual_seed"], default="surface_scan")
    pose_create.set_defaults(handler=_pose_create)

    pose_list = sub.add_parser("pose-list", help="list saved pose sets")
    pose_list.add_argument("--limit", type=int, default=50)
    pose_list.set_defaults(handler=lambda args: (_json(list_pose_sets(args.limit)) or 0))

    pose_show = sub.add_parser("pose-show", help="show one saved pose set")
    pose_show.add_argument("--pose-set-id", required=True)
    pose_show.set_defaults(handler=lambda args: (_json(get_pose_set(args.pose_set_id)) or 0))

    pose_xtb = sub.add_parser("pose-xtb", help="run xTB for all or selected poses in a set")
    pose_xtb.add_argument("--pose-set-id", required=True)
    pose_xtb.add_argument("--poses", nargs="+")
    pose_xtb.add_argument("--solvent", choices=["water", "none"], default="water")
    pose_xtb.add_argument("--uhf", type=int, default=0)
    pose_xtb.add_argument("--optimization-mode", choices=["surface_fixed", "full"], default="surface_fixed")
    pose_xtb.add_argument("--timeout", type=int, default=300)
    pose_xtb.set_defaults(handler=_pose_xtb)

    crest_status = sub.add_parser("crest-status", help="show CREST/GFN-FF backend status")
    crest_status.set_defaults(handler=lambda _: (_json(get_crest_status()) or 0))

    crest_refine = sub.add_parser("crest-refine", help="refine one saved complex pose with constrained CREST")
    crest_refine.add_argument("--model-id", required=True)
    crest_refine.add_argument("--top-k", type=int, default=3)
    crest_refine.add_argument("--solvent", choices=["water", "none"], default="water")
    crest_refine.add_argument("--timeout", type=int, default=600)
    crest_refine.set_defaults(handler=_crest_refine)

    args = parser.parse_args()
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
