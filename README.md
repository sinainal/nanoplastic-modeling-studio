# Nanoplastic Modeling Studio

A local FastAPI workspace for molecular structure generation, adsorption pose preparation, and xTB calculations. NGL provides interactive 3D inspection alongside modeling, simulation, experiment, and analysis pages.

## Quick start

Requires Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
pip install .
nanoplastic-studio
```

Open http://127.0.0.1:8092. The service binds to localhost by default and is intended for trusted local use, not public deployment.

## Workflow

1. Build small molecules or capped polymer oligomers with RDKit; inspect and save their coordinates.
2. Generate surface sites and adsorbate orientations, then inspect selected poses in 3D.
3. Run xTB optimizations on selected structures or prepare sequential experiments.
4. Inspect saved structures, convergence information, energies, and run diagnostics in Analysis.

The command-line interface is available with `python -m modeling.app.cli --help`. Interactive API documentation is at `/docs`.

## Optional engines and storage

xTB and CREST are separate installations. Put their executables on PATH or configure `XTB_BIN` and `CREST_BIN`. External polymer builders are optional integrations and require their own installation; they are not bundled. Missing engines are reported by the application.

Set `MODELING_DATA_DIR` to choose a writable data directory. New standalone installations otherwise use `~/.local/share/nanoplastic-studio`. Keep generated data outside the source tree.

## Scientific scope

This is a research workflow tool, not a validated affinity predictor. A pH label selects a representative protonation state; it does not run constant-pH dynamics. Finite oligomers and surface proxies are not equilibrated bulk nanoplastic particles. Electronic energies must not be interpreted as binding free energies. Comparisons require consistent methods, constraints, solvent settings, and fragment references.

## Development checks

```bash
pip install '.[dev]'
python -m compileall -q .
# With the service running and Google Chrome installed:
python tests/browser_smoke.py
```

This repository contains application code only. Study-specific scripts, datasets, reports, molecular outputs, credentials, and historical research commits are intentionally excluded.

## Third-party software

The bundled NGL browser library is developed by the [NGL project](https://github.com/nglviewer/ngl) and retains its upstream license notices. RDKit, FastAPI, xTB, and CREST are independent projects with their own licenses. No license grant for this application's original code is implied by those dependencies.
