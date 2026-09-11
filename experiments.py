"""Persistent, non-executing experiment plans for the modeling workspace.

Plans deliberately record the intended matrix before any expensive geometry
generation or xTB work is launched.  This keeps batch design reproducible and
makes the later queue implementation able to consume an explicit manifest.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4


from .builder import ROOT
EXPERIMENTS_DIR = ROOT / "experiments"


def _safe_title(value: str) -> str:
    title = " ".join(value.split()).strip()
    if not title:
        raise ValueError("Experiment title is required")
    return title[:120]


def _summary(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": record["id"],
        "title": record["title"],
        "stage": record["stage"],
        "created_at": record["created_at"],
        "status": record.get("status", "planned"),
    }


def create_experiment(*, title: str, stage: str, design: dict[str, Any], notes: str = "") -> dict[str, Any]:
    EXPERIMENTS_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc)
    slug = re.sub(r"[^a-z0-9]+", "-", _safe_title(title).lower()).strip("-")[:42] or "experiment"
    record = {
        "id": f"exp_{timestamp.strftime('%Y%m%d_%H%M%S')}_{slug}_{uuid4().hex[:6]}",
        "title": _safe_title(title),
        "stage": stage,
        "status": "planned",
        "created_at": timestamp.isoformat(),
        "notes": notes.strip()[:1000],
        "design": design,
        "execution": {
            "state": "not_started",
            "note": "Planning manifest only; no structures or xTB jobs were started when this record was created.",
        },
    }
    path = EXPERIMENTS_DIR / f"{record['id']}.json"
    path.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return record


def list_experiments(limit: int = 100) -> list[dict[str, Any]]:
    if not EXPERIMENTS_DIR.exists():
        return []
    records: list[dict[str, Any]] = []
    for path in sorted(EXPERIMENTS_DIR.glob("exp_*.json"), reverse=True):
        try:
            records.append(_summary(json.loads(path.read_text(encoding="utf-8"))))
        except (OSError, json.JSONDecodeError, KeyError):
            continue
        if len(records) >= limit:
            break
    return records


def get_experiment(experiment_id: str) -> dict[str, Any]:
    path = EXPERIMENTS_DIR / f"{experiment_id}.json"
    if not path.is_file():
        raise FileNotFoundError(experiment_id)
    return json.loads(path.read_text(encoding="utf-8"))
