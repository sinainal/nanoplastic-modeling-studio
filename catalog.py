from __future__ import annotations

from typing import Any


DOPAMINE_STATES: dict[str, dict[str, Any]] = {
    "DAH_plus": {
        "label": "Dopaminyum (+1)",
        "smiles": "[NH3+]CCc1ccc(O)c(O)c1",
        "charge": 1,
        "description": "Protonated amine, neutral catechol; pH 4/7 baseline.",
    },
    "DA_neutral": {
        "label": "Nötr dopamine",
        "smiles": "NCCc1ccc(O)c(O)c1",
        "charge": 0,
        "description": "Neutral control microstate; not the dominant pH 7 state.",
    },
    "DA_zwitterion": {
        "label": "Zwitteriyon (0)",
        "smiles": "[NH3+]CCc1ccc([O-])c(O)c1",
        "charge": 0,
        "description": "Ammonium plus one catecholate; pH 9-sensitive state.",
    },
    "DA_anion": {
        "label": "Dopamine anyonu (-1)",
        "smiles": "NCCc1ccc([O-])c(O)c1",
        "charge": -1,
        "description": "Neutral amine plus one catecholate; minor pH 9 control.",
    },
}


POLYMERS: dict[str, dict[str, Any]] = {
    "PE": {
        "label": "Polyethylene",
        "repeat_formula": "C2H4",
        "default_repeats": 6,
        "max_repeats": 18,
        "description": "Capped linear alkane pilot; amorphous density is not implied.",
    },
    "PP": {
        "label": "Polypropylene",
        "repeat_formula": "C3H6",
        "default_repeats": 5,
        "max_repeats": 15,
        "description": "Stereochemistry-unspecified PP pilot; tacticity metadata required later.",
    },
    "PET": {
        "label": "Polyethylene terephthalate",
        "repeat_formula": "C10H8O4",
        "default_repeats": 2,
        "max_repeats": 15,
        "description": "Capped PET oligomer with explicit ester-linked terephthalate units.",
    },
    "PS": {
        "label": "Polystyrene",
        "repeat_formula": "C8H8",
        "default_repeats": 4,
        "max_repeats": 15,
        "description": "Stereochemistry-unspecified PS pilot with explicit phenyl groups.",
    },
}


def resolve_dopamine_state(pH: float, requested: str) -> str:
    if requested != "auto":
        if requested not in DOPAMINE_STATES:
            raise ValueError(f"Unknown dopamine microstate: {requested}")
        return requested
    if pH < 8.4:
        return "DAH_plus"
    return "DA_zwitterion"


def polymer_smiles(polymer: str, repeats: int) -> str:
    polymer = polymer.upper()
    if polymer not in POLYMERS:
        raise ValueError(f"Unknown polymer: {polymer}")
    if not 1 <= repeats <= int(POLYMERS[polymer]["max_repeats"]):
        raise ValueError(f"{polymer} repeat count must be 1-{POLYMERS[polymer]['max_repeats']}")
    if polymer == "PE":
        return "C" * (2 * repeats)
    if polymer == "PP":
        return "CC(C)" * repeats
    if polymer == "PS":
        return "CC(c1ccccc1)" * repeats

    # HO-CH2-CH2-O-[CO-Ph-CO-O-CH2-CH2-O]n-H capped PET oligomer.
    chain = "OCCO"
    for ring_index in range(1, repeats + 1):
        # SMILES uses ``%nn`` for aromatic ring labels above 9.  Without this
        # conversion PET models silently become impossible at the requested
        # 15-mer Pilot 2 length.
        ring = str(ring_index) if ring_index < 10 else f"%{ring_index}"
        chain = f"OCCOC(=O)c{ring}ccc(C(=O){chain})cc{ring}"
    return chain


def public_catalog() -> dict[str, Any]:
    # Imported lazily to avoid a circular import during the standalone catalog
    # module tests.  Availability reflects the local modeling environment.
    from .external_engines import engine_status

    return {
        "small_molecules": {
            "dopamine": {
                "label": "Dopamine",
                "microstates": [
                    {"id": state_id, **state}
                    for state_id, state in DOPAMINE_STATES.items()
                ],
                "pH_values": [4.0, 7.0, 9.0],
            }
        },
        "polymers": [{"id": polymer_id, **definition} for polymer_id, definition in POLYMERS.items()],
        "model_kinds": [
            {"id": "molecule", "label": "Tek molekül / oligomer"},
            {"id": "surface_patch", "label": "2×2 yüzey proxy"},
        ],
        "engines": list(engine_status().values()),
    }
