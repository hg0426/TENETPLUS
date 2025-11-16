"""Enumerations of interactive options for TENET_Plus_for_py.sh."""

from __future__ import annotations

from typing import Dict, List

MODALITY_CHOICES: List[str] = ["rna", "atac", "auto", "none"]
SCREEN_CHOICES: List[str] = [
    "linear",
    "poly",
    "ksg",
    "kernel",
    "gcmi",
    "disc",
    "ordinal",
    "kernel_grid",
]
REFINE_CHOICES: List[str] = ["none", "kernel", "ksg"]
PERM_TOGGLE: List[str] = ["off", "on"]
PERM_FDR_TOGGLE: List[str] = ["off", "on"]
LOCAL_TE_TOGGLE: List[str] = ["off", "on"]
TENET_MODE_DETAILS = {
    "0": {
        "name": "TENET_TF (RNA only)",
        "description": "Runs the original TENET pipeline using only RNA expression data.",
    },
    "1": {
        "name": "TENET_Plus (RNA + ATAC, full matrix)",
        "description": "Generates both TF→gene and TF→peak edges, including peak-source integration.",
    },
    "2": {
        "name": "TENET_Plus rowTF_colGN",
        "description": "Outputs only TF (rows) to gene (columns) edges from TENET_Plus.",
    },
    "3": {
        "name": "TENET_Plus rowTF_colPK",
        "description": "Outputs only TF (rows) to peak (columns) edges from TENET_Plus.",
    },
    "4": {
        "name": "TENET_Plus rowTF_colGN+PK",
        "description": "Combines TF→gene and TF→peak matrices for integrated analysis.",
    },
    "5": {
        "name": "TENET_Plus rowPeak (cis peak-source)",
        "description": "Builds peak→gene cis interactions using peak windows and gene proximity.",
    },
    "6": {
        "name": "TENET_Plus peak→peak (cis)",
        "description": "Builds cis peak→peak interaction network constrained by genomic distance.",
    },
}  # type: Dict[str, Dict[str, str]]
