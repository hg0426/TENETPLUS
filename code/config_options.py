"""Kernel-only interactive options for TENETPLUS_KERNEL."""

from __future__ import annotations

from typing import Dict, List

PERM_TOGGLE: List[str] = ["off", "on"]
PERM_FDR_TOGGLE: List[str] = ["off", "on"]
LOCAL_TE_TOGGLE: List[str] = ["off", "on"]
TENET_MODE_DETAILS = {
    "0": {
        "name": "TENET_TF (RNA only)",
        "description": "Runs the RNA-only TENET TF pair workflow with kernel TE.",
    },
    "1": {
        "name": "TENET_Plus full",
        "description": "Generates TF->gene, TF->peak, and peak->gene outputs.",
    },
    "2": {
        "name": "TF->gene",
        "description": "Outputs only TF to gene edges.",
    },
    "3": {
        "name": "TF->peak",
        "description": "Outputs only TF to peak edges.",
    },
    "4": {
        "name": "TF->gene + TF->peak",
        "description": "Combines TF->gene and TF->peak outputs.",
    },
    "5": {
        "name": "peak->gene (cis)",
        "description": "Builds peak->gene cis interactions.",
    },
    "6": {
        "name": "peak->peak (cis)",
        "description": "Builds cis peak->peak interactions.",
    },
}  # type: Dict[str, Dict[str, str]]
