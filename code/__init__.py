"""
TENETPLUS code package.

Provides a central place to locate project directories so that scripts
can reliably reference data after the repository restructure.
"""

from pathlib import Path


def project_root() -> Path:
    """Return the repository root directory."""
    return Path(__file__).resolve().parents[1]


def code_dir() -> Path:
    """Return the path to the code directory."""
    return Path(__file__).resolve().parent


def input_dir() -> Path:
    """Return the path to the input data directory."""
    return project_root() / "input"


def output_dir() -> Path:
    """Return the path to the output directory."""
    return project_root() / "output"


__all__ = ["project_root", "code_dir", "input_dir", "output_dir"]
