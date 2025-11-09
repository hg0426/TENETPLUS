"""
Utility helpers for accessing the reorganised project directories.
"""
"""
Utility helpers for accessing the reorganised project directories.
"""

from pathlib import Path
from typing import Union

from . import input_dir, output_dir, project_root

Pathish = Union[str, Path]


def resolve_input(*parts: Pathish) -> Path:
    """Return a path inside the input directory."""
    return input_dir().joinpath(*map(str, parts))


def resolve_output(*parts: Pathish) -> Path:
    """Return a path inside the output directory."""
    return output_dir().joinpath(*map(str, parts))


def coerce_input_path(path: Pathish) -> Path:
    """
    Return a path for reading, preferring the provided location and falling
    back to the input directory when necessary.
    """
    candidate = Path(path)
    if candidate.exists() or candidate.is_absolute():
        return candidate
    return resolve_input(candidate)


def coerce_output_path(path: Pathish) -> Path:
    """
    Return a path for writing into the output directory when a bare filename
    is provided.
    """
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return resolve_output(candidate)


def locate_file(path: Pathish) -> Path:
    """
    Locate an existing file by checking the provided path first, followed by
    the input and then the output directory.
    """
    candidate = Path(path)
    if candidate.exists():
        return candidate
    in_path = resolve_input(candidate)
    if in_path.exists():
        return in_path
    out_path = resolve_output(candidate)
    if out_path.exists():
        return out_path
    return candidate


def ensure_output_subdir(*parts: Pathish) -> Path:
    """
    Ensure that an output subdirectory exists and return its path.
    """
    path = resolve_output(*parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


__all__ = [
    "project_root",
    "resolve_input",
    "resolve_output",
    "coerce_input_path",
    "coerce_output_path",
    "locate_file",
    "ensure_output_subdir",
]
