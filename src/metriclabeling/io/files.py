"""Transparent opening of plain / gzip / bzip2 / zip text files."""

from __future__ import annotations

import bz2
import gzip
import io
import zipfile
from pathlib import Path
from typing import IO, Union

PathLike = Union[str, Path]


def open_text(path: PathLike, member: str | None = None) -> IO[str]:
    """Open ``path`` for text reading, decompressing on the fly.

    Supports ``.gz``, ``.bz2`` and ``.zip`` (first member, or ``member`` if
    given).  Decoding errors are replaced so that odd bytes in old database
    dumps never abort a run.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    if suffix == ".bz2":
        return io.TextIOWrapper(bz2.open(path, "rb"), encoding="utf-8", errors="replace")
    if suffix == ".zip":
        zf = zipfile.ZipFile(path)
        names = [n for n in zf.namelist() if not n.endswith("/")]
        if member is None:
            if not names:
                raise ValueError(f"empty zip archive: {path}")
            member = names[0]
        else:
            matches = [n for n in names if member in n]
            if not matches:
                raise KeyError(f"{member!r} not found in {path}")
            member = matches[0]
        return io.TextIOWrapper(zf.open(member), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")
