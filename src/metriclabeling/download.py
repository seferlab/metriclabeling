"""Download the raw data (see ``configs/species.yaml``)."""

from __future__ import annotations

import logging
import shutil
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional

from .dataset import REPO_ROOT, load_species_config, raw_dir

log = logging.getLogger(__name__)


def _fetch(url: str, dest: Path, force: bool = False) -> Path:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        log.info("exists: %s", dest)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    log.info("downloading %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": "metriclabeling-replication"})
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 20)
    tmp.replace(dest)
    return dest


def _extract(archive: Path, target: Path) -> None:
    with tarfile.open(archive) as tf:
        members = [m for m in tf.getmembers() if not Path(m.name).name.startswith("._") and ".DS_Store" not in m.name]
        for m in members:
            if m.name.startswith("/") or ".." in Path(m.name).parts:
                raise ValueError(f"unsafe path in archive: {m.name}")
        tf.extractall(target, members=members)


def download(source: str = "original", data_dir: Path = REPO_ROOT / "data", force: bool = False,
             only: Optional[list] = None) -> Path:
    """Download (and unpack) all files of a data snapshot into ``data/raw/<source>``."""
    cfg = load_species_config()
    src = cfg["sources"][source]
    base = raw_dir(data_dir, source)
    base.mkdir(parents=True, exist_ok=True)
    for name, url in src["urls"].items():
        if only and name not in only:
            continue
        path = _fetch(url, base / name, force)
        if name in src.get("extract", []):
            _extract(path, base)
    return base
