"""Checksum verification and run metadata for CAP."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

from cap.paths import ProjectPaths


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_frozen_inputs(paths: ProjectPaths) -> None:
    manifest_path = paths.manifest_dir / "frozen_inputs.json"
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)

    mismatches: list[str] = []
    for filename, record in manifest.items():
        file_path = paths.frozen_dir / filename
        paths.require_file(file_path, f"frozen input {filename}")
        observed = sha256_file(file_path)
        expected = record["sha256"]
        if observed != expected:
            mismatches.append(
                f"{filename}: expected {expected}, observed {observed}"
            )

    if mismatches:
        raise ValueError("Frozen input checksum mismatch:\n" + "\n".join(mismatches))


def write_run_provenance(paths: ProjectPaths, stage_dir: Path) -> Path:
    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=paths.root, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        git_commit = None

    payload = {
        "run_id": paths.run_id,
        "created_at_utc": datetime.now(UTC).isoformat(),
        "repository_root": str(paths.root),
        "data_root": str(paths.data_root),
        "git_commit": git_commit,
        "python_version": sys.version,
        "platform": platform.platform(),
    }
    destination = stage_dir / "run_provenance.json"
    destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return destination