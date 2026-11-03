"""BED and JSON provenance exports for audit-friendly, reproducible reruns.

The manifest separates *what was asked* (configuration), *what it ran with*
(tool/model versions, DB and file fingerprints) and *what came out* (a
content hash of the result rows). The content hash excludes the timestamp so
reruns on pinned references can be compared byte-for-byte.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from typing import Dict, Mapping, Optional, Sequence

from . import __version__

BED_COLUMNS = ("参照染色体", "参照産物開始(1-based)", "参照産物終了(1-based)")


def bed_bytes(rows: Sequence[dict]) -> bytes:
    """BED6 of reference-anchored amplicons (0-based half-open).

    Rows without a reference anchor are skipped rather than given invented
    coordinates. The score column is the design score scaled to 0-1000.
    """
    lines = []
    seen = set()
    for row in rows:
        if not all(row.get(column) not in (None, "") for column in BED_COLUMNS):
            continue
        chrom = str(row["参照染色体"])
        start = int(row["参照産物開始(1-based)"]) - 1
        end = int(row["参照産物終了(1-based)"])
        name = "|".join(str(part) for part in (
            row.get("酵素", ""), row.get("設計方式", ""),
            row.get("特異性判定", ""),
        ))
        key = (chrom, start, end, name)
        if key in seen:
            continue
        seen.add(key)
        score = int(round(10 * max(0.0, min(100.0, float(row.get("総合設計スコア") or 0)))))
        lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t+")
    lines.sort(key=lambda line: (line.split("\t")[0], int(line.split("\t")[1]),
                                 line.split("\t")[3]))
    header = 'track name="CAPS-Maker" description="reference-anchored CAPS/dCAPS amplicons"'
    return ("\n".join([header, *lines]) + "\n").encode("utf-8") if lines else b""


def rows_digest(rows: Sequence[dict]) -> str:
    payload = json.dumps(list(rows), ensure_ascii=False, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_manifest(
    config: Mapping,
    rows: Sequence[dict],
    databases: Sequence[str] = (),
    reference_files: Optional[Mapping[str, Optional[str]]] = None,
    summary: Optional[Mapping] = None,
    now: Optional[datetime] = None,
    strong_hashes: bool = False,
) -> Dict:
    """Configuration, versions, input fingerprints and a result content hash."""
    from primerblast_oss.provenance import (
        db_fingerprint, file_fingerprint, model_versions, tool_versions,
    )

    versions = tool_versions()
    versions["caps_maker"] = __version__
    return {
        "schema": "caps-maker-provenance/1",
        "generated_at": (now or datetime.now(timezone.utc)).isoformat(),
        "config": dict(config),
        "tool_versions": versions,
        "model_versions": model_versions(),
        "databases": [db_fingerprint(db, strong=strong_hashes) for db in databases],
        "reference_files": {
            label: file_fingerprint(path, strong=strong_hashes)
            for label, path in sorted((reference_files or {}).items())
        },
        "summary": dict(summary or {}),
        "result_rows": len(rows),
        "result_sha256": rows_digest(rows),
    }


def manifest_bytes(manifest: Mapping) -> bytes:
    return json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True,
                      default=str).encode("utf-8")
