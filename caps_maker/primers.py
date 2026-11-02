"""Primer-BLAST design through the upstream primerblast-oss CLI.

A single shared primer pair must bind both alleles exactly. Genomic off-target
claims come exclusively from the upstream BLAST result, not a local sequence
match. A pair is never silently accepted when the B allele cannot be mapped.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from primerblast_oss.caps import revcomp
from .core import AlleleAlignment, Variant, mapped_interval, scan_caps


def _amplicons_for_pair(pair: dict, a: str, b: str,
                        alignment: AlleleAlignment) -> Tuple[str, str] | None:
    left = int(pair["left_start"])
    left_len = int(pair["left_len"])
    right = int(pair["right_start"])
    right_len = int(pair["right_len"])
    right_left = right - right_len + 1
    if left < 0 or right >= len(a) or right_left <= left:
        return None
    # Primer3's right_start denotes the 5'-most coordinate of the R primer
    # on the forward template (rightmost in the reference).
    forward = str(pair["forward"]).upper()
    reverse_binding = revcomp(str(pair["reverse"]))
    if a[left:left + left_len] != forward:
        return None
    if a[right_left:right + 1] != reverse_binding:
        return None

    left_b = mapped_interval(alignment.a_to_b, left, left + left_len)
    right_b = mapped_interval(alignment.a_to_b, right_left, right + 1)
    if not left_b or not right_b or left_b[1] > right_b[0]:
        return None
    if b[left_b[0]:left_b[1]] != forward:
        return None
    if b[right_b[0]:right_b[1]] != reverse_binding:
        return None
    return a[left:right + 1], b[left_b[0]:right_b[1]]


def design_caps(
    a: str,
    b: str,
    alignment: AlleleAlignment,
    variant: Variant,
    enzymes: Sequence[str],
    databases: Sequence[str],
    product_min: int = 200,
    product_max: int = 800,
    num_pairs: int = 15,
    min_gap: int = 25,
    timeout: int = 600,
) -> Tuple[List[dict], dict]:
    if not databases:
        raise ValueError("Primer-BLASTにはローカルBLAST DBが必要です。")
    if product_min < 70 or product_max < product_min:
        raise ValueError("PCR産物長の範囲が不正です。")
    if num_pairs < 1 or num_pairs > 100:
        raise ValueError("プライマー候補数は1～100にしてください。")
    if not 0 <= variant.position_a < len(a):
        raise ValueError("標的多型の座標が不正です。")

    with tempfile.TemporaryDirectory(prefix="caps-maker-") as directory:
        query = Path(directory) / "allele_a.fa"
        query.write_text(">allele_A\n" + a + "\n", encoding="utf-8")
        command = [
            sys.executable, "-m", "primerblast_oss", "design",
            "--template-fasta", str(query),
            "--target", f"{variant.position_a},1",
            "--product-size", f"{product_min}-{product_max}",
            "--num-return", str(num_pairs),
            "--format", "json",
        ]
        for db in databases:
            command.extend(["--db", db])
        try:
            run = subprocess.run(command, capture_output=True, text=True,
                                 timeout=timeout, check=False)
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError("Primer-BLASTが制限時間内に完了しませんでした。") from exc
        if run.returncode:
            raise RuntimeError("primerblast-oss: " + (run.stderr or run.stdout).strip()[-1500:])
        try:
            data = json.loads(run.stdout)
        except json.JSONDecodeError as exc:
            raise RuntimeError("PrimerBlastOssからJSON形式の結果を取得できませんでした。") from exc

    rows = []
    mapped = 0
    for rank, pair in enumerate(data.get("pairs", []), 1):
        amplicons = _amplicons_for_pair(pair, a, b, alignment)
        if amplicons is None:
            continue
        mapped += 1
        amplicon_a, amplicon_b = amplicons
        evidence = pair.get("specificity") or {}
        digests = scan_caps(amplicon_a, amplicon_b, enzymes, min_gap)
        for digest in digests:
            rows.append({
                "順位": rank,
                "標的多型": variant.label,
                "Forward(5'-3')": pair["forward"],
                "Reverse(5'-3')": pair["reverse"],
                "Tm F(°C)": round(float(pair["tm_f"]), 1),
                "Tm R(°C)": round(float(pair["tm_r"]), 1),
                "A産物(bp)": len(amplicon_a),
                "B産物(bp)": len(amplicon_b),
                "BLAST判定": evidence.get("specificity_status", "未確認"),
                "BLAST探索完全性": evidence.get("search_completeness", "未確認"),
                "オフターゲット数": evidence.get("total_off_target", ""),
                **digest,
            })
    summary = {
        "primer3_pairs": len(data.get("pairs", [])),
        "both_alleles_bind": mapped,
        "caps_candidates": len(rows),
        "primer3_explain": data.get("primer3_explain", ""),
        "databases": list(databases),
    }
    return rows, summary
