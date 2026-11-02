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


def _pair_evidence(per_database: Sequence[dict]) -> List[dict]:
    """Adapt the upstream design JSON into locus-anchored PCR products."""
    return [{
        "db": entry.get("db"),
        "search_completeness": entry.get("search_completeness"),
        "search_complete": entry.get("search_complete"),
        "products": (
            list(entry.get("on_target") or []) + list(entry.get("off_target") or [])
        ),
    } for entry in per_database]


def _qc_row(
    row: dict, *,
    evidence: Sequence[dict],
    databases: Sequence[str],
    window,
    left: int, right: int,
    gel_profile,
    specificity_policy: str,
    fasta_by_db: dict | None,
) -> dict:
    from .validation import gel_specificity, parse_fragment_column

    verdict = gel_specificity(
        evidence=evidence, databases=databases, reference=window,
        left=left, right=right, forward=row["Forward(5'-3')"],
        reverse=row["Reverse(5'-3')"], enzyme=row["酵素"],
        fragments_a=parse_fragment_column(row["A断片(bp)"]),
        fragments_b=parse_fragment_column(row["B断片(bp)"]),
        gel_profile=gel_profile, policy=specificity_policy,
        genome_fasta_by_db=fasta_by_db,
    )
    accepted = bool(verdict["accepted"]) and row["ゲル判別"] == "可能"
    row["特異性判定"] = verdict["status"]
    row["参照検証採用可"] = accepted
    row["特異性詳細"] = verdict.get("note", "") or "; ".join(
        str(e.get("db")) + ":off=" + str(e.get("n_off_target", "?"))
        for e in verdict.get("per_db", [])
    )
    score = float(row.get("実験適合スコア", 0))
    # Not a probability. A high gel score must not override unsafe specificity.
    row["総合設計スコア"] = round(
        0.8 * score + (20 if accepted else 0), 1)
    return row


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
    gel_profile=None,
    dcaps_mode: str = "none",
    dcaps_max_mismatches: int = 1,
    specificity_policy: str = "strict",
    window=None,
    genome_fasta_by_db: dict | None = None,
    max_dcaps_per_pair: int = 4,
    max_dcaps_parent_pairs: int = 2,
) -> Tuple[List[dict], dict]:
    """Design allele-shared CAPS; optionally generate dCAPS with new BLAST checks.

    Acceptance is always conservative: without a mapped reference locus the
    design may be viewed/exported, but it is NOT labeled specificity-validated.
    """
    from .experiments import GelProfile
    from .validation import run_primer_check
    from .dcaps import design_dcaps_for_pair

    if not databases:
        raise ValueError("Primer-BLASTにはローカルBLAST DBが必要です。")
    if product_min < 70 or product_max < product_min:
        raise ValueError("PCR産物長の範囲が不正です。")
    if num_pairs < 1 or num_pairs > 100:
        raise ValueError("プライマー候補数は1～100にしてください。")
    if not 0 <= variant.position_a < len(a):
        raise ValueError("標的多型の座標が不正です。")
    if dcaps_mode not in ("none", "fallback", "always"):
        raise ValueError("dCAPSモードが不正です。")
    if dcaps_max_mismatches not in (1, 2):
        raise ValueError("dCAPSの人工ミスマッチ数は1または2にしてください。")
    if specificity_policy not in ("strict", "gel-tolerant"):
        raise ValueError("特異性モードが不正です。")
    gel_profile = gel_profile or GelProfile()
    gel_profile.validate()

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

    rows: List[dict] = []
    mapped = dcaps_count = 0
    checked: Dict[Tuple[str, str], Sequence[dict]] = {}
    for rank, pair in enumerate(data.get("pairs", []), 1):
        amplicons = _amplicons_for_pair(pair, a, b, alignment)
        if amplicons is None:
            continue
        mapped += 1
        amplicon_a, amplicon_b = amplicons
        evidence = pair.get("specificity") or {}
        evidence_records = _pair_evidence(evidence.get("per_db") or [])
        digests = scan_caps(amplicon_a, amplicon_b, enzymes, min_gap,
                            profile=gel_profile)
        for digest in digests:
            row = {
                "順位": rank,
                "設計方式": "CAPS",
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
            }
            rows.append(_qc_row(
                row, evidence=evidence_records, databases=databases,
                window=window, left=int(pair["left_start"]),
                right=int(pair["right_start"]),
                gel_profile=gel_profile, specificity_policy=specificity_policy,
                fasta_by_db=genome_fasta_by_db,
            ))

        useful_caps = any(d.get("ゲル判別") == "可能" for d in digests)
        if mapped <= max_dcaps_parent_pairs and (
                dcaps_mode == "always" or
                (dcaps_mode == "fallback" and not useful_caps)):
            candidates = design_dcaps_for_pair(
                a, b, alignment, variant, pair, enzymes, min_gap=min_gap,
                max_mismatches=dcaps_max_mismatches,
                max_candidates=max_dcaps_per_pair,
                profile=gel_profile,
            )
            for row in candidates:
                f, rev = row["Forward(5'-3')"], row["Reverse(5'-3')"]
                key = (f, rev)
                if key not in checked:
                    checked[key] = run_primer_check(
                        f, rev, databases, timeout=timeout,
                        max_product=max(10000, product_max + 3000),
                    )
                left = row.pop("_left_start")
                right = row.pop("_right_start")
                row["順位"] = rank
                row["BLAST判定"] = "dCAPS: independent primer check"
                row["BLAST探索完全性"] = ";".join(
                    str(entry.get("search_completeness", "unknown"))
                    for entry in checked[key])
                row["オフターゲット数"] = "see specificity details"
                rows.append(_qc_row(
                    row, evidence=checked[key], databases=databases,
                    window=window, left=left, right=right,
                    gel_profile=gel_profile, specificity_policy=specificity_policy,
                    fasta_by_db=genome_fasta_by_db,
                ))
                dcaps_count += 1

    rows.sort(key=lambda row: (
        not row.get("参照検証採用可", False),
        -float(row.get("総合設計スコア", 0)),
        row.get("順位", 999),
        row.get("酵素", ""),
    ))
    summary = {
        "primer3_pairs": len(data.get("pairs", [])),
        "both_alleles_bind": mapped,
        "caps_candidates": sum(r.get("設計方式") == "CAPS" for r in rows),
        "dcaps_candidates": dcaps_count,
        "total_candidates": len(rows),
        "validated_reference_candidates": sum(
            bool(r.get("参照検証採用可")) for r in rows),
        "primer3_explain": data.get("primer3_explain", ""),
        "databases": list(databases),
        "interpretation": (
            "BLAST is reference-scoped and cannot guarantee wet-PCR outcomes "
            "or the absence of alleles not represented in the selected DB."
        ),
    }
    return rows, summary
