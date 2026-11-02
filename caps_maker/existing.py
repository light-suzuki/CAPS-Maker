"""Analyze existing PCR primers: use their amplified span as the CAPS interval."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from primerblast_oss.caps import revcomp
from primerblast_oss.dcaps_workflow import _tm

from .core import AlleleAlignment, scan_caps
from .experiments import GelProfile
from .validation import gel_specificity, run_primer_check


def _find_all(sequence: str, needle: str) -> List[int]:
    results = []
    index = -1
    while True:
        index = sequence.find(needle, index + 1)
        if index < 0:
            break
        results.append(index)
    return results


def primer_intervals(template: str, forward: str, reverse: str,
                     min_bp: int = 70, max_bp: int = 4000
                     ) -> List[Tuple[int, int]]:
    """Find exact F/R matches to the same plus-strand amplicon."""
    forward = forward.strip().upper()
    reverse = reverse.strip().upper()
    if not (18 <= len(forward) <= 100 and 18 <= len(reverse) <= 100):
        raise ValueError("プライマー長は18～100塩基にしてください。")
    if set(forward + reverse) - set("ACGT"):
        raise ValueError("プライマーにはACGTの塩基を指定してください。")
    right_binding = revcomp(reverse)
    intervals = []
    for left in _find_all(template, forward):
        for right_start in _find_all(template, right_binding):
            end = right_start + len(reverse) - 1
            size = end - left + 1
            if right_start < left + len(forward) or not min_bp <= size <= max_bp:
                continue
            intervals.append((left, end))
    return sorted(set(intervals))



def locate_reference_product(
    forward: str, reverse: str, databases: Sequence[str],
    chrom: str, interval_start: int, interval_end: int,
    max_product: int = 10000,
) -> Tuple[Tuple[int, int], Sequence[dict]]:
    """Find one F/R product from a complete BLAST search in a reference interval."""
    if not databases:
        raise ValueError("標的座位の検索にはBLAST DBが必要です。")
    evidence = run_primer_check(forward.upper(), reverse.upper(), databases,
                                max_product=max_product)
    primary = evidence[0]
    if primary.get("search_completeness") != "complete":
        raise ValueError("BLAST探索が不完全です。ヒット数上限を上げてください。")
    matches = []
    for product in primary.get("products") or []:
        start, end = sorted((int(product["start"]), int(product["end"])))
        if (product.get("subject") == chrom
                and product.get("orientation") == "F/R"
                and start >= interval_start and end <= interval_end):
            matches.append((start, end))
    matches = sorted(set(matches))
    if len(matches) != 1:
        raise ValueError(
            f"指定区間内でプライマーの予測PCR産物は{len(matches)}件です。"
            "一意の区間に絞るか、プライマーの向きと配列を確認してください。")
    return matches[0], evidence


def analyze_primers(
    a: str, b: str, alignment: AlleleAlignment,
    forward: str, reverse: str,
    enzymes: Sequence[str],
    databases: Sequence[str],
    min_bp: int = 70, max_bp: int = 4000,
    min_gap: int = 25,
    gel_profile: Optional[GelProfile] = None,
    specificity_policy: str = "strict",
    window=None,
    check_evidence: Optional[Sequence[dict]] = None,
) -> Tuple[List[dict], dict]:
    from .primers import _amplicons_for_pair

    if not databases:
        raise ValueError("既存プライマーの特異性検証にはBLAST DBが必要です。")
    gel_profile = gel_profile or GelProfile()
    intervals = primer_intervals(a, forward, reverse, min_bp=min_bp, max_bp=max_bp)
    if not intervals:
        raise ValueError("配列Aに指定したF/Rプライマーが同方向で見つからないか、PCR産物長の条件外です。")
    # Ambiguous products are reported but not auto-endorsed.
    if len(intervals) > 1:
        raise ValueError(f"配列AにPCR候補が{len(intervals)}箇所あります。"
                         "解析する標的を狭い範囲に限定してください。")
    left, right = intervals[0]
    pair = {
        "left_start": left, "left_len": len(forward),
        "right_start": right, "right_len": len(reverse),
        "forward": forward.upper(), "reverse": reverse.upper(),
    }
    products = _amplicons_for_pair(pair, a, b, alignment)
    if products is None:
        raise ValueError("配列Bで両方のプライマーの完全一致結合が確認できません。"
                         "少なくとも片方のプライマーが多型に重なっています。")
    seq_a, seq_b = products
    results = check_evidence if check_evidence is not None else run_primer_check(
        forward.upper(), reverse.upper(), databases,
        max_product=max(10000, max_bp + 3000),
    )
    rows = scan_caps(seq_a, seq_b, enzymes, min_gap, profile=gel_profile)
    for row in rows:
        from .validation import parse_fragment_column
        verdict = gel_specificity(
            evidence=results, databases=databases, reference=window,
            left=left, right=right, forward=forward.upper(),
            reverse=reverse.upper(), enzyme=row["酵素"],
            fragments_a=parse_fragment_column(row["A断片(bp)"]),
            fragments_b=parse_fragment_column(row["B断片(bp)"]),
            gel_profile=gel_profile, policy=specificity_policy,
        )
        good = bool(verdict["accepted"]) and row["ゲル判別"] == "可能"
        row.update({
            "設計方式": "既存プライマー",
            "Forward(5'-3')": forward.upper(),
            "Reverse(5'-3')": reverse.upper(),
            "A産物(bp)": len(seq_a),
            "B産物(bp)": len(seq_b),
            "Tm F(°C)": _tm(forward),
            "Tm R(°C)": _tm(reverse),
            "特異性判定": verdict["status"],
            "参照検証採用可": good,
            "総合設計スコア": round(
                0.8 * float(row["実験適合スコア"]) + (20 if good else 0), 1),
        })
    rows.sort(key=lambda x: (not x["参照検証採用可"], -x["総合設計スコア"]))
    return rows, {
        "n_existing_pairs": 1,
        "caps_candidates": len(rows),
        "amplicon_start_a": left + 1,
        "amplicon_end_a": right + 1,
    }
