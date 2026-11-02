"""Experimental gel profiles, visibility-aware classification, and ranking.

Scores are ordinal design heuristics, not probabilities of wet-lab success.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Dict, Optional, Sequence, Tuple

from primerblast_oss.gel import (
    GEL_RANGES, analyze_digest_patterns, genotype_pattern_discrimination,
)


@dataclass(frozen=True)
class GelProfile:
    ladder: str = "auto"  # auto / 100bp / 1kb / custom
    custom_bands: Tuple[int, ...] = ()
    gel_percent: Optional[float] = None
    min_visible_bp: int = 50
    min_log10_gap: float = 0.04
    max_visible_bp: Optional[int] = None

    def validate(self) -> None:
        if self.ladder not in ("auto", "100bp", "1kb", "custom"):
            raise ValueError("対応していないラダープリセットです。")
        if self.ladder == "custom" and not self.custom_bands:
            raise ValueError("カスタムラダーの断片長を指定してください。")
        if self.gel_percent is not None and self.gel_percent not in GEL_RANGES:
            raise ValueError("対応していないゲル濃度です。")
        if self.min_visible_bp < 1:
            raise ValueError("最小可視断片長は1 bp以上にしてください。")
        if not isfinite(self.min_log10_gap) or not 0 < self.min_log10_gap < 1:
            raise ValueError("分離閾値は0～1の範囲で指定してください。")
        if self.max_visible_bp is not None and self.max_visible_bp < self.min_visible_bp:
            raise ValueError("最大可視断片長は最小値以上にしてください。")


PROFILE_PRESETS = {
    "標準（2%アガロース）": GelProfile(
        ladder="100bp", gel_percent=2.0, min_visible_bp=50, min_log10_gap=0.04),
    "小断片重視（3%アガロース）": GelProfile(
        ladder="100bp", gel_percent=3.0, min_visible_bp=25, min_log10_gap=0.04),
    "長断片重視（1%アガロース）": GelProfile(
        ladder="1kb", gel_percent=1.0, min_visible_bp=150, min_log10_gap=0.04),
    "自動推奨": GelProfile(),
}


def parse_custom_ladder(text: str) -> Tuple[int, ...]:
    try:
        values = sorted({int(part.strip()) for part in text.replace("\n", ",").split(",") if part.strip()})
    except ValueError as exc:
        raise ValueError("ラダーの各断片長は整数で入力してください。") from exc
    if not values or values[0] < 1:
        raise ValueError("ラダーの断片長を正の整数で指定してください。")
    return tuple(values)


def assess_gel(fragments_a: Sequence[int], fragments_b: Sequence[int],
               profile: Optional[GelProfile] = None) -> Dict:
    profile = profile or GelProfile()
    profile.validate()
    standard = analyze_digest_patterns(
        fragments_a, fragments_b,
        ladder=profile.ladder,
        custom_ladder_bands=profile.custom_bands or None,
        gel_percent=profile.gel_percent,
        min_visible_bp=profile.min_visible_bp,
    )
    # The library's basic genotype calls do not filter invisible fragments.
    # Model gel-visible fragments separately before declaring co-dominance.
    percent = standard["gel_percent"]
    gel_min, gel_max = GEL_RANGES[percent]
    lower = max(gel_min, profile.min_visible_bp)
    upper = min(gel_max, profile.max_visible_bp or gel_max)
    visible_a = [int(f) for f in fragments_a if lower <= int(f) <= upper]
    visible_b = [int(f) for f in fragments_b if lower <= int(f) <= upper]
    pattern = genotype_pattern_discrimination(
        visible_a, visible_b, min_log10_gap=profile.min_log10_gap)
    comparisons = pattern["comparisons"]
    aa_ab = comparisons["AA_vs_AB"]["distinguishable"]
    aa_bb = comparisons["AA_vs_BB"]["distinguishable"]
    ab_bb = comparisons["AB_vs_BB"]["distinguishable"]
    if aa_ab and aa_bb and ab_bb:
        kind = "共優性（3遺伝子型を識別）"
    elif aa_bb:
        if aa_ab and not ab_bb:
            kind = "優性型（ABとBBが同じ）"
        if not aa_ab and ab_bb:
            kind = "優性型（AAとABが同じ）"
        kind = "2群識別（ヘテロ判別不可）"
    else:
        kind = "判別困難"

    all_visible = visible_a + visible_b
    visible_fraction = len(all_visible) / max(1, len(fragments_a) + len(fragments_b))
    diagnostic = list((set(visible_a) ^ set(visible_b)))
    min_diag = min(diagnostic) if diagnostic else None
    # Intentionally a transparent quality *score*, never a probability.
    score = (
        0.45 * float(standard["score"])
        + (25 if pattern["distinguishable"] else 8 if aa_bb else 0)
        + 15 * visible_fraction
        + (10 if diagnostic else 0)
        + (5 if len(set(all_visible)) <= 5 else 0)
    )
    score = round(max(0, min(100, score)), 1)
    standard.update({
        "genotype_discrimination": pattern,
        "observed_lanes": pattern["lanes"],
        "visible_range_bp": (lower, upper),
        "visible_fraction": round(visible_fraction, 3),
        "min_diagnostic_visible_bp": min_diag,
        "marker_type": kind,
        "utility_score": score,
        "score_kind": "heuristic_not_success_probability",
    })
    return standard


def off_target_verdict(
    *,
    rows_by_db: Sequence[Dict],
    genomes_by_db: Dict[str, object],
    primers: Dict[str, str],
    enzyme: str,
    a_fragments: Sequence[int],
    b_fragments: Sequence[int],
    profile: GelProfile,
    target_chrom: Optional[str] = None,
    target_interval: Optional[Tuple[int, int]] = None,
    policy: str = "strict",
) -> Dict:
    """Classify specificity using exact intended-locus coordinates when available.

    Gel-tolerant requires complete BLAST search, exhaustive retrieval of every
    off-target product, complete enzyme digestion, and distinguishable AA/AB/BB
    after adding the off-target bands to all lanes. Unmapped or incomplete
    evidence must never be accepted.
    """
    if policy not in ("strict", "gel-tolerant"):
        raise ValueError("特異性ポリシーが不正です。")
    from primerblast_oss.amplicon_digest import digest_offtarget_products

    profile.validate()
    details = []
    for record in rows_by_db:
        products = list(record.get("products") or [])
        observed = []
        background = []
        for p in products:
            start, end = sorted((int(p["start"]), int(p["end"])))
            same_locus = (
                target_chrom is not None and target_interval is not None
                and p.get("subject") == target_chrom
                and abs(start - target_interval[0]) <= 10
                and abs(end - target_interval[1]) <= 10
                and p.get("orientation") in ("F/R", "R/F")
            )
            if same_locus:
                observed.append(p)
            else:
                background.append({**p, "on_target": False})
        located = target_chrom is not None and target_interval is not None
        unique_locus = located and len(observed) == 1
        complete = record.get("search_completeness") == "complete" and record.get("search_complete") is True
        entry = {
            "db": record.get("db"), "n_off_target": len(background),
            "expected_product_found": len(observed),
            "complete": complete, "unique_locus": unique_locus,
            "products": background,
        }
        details.append(entry)

    clean = bool(details) and all(
        item["complete"] and item["unique_locus"] and item["n_off_target"] == 0
        for item in details)
    if clean:
        return {"status": "strict_pass", "accepted": True, "background_bands": 0,
                "per_db": details}
    if policy == "strict":
        return {"status": "not_strict", "accepted": False,
                "background_bands": None, "per_db": details}
    if not details or any(not i["complete"] or not i["unique_locus"] for i in details):
        return {"status": "indeterminate", "accepted": False,
                "background_bands": None, "per_db": details}

    annotated = digest_offtarget_products(
        details, genomes_by_db, primers, enzyme)
    all_pass = True
    for result in annotated:
        dig = result["offtarget_digest"]
        if not dig["complete"]:
            result["gel_result"] = "digest_unresolved"
            all_pass = False
            continue
        background = [int(f) for f in dig["background_fragments"]]
        # Match the user-selected visibility and resolution settings.
        gel = assess_gel(a_fragments, b_fragments, profile)
        minimum, maximum = gel["visible_range_bp"]
        ga = [int(v) for v in a_fragments if minimum <= int(v) <= maximum]
        gb = [int(v) for v in b_fragments if minimum <= int(v) <= maximum]
        bg = [v for v in background if minimum <= v <= maximum]
        pattern = genotype_pattern_discrimination(
            ga, gb, background_sizes=bg,
            min_log10_gap=profile.min_log10_gap,
        )
        result["gel_result"] = "distinguishable" if pattern["distinguishable"] else "ambiguous"
        result["background_bands"] = bg
        if not pattern["distinguishable"]:
            all_pass = False
    return {
        "status": "gel_tolerant_pass" if all_pass else "gel_not_resolvable",
        "accepted": all_pass, "per_db": annotated,
        "background_bands": sum(len(i.get("background_bands", [])) for i in annotated),
    }
