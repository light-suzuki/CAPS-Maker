"""Planning a sparse or enzyme-site-dense CAPS marker panel.

Discover restriction polymorphisms in short local windows first. This avoids
running full-genome BLAST for every variant in a large interval. The resulting
table is a *planning* view; final bands must be computed on designed amplicons.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from primerblast_oss.caps import (
    ENZYME_METADATA, caps_scan, dcaps_candidates, enzyme_records,
)
from .sources import SampleVariant, source_window


@dataclass(frozen=True)
class MarkerSite:
    chrom: str
    position: int  # reference genomic coordinate, 1-based
    enzyme: str
    site_position: int
    kind: str  # CAPS / dCAPS
    source_variant: str

    @property
    def identifier(self) -> str:
        return f"{self.chrom}:{self.site_position}:{self.enzyme}:{self.kind}"


def _event_refpos(window, event, allele: str) -> Optional[int]:
    inverse = window.inverse(allele)
    if event.site_pos in inverse:
        return inverse[event.site_pos]
    # Recognition may start in an inserted allele; return next anchored base,
    # explicitly treated as a reference-neighborhood coordinate.
    for i in range(event.site_pos, event.site_pos + len(event.recognition)):
        if i in inverse:
            return inverse[i]
    return None


def discover_sites(source: str, locator: str, chrom: str,
                   region_start: int, region_end: int,
                   variants: Sequence[SampleVariant],
                   enzymes: Sequence[str],
                   include_dcaps: bool = True,
                   max_snps: int = 1000) -> Tuple[List[MarkerSite], Dict]:
    if not enzymes:
        raise ValueError("解析する制限酵素を1つ以上選んでください。")
    if region_end < region_start:
        raise ValueError("区間の終了座標が開始座標より小さくなっています。")
    diffs = [v for v in variants if v.different]
    if len(diffs) > max_snps:
        raise ValueError(
            f"差のある変異が{len(diffs)}件あります。上限{max_snps}件です。"
            "対象区間を狭くするか、詳細設定で探索上限を変更してください。"
        )
    enz = {name: ENZYME_METADATA[name] for name in enzymes}
    max_recognition = max(len(ENZYME_METADATA[name].recognition) for name in enzymes)
    flank = max(35, max_recognition + 25)
    sites: Dict[str, MarkerSite] = {}
    covered_dcaps = 0
    for var in diffs:
        low = max(region_start, var.pos - flank)
        high = min(region_end, var.pos + max(len(var.reference), flank))
        if low > var.pos or var.pos > high:
            continue
        window = source_window(source, locator, chrom, low, high, variants)
        # Critical: an indel can change amplicon *length* without changing
        # restriction sites. Only site presence/absence is a true CAPS event.
        digests = caps_scan(window.a, window.b, enzymes=enz, gel_min_gap=1)
        natural = False
        for digest in digests:
            a_events = {
                (_event_refpos(window, event, "A"), event.site_strand)
                for event in digest.allele_a_cuts
            }
            b_events = {
                (_event_refpos(window, event, "B"), event.site_strand)
                for event in digest.allele_b_cuts
            }
            a_events.discard((None, "+"))
            a_events.discard((None, "-"))
            b_events.discard((None, "+"))
            b_events.discard((None, "-"))
            delta = a_events ^ b_events
            if not delta:
                continue
            natural = True
            for position, strand in delta:
                if position is None:
                    continue
                site = MarkerSite(
                    chrom, var.pos, digest.enzyme, position, "CAPS",
                    f"{var.pos}:{var.allele_a}>{var.allele_b}",
                )
                sites.setdefault(site.identifier, site)

        # Fall back to SNP-specific dCAPS frames, not arbitrary sequence length
        # differences. Actual primer geometry is checked during design.
        if include_dcaps and not natural and len(var.allele_a) == len(var.allele_b) == 1:
            local = window.ref_to_a.get(var.pos - window.start)
            if local is None:
                continue
            # Need corresponding B base at a single mapped SNP.
            if window.a[local] != var.allele_a:
                continue
            choices = dcaps_candidates(
                window.a, local, var.allele_a, var.allele_b,
                enzymes=enz, max_primer_mismatch=1,
            )
            if choices:
                covered_dcaps += 1
            for item in choices:
                if not item["engineered_changes"]:
                    continue
                site = MarkerSite(
                    chrom, var.pos, item["enzyme"], var.pos, "dCAPS",
                    f"{var.pos}:{var.allele_a}>{var.allele_b}",
                )
                sites.setdefault(site.identifier, site)
    result = sorted(sites.values(), key=lambda item: (
        item.position, 0 if item.kind == "CAPS" else 1, item.enzyme,
        item.site_position))
    return result, {
        "n_differential_variants": len(diffs),
        "n_sites": len(result),
        "n_dcaps_fallback_variants": covered_dcaps,
        "purpose": "discovery_only_not_validated_amplicons",
    }


def choose_sites(sites: Sequence[MarkerSite], mode: str,
                 region_start: int, region_end: int,
                 spacing_bp: int = 10000,
                 target_count: int = 20,
                 limit: int = 100) -> Tuple[List[MarkerSite], Dict]:
    if region_end < region_start or limit < 1:
        raise ValueError("区間または設計上限が不正です。")
    if mode not in ("全酵素サイト", "一定間隔", "指定本数"):
        raise ValueError("マーカーパネルの設計方式が不明です。")
    if spacing_bp < 1 or target_count < 1:
        raise ValueError("間隔と目標本数は正の整数にしてください。")
    total = region_end - region_start + 1
    width = spacing_bp if mode == "一定間隔" else max(1, (total + target_count - 1) // target_count)
    if mode == "全酵素サイト":
        chosen = list(sites)
    else:
        grouped: Dict[int, List[MarkerSite]] = {}
        for site in sites:
            bucket = (site.site_position - region_start) // width
            if bucket < 0:
                continue
            grouped.setdefault(bucket, []).append(site)
        chosen = []
        for bucket, options in sorted(grouped.items()):
            middle = region_start + bucket * width + width // 2
            # CAPS generally avoids engineered mismatches and is preferred.
            selected = min(
                options, key=lambda site: (
                    0 if site.kind == "CAPS" else 1,
                    abs(site.site_position - middle),
                    site.enzyme,
                )
            )
            chosen.append(selected)
    chosen.sort(key=lambda site: (site.site_position, site.enzyme, site.kind))
    return chosen[:limit], {
        "discovered": len(sites),
        "panel_candidates": len(chosen),
        "selected": min(len(chosen), limit),
        "truncated": len(chosen) > limit,
        "spacing_bp": width if mode != "全酵素サイト" else None,
        "notice": "No extra markers are invented for bins with zero assayable sites.",
    }
