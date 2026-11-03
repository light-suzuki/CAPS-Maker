"""Planning a sparse or enzyme-site-dense CAPS marker panel.

Discover restriction polymorphisms in short local windows first. This avoids
running full-genome BLAST for every variant in a large interval. The resulting
table is a *planning* view; final bands must be computed on designed amplicons.
"""
from __future__ import annotations

from dataclasses import dataclass
from html import escape
from typing import Dict, List, Optional, Sequence, Tuple

from primerblast_oss.caps import (
    ENZYME_METADATA, caps_scan, dcaps_candidates,
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


PANEL_MODES = ("全酵素サイト", "一定間隔", "指定本数")


def _locus_key(site: MarkerSite) -> Tuple[str, str]:
    """One polymorphism is one map marker, whatever the enzyme options."""
    return site.chrom, site.source_variant


def group_by_locus(sites: Sequence[MarkerSite]) -> List[Dict]:
    """Collapse enzyme/CAPS/dCAPS options at one polymorphism into one locus."""
    groups: Dict[Tuple[str, str], List[MarkerSite]] = {}
    for site in sites:
        groups.setdefault(_locus_key(site), []).append(site)
    loci = []
    for (chrom, variant), options in groups.items():
        options = sorted(options, key=lambda s: (
            0 if s.kind == "CAPS" else 1, s.enzyme, s.site_position))
        loci.append({
            "chrom": chrom,
            "variant": variant,
            "position": options[0].position,
            "options": options,
            "has_caps": any(s.kind == "CAPS" for s in options),
        })
    loci.sort(key=lambda item: (item["chrom"], item["position"], item["variant"]))
    return loci


def choose_sites(sites: Sequence[MarkerSite], mode: str,
                 region_start: int, region_end: int,
                 spacing_bp: int = 10000,
                 target_count: int = 20,
                 limit: int = 100,
                 min_spacing_bp: int = 0,
                 search_radius_bp: Optional[int] = None,
                 scores: Optional[Dict[str, float]] = None,
                 ) -> Tuple[List[MarkerSite], Dict]:
    """Select a marker panel without inventing markers for empty bins.

    Grid modes pick at most one *polymorphic locus* per bin (enzyme options at
    the same variant are alternatives, not independent map markers), honour a
    minimum physical spacing between selected loci and a search radius around
    each bin centre, and report every uncovered bin with its reason.
    """
    if region_end < region_start or limit < 1:
        raise ValueError("区間または設計上限が不正です。")
    if mode not in PANEL_MODES:
        raise ValueError("マーカーパネルの設計方式が不明です。")
    if spacing_bp < 1 or target_count < 1:
        raise ValueError("間隔と目標本数は正の整数にしてください。")
    if min_spacing_bp < 0 or (search_radius_bp is not None and search_radius_bp < 0):
        raise ValueError("最小間隔と探索半径は0以上にしてください。")
    scores = scores or {}
    total = region_end - region_start + 1
    inside = [s for s in sites if region_start <= s.position <= region_end]
    loci = group_by_locus(inside)

    def rank(site: MarkerSite, centre: float) -> Tuple:
        # CAPS avoids engineered mismatches; then assay score; then centrality.
        return (0 if site.kind == "CAPS" else 1,
                -float(scores.get(site.identifier, 0.0)),
                abs(site.position - centre), site.enzyme, site.site_position)

    if mode == "全酵素サイト":
        kept: List[MarkerSite] = []
        last = None
        for locus in loci:
            if (min_spacing_bp and last is not None
                    and locus["position"] - last < min_spacing_bp):
                continue
            last = locus["position"]
            kept.extend(locus["options"])
        chosen = sorted(kept, key=lambda s: (s.position, s.site_position,
                                             s.enzyme, s.kind))
        selected = chosen[:limit]
        return selected, {
            "discovered": len(sites),
            "loci": len(loci),
            "panel_candidates": len(chosen),
            "selected": len(selected),
            "selected_loci": len({_locus_key(s) for s in selected}),
            "truncated": len(chosen) > limit,
            "spacing_bp": None,
            "min_spacing_bp": min_spacing_bp,
            "uncovered_bins": [],
            "notice": "Exhaustive mode lists every enzyme option; "
                      "count loci, not rows, as map markers.",
        }

    width = spacing_bp if mode == "一定間隔" else max(1, -(-total // target_count))
    n_bins = -(-total // width)
    by_bin: Dict[int, List[Dict]] = {}
    for locus in loci:
        by_bin.setdefault((locus["position"] - region_start) // width, []).append(locus)

    chosen: List[MarkerSite] = []
    taken: List[int] = []
    uncovered = []
    for index in range(n_bins):
        bin_start = region_start + index * width
        bin_end = min(region_end, bin_start + width - 1)
        centre = (bin_start + bin_end) / 2
        options = by_bin.get(index, [])
        reason = None
        if not options:
            reason = "多型サイトなし"
        elif search_radius_bp is not None:
            options = [l for l in options
                       if abs(l["position"] - centre) <= search_radius_bp]
            if not options:
                reason = "探索半径内にサイトなし"
        if reason is None:
            options = [l for l in options
                       if all(abs(l["position"] - p) >= min_spacing_bp for p in taken)]
            if not options:
                reason = "最小間隔の制約で選択不可"
        if reason is not None:
            uncovered.append({"bin": index + 1, "start": bin_start,
                              "end": bin_end, "reason": reason})
            continue
        best = min((min(l["options"], key=lambda s: rank(s, centre)) for l in options),
                   key=lambda s: rank(s, centre))
        chosen.append(best)
        taken.append(best.position)
    chosen.sort(key=lambda site: (site.position, site.enzyme, site.kind))
    selected = chosen[:limit]
    return selected, {
        "discovered": len(sites),
        "loci": len(loci),
        "panel_candidates": len(chosen),
        "selected": len(selected),
        "selected_loci": len(selected),
        "truncated": len(chosen) > limit,
        "spacing_bp": width,
        "bins": n_bins,
        "min_spacing_bp": min_spacing_bp,
        "search_radius_bp": search_radius_bp,
        "uncovered_bins": uncovered,
        "notice": "No extra markers are invented for bins with zero assayable sites.",
    }


def marker_track_svg(sites: Sequence[MarkerSite], region_start: int,
                     region_end: int, uncovered: Sequence[Dict] = (),
                     width: int = 900) -> str:
    """Physical marker track: selected loci as ticks, uncovered bins shaded."""
    span = max(1, region_end - region_start)
    left, right = 20, width - 20

    def x(pos: float) -> float:
        return left + (pos - region_start) / span * (right - left)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="90" '
        f'viewBox="0 0 {width} 90" font-family="sans-serif" font-size="11">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    for item in uncovered:
        x0, x1 = x(item["start"]), x(item["end"])
        parts.append(
            f'<rect x="{x0:.1f}" y="30" width="{max(1.0, x1 - x0):.1f}" height="20" '
            f'fill="#f3d6d6"><title>{escape(str(item["reason"]))}</title></rect>')
    parts.append(f'<line x1="{left}" y1="40" x2="{right}" y2="40" '
                 'stroke="#333" stroke-width="2"/>')
    for site in sites:
        colour = "#2f6f4f" if site.kind == "CAPS" else "#b5651d"
        px = x(site.position)
        parts.append(
            f'<line x1="{px:.1f}" y1="26" x2="{px:.1f}" y2="54" stroke="{colour}" '
            f'stroke-width="2"><title>{escape(site.identifier)}</title></line>')
    parts.append(f'<text x="{left}" y="72">{region_start:,}</text>')
    parts.append(f'<text x="{right}" y="72" text-anchor="end">{region_end:,}</text>')
    parts.append(f'<text x="{left}" y="16">緑=CAPS 橙=dCAPS 赤帯=未カバー区画'
                 f'（{len(sites)}マーカー）</text>')
    parts.append("</svg>")
    return "".join(parts)
