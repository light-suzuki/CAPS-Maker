"""Pure analysis and export helpers for CAPS-Maker.

Restriction-enzyme cleavage geometry is delegated to primerblast-oss. Sequence
alignment maps the reference allele to the second allele without fabricating
coordinates through insertions or deletions.
"""
from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from Bio.Align import PairwiseAligner
from primerblast_oss.caps import ENZYME_METADATA, caps_scan
from .experiments import GelProfile, assess_gel

MAX_SEQUENCE_LENGTH = 8000
IUPAC = frozenset("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class SequenceRecord:
    name: str
    sequence: str


@dataclass(frozen=True)
class Variant:
    position_a: int   # 0-based; insertion uses the adjacent A anchor
    position_b: int
    allele_a: str
    allele_b: str
    kind: str

    @property
    def label(self) -> str:
        a = self.allele_a or "-"
        b = self.allele_b or "-"
        return f"{self.kind}: A:{self.position_a + 1} {a} → {b}"


@dataclass
class AlleleAlignment:
    variants: List[Variant]
    a_to_b: Dict[int, int]


def parse_sequence(text: str) -> SequenceRecord:
    """Accept exactly one FASTA record or one bare DNA sequence."""
    text = text.strip().lstrip("\ufeff")
    if not text:
        raise ValueError("DNA配列を入力してください。")
    if text.startswith(">"):
        records = []
        name, parts = None, []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                header = line[1:].strip().split()
                name = header[0] if header else "sequence"
                parts = []
            elif name is None:
                raise ValueError("FASTAヘッダーの前に配列があります。")
            else:
                parts.append(line)
        records.append((name, "".join(parts)))
        if len(records) != 1:
            raise ValueError("A/Bにそれぞれ1レコードのFASTAを指定してください。")
        name, sequence = records[0]
    else:
        if ">" in text:
            raise ValueError("FASTAヘッダーの位置が不正です。")
        name, sequence = "sequence", text
    sequence = re.sub(r"\s+", "", sequence).upper()
    if not sequence:
        raise ValueError("FASTAに塩基配列がありません。")
    invalid = sorted(set(sequence) - IUPAC)
    if invalid:
        raise ValueError("DNA以外の文字が含まれています: " + ", ".join(invalid))
    if len(sequence) > MAX_SEQUENCE_LENGTH:
        raise ValueError(f"現在のGUIは各配列{MAX_SEQUENCE_LENGTH} bp以下に対応します。")
    return SequenceRecord(name, sequence)


def align_alleles(a: str, b: str) -> AlleleAlignment:
    """Global align with affine-gap scoring, returning A→B coordinate mapping."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(a, b)[0]
    blocks_a, blocks_b = alignment.aligned
    mapping: Dict[int, int] = {}
    variants: List[Variant] = []
    prev_a = prev_b = 0

    for (start_a, end_a), (start_b, end_b) in zip(blocks_a, blocks_b):
        start_a, end_a = int(start_a), int(end_a)
        start_b, end_b = int(start_b), int(end_b)
        if start_a > prev_a or start_b > prev_b:
            ref, alt = a[prev_a:start_a], b[prev_b:start_b]
            kind = "挿入" if not ref else "欠失" if not alt else "複合変異"
            anchor = min(prev_a, len(a) - 1)
            variants.append(Variant(anchor, prev_b, ref, alt, kind))
        for pos_a, pos_b in zip(range(start_a, end_a), range(start_b, end_b)):
            mapping[pos_a] = pos_b
            if a[pos_a] != b[pos_b]:
                variants.append(Variant(pos_a, pos_b, a[pos_a], b[pos_b], "SNP"))
        prev_a, prev_b = end_a, end_b

    if prev_a < len(a) or prev_b < len(b):
        ref, alt = a[prev_a:], b[prev_b:]
        kind = "挿入" if not ref else "欠失" if not alt else "複合変異"
        variants.append(Variant(min(prev_a, len(a) - 1), prev_b, ref, alt, kind))
    return AlleleAlignment(variants, mapping)


def available_enzymes() -> List[str]:
    """Exclude non-PCR and unverified restriction enzymes from recommendations."""
    return sorted(
        name for name, enzyme in ENZYME_METADATA.items()
        if enzyme.pcr_compatible and enzyme.recommendable
        and enzyme.cut_model == "verified"
    )


def _format_fragments(sizes: Sequence[int]) -> str:
    """Preserve co-migrating fragments as multiplicities (e.g. 200×2)."""
    from collections import Counter

    counts = Counter(int(size) for size in sizes)
    return " + ".join(
        f"{size}×{count}" if count > 1 else str(size)
        for size, count in sorted(counts.items(), reverse=True)
    )


def marker_classification(gel_analysis: dict) -> str:
    """Classify visible genotyping patterns; not biological dominance."""
    comparisons = (gel_analysis.get("genotype_discrimination") or {}).get("comparisons") or {}
    aa_ab = bool((comparisons.get("AA_vs_AB") or {}).get("distinguishable"))
    aa_bb = bool((comparisons.get("AA_vs_BB") or {}).get("distinguishable"))
    ab_bb = bool((comparisons.get("AB_vs_BB") or {}).get("distinguishable"))
    if aa_ab and ab_bb and aa_bb:
        return "共優性（3遺伝子型を識別）"
    if aa_bb:
        if not aa_ab and ab_bb:
            return "優性型（AAとABが同じ）"
        if aa_ab and not ab_bb:
            return "優性型（ABとBBが同じ）"
        return "2群識別（ヘテロ判別不可）"
    return "判別困難"


def scan_caps(
    a: str,
    b: str,
    enzymes: Sequence[str],
    min_gap: int = 25,
    profile: GelProfile | None = None,
) -> List[dict]:
    """Simulate restriction of the *supplied* molecules, not hypothetical PCR.

    Use on actual amplicons for marker validation. Applying to longer templates
    is an exploratory scan and must not be interpreted as an amplicon prediction.
    """
    if not enzymes:
        return []
    unknown = set(enzymes) - set(available_enzymes())
    if unknown:
        raise ValueError("未対応またはPCRに適さない酵素: " + ", ".join(sorted(unknown)))
    selected = {name: ENZYME_METADATA[name] for name in enzymes}
    profile = profile or GelProfile()
    profile.validate()
    candidates = caps_scan(
        a, b, enzymes=selected, gel_min_gap=min_gap,
        ladder=profile.ladder,
        custom_ladder_bands=profile.custom_bands or None,
        gel_percent=profile.gel_percent,
    )
    rows = []
    for result in candidates:
        if not result.allele_a_cuts and not result.allele_b_cuts:
            continue  # length-only InDel, not a restriction-based CAPS marker
        gap = result.min_gel_gap
        gel = assess_gel(
            result.allele_a_fragments, result.allele_b_fragments, profile,
        )
        hetero = gel.get("heterozygote_bands") or []
        hetero_sizes = [
            int(band["size"])
            for band in hetero
            for _ in range(int(band.get("copies", 1)))
        ]
        cut_a = bool(result.allele_a_cuts)
        cut_b = bool(result.allele_b_cuts)
        cut_allele = (
            "Aのみ" if cut_a and not cut_b
            else "Bのみ" if cut_b and not cut_a
            else "A・B両方" if cut_a and cut_b
            else "切断なし"
        )
        rows.append({
            "酵素": result.enzyme,
            "認識配列": result.recognition,
            "切断アレル": cut_allele,
            "A切断部位数": len(result.allele_a_cuts),
            "B切断部位数": len(result.allele_b_cuts),
            "A断片(bp)": _format_fragments(result.allele_a_fragments),
            "B断片(bp)": _format_fragments(result.allele_b_fragments),
            "AB断片(bp)": _format_fragments(hetero_sizes),
            "マーカー型": gel["marker_type"],
            "実験適合スコア": gel["utility_score"],
            "最大バンド差(bp)": gap,
            "最小断片(bp)": min(result.allele_a_fragments + result.allele_b_fragments),
            "ゲル判別": (
                "可能" if result.recommendation_eligible
                and gel["marker_type"] != "判別困難" else "要検討"
            ),
            "認識サイト数の差": len(result.allele_a_cuts) != len(result.allele_b_cuts),
            "ゲル評価": round(float(gel.get("score") or 0), 2),
            "推奨ゲル濃度(%)": gel.get("gel_percent", ""),
            "ラダー": gel.get("ladder", ""),
            "可視断片範囲(bp)": "-".join(map(str, gel["visible_range_bp"])),
            "可視識別断片下限(bp)": gel.get("min_diagnostic_visible_bp") or "",
            "ゲル注意": "; ".join(
                reason for reason in gel.get("reasons", [])
                if any(word in reason for word in ("below", "outside", "penalty"))
            ),
            "備考": result.note,
        })
    return rows


def mapped_interval(mapping: Dict[int, int], start: int, end: int) -> Tuple[int, int] | None:
    """Map an A interval to B only when every base is aligned and consecutive."""
    if start < 0 or end <= start:
        return None
    positions = [mapping.get(position) for position in range(start, end)]
    if any(position is None for position in positions):
        return None
    if positions != list(range(positions[0], positions[0] + len(positions))):
        return None
    return positions[0], positions[-1] + 1


def table_bytes(rows: Sequence[dict], delimiter: str = ",") -> bytes:
    """Excel-compatible UTF-8 BOM for CSV; standard UTF-8 for TSV."""
    if delimiter not in (",", "\t"):
        raise ValueError("delimiter must be comma or tab")
    if not rows:
        return b""
    output = io.StringIO(newline="")
    columns = list(dict.fromkeys(key for row in rows for key in row))
    writer = csv.DictWriter(output, fieldnames=columns, delimiter=delimiter,
                            extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    encoding = "utf-8-sig" if delimiter == "," else "utf-8"
    return output.getvalue().encode(encoding)
