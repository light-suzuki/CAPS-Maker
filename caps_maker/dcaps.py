"""One-SNP dCAPS candidate construction on both actual allele amplicons.

The altered primer bases are copied into both PCR products. The non-engineered
partner primer is checked against both allele sequences, and no SNP is silently
masked by the primer footprint. BLAST evaluation is performed separately.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from primerblast_oss.caps import (
    ENZYME_METADATA, dcaps_candidates, materialize_dcaps_primers, revcomp,
)
from primerblast_oss.dcaps_workflow import _tm
from .core import AlleleAlignment, Variant, scan_caps
from .experiments import GelProfile


def _primer_product(pair: dict, a: str, b: str, alignment: AlleleAlignment
                    ) -> Tuple[str, str] | None:
    # Late import avoids recursive dependency when used by design_caps().
    from .primers import _amplicons_for_pair
    return _amplicons_for_pair(pair, a, b, alignment)


def _gc(seq: str) -> float:
    return round(100 * sum(base in "GC" for base in seq) / len(seq), 1)


def design_dcaps_for_pair(
    a: str,
    b: str,
    alignment: AlleleAlignment,
    variant: Variant,
    parent_pair: Dict,
    enzymes: Sequence[str],
    min_gap: int = 25,
    max_mismatches: int = 1,
    max_candidates: int = 25,
    profile: GelProfile | None = None,
) -> List[Dict]:
    if (variant.kind != "SNP"
            or len(variant.allele_a) != 1
            or len(variant.allele_b) != 1
            or not 0 <= variant.position_a < len(a)):
        return []
    parent_start = int(parent_pair["left_start"])
    parent_end = int(parent_pair["right_start"])
    if not parent_start < variant.position_a < parent_end:
        return []
    amplicons = _primer_product(parent_pair, a, b, alignment)
    if not amplicons:
        return []
    ref_product, alt_product = amplicons
    snp_index = variant.position_a - parent_start
    if not 0 < snp_index < len(ref_product) - 1:
        return []
    if ref_product[snp_index] != variant.allele_a:
        return []

    selected = {name: ENZYME_METADATA[name] for name in enzymes}
    frames = dcaps_candidates(
        ref_product, snp_index, variant.allele_a, variant.allele_b,
        enzymes=selected, max_primer_mismatch=max_mismatches,
    )
    choices = materialize_dcaps_primers(
        ref_product, snp_index, frames, min_size=18, opt_size=20, max_size=25,
    )
    results: List[Dict] = []
    seen = set()
    for option in choices:
        if len(results) >= max_candidates:
            break
        role = option["primer_role"]
        modified = option["primer_sequence"]
        if role == "F":
            left = parent_start + option["primer_start"]
            right = parent_end
            unmodified_forward = a[left:left + len(modified)]
            primer_for_matching = {
                **parent_pair, "forward": unmodified_forward,
                "left_start": left, "left_len": len(modified),
            }
            actual_f, actual_r = modified, parent_pair["reverse"]
        else:
            left = parent_start
            right = parent_start + option["primer_end"]
            original_reverse = revcomp(a[right - len(modified) + 1:right + 1])
            primer_for_matching = {
                **parent_pair, "reverse": original_reverse,
                "right_start": right, "right_len": len(modified),
            }
            actual_f, actual_r = parent_pair["forward"], modified
        if left < 0 or right >= len(a):
            continue
        # Checks that the *unmodified* recognition-region primer binding bases
        # match BOTH allele templates. The intended 1 bp mismatch is explicit.
        raw_products = _primer_product(primer_for_matching, a, b, alignment)
        if not raw_products:
            continue
        original_a, original_b = raw_products
        if role == "F":
            product_a = modified + original_a[len(modified):]
            product_b = modified + original_b[len(modified):]
        else:
            encoded = revcomp(modified)
            product_a = original_a[:-len(modified)] + encoded
            product_b = original_b[:-len(modified)] + encoded

        key = (actual_f, actual_r, option["enzyme"])
        if key in seen:
            continue
        seen.add(key)
        digests = scan_caps(
            product_a, product_b, [option["enzyme"]], min_gap, profile=profile)
        if not digests:
            continue
        difference = ",".join(
            f"{parent_start + int(change['position']) + 1}:{change['from']}>{change['to']}"
            for change in option["engineered_changes"]
        )
        for digest in digests:
            results.append({
                "設計方式": "dCAPS",
                "_left_start": left,
                "_right_start": right,
                "標的多型": variant.label,
                "Forward(5'-3')": actual_f,
                "Reverse(5'-3')": actual_r,
                "Tm F(°C)": _tm(actual_f) if role == "F" else parent_pair["tm_f"],
                "Tm R(°C)": _tm(actual_r) if role == "R" else parent_pair["tm_r"],
                "GC F(%)": _gc(actual_f),
                "GC R(%)": _gc(actual_r),
                "A産物(bp)": len(product_a),
                "B産物(bp)": len(product_b),
                "人工ミスマッチ数": len(option["engineered_changes"]),
                "人工ミスマッチ位置(A配列1-based)": difference,
                "改変プライマー": role,
                "SNP距離(3'末端)": option["three_prime_distance_to_snp"],
                "備考": digest.get("備考", "") + "; engineered primer mismatch; validate experimentally",
                **{key: value for key, value in digest.items() if key != "備考"},
            })
    return results
