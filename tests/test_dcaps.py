"""Ensure dCAPS substitutes primer bases in both actual amplicons."""
from primerblast_oss.caps import revcomp

from caps_maker.core import align_alleles
from caps_maker.dcaps import design_dcaps_for_pair


def test_engineered_one_mismatch_dcaps_not_natural_caps():
    # EcoRI: GAATTC. A engineered C→A at third position gives GAATTC;
    # B (SNP at last position) gives GAATTA. Neither original is cut.
    a = "A" * 150 + "GACTTC" + "A" * 150
    b = "A" * 150 + "GACTTA" + "A" * 150
    alignment = align_alleles(a, b)
    assert len(alignment.variants) == 1
    variant = alignment.variants[0]
    assert variant.kind == "SNP"
    parent = {
        "left_start": 0, "left_len": 20,
        "right_start": len(a) - 1, "right_len": 20,
        "forward": a[:20],
        "reverse": revcomp(a[-20:]),
        "tm_f": 50, "tm_r": 50,
    }
    rows = design_dcaps_for_pair(
        a, b, alignment, variant, parent,
        ["EcoRI"], max_mismatches=1, max_candidates=20,
    )
    assert rows, "dCAPS should yield an engineered F primer"
    assert any(row["設計方式"] == "dCAPS"
               and row["改変プライマー"] == "F"
               and row["人工ミスマッチ数"] == 1
               and row["切断アレル"] == "Aのみ"
               and row["_left_start"] > 0
               for row in rows)
