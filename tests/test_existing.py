from caps_maker.existing import primer_intervals
from caps_maker.core import align_alleles
from caps_maker.primers import _amplicons_for_pair


def test_primers_define_target_interval_and_both_alleles_bind():
    forward = "G" * 20
    reverse = "C" * 20
    a = forward + "ATGC" * 50 + "G" * 20
    b = forward + "ATGC" * 25 + "ATGA" + "ATGC" * 24 + "G" * 20
    pairs = primer_intervals(a, forward, reverse, 70, 300)
    assert (0, len(a) - 1) in pairs
    pair = {
        "left_start": 0, "left_len": 20,
        "right_start": len(a) - 1, "right_len": 20,
        "forward": forward, "reverse": reverse,
    }
    assert _amplicons_for_pair(pair, a, b, align_alleles(a, b)) == (a, b)
