from unittest.mock import patch
import json
from types import SimpleNamespace

from caps_maker.core import align_alleles
from caps_maker.primers import _amplicons_for_pair, design_caps


def example_pair(a):
    return {
        "forward": "G" * 20,
        "reverse": "G" * 20,
        "left_start": 0,
        "left_len": 20,
        "right_start": len(a) - 1,
        "right_len": 20,
        "product_size": len(a),
        "tm_f": 60,
        "tm_r": 60,
        "specificity": {
            "specificity_status": "specific",
            "search_completeness": "complete",
            "total_off_target": 0,
        }
    }


def test_allele_binding_and_digest_on_real_products():
    a = "G" * 80 + "GAATTC" + "C" * 80
    b = "G" * 80 + "GACTTC" + "C" * 80
    alignment = align_alleles(a, b)
    pair = example_pair(a)
    assert _amplicons_for_pair(pair, a, b, alignment) == (a, b)
    mock_response = SimpleNamespace(
        returncode=0,
        stdout=json.dumps({"pairs": [pair], "primer3_explain": ""}),
        stderr="",
    )
    with patch("caps_maker.primers.subprocess.run", return_value=mock_response):
        rows, summary = design_caps(
            a, b, alignment, alignment.variants[0], ["EcoRI"], ["mock_db"]
        )
    assert summary["both_alleles_bind"] == 1
    assert rows and rows[0]["酵素"] == "EcoRI"
    assert rows[0]["BLAST判定"] == "specific"
    assert rows[0]["Forward(5'-3')"] == "G" * 20


def test_allele_specific_primer_is_not_labeled_shared():
    a = "G" * 80 + "GAATTC" + "C" * 80
    b = "A" + a[1:]
    alignment = align_alleles(a, b)
    assert _amplicons_for_pair(example_pair(a), a, b, alignment) is None
