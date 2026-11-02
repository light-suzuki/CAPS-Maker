import csv
import io

import pytest

from caps_maker.core import (
    align_alleles, available_enzymes, mapped_interval,
    parse_sequence, scan_caps, table_bytes,
)


def test_one_record_fasta_and_bare_sequence():
    assert parse_sequence(">A description\nacgt\nTT\n").sequence == "ACGTTT"
    assert parse_sequence("acgt\n").sequence == "ACGT"
    assert parse_sequence(">\nACGT\n").name == "sequence"


def test_multiple_records_rejected():
    with pytest.raises(ValueError):
        parse_sequence(">a\nACGT\n>b\nACGT\n")


def test_bad_characters_rejected():
    with pytest.raises(ValueError):
        parse_sequence(">a\nACGT1\n")


def test_snp_alignment_and_mapping():
    result = align_alleles("ACGTAACG", "ACGTCACG")
    assert len(result.variants) == 1
    assert result.variants[0].position_a == 4
    assert mapped_interval(result.a_to_b, 0, 8) == (0, 8)


def test_insertion_preserves_contiguous_binding_mapping():
    result = align_alleles("AACCTTTGG", "AACCTTTAGG")
    assert result.variants
    assert mapped_interval(result.a_to_b, 0, 4) == (0, 4)
    # The inserted base makes the middle cross-indel interval nonconsecutive.
    assert mapped_interval(result.a_to_b, 0, 9) is None


def test_caps_eco_ri_and_enzyme_selection():
    a = "C" * 60 + "GAATTC" + "C" * 60
    b = "C" * 60 + "GACTTC" + "C" * 60
    assert "EcoRI" in available_enzymes()
    results = scan_caps(a, b, ["EcoRI"], min_gap=25)
    assert len(results) == 1
    row = results[0]
    assert row["酵素"] == "EcoRI"
    assert row["認識サイト数の差"]
    assert row["A切断部位数"] == 1
    assert row["B切断部位数"] == 0
    assert row["切断アレル"] == "Aのみ"
    assert row["AB断片(bp)"]
    assert row["マーカー型"].startswith("共優性")
    assert row["最小断片(bp)"] > 0


def test_csv_excel_and_tsv_roundtrip():
    rows = [{"酵素": "EcoRI", "A断片(bp)": "60 + 66"},
            {"酵素": "XbaI", "A断片(bp)": "126"}]
    csv_bytes = table_bytes(rows)
    assert csv_bytes.startswith(b"\xef\xbb\xbf")
    parsed = list(csv.DictReader(io.StringIO(csv_bytes.decode("utf-8-sig"))))
    assert parsed[1]["酵素"] == "XbaI"
    assert b"\t" in table_bytes(rows, "\t")
