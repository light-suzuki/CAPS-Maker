"""Population-aware marker gate, band calls, scoring-sheet QC and exports (synthetic)."""
import csv
import io

import pytest

from caps_maker.experiments import GelProfile
from caps_maker.genotyping import (
    GOALS, bench_sheet, call_from_bands, expected_lanes, genotype_table_bytes, normalize_call,
    qc_report, read_scoring_sheet, recombination_report, rqtl_csv_bytes, scoring_sheet,
    segregation_test, select_by_density, unique_marker_names, usable_marker,
)

F2, RIL, BC_A, BC_B = GOALS


def marker(kind, a="316 + 141", b="457", chrom="chr1", start=1001, end=1457, enzyme="EcoRI"):
    return {
        "Forward(5'-3')": "ACGTACGTACGTACGTAC", "Reverse(5'-3')": "TGCATGCATGCATGCATG",
        "A産物(bp)": 457, "B産物(bp)": 457, "参照検証採用可": True, "ゲル判別": "可能",
        "特異性判定": "strict_pass", "特異性区分": "A 厳密合格", "マーカー型": kind,
        "A断片(bp)": a, "B断片(bp)": b, "参照染色体": chrom,
        "参照産物開始(1-based)": start, "参照産物終了(1-based)": end, "酵素": enzyme,
        "設計方式": "CAPS",
    }


@pytest.mark.parametrize("kind, accepted", [
    ("共優性（3遺伝子型を識別）", {F2, RIL, BC_A, BC_B}),
    ("優性型（ABとBBが同じ）", {RIL, BC_A}),
    ("優性型（AAとABが同じ）", {RIL, BC_B}),
    ("2群識別（ヘテロ判別不可）", {RIL}),
    ("判別困難", set()),
])
def test_population_goal_requires_the_needed_genotype_pair(kind, accepted):
    row = marker(kind)
    assert {g for g in GOALS if usable_marker(row, goal=g)} == accepted
    # Legacy call sites keep their meaning.
    assert usable_marker(row) == (F2 in accepted)
    assert usable_marker(row, three_genotypes=False) == (RIL in accepted)


def test_backcross_rejects_marker_that_cannot_split_aa_from_ab():
    # The old "AA vs BB" setting would have offered this to a BC population.
    assert usable_marker(marker("2群識別（ヘテロ判別不可）"), three_genotypes=False)
    assert not usable_marker(marker("2群識別（ヘテロ判別不可）"), goal=BC_A)


def test_band_reading_calls_the_matching_genotype():
    lanes = expected_lanes(marker("共優性（3遺伝子型を識別）"), GelProfile())
    assert call_from_bands([320, 140], lanes) == ["AA"]
    assert call_from_bands([455], lanes) == ["BB"]
    assert call_from_bands([460, 318, 139], lanes) == ["AB"]
    assert call_from_bands([200], lanes) == []
    with pytest.raises(ValueError):
        call_from_bands([200], lanes, tolerance=0)


def test_marker_names_are_positional_and_unique():
    rows = [marker("共優性（3遺伝子型を識別）"), marker("共優性（3遺伝子型を識別）")]
    assert unique_marker_names(rows) == ["chr1_1229_EcoRI", "chr1_1229_EcoRI_2"]


def test_bench_sheet_lists_visible_bands_per_genotype():
    sheet = bench_sheet([marker("共優性（3遺伝子型を識別）")], GelProfile(), F2)[0]
    assert sheet["AA（親A型）"] == "316 + 141"
    assert sheet["AB（ヘテロ）"] == "457 + 316 + 141"
    assert sheet["BB（親B型）"] == "457"
    assert sheet["記録できる型"] == "AA・AB・BB" and sheet["目的に適合"] is True


@pytest.mark.parametrize("value, code", [
    ("AA", "A"), ("a", "A"), ("0", "A"), ("AB", "H"), ("1.0", "H"), ("BB", "B"),
    ("A-", "C"), ("D", "D"), ("", "-"), ("NA", "-"),
])
def test_call_codes(value, code):
    assert normalize_call(value) == code


def test_segregation_tests():
    assert segregation_test({"A": 25, "H": 50, "B": 25}, "F2（自殖・F1×F1）")["p値"] == pytest.approx(1.0)
    distorted = segregation_test({"A": 60, "H": 30, "B": 10}, "F2（自殖・F1×F1）")
    assert distorted["自由度"] == 2 and distorted["p値"] < 0.001
    assert segregation_test({"A": 30, "D": 90}, "F2（自殖・F1×F1）")["期待比"].endswith("1:3")
    assert segregation_test({"A": 40, "H": 40}, "戻し交雑（反復親A）")["p値"] == pytest.approx(1.0)
    assert segregation_test({"A": 2, "B": 1}, "RIL・DH（ホモ固定）") is None


def _filled_sheet(rows, calls):
    raw = scoring_sheet(rows, GelProfile(), [name for name, _ in calls[3:]])
    table = list(csv.reader(io.StringIO(raw.decode("utf-8-sig"))))
    by_line = {r[1]: r for r in table if r[0] in ("Parent_A", "Parent_B", "F1", "Progeny")}
    for line, values in calls:
        by_line[line][2:] = values
    out = io.StringIO()
    csv.writer(out, lineterminator="\n").writerows(table)
    return out.getvalue().encode("utf-8-sig")


def test_scoring_sheet_round_trip_qc_and_exports():
    rows = [marker("共優性（3遺伝子型を識別）"),
            marker("優性型（ABとBBが同じ）", a="316 + 141", b="457", start=5001, end=5457, enzyme="HinfI")]
    progeny = [(f"F2_{i:03d}", v) for i, v in enumerate(
        [["AA", "A"], ["AB", "B-"], ["BB", "B-"], ["AB", "AB"], ["AA", "AA"], ["", "B-"]], 1)]
    data = _filled_sheet(rows, [("Parent_A", ["AA", "AA"]), ("Parent_B", ["AA", "BB"]),
                                ("F1", ["AB", "B-"]), *progeny])
    scored = read_scoring_sheet(data)
    assert scored.markers == ["chr1_1229_EcoRI", "chr1_5229_HinfI"]
    assert scored.info["Expected_AA"]["chr1_1229_EcoRI"] == "316 + 141"
    table, problems = qc_report(scored, "F2（自殖・F1×F1）")
    first = table[0]
    assert (first["AA"], first["AB"], first["BB"], first["欠測"]) == (2, 2, 1, 1)
    messages = " ".join(p["内容"] for p in problems)
    assert "コントロール不一致" in messages        # Parent_B scored AA on marker 1
    assert "判別できません" in messages            # AB / BB on a dominant marker

    wide = list(csv.reader(io.StringIO(genotype_table_bytes(scored).decode("utf-8-sig"))))
    assert wide[0][:3] == ["Role", "Line", "Phenotype"]
    assert wide[1][3:] == ["chr1", "chr1"] and wide[2][3:] == ["1229", "5229"]
    assert wide[-1] == ["Progeny", "F2_006", "", "", "B-"]

    rqtl = list(csv.reader(io.StringIO(rqtl_csv_bytes(scored).decode())))
    assert rqtl[0] == ["id", "chr1_1229_EcoRI", "chr1_5229_HinfI"]
    assert rqtl[2][1:] == ["0.001229", "0.005229"]
    assert [r[0] for r in rqtl[3:]] == [name for name, _ in progeny]
    assert rqtl[4][1:] == ["H", "D"] and rqtl[-1][1:] == ["-", "D"]


def test_scoring_sheet_rejects_foreign_layout():
    with pytest.raises(ValueError, match="記録シート形式"):
        read_scoring_sheet("sample,M1\nx,AA\n".encode())


def _at(pos, kind="共優性（3遺伝子型を識別）", chrom="chr1"):
    row = marker(kind, chrom=chrom, start=pos - 200, end=pos + 200)
    row["標的座標(1-based)"] = pos
    return row


def test_density_picks_one_locus_per_bin_and_reports_empty_bins():
    rows = [_at(10_000), _at(10_000),                       # alternative pair, same locus
            _at(40_000, "優性型（ABとBBが同じ）"), _at(60_000),  # same bin: co-dominant wins
            _at(250_000), _at(55_000, chrom="chr2")]
    chosen, empty = select_by_density(rows, "一定間隔", "chr1", 1, 300_000, spacing_bp=100_000)
    assert chosen == [3, 4]          # 60 kb is co-dominant and nearest the 50 kb centre
    assert [e["区画"] for e in empty] == [2]

    # Narrowing: same candidates, denser inside a sub-interval.
    chosen, empty = select_by_density(rows, "一定間隔", "chr1", 1, 100_000, spacing_bp=25_000)
    assert chosen == [0, 2, 3] and len(empty) == 1
    chosen, _ = select_by_density(rows, "指定本数", "chr1", 1, 300_000, count=1)
    assert len(chosen) == 1
    chosen, _ = select_by_density(rows, "全座位", "chr1", 1, 300_000)
    assert chosen == [0, 2, 3, 4]
    with pytest.raises(ValueError):
        select_by_density(rows, "一定間隔", "chr1", 10, 1)


def test_recombination_report_locates_breakpoints():
    rows = [_at(p) for p in (100_000, 200_000, 300_000, 400_000)]
    progeny = [("F2_001", ["AA", "AA", "AB", "AB"]),   # break 200k-300k
               ("F2_002", ["AA", "", "AB", "AB"]),     # break spans 100k-300k (missing)
               ("F2_003", ["A-", "AB", "AB", "BB"]),   # A- vs AB compatible; break 300k-400k
               ("F2_004", ["BB", "BB", "BB", "BB"])]
    data = _filled_sheet(rows, [("Parent_A", ["AA"] * 4), ("Parent_B", ["BB"] * 4),
                                ("F1", ["AB"] * 4), *progeny])
    breaks, intervals = recombination_report(read_scoring_sheet(data))
    assert [(b["個体"], b["左マーカー"].split("_")[1], b["右マーカー"].split("_")[1]) for b in breaks] == [
        ("F2_001", "200000", "300000"), ("F2_002", "100000", "300000"), ("F2_003", "300000", "400000")]
    counts = [(r["組換え個体"], r["うち区間が確定"]) for r in intervals]
    assert counts == [(1, 0), (2, 1), (1, 1)]
