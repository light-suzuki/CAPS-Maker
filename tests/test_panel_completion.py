"""Issue #2 acceptance: spacing/uncovered bins, product reporting, A/B/C/D
specificity classes, BED and deterministic provenance."""
import json
import random
from datetime import datetime, timezone

import pytest

from caps_maker import existing
from caps_maker.batch import (
    MarkerSite, choose_sites, group_by_locus, marker_track_svg,
)
from caps_maker.experiments import GelProfile, off_target_verdict, specificity_class
from caps_maker.export import bed_bytes, build_manifest, manifest_bytes, rows_digest
from caps_maker.primers import apply_verdict


def synthetic_1mb_sites(seed=7):
    """Known SNP loci across 1 Mb, with an empty 300-500 kb stretch."""
    rng = random.Random(seed)
    sites = []
    for pos in sorted(rng.sample(range(1, 1_000_001), 120)):
        if 300_000 <= pos < 500_000:
            continue
        var = f"{pos}:A>G"
        sites.append(MarkerSite("Chr1", pos, "EcoRI", pos, "CAPS", var))
        # A second enzyme at the same polymorphism is an alternative assay,
        # not an independent map marker.
        sites.append(MarkerSite("Chr1", pos, "HinfI", pos + 2, "CAPS", var))
    return sites


def test_grid_respects_min_spacing_and_reports_uncovered_bins():
    sites = synthetic_1mb_sites()
    chosen, stats = choose_sites(sites, "一定間隔", 1, 1_000_000,
                                 spacing_bp=100_000, min_spacing_bp=60_000)
    positions = [s.position for s in chosen]
    assert all(b - a >= 60_000 for a, b in zip(positions, positions[1:]))
    assert len({s.source_variant for s in chosen}) == len(chosen)  # one per locus
    assert stats["bins"] == 10
    empty = [b for b in stats["uncovered_bins"] if b["reason"] == "多型サイトなし"]
    assert {b["bin"] for b in empty} == {4, 5}
    assert len(chosen) + len(stats["uncovered_bins"]) == 10
    assert stats["loci"] == len(sites) // 2


def test_target_count_never_fills_empty_bins():
    sites = synthetic_1mb_sites()
    chosen, stats = choose_sites(sites, "指定本数", 1, 1_000_000, target_count=10)
    assert len(chosen) == 8
    assert len(stats["uncovered_bins"]) == 2


def test_search_radius_limits_candidates_to_bin_centre():
    sites = [MarkerSite("Chr1", 10, "EcoRI", 10, "CAPS", "10:A>G"),
             MarkerSite("Chr1", 480, "EcoRI", 480, "CAPS", "480:A>G")]
    chosen, stats = choose_sites(sites, "一定間隔", 1, 1000, spacing_bp=1000,
                                 search_radius_bp=100)
    assert [s.position for s in chosen] == [480]
    chosen, stats = choose_sites(sites, "一定間隔", 1, 1000, spacing_bp=1000,
                                 search_radius_bp=5)
    assert chosen == [] and stats["uncovered_bins"][0]["reason"] == "探索半径内にサイトなし"


def test_grid_prefers_caps_then_score():
    sites = [MarkerSite("Chr1", 500, "MseI", 500, "dCAPS", "500:A>G"),
             MarkerSite("Chr1", 900, "EcoRI", 900, "CAPS", "900:C>T"),
             MarkerSite("Chr1", 910, "HinfI", 910, "CAPS", "910:C>T")]
    chosen, _ = choose_sites(sites, "一定間隔", 1, 1000, spacing_bp=1000)
    assert chosen[0].kind == "CAPS"
    scores = {"Chr1:910:HinfI:CAPS": 90.0}
    chosen, _ = choose_sites(sites, "一定間隔", 1, 1000, spacing_bp=1000, scores=scores)
    assert chosen[0].enzyme == "HinfI"


def test_exhaustive_lists_options_but_counts_loci():
    sites = synthetic_1mb_sites()
    chosen, stats = choose_sites(sites, "全酵素サイト", 1, 1_000_000, limit=10_000)
    assert len(chosen) == len(sites)
    assert stats["selected_loci"] == len(group_by_locus(sites)) == len(sites) // 2
    spaced, s_stats = choose_sites(sites, "全酵素サイト", 1, 1_000_000,
                                   limit=10_000, min_spacing_bp=50_000)
    loci = sorted({s.position for s in spaced})
    assert all(b - a >= 50_000 for a, b in zip(loci, loci[1:]))


def test_marker_track_svg_shows_markers_and_gaps():
    sites = synthetic_1mb_sites()
    chosen, stats = choose_sites(sites, "一定間隔", 1, 1_000_000, spacing_bp=100_000)
    svg = marker_track_svg(chosen, 1, 1_000_000, stats["uncovered_bins"])
    assert svg.startswith("<svg") and svg.count("<line") == len(chosen) + 1
    assert "#f3d6d6" in svg


def _search(complete=True, target_present=True, off=False):
    products = []
    if target_present:
        products.append({"subject": "Chr4", "start": 100, "end": 599, "size": 500,
                         "orientation": "F/R"})
    if off:
        products.append({"subject": "Chr2", "start": 1, "end": 900, "size": 900,
                         "orientation": "F/R"})
    return [{"db": "ref", "search_completeness": "complete" if complete else "possibly_truncated",
             "search_complete": complete, "products": products}]


def _verdict(rows, policy="strict", locate=True):
    return off_target_verdict(
        rows_by_db=rows, genomes_by_db={}, primers={"F": "G" * 20, "R": "C" * 20},
        enzyme="EcoRI", a_fragments=[300, 200], b_fragments=[500],
        profile=GelProfile(ladder="100bp", gel_percent=2.0),
        target_chrom="Chr4" if locate else None,
        target_interval=(100, 599) if locate else None, policy=policy)


def test_specificity_ladder_never_collapses():
    assert _verdict(_search())["status"] == "strict_pass"
    # Incomplete search is indeterminate even under the strict policy.
    assert _verdict(_search(complete=False, off=True))["status"] == "indeterminate"
    assert _verdict(_search(complete=False))["status"] == "indeterminate"
    assert _verdict(_search(), locate=False)["status"] == "indeterminate"
    assert _verdict(_search(target_present=False, off=True))["status"] == "target_not_unique"
    assert _verdict(_search(off=True))["status"] == "not_strict"
    classes = {specificity_class(s)[0] for s in (
        "strict_pass", "gel_tolerant_pass", "not_strict", "indeterminate")}
    assert classes == {"A", "B", "C", "D"}


def test_gel_tolerant_pass_ranks_below_strict_pass():
    base = {"ゲル判別": "可能", "実験適合スコア": 80, "酵素": "EcoRI"}
    strict = apply_verdict(dict(base), {"status": "strict_pass", "accepted": True},
                           None, 0, 10)
    tolerant = apply_verdict(dict(base), {"status": "gel_tolerant_pass", "accepted": True},
                             None, 0, 10)
    assert strict["総合設計スコア"] > tolerant["総合設計スコア"]
    assert strict["特異性区分"].startswith("A") and tolerant["特異性区分"].startswith("B")


def test_existing_primers_report_all_products(monkeypatch):
    evidence = [
        {"db": "ref", "search_completeness": "complete", "products": [
            {"subject": "Chr4", "start": 1000, "end": 1499, "size": 500, "orientation": "F/R"},
            {"subject": "Chr4", "start": 1200, "end": 1899, "size": 700, "orientation": "F/R"},
            {"subject": "Chr7", "start": 50, "end": 349, "size": 300, "orientation": "F/F"},
        ]},
        {"db": "other", "search_completeness": "complete", "products": [
            {"subject": "ctg9", "start": 10, "end": 509, "size": 500, "orientation": "F/R"},
        ]},
    ]
    monkeypatch.setattr(existing, "run_primer_check", lambda *a, **k: evidence)
    with pytest.raises(existing.AmbiguousProductError) as info:
        existing.locate_reference_product("G" * 20, "C" * 20, ["ref", "other"],
                                          "Chr4", 900, 2000)
    assert len(info.value.products) == 4
    exact, _ = existing.locate_reference_product(
        "G" * 20, "C" * 20, ["ref", "other"], "Chr4", 900, 1600)
    assert exact == (1000, 1499)
    table = existing.product_table(evidence, ("Chr4", 1000, 1499))
    assert table[0]["区分"] == "標的産物"
    assert sum(r["区分"] == "標的産物" for r in table) == 1
    # Same coordinates in a second assembly are not the intended design locus.
    assert all(r["区分"] == "その他の産物" for r in table if r["DB"] == "other")


def test_bed_and_deterministic_provenance():
    rows = [
        {"酵素": "EcoRI", "設計方式": "CAPS", "特異性判定": "strict_pass",
         "参照染色体": "Chr4", "参照産物開始(1-based)": 100,
         "参照産物終了(1-based)": 599, "総合設計スコア": 88.4},
        {"酵素": "MseI", "設計方式": "CAPS", "特異性判定": "no_reference_anchor",
         "総合設計スコア": 40},
    ]
    bed = bed_bytes(rows).decode().splitlines()
    assert bed[0].startswith("track") and bed[1] == "Chr4\t99\t599\tEcoRI|CAPS|strict_pass\t884\t+"
    assert len(bed) == 2  # unanchored row is not given invented coordinates
    assert bed_bytes(rows[1:]) == b""
    fixed = datetime(2026, 1, 1, tzinfo=timezone.utc)
    m1 = build_manifest({"enzymes": ["EcoRI"]}, rows, now=fixed)
    m2 = build_manifest({"enzymes": ["EcoRI"]}, list(rows), now=fixed)
    assert manifest_bytes(m1) == manifest_bytes(m2)
    assert m1["result_sha256"] == rows_digest(rows)
    assert m1["tool_versions"]["caps_maker"]
    assert rows_digest(rows) != rows_digest(rows[:1])
    json.loads(manifest_bytes(m1))
