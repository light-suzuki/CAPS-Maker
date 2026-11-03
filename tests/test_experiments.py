"""Experimental profiles, phenotype interpretation and strict specificity."""
import pytest

from caps_maker.experiments import (
    GelProfile, assess_gel, parse_custom_ladder, off_target_verdict,
)
from caps_maker.validation import parse_fragment_column


def test_profile_and_three_genotypes():
    profile = GelProfile(ladder="100bp", gel_percent=2.0, min_visible_bp=50)
    result = assess_gel([300, 200], [500], profile)
    assert result["marker_type"].startswith("共優性")
    assert result["visible_range_bp"] == (100, 2000)
    assert 0 <= result["utility_score"] <= 100
    assert result["score_kind"] == "heuristic_not_success_probability"


def test_invisible_diagnostic_bands_are_not_claimed_codom():
    profile = GelProfile(ladder="100bp", gel_percent=2.0, min_visible_bp=250)
    result = assess_gel([300, 200], [500], profile)
    assert result["marker_type"].startswith("共優性")  # 300 vs 500 remains visible
    profile = GelProfile(ladder="1kb", gel_percent=1.0, min_visible_bp=600)
    result = assess_gel([300, 200], [500], profile)
    assert result["marker_type"] == "判別困難"


def test_custom_ladder_validation():
    assert parse_custom_ladder("500,100,100,300") == (100, 300, 500)
    with pytest.raises(ValueError):
        parse_custom_ladder("not-a-size")


class FakeGenome:
    def fetch(self, name, start, end, strand="+"):
        assert name == "off"
        return "G" * (end - start + 1)


def fake_search():
    return [{
        "db": "mockDB", "search_completeness": "complete", "search_complete": True,
        "products": [
            {"subject": "Chr4", "start": 100, "end": 599, "size": 500,
             "orientation": "F/R", "fwd_primer": "F", "rev_primer": "R"},
            {"subject": "off", "start": 1, "end": 1500, "size": 1500,
             "orientation": "F/R", "fwd_primer": "F", "rev_primer": "R"},
        ],
    }]


def test_strict_rejects_offtarget_and_gel_tolerant_reconstructs_it():
    kwargs = {
        "rows_by_db": fake_search(), "genomes_by_db": {"mockDB": FakeGenome()},
        "primers": {"F": "G" * 20, "R": "C" * 20}, "enzyme": "EcoRI",
        "a_fragments": [300, 200], "b_fragments": [500],
        "profile": GelProfile(ladder="100bp", gel_percent=2.0),
        "target_chrom": "Chr4", "target_interval": (100, 599),
    }
    assert off_target_verdict(**kwargs, policy="strict")["accepted"] is False
    positive = off_target_verdict(**kwargs, policy="gel-tolerant")
    assert positive["accepted"] is True
    assert positive["status"] == "gel_tolerant_pass"
    missing = off_target_verdict(**{**kwargs, "genomes_by_db": {}},
                                 policy="gel-tolerant")
    assert missing["accepted"] is False


def test_specificity_unknown_if_search_truncated():
    rows = fake_search()
    rows[0]["search_completeness"] = "possibly_truncated"
    outcome = off_target_verdict(
        rows_by_db=rows, genomes_by_db={"mockDB": FakeGenome()},
        primers={"F": "G" * 20, "R": "C" * 20}, enzyme="EcoRI",
        a_fragments=[300, 200], b_fragments=[500],
        profile=GelProfile(), target_chrom="Chr4",
        target_interval=(100, 599), policy="gel-tolerant",
    )
    assert not outcome["accepted"]


def test_parse_band_multiplicity():
    assert parse_fragment_column("300 + 200×2 + 100") == [300, 200, 200, 100]
