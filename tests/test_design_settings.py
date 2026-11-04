import json

import pytest

from caps_maker.catalog import load_catalog, upsert_entry
from caps_maker.design_settings import BUILTINS, apply_settings, snapshot_settings
from caps_maker import existing


def test_template_roundtrip_keeps_conditions_and_excludes_inputs(tmp_path):
    state = {"primary_db": "/private/reference", "chrom": "Chr4", "spacing": 50000,
             "run_mode": "区間マーカーを一括設計", "active_enzymes": ["EcoRI"],
             "preset": "詳細設定", "gel_percentage": 3.0, "min_visible": 30,
             "max_visible": 800, "min_log_gap": 0.05, "gel_min_gap": 20}
    settings = snapshot_settings(state)
    path = tmp_path / "catalog.json"
    upsert_entry("design_templates", {"name": "50 kb", "settings": settings}, path)
    restored = {}
    apply_settings(load_catalog(path)["design_templates"][0]["settings"], restored)
    assert restored["spacing"] == 50000 and restored["run_mode"] == "区間マーカーを一括設計"
    assert restored["enzyme_selection_manual"] == ["EcoRI"]
    assert restored["gel_percentage"] == 3.0 and restored["max_visible"] == 800
    assert "primary_db" not in settings and "chrom" not in settings


def test_old_catalog_keeps_resources_when_templates_are_added(tmp_path):
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps({"schema": 1, "databases": [{"id": "old", "name": "Reference"}], "enzyme_sets": []}))
    upsert_entry("design_templates", {"name": "test", "settings": {}}, path)
    assert load_catalog(path)["databases"][0]["id"] == "old"


@pytest.mark.parametrize("settings", [
    {"spacing": 0}, {"product_min": 800, "product_max": 200}, {"run_mode": "unsupported"},
    {"enzyme_selection_manual": ["not-an-enzyme"]},
    {"gel_options": {"min_visible": 500, "max_visible": 200}},
])
def test_invalid_template_does_not_partly_change_current_settings(settings):
    state = {"spacing": 75000, "gel_percentage": 2.0}
    with pytest.raises(ValueError):
        apply_settings(settings, state)
    assert state == {"spacing": 75000, "gel_percentage": 2.0}


def test_builtin_changes_route_and_resets_previous_custom_gel():
    state = {"gel_percentage": 3.0, "min_visible": 10}
    apply_settings(BUILTINS["spaced"]["settings"], state)
    assert state["spacing"] == 100000 and state["min_spacing"] == 50000
    assert state["run_mode"] == "区間マーカーを一括設計"
    assert "gel_percentage" not in state


@pytest.mark.parametrize("count,complete", [(0, True), (2, True), (1, False)])
def test_primer_only_location_rejects_ambiguous_or_incomplete_search(monkeypatch, count, complete):
    products = [{"subject": "Chr4", "start": 100 + i * 1000, "end": 599 + i * 1000,
                 "orientation": "F/R"} for i in range(count)]
    monkeypatch.setattr(existing, "run_primer_check", lambda *a, **k: [
        {"db": "ref", "search_completeness": "complete" if complete else "possibly_truncated", "products": products}])
    with pytest.raises(existing.AmbiguousProductError):
        existing.locate_primer_product("G" * 20, "C" * 20, ["ref"])
