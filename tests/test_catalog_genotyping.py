"""Persistent named resources and conservative PCR-marker display gates."""
import os
import shutil
import subprocess

import pytest

from caps_maker.catalog import load_catalog, local_path, upsert_entry, validate_resource
from caps_maker.genotyping import usable_marker


def test_named_entries_update_without_losing_other_resources(tmp_path):
    path = tmp_path / "settings" / "catalog.json"
    first = upsert_entry("databases", {"name": "Reference", "kind": "BLAST DB", "path": "/ref"}, path)
    upsert_entry("enzyme_sets", {"name": "Lab", "enzymes": ["EcoRI"]}, path)
    with pytest.raises(ValueError, match="同じ名前"):
        upsert_entry("databases", {"name": "reference", "path": "/different"}, path)
    upsert_entry("databases", {**first, "path": "/new-reference"}, path)
    reopened = load_catalog(path)
    assert len(reopened["databases"]) == 1
    assert reopened["databases"][0]["path"] == "/new-reference"
    assert reopened["enzyme_sets"][0]["enzymes"] == ["EcoRI"]


def test_corrupt_catalog_is_never_silently_replaced(tmp_path):
    path = tmp_path / "catalog.json"
    path.write_text("broken-json")
    with pytest.raises(ValueError, match="上書きしていません"):
        upsert_entry("databases", {"name": "test"}, path)
    assert path.read_text() == "broken-json"


def test_windows_paths_are_usable_on_wsl():
    expected = "D:\\genome\\reference" if os.name == "nt" else "/mnt/d/genome/reference"
    assert local_path('"D:\\genome\\reference"') == expected
    with pytest.raises(ValueError):
        local_path("")


@pytest.mark.integration
def test_registration_checks_real_blast_database(tmp_path):
    if not (shutil.which("blastdbcmd") and shutil.which("makeblastdb")):
        pytest.skip("requires BLAST+")
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">Chr1\n" + "ACGT" * 40 + "\n")
    db = tmp_path / "reference"
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(db)],
                   check=True, capture_output=True)
    checked, _ = validate_resource("BLAST DB", str(db))
    assert checked == str(db)
    with pytest.raises(ValueError, match="DBを開けません"):
        validate_resource("BLAST DB", str(tmp_path / "missing"))


def valid_row():
    return {"Forward(5'-3')": "A" * 20, "Reverse(5'-3')": "C" * 20,
            "A産物(bp)": 500, "B産物(bp)": 500, "参照検証採用可": True,
            "ゲル判別": "可能", "特異性判定": "strict_pass",
            "マーカー型": "共優性（3遺伝子型を識別）"}


@pytest.mark.parametrize("change", [
    {"Forward(5'-3')": ""}, {"B産物(bp)": 0}, {"参照検証採用可": False},
    {"特異性判定": "indeterminate"}, {"特異性判定": "no_reference_anchor"},
    {"特異性判定": "not_strict"}, {"ゲル判別": "要検討"}, {"マーカー型": "判別困難"},
])
def test_unusable_or_unchecked_markers_are_hidden(change):
    assert usable_marker(valid_row())
    assert not usable_marker({**valid_row(), **change})


def test_two_group_marker_requires_explicit_homozygote_only_goal():
    row = {**valid_row(), "マーカー型": "優性型（AAとABが同じ）"}
    assert not usable_marker(row)
    assert usable_marker(row, three_genotypes=False)


def test_gel_tolerant_requires_explicit_policy():
    row = {**valid_row(), "特異性判定": "gel_tolerant_pass"}
    assert not usable_marker(row)
    assert usable_marker(row, policy="gel-tolerant")
