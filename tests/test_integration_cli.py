"""Optional real-tool smoke test (enabled on CI with Primer3 and BLAST+)."""
import random
import shutil
import subprocess

import pytest

from caps_maker.core import align_alleles
from caps_maker.primers import design_caps


@pytest.mark.integration
def test_real_primer3_blast_caps(tmp_path):
    if not (shutil.which("primer3_core") and shutil.which("blastn")
            and shutil.which("makeblastdb")):
        pytest.skip("Primer3 and BLAST+ executables are required")
    rng = random.Random(1841)
    # Avoid extra EcoRI sites so the diagnostic digest is predictable.
    while True:
        left = "".join(rng.choice("ACGT") for _ in range(360))
        right = "".join(rng.choice("ACGT") for _ in range(360))
        a = left + "GAATTC" + right
        if a.count("GAATTC") == 1:
            break
    b = left + "GACTTC" + right
    fasta = tmp_path / "synthetic.fa"
    fasta.write_text(">target\n" + a + "\n", encoding="utf-8")
    db = tmp_path / "synthetic-db"
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
                    "-out", str(db)], check=True, capture_output=True, text=True)
    alignment = align_alleles(a, b)
    assert len(alignment.variants) == 1
    rows, summary = design_caps(
        a, b, alignment, alignment.variants[0],
        ["EcoRI"], [str(db)], product_min=200, product_max=500,
        num_pairs=20, timeout=180,
    )
    assert summary["primer3_pairs"] > 0
    assert summary["both_alleles_bind"] > 0
    assert rows
    for row in rows:
        assert row["切断アレル"] == "Aのみ"
        assert row["A産物(bp)"] == sum(
            int(part.split("×")[0]) * int(part.split("×")[1] if "×" in part else 1)
            for part in row["A断片(bp)"].split(" + ")
        )
        assert 200 <= row["A産物(bp)"] <= 500
        assert row["BLAST判定"] in ("specific", "indeterminate", "nonspecific",
                                     "ambiguous")
