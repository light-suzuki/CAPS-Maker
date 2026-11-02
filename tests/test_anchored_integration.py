"""Real reference+sample-VCF → Primer3 → BLAST → anchored CAPS integration."""
import random
import shutil
import subprocess

import pysam
import pytest

from caps_maker.core import align_alleles
from caps_maker.primers import design_caps
from caps_maker.sources import read_sample_variants, source_window


@pytest.mark.integration
def test_indexed_vcf_and_reference_anchor_end_to_end(tmp_path):
    if not all(shutil.which(p) for p in ("makeblastdb", "primer3_core", "blastn", "blastdbcmd")):
        pytest.skip("requires real external tools")
    rng = random.Random(1841)
    while True:
        bases = [rng.choice("ACGT") for _ in range(850)]
        bases[418:424] = list("GAATTC")
        a = "".join(bases)
        if a.count("GAATTC") == 1:
            break
    b = a[:420] + "C" + a[421:]
    pos = 421
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">Chr4\n" + a + "\n", encoding="utf-8")
    pysam.faidx(str(fasta))
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Chr4,length=850>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        f"Chr4\t{pos}\t.\tA\tC\t60\tPASS\t.\tGT\t0/0\t1/1\n",
        encoding="utf-8",
    )
    vcfgz = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(vcf), str(vcfgz))
    pysam.tabix_index(str(vcfgz), preset="vcf")
    db = tmp_path / "reference-db"
    subprocess.run(
        ["makeblastdb", "-dbtype", "nucl", "-parse_seqids",
         "-in", str(fasta), "-out", str(db)],
        check=True, capture_output=True, text=True,
    )
    variants = read_sample_variants(str(vcfgz), "Chr4", 1, 850, "A", "B")
    window = source_window("FASTA", str(fasta), "Chr4", 1, 850, variants)
    assert window.a == a
    assert window.b == b
    alignment = align_alleles(window.a, window.b)
    assert len(alignment.variants) == 1
    rows, summary = design_caps(
        window.a, window.b, alignment, alignment.variants[0],
        ["EcoRI"], [str(db)], product_min=200, product_max=500,
        num_pairs=20, dcaps_mode="none", window=window, timeout=300,
    )
    assert summary["both_alleles_bind"] > 0
    assert rows
    assert any(row["参照検証採用可"] for row in rows), (
        "expected single-amplicon locus to pass strict anchored specificity"
    )
    assert all(row["特異性判定"] != "primary_reference_mismatch" for row in rows)
