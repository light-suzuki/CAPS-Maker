"""Reference-interval, two-sample genotype and panel-spacing regression tests."""
import pytest
import pysam

from caps_maker.batch import MarkerSite, choose_sites, discover_sites
from caps_maker.sources import (
    construct_window, fetch_reference, read_sample_variants, sample_names,
    source_window,
)


def demo_reference(tmp_path, prefix="reference"):
    fasta = tmp_path / (prefix + ".fa")
    # One A allele recognition site at reference genomic positions 61-66.
    a = "C" * 60 + "GAATTC" + "C" * 100
    b = "C" * 60 + "GACTTC" + "C" * 100
    fasta.write_text(">Chr4\n" + a + "\n")
    pysam.faidx(str(fasta))
    variant_pos = 63  # A allele A -> B allele C
    vcf = tmp_path / "strains.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Chr4,length=166>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        f"Chr4\t{variant_pos}\t.\tA\tC\t60\tPASS\t.\tGT\t0/0\t1/1\n"
    )
    vcfgz = tmp_path / "strains.vcf.gz"
    pysam.tabix_compress(str(vcf), str(vcfgz))
    pysam.tabix_index(str(vcfgz), preset="vcf")
    return str(fasta), str(vcfgz), a, b


def test_two_samples_and_reference_interval(tmp_path):
    fasta, vcf, a, b = demo_reference(tmp_path)
    assert sample_names(vcf) == ["A", "B"]
    assert fetch_reference("FASTA", fasta, "Chr4", 1, len(a)) == a
    variants = read_sample_variants(vcf, "Chr4", 1, len(a), "A", "B")
    assert len(variants) == 1
    window = source_window("FASTA", fasta, "Chr4", 1, len(a), variants)
    assert (window.a, window.b) == (a, b)
    assert window.inverse("A")[62] == 63
    assert window.inverse("B")[62] == 63


def test_reject_heterozygote_in_synthetic_allele(tmp_path):
    fasta, vcf, *_ = demo_reference(tmp_path)
    # A direct genotype change must be explicitly rejected, never assumed REF.
    altered = tmp_path / "het.vcf"
    altered.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=Chr4,length=166>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        "Chr4\t63\t.\tA\tC\t60\tPASS\t.\tGT\t0/1\t1/1\n"
    )
    bgz = tmp_path / "het.vcf.gz"
    pysam.tabix_compress(str(altered), str(bgz))
    pysam.tabix_index(str(bgz), preset="vcf")
    with pytest.raises(ValueError, match="ヘテロ接合"):
        read_sample_variants(str(bgz), "Chr4", 1, 166, "A", "B")


def test_dense_scan_detects_real_restriction_site(tmp_path):
    fasta, vcf, a, b = demo_reference(tmp_path)
    variants = read_sample_variants(vcf, "Chr4", 1, len(a), "A", "B")
    sites, stats = discover_sites(
        "FASTA", fasta, "Chr4", 1, len(a), variants, ["EcoRI"],
        include_dcaps=False,
    )
    assert stats["n_differential_variants"] == 1
    assert any(site.enzyme == "EcoRI" and site.site_position == 61 for site in sites)
    assert all(site.kind == "CAPS" for site in sites)


def test_fixed_spacing_and_dense_modes():
    candidates = [
        MarkerSite("Chr4", 120, "EcoRI", 120, "CAPS", "v1"),
        MarkerSite("Chr4", 150, "MseI", 150, "CAPS", "v1"),
        MarkerSite("Chr4", 510, "XbaI", 510, "CAPS", "v2"),
        MarkerSite("Chr4", 900, "EcoRI", 900, "CAPS", "v3"),
    ]
    dense, d_stats = choose_sites(candidates, "全酵素サイト", 1, 1000)
    assert len(dense) == 4
    sparse, s_stats = choose_sites(candidates, "一定間隔", 1, 1000,
                                   spacing_bp=400)
    assert len(sparse) == 3
    fixed, f_stats = choose_sites(candidates, "指定本数", 1, 1000,
                                  target_count=2)
    assert len(fixed) == 2
    limited, l_stats = choose_sites(candidates, "全酵素サイト", 1, 1000, limit=2)
    assert l_stats["truncated"] and len(limited) == 2
