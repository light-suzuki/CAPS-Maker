"""Reference-interval and two-sample VCF input for CAPS-Maker.

Input coordinates are ALWAYS reference-based, 1-based inclusive. Variant sites
must have confident, homozygous GT calls in both samples. Missing or mixed
genotypes are errors, never silently treated as the reference allele.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from primerblast_oss.caps import revcomp

IUPAC_BASES = frozenset("ACGT")
MAX_RECORDS = 250_000


@dataclass(frozen=True)
class SampleVariant:
    chrom: str
    pos: int
    reference: str
    allele_a: str
    allele_b: str

    @property
    def different(self) -> bool:
        return self.allele_a != self.allele_b


@dataclass
class LocusWindow:
    chrom: str
    start: int
    end: int
    reference: str
    a: str
    b: str
    ref_to_a: Dict[int, Optional[int]]
    ref_to_b: Dict[int, Optional[int]]
    variants: List[SampleVariant]
    source_kind: str = ""
    source_locator: str = ""

    def inverse(self, allele: str) -> Dict[int, int]:
        mapping = self.ref_to_a if allele.upper() == "A" else self.ref_to_b
        return {value: self.start + key for key, value in mapping.items() if value is not None}


class BlastGenome:
    """Use blastdbcmd as an indexed reference/secondary genome sequence provider."""
    def __init__(self, database: str):
        self.database = str(database)

    def fetch(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        return fetch_reference("BLAST DB", self.database, chrom, start, end, strand)


def fetch_reference(source: str, locator: str, chrom: str,
                    start: int, end: int, strand: str = "+") -> str:
    if not chrom or not locator or not 1 <= start <= end:
        raise ValueError("染色体・リファレンスと1-based座標を指定してください。")
    if end - start + 1 > 1_000_000:
        raise ValueError("1回に取得する参照配列は1 Mb以内にしてください。")
    if source == "FASTA":
        from primerblast_oss.genome import Genome
        fasta = str(Path(locator).expanduser())
        if not Path(fasta + ".fai").exists():
            raise ValueError("FASTAに.faiがありません。samtools faidxで索引を作成してください。")
        bases = Genome(fasta).fetch(chrom, start, end, strand)
    elif source == "BLAST DB":
        executable = shutil.which("blastdbcmd")
        if not executable:
            raise RuntimeError("BLAST+のblastdbcmdが見つかりません。")
        command = [
            executable, "-db", locator, "-entry", chrom,
            "-range", f"{start}-{end}", "-strand",
            "plus" if strand == "+" else "minus", "-outfmt", "%s",
        ]
        completed = subprocess.run(command, capture_output=True, text=True,
                                   timeout=60, check=False)
        if completed.returncode:
            raise RuntimeError("blastdbcmd: " + completed.stderr[-800:])
        bases = "".join(completed.stdout.split()).upper()
    else:
        raise ValueError("参照ソースはFASTAまたはBLAST DBを指定してください。")
    if len(bases) != end - start + 1:
        raise ValueError("リファレンスの配列長と指定区間が一致しません。"
                         "配列ID、索引、座標系を確認してください。")
    if set(bases) - IUPAC_BASES:
        raise ValueError("参照区間にN/曖昧塩基が存在します。別の領域を選択してください。")
    return bases


def sample_names(vcf_path: str) -> List[str]:
    import pysam

    with pysam.VariantFile(vcf_path) as handle:
        return list(handle.header.samples)


def read_sample_variants(vcf_path: str, chrom: str, start: int, end: int,
                         sample_a: str, sample_b: str,
                         max_records: int = MAX_RECORDS) -> List[SampleVariant]:
    import pysam

    if sample_a == sample_b:
        raise ValueError("2つの異なる系統（VCFのサンプル）を選んでください。")
    if not 1 <= start <= end:
        raise ValueError("開始と終了座標が不正です。")
    result = []
    with pysam.VariantFile(vcf_path) as handle:
        for name in (sample_a, sample_b):
            if name not in handle.header.samples:
                raise ValueError("VCFにサンプルがありません: " + name)
        try:
            records = handle.fetch(chrom, start - 1, end)
        except (ValueError, OSError) as exc:
            raise ValueError("VCFの区間検索に失敗しました。bgzip+tabix/CSI索引と配列IDを確認してください。") from exc
        for count, record in enumerate(records, 1):
            if count > max_records:
                raise ValueError("VCFレコード数が上限を超えました。区間を小さくしてください。")
            pos = int(record.pos)
            ref = str(record.ref).upper()
            if pos < start or pos + len(ref) - 1 > end:
                # Overlapping records at window edges must be handled explicitly.
                raise ValueError(
                    f"VCFの変異{chrom}:{pos}が解析ウィンドウ境界をまたいでいます。"
                    "参照区間を広げるか境界を調整してください。"
                )
            if not ref or set(ref) - IUPAC_BASES:
                raise ValueError(f"{chrom}:{pos}のREFが通常のDNA配列ではありません。")
            filters = set(record.filter.keys())
            if filters and filters != {"PASS"}:
                raise ValueError(
                    f"{chrom}:{pos}に未通過FILTER（{','.join(sorted(filters))}）があります。"
                    "対象を確認・正規化してください。"
                )
            alleles = tuple(str(v).upper() for v in record.alleles)
            calls = []
            for sample in (sample_a, sample_b):
                gt = record.samples[sample].get("GT")
                if not gt or any(v is None for v in gt) or len(set(gt)) != 1:
                    raise ValueError(
                        f"{sample}の{chrom}:{pos}はヘテロ接合/欠損のため単一配列に"
                        "再構成できません。ホモ接合サンプルか位相済みハプロタイプを使用してください。"
                    )
                index = int(gt[0])
                if not 0 <= index < len(alleles):
                    raise ValueError(f"{sample}の{chrom}:{pos}のGTが不正です。")
                allele = alleles[index]
                if not allele or set(allele) - IUPAC_BASES or len(allele) > 200:
                    raise ValueError(
                        f"{sample}の{chrom}:{pos}は記号アレル/長い変異です。"
                        "現行の小変異モードは対応していません。"
                    )
                calls.append(allele)
            if calls[0] != ref or calls[1] != ref:
                result.append(SampleVariant(chrom, pos, ref, calls[0], calls[1]))
    result.sort(key=lambda x: x.pos)
    for left, right in zip(result, result[1:]):
        if left.pos + len(left.reference) > right.pos:
            raise ValueError("VCFに重複する変異があります。bcftools norm等で整形してください。")
    return result


def _apply_sample(reference: str, window_start: int,
                  variants: Sequence[SampleVariant], sample: str
                  ) -> Tuple[str, Dict[int, Optional[int]]]:
    output: List[str] = []
    mapping: Dict[int, Optional[int]] = {}
    ref_cursor = allele_cursor = 0
    for variant in variants:
        start = variant.pos - window_start
        if start < ref_cursor or start + len(variant.reference) > len(reference):
            raise ValueError("変異が区間外または互いに重複しています。")
        if reference[start:start + len(variant.reference)] != variant.reference:
            raise ValueError(f"{variant.chrom}:{variant.pos}でFASTAとVCF REFが不一致です。")
        prefix = reference[ref_cursor:start]
        output.append(prefix)
        for index in range(ref_cursor, start):
            mapping[index] = allele_cursor + index - ref_cursor
        allele_cursor += len(prefix)
        allele = variant.allele_a if sample == "A" else variant.allele_b
        output.append(allele)
        for i in range(len(variant.reference)):
            mapping[start + i] = allele_cursor + i if i < len(allele) else None
        allele_cursor += len(allele)
        ref_cursor = start + len(variant.reference)
    output.append(reference[ref_cursor:])
    for index in range(ref_cursor, len(reference)):
        mapping[index] = allele_cursor + index - ref_cursor
    return "".join(output), mapping


def construct_window(chrom: str, start: int, reference: str,
                     variants: Sequence[SampleVariant]) -> LocusWindow:
    end = start + len(reference) - 1
    subset = [
        v for v in variants if start <= v.pos and
        v.pos + len(v.reference) - 1 <= end
    ]
    # If a record overlaps a boundary, do not silently replace its contribution.
    for v in variants:
        if v.pos < start <= v.pos + len(v.reference) - 1 or (
                v.pos <= end < v.pos + len(v.reference) - 1):
            raise ValueError(f"変異{chrom}:{v.pos}が配列ウィンドウ境界をまたぎます。")
    a, ref_to_a = _apply_sample(reference, start, subset, "A")
    b, ref_to_b = _apply_sample(reference, start, subset, "B")
    return LocusWindow(chrom, start, end, reference, a, b, ref_to_a, ref_to_b, subset)


def source_window(source: str, locator: str, chrom: str, start: int, end: int,
                  variants: Sequence[SampleVariant]) -> LocusWindow:
    window = construct_window(
        chrom, start, fetch_reference(source, locator, chrom, start, end),
        variants)
    window.source_kind = source
    window.source_locator = locator
    return window
