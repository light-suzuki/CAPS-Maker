"""Anchor-aware BLAST interpretation for CAPS and dCAPS primers."""
from __future__ import annotations

import json
import subprocess
import sys
from typing import Dict, Optional, Sequence, Tuple

from .experiments import GelProfile, off_target_verdict
from .sources import BlastGenome, LocusWindow


def reference_anchor(window: Optional[LocusWindow],
                     left: int, right: int) -> Optional[Tuple[str, Tuple[int, int]]]:
    if window is None:
        return None
    mapping = window.inverse("A")
    if left not in mapping or right not in mapping:
        return None
    a, b = mapping[left], mapping[right]
    if b <= a:
        return None
    return window.chrom, (a, b)


def run_primer_check(forward: str, reverse: str, databases: Sequence[str],
                     max_product: int = 10000, timeout: int = 600,
                     max_target_seqs: int = 5000) -> Sequence[dict]:
    if not databases:
        raise ValueError("BLAST DBが指定されていません。")
    cmd = [
        sys.executable, "-m", "primerblast_oss", "check",
        "--forward", forward, "--reverse", reverse,
        "--max-product", str(max_product),
        "--max-target-seqs", str(max_target_seqs),
        "--format", "json",
    ]
    for database in databases:
        cmd.extend(["--db", database])
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                            check=False)
    if result.returncode:
        raise RuntimeError("primerblast-oss check: " + (result.stderr or result.stdout)[-1500:])
    try:
        return json.loads(result.stdout)["results"]
    except (ValueError, KeyError, TypeError) as exc:
        raise RuntimeError("PrimerBlastOssの特異性レポートが不正です。") from exc


def gel_specificity(
    *,
    evidence: Sequence[dict],
    databases: Sequence[str],
    reference: Optional[LocusWindow],
    left: int, right: int,
    forward: str, reverse: str,
    enzyme: str,
    fragments_a: Sequence[int],
    fragments_b: Sequence[int],
    gel_profile: GelProfile,
    policy: str = "strict",
    genome_fasta_by_db: Optional[Dict[str, str]] = None,
) -> Dict:
    anchor = reference_anchor(reference, left, right)
    if not anchor:
        return {
            "status": "no_reference_anchor",
            "accepted": False,
            "note": "参照座標で標的座位を確認できないため、自動採用不可",
        }
    from primerblast_oss.genome import Genome

    genome_fasta_by_db = genome_fasta_by_db or {}
    genomes = {}
    # For gel-tolerant, every off-target product must be reconstructed from
    # its actual DB. Never substitute a different assembly's reference.
    if policy == "gel-tolerant":
        for db in databases:
            fasta = genome_fasta_by_db.get(db)
            try:
                genomes[db] = Genome(fasta) if fasta else BlastGenome(db)
            except (OSError, RuntimeError, ValueError):
                # Underlying off-target digest then becomes unresolved.
                pass
    return off_target_verdict(
        rows_by_db=evidence,
        genomes_by_db=genomes,
        primers={"F": forward, "R": reverse},
        enzyme=enzyme,
        a_fragments=fragments_a,
        b_fragments=fragments_b,
        profile=gel_profile,
        target_chrom=anchor[0],
        target_interval=anchor[1],
        policy=policy,
    )


def parse_fragment_column(text: str) -> Sequence[int]:
    lengths = []
    for token in text.split(" + "):
        if "×" in token:
            size, copies = token.split("×")
            lengths.extend([int(size)] * int(copies))
        elif token:
            lengths.append(int(token))
    return lengths
