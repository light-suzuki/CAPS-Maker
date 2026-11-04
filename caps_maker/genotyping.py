"""Display only reference-validated, genotype-informative PCR assays."""
from __future__ import annotations


def usable_marker(row: dict, three_genotypes: bool = True,
                  policy: str = "strict") -> bool:
    # Rows from design_caps/analyze_primers already require binding in A and B.
    if not (row.get("Forward(5'-3')") and row.get("Reverse(5'-3')")
            and row.get("A産物(bp)", 0) > 0 and row.get("B産物(bp)", 0) > 0):
        return False
    if row.get("参照検証採用可") is not True or row.get("ゲル判別") != "可能":
        return False
    status = row.get("特異性判定")
    if status != "strict_pass" and not (policy == "gel-tolerant" and status == "gel_tolerant_pass"):
        return False
    kind = row.get("マーカー型", "")
    if three_genotypes:
        return "共優性" in kind
    return "共優性" in kind or kind.startswith(("優性型", "2群識別"))
