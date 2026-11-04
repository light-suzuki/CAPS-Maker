"""Genotyping use of validated markers.

Selection gate, bench sheet, scoring sheet, call QC and genotype-table export.
Everything here predicts or checks *gel patterns*; it never infers a genotype
that the gel cannot show, and segregation tests are descriptive QC only.
"""
from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

PAIRS = ("AA_vs_AB", "AA_vs_BB", "AB_vs_BB")

# What each population needs the gel to separate. The first two labels are the
# historical settings values and must stay unchanged for saved templates.
GOALS: Dict[str, Tuple[str, ...]] = {
    "AA・AB・BBを区別（ヘテロも判別）": PAIRS,
    "AAとBBを区別（ヘテロ判別は不要）": ("AA_vs_BB",),
    "AAとABを区別（Aへの戻し交雑）": ("AA_vs_AB",),
    "ABとBBを区別（Bへの戻し交雑）": ("AB_vs_BB",),
}
GOAL_CAPTIONS = (
    "F2・F2:3など分離集団。共優性マーカーのみ",
    "RIL・DH・固定系統、または優性型でもよいF2",
    "BC集団（反復親A）。AAとABが別パターンのもの",
    "BC集団（反復親B）。ABとBBが別パターンのもの",
)
GOAL_SHORT = {
    "AA・AB・BBを区別（ヘテロも判別）": "F2：AA・AB・BB",
    "AAとBBを区別（ヘテロ判別は不要）": "RIL・DH：AA・BB",
    "AAとABを区別（Aへの戻し交雑）": "BC→A：AA・AB",
    "ABとBBを区別（Bへの戻し交雑）": "BC→B：AB・BB",
}


def separable_pairs(marker_type: str) -> frozenset:
    """Genotype pairs a marker type separates on the gel (see experiments.assess_gel)."""
    kind = str(marker_type or "")
    if "共優性" in kind:
        return frozenset(PAIRS)
    if kind.startswith("優性型（ABとBBが同じ）"):
        return frozenset({"AA_vs_AB", "AA_vs_BB"})
    if kind.startswith("優性型（AAとABが同じ）"):
        return frozenset({"AB_vs_BB", "AA_vs_BB"})
    if kind.startswith("2群識別"):
        return frozenset({"AA_vs_BB"})
    return frozenset()


def usable_marker(row: dict, three_genotypes: bool = True,
                  policy: str = "strict", goal: Optional[str] = None) -> bool:
    # Rows from design_caps/analyze_primers already require binding in A and B.
    if not (row.get("Forward(5'-3')") and row.get("Reverse(5'-3')")
            and row.get("A産物(bp)", 0) > 0 and row.get("B産物(bp)", 0) > 0):
        return False
    if row.get("参照検証採用可") is not True or row.get("ゲル判別") != "可能":
        return False
    status = row.get("特異性判定")
    if status != "strict_pass" and not (policy == "gel-tolerant" and status == "gel_tolerant_pass"):
        return False
    if goal is None:
        goal = next(iter(GOALS)) if three_genotypes else "AAとBBを区別（ヘテロ判別は不要）"
    if goal not in GOALS:
        raise ValueError("未対応の判別目的です。")
    return set(GOALS[goal]) <= separable_pairs(row.get("マーカー型", ""))


# ---------------------------------------------------------------------------
# Marker identity and bench sheet

def marker_position(row: Mapping) -> Optional[int]:
    """Target coordinate when known, else the amplicon midpoint on the reference."""
    for key in ("標的座標(1-based)",):
        if row.get(key) not in (None, ""):
            return int(row[key])
    start, end = row.get("参照産物開始(1-based)"), row.get("参照産物終了(1-based)")
    if start not in (None, "") and end not in (None, ""):
        return (int(start) + int(end)) // 2
    return None


def marker_name(row: Mapping, index: int = 0) -> str:
    chrom = row.get("参照染色体") or row.get("染色体")
    position = marker_position(row)
    enzyme = row.get("酵素", "")
    kind = "d" if row.get("設計方式") == "dCAPS" else ""
    if chrom and position:
        return f"{chrom}_{position}_{kind}{enzyme}"
    return f"M{index + 1:03d}_{kind}{enzyme}"


def unique_marker_names(rows: Sequence[Mapping]) -> List[str]:
    names, seen = [], {}
    for i, row in enumerate(rows):
        base = marker_name(row, i)
        seen[base] = seen.get(base, 0) + 1
        names.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return names


def _bands(sizes: Iterable[int]) -> str:
    values = sorted({int(v) for v in sizes}, reverse=True)
    return " + ".join(str(v) for v in values) if values else "（可視バンドなし）"


def expected_lanes(row: Mapping, profile) -> Dict[str, List[int]]:
    """Visible band sizes per genotype under the user's gel profile."""
    from .experiments import assess_gel
    from .validation import parse_fragment_column

    gel = assess_gel(parse_fragment_column(row["A断片(bp)"]),
                     parse_fragment_column(row["B断片(bp)"]), profile)
    lanes = gel["observed_lanes"]
    return {"AA": list(lanes["AA"]), "AB": list(lanes["AB"]), "BB": list(lanes["BB"]),
            "gel_percent": gel["gel_percent"], "ladder": gel["ladder"]}


def bench_sheet(rows: Sequence[Mapping], profile, goal: str) -> List[dict]:
    """One row per marker: what to run and what each genotype looks like."""
    out = []
    for name, row in zip(unique_marker_names(rows), rows):
        lanes = expected_lanes(row, profile)
        pairs = separable_pairs(row.get("マーカー型", ""))
        readable = calls_readable(pairs)
        out.append({
            "マーカー名": name,
            "染色体": row.get("参照染色体") or row.get("染色体", ""),
            "位置(bp)": marker_position(row) or "",
            "設計方式": row.get("設計方式", ""),
            "酵素": row.get("酵素", ""),
            "Forward(5'-3')": row.get("Forward(5'-3')", ""),
            "Reverse(5'-3')": row.get("Reverse(5'-3')", ""),
            "Tm F(°C)": row.get("Tm F(°C)", ""),
            "Tm R(°C)": row.get("Tm R(°C)", ""),
            "PCR産物 A/B(bp)": f"{row.get('A産物(bp)', '')} / {row.get('B産物(bp)', '')}",
            "AA（親A型）": _bands(lanes["AA"]),
            "AB（ヘテロ）": _bands(lanes["AB"]),
            "BB（親B型）": _bands(lanes["BB"]),
            "ゲル": f"{lanes['gel_percent']}%・{lanes['ladder']}ラダー",
            "マーカー型": row.get("マーカー型", ""),
            "記録できる型": "・".join(CALL_LABEL[c] for c in readable),
            "目的に適合": set(GOALS[goal]) <= pairs,
            "特異性区分": row.get("特異性区分", ""),
        })
    return out


# ---------------------------------------------------------------------------
# Calls

CALL_ALIASES = {
    "A": "A", "AA": "A", "0": "A", "P1": "A",
    "H": "H", "AB": "H", "BA": "H", "1": "H", "HET": "H", "F1": "H",
    "B": "B", "BB": "B", "2": "B", "P2": "B",
    # R/qtl dominant codes: C = not BB (AA or AB), D = not AA (AB or BB)
    "C": "C", "A-": "C", "NOTB": "C", "D": "D", "B-": "D", "NOTA": "D",
    "": "-", "-": "-", "NA": "-", "N/A": "-", "N": "-", ".": "-", "?": "-", "NAN": "-",
}
CALL_LABEL = {"A": "AA", "H": "AB", "B": "BB", "C": "A-", "D": "B-", "-": ""}


def normalize_call(value) -> str:
    key = str(value if value is not None else "").strip().upper().replace(" ", "")
    if key.endswith(".0") and key[:-2].isdigit():
        key = key[:-2]
    if key not in CALL_ALIASES:
        raise ValueError(f"解釈できない遺伝子型「{value}」")
    return CALL_ALIASES[key]


def calls_readable(pairs: frozenset) -> List[str]:
    """Calls a gel can actually support for a marker separating these pairs."""
    if pairs >= set(PAIRS):
        return ["A", "H", "B"]
    if pairs == {"AA_vs_AB", "AA_vs_BB"}:
        return ["A", "D"]
    if pairs == {"AB_vs_BB", "AA_vs_BB"}:
        return ["B", "C"]
    if pairs == {"AA_vs_BB"}:
        return ["A", "B"]
    return []


def call_from_bands(observed: Sequence[float], lanes: Mapping[str, Sequence[int]],
                    tolerance: float = 0.10, min_log10_gap: float = 0.04) -> List[str]:
    """Genotypes whose predicted visible bands match the bands read off a gel.

    Each predicted band (bands closer than the gel's resolving gap count as
    one) must have an observed band within ``tolerance`` (relative size
    error) and every observed band must be explained by a predicted band.
    """
    if not 0 < tolerance < 1:
        raise ValueError("許容誤差は0～100%未満で指定してください。")
    seen = [float(v) for v in observed if float(v) > 0]
    if not seen:
        return []
    limit = math.log10(1 + tolerance)

    def merged(sizes):
        groups: List[List[int]] = []
        for size in sorted(int(s) for s in sizes):
            if groups and math.log10(size / groups[-1][-1]) < min_log10_gap:
                groups[-1].append(size)
            else:
                groups.append([size])
        return groups

    def near(size, group):
        return any(abs(math.log10(size / g)) <= limit for g in group)

    matches = []
    for genotype in ("AA", "AB", "BB"):
        groups = merged(lanes.get(genotype, []))
        if not groups:
            continue
        if all(any(near(o, g) for o in seen) for g in groups) and \
                all(any(near(o, g) for g in groups) for o in seen):
            matches.append(genotype)
    return matches


# ---------------------------------------------------------------------------
# Scoring sheet: one CSV that is both the gel-reading key and the record

INFO_ROLES = ("Chromosome", "Position", "Enzyme", "Expected_AA", "Expected_AB",
              "Expected_BB", "Readable")
CONTROL_ROLES = {"Parent_A": "A", "Parent_B": "B", "F1": "H"}


def scoring_sheet(rows: Sequence[Mapping], profile, samples: Sequence[str],
                  include_f1: bool = True) -> bytes:
    names = unique_marker_names(rows)
    sheet = bench_sheet(rows, profile, next(iter(GOALS)))
    table = [["Role", "Line", *names]]
    table.append(["Chromosome", "", *[s["染色体"] for s in sheet]])
    table.append(["Position", "", *[s["位置(bp)"] for s in sheet]])
    table.append(["Enzyme", "", *[s["酵素"] for s in sheet]])
    table.append(["Expected_AA", "", *[s["AA（親A型）"] for s in sheet]])
    table.append(["Expected_AB", "", *[s["AB（ヘテロ）"] for s in sheet]])
    table.append(["Expected_BB", "", *[s["BB（親B型）"] for s in sheet]])
    table.append(["Readable", "", *["/".join(calls_readable(separable_pairs(r.get("マーカー型", ""))))
                                    for r in rows]])
    table.append(["Parent_A", "Parent_A", *[""] * len(names)])
    table.append(["Parent_B", "Parent_B", *[""] * len(names)])
    if include_f1:
        table.append(["F1", "F1", *[""] * len(names)])
    for sample in samples:
        table.append(["Progeny", sample, *[""] * len(names)])
    buffer = io.StringIO(newline="")
    csv.writer(buffer, lineterminator="\n").writerows(table)
    return buffer.getvalue().encode("utf-8-sig")


@dataclass
class ScoredSheet:
    markers: List[str]
    info: Dict[str, Dict[str, str]]        # role -> marker -> value
    individuals: List[Tuple[str, str]]     # (role, line)
    calls: Dict[Tuple[str, str], Dict[str, str]]  # (role, line) -> marker -> code
    problems: List[dict]


def read_scoring_sheet(data: bytes) -> ScoredSheet:
    text = data.decode("utf-8-sig")
    dialect = "\t" if text.split("\n", 1)[0].count("\t") > text.split("\n", 1)[0].count(",") else ","
    rows = [r for r in csv.reader(io.StringIO(text), delimiter=dialect) if any(c.strip() for c in r)]
    if not rows or [c.strip() for c in rows[0][:2]] != ["Role", "Line"]:
        raise ValueError("1行目は Role, Line, マーカー名… の記録シート形式にしてください。")
    header = [c.strip() for c in rows[0]]
    markers = [m for m in header[2:] if m]
    if len(set(markers)) != len(markers):
        raise ValueError("マーカー名が重複しています。")
    info: Dict[str, Dict[str, str]] = {}
    individuals, calls, problems = [], {}, []
    for line_no, raw in enumerate(rows[1:], start=2):
        raw = raw + [""] * (len(header) - len(raw))
        role, line = raw[0].strip(), raw[1].strip()
        values = dict(zip(header[2:], (v.strip() for v in raw[2:])))
        if role in INFO_ROLES:
            info[role] = {m: values.get(m, "") for m in markers}
            continue
        if role not in (*CONTROL_ROLES, "Progeny"):
            problems.append({"行": line_no, "個体": line, "マーカー": "", "内容": f"不明なRole「{role}」の行を無視"})
            continue
        key = (role, line or role)
        if key in calls:
            raise ValueError(f"{line_no}行目: 個体「{key[1]}」が重複しています。")
        individuals.append(key)
        calls[key] = {}
        for marker in markers:
            try:
                calls[key][marker] = normalize_call(values.get(marker, ""))
            except ValueError as exc:
                calls[key][marker] = "-"
                problems.append({"行": line_no, "個体": key[1], "マーカー": marker,
                                 "内容": str(exc) + "（欠測として扱う）"})
    return ScoredSheet(markers, info, individuals, calls, problems)


# ---------------------------------------------------------------------------
# QC

POPULATIONS = {
    "F2（自殖・F1×F1）": "F2",
    "戻し交雑（反復親A）": "BC_A",
    "戻し交雑（反復親B）": "BC_B",
    "RIL・DH（ホモ固定）": "RIL",
}
POPULATION_FOR_GOAL = {
    "AA・AB・BBを区別（ヘテロも判別）": "F2（自殖・F1×F1）",
    "AAとBBを区別（ヘテロ判別は不要）": "RIL・DH（ホモ固定）",
    "AAとABを区別（Aへの戻し交雑）": "戻し交雑（反復親A）",
    "ABとBBを区別（Bへの戻し交雑）": "戻し交雑（反復親B）",
}
POSSIBLE = {"F2": {"A", "H", "B", "C", "D"}, "BC_A": {"A", "H", "C"},
            "BC_B": {"H", "B", "D"}, "RIL": {"A", "B", "H"}}


def _chi2_p(stat: float, df: int) -> float:
    if df == 1:
        return math.erfc(math.sqrt(stat / 2))
    if df == 2:
        return math.exp(-stat / 2)
    raise ValueError("df must be 1 or 2")


def segregation_test(counts: Mapping[str, int], population: str) -> Optional[dict]:
    """Chi-square goodness of fit for the classes the calls actually contain."""
    pop = POPULATIONS.get(population, population)
    observed = {k: v for k, v in counts.items() if v and k in "AHBCD"}
    classes = set(observed)
    if pop == "F2":
        if classes <= {"A", "H", "B"} and len(classes) >= 2:
            expected = {"A": 1, "H": 2, "B": 1}
        elif classes <= {"A", "D"}:
            expected = {"A": 1, "D": 3}
        elif classes <= {"B", "C"}:
            expected = {"C": 3, "B": 1}
        else:
            return None
    elif pop == "BC_A" and classes <= {"A", "H"}:
        expected = {"A": 1, "H": 1}
    elif pop == "BC_B" and classes <= {"H", "B"}:
        expected = {"H": 1, "B": 1}
    elif pop == "RIL" and classes <= {"A", "B"}:
        expected = {"A": 1, "B": 1}
    else:
        return None
    n = sum(observed.get(k, 0) for k in expected)
    if n < 5:
        return None
    total = sum(expected.values())
    stat = sum((observed.get(k, 0) - n * w / total) ** 2 / (n * w / total)
               for k, w in expected.items())
    df = len(expected) - 1
    ratio = ":".join(str(w) for w in expected.values())
    return {"期待比": f"{'/'.join(CALL_LABEL[k] for k in expected)} = {ratio}",
            "χ²": round(stat, 3), "自由度": df, "p値": _chi2_p(stat, df)}


def qc_report(sheet: ScoredSheet, population: str) -> Tuple[List[dict], List[dict]]:
    """Per-marker QC table and a list of individual problems to look at."""
    pop = POPULATIONS.get(population, population)
    problems = list(sheet.problems)
    readable = {}
    for marker in sheet.markers:
        text = sheet.info.get("Readable", {}).get(marker, "")
        readable[marker] = {normalize_call(c) for c in text.split("/") if c.strip()} or None
    progeny = [k for k in sheet.individuals if k[0] == "Progeny"]
    table = []
    for marker in sheet.markers:
        for role, want in CONTROL_ROLES.items():
            for key in (k for k in sheet.individuals if k[0] == role):
                got = sheet.calls[key][marker]
                if got == "-":
                    continue
                allowed = {want} | ({"C"} if want in "AH" else set()) | ({"D"} if want in "HB" else set())
                if got not in allowed:
                    problems.append({"行": "", "個体": key[1], "マーカー": marker,
                                     "内容": f"コントロール不一致：{role}は{CALL_LABEL[want]}の想定、記録は{CALL_LABEL[got]}"})
        counts = {c: 0 for c in "AHBCD-"}
        for key in progeny:
            code = sheet.calls[key][marker]
            counts[code] += 1
            ok = readable[marker]
            if code != "-" and ok is not None and code not in ok:
                problems.append({"行": "", "個体": key[1], "マーカー": marker,
                                 "内容": f"このマーカーのゲルでは{CALL_LABEL[code]}を判別できません（記録可：{'/'.join(CALL_LABEL[c] for c in sorted(ok))}）"})
            elif code != "-" and code not in POSSIBLE[pop]:
                problems.append({"行": "", "個体": key[1], "マーカー": marker,
                                 "内容": f"{population}では通常現れない{CALL_LABEL[code]}"})
        n = len(progeny)
        test = segregation_test(counts, population)
        table.append({
            "マーカー名": marker,
            "染色体": sheet.info.get("Chromosome", {}).get(marker, ""),
            "位置(bp)": sheet.info.get("Position", {}).get(marker, ""),
            "AA": counts["A"], "AB": counts["H"], "BB": counts["B"],
            "A-（BB以外）": counts["C"], "B-（AA以外）": counts["D"],
            "欠測": counts["-"],
            "欠測率": round(counts["-"] / n, 3) if n else None,
            "期待比": test["期待比"] if test else "",
            "χ²": test["χ²"] if test else None,
            "p値": round(test["p値"], 4) if test else None,
            "分離のずれ": ("p<0.05" if test and test["p値"] < 0.05 else "") if test else "検定対象外",
        })
    return table, problems


# ---------------------------------------------------------------------------
# Export

def genotype_table_bytes(sheet: ScoredSheet) -> bytes:
    """Wide table: Role, Line, Phenotype, markers; AA/AB/BB, A-/B- for dominant calls."""
    header = ["Role", "Line", "Phenotype", *sheet.markers]
    lines = [header,
             ["Chromosome", "", "", *[sheet.info.get("Chromosome", {}).get(m, "") for m in sheet.markers]],
             ["Position", "", "", *[sheet.info.get("Position", {}).get(m, "") for m in sheet.markers]]]
    for role, line in sheet.individuals:
        lines.append([role, line, "", *[CALL_LABEL[sheet.calls[(role, line)][m]] for m in sheet.markers]])
    buffer = io.StringIO(newline="")
    csv.writer(buffer, lineterminator="\n").writerows(lines)
    return buffer.getvalue().encode("utf-8-sig")


def rqtl_csv_bytes(sheet: ScoredSheet) -> bytes:
    """R/qtl "csv" layout (progeny only): A/H/B/C/D, '-' missing, positions in Mb."""
    def mb(value):
        try:
            return f"{int(float(value)) / 1e6:.6f}"
        except (TypeError, ValueError):
            return ""
    lines = [["id", *sheet.markers],
             ["", *[sheet.info.get("Chromosome", {}).get(m, "") for m in sheet.markers]],
             ["", *[mb(sheet.info.get("Position", {}).get(m, "")) for m in sheet.markers]]]
    for role, line in sheet.individuals:
        if role == "Progeny":
            lines.append([line, *[sheet.calls[(role, line)][m] for m in sheet.markers]])
    buffer = io.StringIO(newline="")
    csv.writer(buffer, lineterminator="\n").writerows(lines)
    return buffer.getvalue().encode("utf-8")


# ---- marker density and narrowing ---------------------------------------------

DENSITY_MODES = ("一定間隔", "指定本数", "全座位")


def _chrom(row: Mapping) -> str:
    return str(row.get("参照染色体") or row.get("染色体") or "")


def locus_span(rows: Sequence[Mapping]) -> Dict[str, Tuple[int, int]]:
    """Chromosome -> (min, max) position of the candidate loci."""
    spans: Dict[str, Tuple[int, int]] = {}
    for row in rows:
        pos = marker_position(row)
        if pos is None:
            continue
        lo, hi = spans.get(_chrom(row), (pos, pos))
        spans[_chrom(row)] = (min(lo, pos), max(hi, pos))
    return spans


def select_by_density(rows: Sequence[Mapping], mode: str, chrom: str, start: int, end: int,
                      spacing_bp: int = 100_000, count: int = 10) -> Tuple[List[int], List[dict]]:
    """Pick at most one assay per bin inside chrom:start-end.

    Candidates are ranked; alternative pairs at one locus are alternatives, not
    extra markers, so only the best-ranked row per locus competes. Within a bin
    a co-dominant marker beats a dominant one, then the locus nearest the bin
    centre wins. Empty bins are reported, never filled.
    """
    if mode not in DENSITY_MODES:
        raise ValueError("マーカー密度の指定方法が不明です。")
    if end < start:
        raise ValueError("区間の終了が開始より前です。")
    if spacing_bp < 1 or count < 1:
        raise ValueError("間隔と本数は1以上にしてください。")
    best: Dict[int, int] = {}
    for i, row in enumerate(rows):
        pos = marker_position(row)
        if pos is None or _chrom(row) != chrom or not start <= pos <= end:
            continue
        best.setdefault(pos, i)          # rows are ranked: first one wins
    loci = sorted(best.items())
    if mode == "全座位":
        return [i for _, i in loci], []
    total = end - start + 1
    width = spacing_bp if mode == "一定間隔" else max(1, -(-total // count))
    chosen, empty = [], []
    for b in range(-(-total // width)):
        lo = start + b * width
        hi = min(end, lo + width - 1)
        centre = (lo + hi) / 2
        inside = [(pos, i) for pos, i in loci if lo <= pos <= hi]
        if not inside:
            empty.append({"区画": b + 1, "開始(bp)": lo, "終了(bp)": hi})
            continue
        pos, i = min(inside, key=lambda item: (
            not str(rows[item[1]].get("マーカー型", "")).startswith("共優性"),
            abs(item[0] - centre), item[1]))
        chosen.append(i)
    return chosen, empty


_CALL_SETS = {"A": {"A"}, "H": {"H"}, "B": {"B"}, "C": {"A", "H"}, "D": {"H", "B"}}


def _ordered_markers(sheet: ScoredSheet) -> List[Tuple[str, int, str]]:
    out = []
    for m in sheet.markers:
        try:
            pos = int(float(sheet.info.get("Position", {}).get(m, "")))
        except ValueError:
            continue
        out.append((sheet.info.get("Chromosome", {}).get(m, ""), pos, m))
    return sorted(out)


def recombination_report(sheet: ScoredSheet) -> Tuple[List[dict], List[dict]]:
    """Breakpoints per progeny and recombinant counts per adjacent marker interval.

    A breakpoint is placed between consecutive scored markers whose calls cannot
    be the same genotype (A- vs AB is compatible, A- vs BB is not). Missing
    calls are skipped, so a breakpoint interval can span several markers; a
    change hidden behind a dominant call is not counted.
    """
    order = _ordered_markers(sheet)
    progeny = [k for k in sheet.individuals if k[0] == "Progeny"]
    adjacent = {}
    for (c1, p1, m1), (c2, p2, m2) in zip(order, order[1:]):
        if c1 == c2:
            adjacent[(m1, m2)] = {"染色体": c1, "左マーカー": m1, "右マーカー": m2,
                                  "区間(bp)": f"{p1:,}–{p2:,}", "幅(kb)": round((p2 - p1) / 1000, 1),
                                  "組換え個体": 0, "うち区間が確定": 0}
    names = [o[2] for o in order]
    breaks = []
    for key in progeny:
        last = None
        for chrom, pos, m in order:
            code = sheet.calls[key][m]
            if code == "-":
                continue
            if last and last[0] == chrom and not (_CALL_SETS[last[3]] & _CALL_SETS[code]):
                breaks.append({"個体": key[1], "染色体": chrom, "左マーカー": last[2],
                               "右マーカー": m, "区間(bp)": f"{last[1]:,}–{pos:,}",
                               "型の変化": f"{CALL_LABEL[last[3]]}→{CALL_LABEL[code]}"})
                span = names[names.index(last[2]):names.index(m) + 1]
                for pair in zip(span, span[1:]):
                    adjacent[pair]["組換え個体"] += 1
                    if len(span) == 2:
                        adjacent[pair]["うち区間が確定"] += 1
            last = (chrom, pos, m, code)
    return breaks, list(adjacent.values())
