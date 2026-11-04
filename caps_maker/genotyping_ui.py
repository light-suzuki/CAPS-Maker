"""Genotyping tab: bench sheet, scoring sheet, call QC and genotype-table export."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from . import ui
from .core import table_bytes
from .genotyping import (
    CALL_LABEL, DENSITY_MODES, POPULATION_FOR_GOAL, POPULATIONS, bench_sheet, genotype_table_bytes,
    locus_span, qc_report, read_scoring_sheet, recombination_report, rqtl_csv_bytes, scoring_sheet,
    marker_position, select_by_density, unique_marker_names,
)


def _sample_names() -> list[str]:
    with st.container(key="seg_gt_sample_source"):
        source = st.radio("個体名の指定", ["連番で作る", "一覧を貼り付け"], horizontal=True,
                          key="gt_sample_source")
    if source == "連番で作る":
        c1, c2 = st.columns(2)
        with c1:
            prefix = st.text_input("個体名の接頭辞", value="F2_", key="gt_prefix")
        with c2:
            count = st.number_input("個体数", 1, 5000, 96, key="gt_count")
        width = max(3, len(str(int(count))))
        return [f"{prefix}{i:0{width}d}" for i in range(1, int(count) + 1)]
    text = st.text_area("個体名（1行に1つ）", height=140, key="gt_sample_text",
                        placeholder="F2_001\nF2_002\n…")
    names = [line.strip() for line in text.splitlines() if line.strip()]
    if len(set(names)) != len(names):
        st.error("個体名が重複しています。")
        return []
    return names


def _density_picker(markers: list[dict], names: list[str]) -> list[int]:
    """Region + density controls, then an editable multiselect seeded from them."""
    spans = locus_span(markers)
    if not spans:
        return st.multiselect("記録するマーカー", list(range(len(markers))), key="gt_markers",
                              format_func=lambda i: names[i])
    ui.subheading("範囲と密度")
    chroms = sorted(spans)
    c1, c2, c3 = st.columns([1.2, 1, 1])
    with c1:
        chrom = st.selectbox("染色体", chroms, key="gt_chrom")
    lo, hi = spans[chrom]
    with c2:
        start = st.number_input("開始(bp)", 1, None, lo, 1000, key=f"gt_start_{chrom}")
    with c3:
        end = st.number_input("終了(bp)", 1, None, hi, 1000, key=f"gt_end_{chrom}")
    d1, d2 = st.columns([1.6, 1])
    with d1, st.container(key="seg_gt_density"):
        mode = st.radio("密度の決め方", DENSITY_MODES, horizontal=True, key="gt_density")
    with d2:
        if mode == "一定間隔":
            spacing_kb = st.number_input("マーカー間隔(kb)", 0.1, 100_000.0, 100.0, 10.0,
                                         key="gt_spacing_kb", format="%g")
            count = 1
        elif mode == "指定本数":
            count = st.number_input("本数", 1, 1000, 10, key="gt_count_markers")
            spacing_kb = 1.0
        else:
            spacing_kb, count = 1.0, 1
            st.caption("範囲内の全座位から各座位1組ずつ選びます。")
    try:
        seeded, empty = select_by_density(markers, mode, chrom, int(start), int(end),
                                          spacing_bp=max(1, round(spacing_kb * 1000)), count=int(count))
    except ValueError as exc:
        st.error(str(exc))
        seeded, empty = [], []
    n_loci = len({(m.get("参照染色体"), marker_position(m)) for m in markers})
    width_kb = (int(end) - int(start) + 1) / 1000
    msg = (f"{chrom}:{int(start):,}–{int(end):,}（{width_kb:,.1f} kb）で {len(seeded)} マーカー"
           f"（候補 {len(markers)} 組・全 {n_loci} 座位）。")
    if empty:
        msg += f" 候補のない区画が {len(empty)} 個あります。"
    st.caption(msg + "区画内では共優性→区画中心に近い座位の順に選びます。下の一覧で追加・削除できます。")
    if empty:
        with st.expander(f"候補のない区画 {len(empty)} 個"):
            st.dataframe(pd.DataFrame(empty), hide_index=True, width="stretch")
    # Re-seed the editable list whenever the region or density changes.
    signature = (chrom, int(start), int(end), mode, spacing_kb, int(count), len(markers))
    if st.session_state.get("gt_density_sig") != signature:
        st.session_state["gt_density_sig"] = signature
        st.session_state["gt_markers"] = seeded
    return st.multiselect(
        "記録するマーカー", list(range(len(markers))), key="gt_markers",
        format_func=lambda i: (f"{names[i]}（{markers[i].get('A産物(bp)', '')} bp・"
                               f"{str(markers[i].get('マーカー型', '')).split('（')[0]}）"))


def genotyping_panel(markers: list[dict], goal: str, profile):
    """Population, marker panel, scoring sheet and QC. Returns the read sheet or None."""
    step = ui.StepCounter()

    with st.container(key="panel_gt_population"):
        step("集団と個体", "集団の種類で、親・F1コントロールと分離比の検定内容が決まります。")
        options = list(POPULATIONS)
        default = POPULATION_FOR_GOAL.get(goal, options[0])
        c1, c2 = st.columns([1, 2], gap="large")
        with c1:
            # Keyed per design goal so that changing the goal resets this default.
            population = st.selectbox("集団の種類", options, key="gt_population_" + str(options.index(default)),
                                      index=options.index(default))
            include_f1 = st.checkbox("F1（ヘテロ）のコントロール行を入れる",
                                     value=POPULATIONS[population] != "RIL", key="gt_include_f1")
        with c2:
            samples = _sample_names()

    with st.container(key="panel_gt_markers"):
        step("マーカーと記録シート",
             "使える候補から範囲と密度を決めて選びます。粗い間隔で全体を見て、原因区間が絞れたらその区間だけを"
             "高密度にします。記録シートには遺伝子型ごとの予測バンドが入り、ゲル判読の早見表になります。")
        if not markers:
            st.info("使えるマーカー候補がまだありません。「マーカー設計」で解析すると、ここで選べます。"
                    "記入済みシートの読み込み（下）は候補なしでも使えます。")
        else:
            names = unique_marker_names(markers)
            chosen = _density_picker(markers, names)
            chosen = sorted(chosen, key=lambda i: (str(markers[i].get("参照染色体", "")),
                                                   marker_position(markers[i]) or 0, i))
            picked = [markers[i] for i in chosen]
            if picked:
                sheet = bench_sheet(picked, profile, goal)
                st.dataframe(pd.DataFrame(sheet).drop(columns=["Forward(5'-3')", "Reverse(5'-3')"]),
                             hide_index=True, width="stretch")
                unfit = [s["マーカー名"] for s in sheet if not s["目的に適合"]]
                if unfit:
                    st.warning("選んだ集団の判別に足りないマーカーがあります: " + "、".join(unfit))
                d1, d2 = st.columns(2)
                with d1:
                    st.download_button("ジェノタイピング実験表CSV", table_bytes(sheet),
                                       "caps_maker_genotyping_plan.csv", "text/csv", width="stretch")
                with d2:
                    st.download_button(
                        "記録シートCSV（判読キー付き）",
                        scoring_sheet(picked, profile, samples, include_f1=include_f1),
                        "caps_maker_scoring_sheet.csv", "text/csv", width="stretch",
                        disabled=not samples)
                st.caption("記録はAA/AB/BB（A/H/B、0/1/2も可）。優性型マーカーでは A-（BB以外）・B-（AA以外）、"
                           "欠測は空欄か「-」。Excelで編集したらCSVで保存してください。")

    with st.container(key="panel_gt_qc"):
        step("記録した遺伝子型の確認", "記入済みの記録シートを読み込み、コントロール・判別できない記録・欠測・分離比を確認します。")
        upload = st.file_uploader("記入済みの記録シート（CSV / TSV）", type=["csv", "tsv", "txt"], key="gt_upload")
        if not upload:
            return None
        try:
            scored = read_scoring_sheet(upload.getvalue())
        except (UnicodeDecodeError, ValueError) as exc:
            st.error(str(exc))
            return None
        table, problems = qc_report(scored, population)
        progeny = sum(1 for k in scored.individuals if k[0] == "Progeny")
        distorted = sum(1 for r in table if r["分離のずれ"] == "p<0.05")
        ui.stats([("個体（後代）", progeny, "key"), ("マーカー", len(scored.markers), ""),
                  ("要確認の記録", len(problems), "" if problems else "ok"),
                  ("分離比のずれ p<0.05", distorted, "")])
        st.dataframe(pd.DataFrame(table), hide_index=True, width="stretch")
        st.caption("分離比はχ²適合度検定（記録された型のクラスで期待比を選択）。多重比較の補正はしていません。"
                   "ずれは選抜・誤判読・連鎖の手がかりで、それ自体は誤りの証拠ではありません。")
        if problems:
            with st.expander(f"要確認の記録 {len(problems)} 件", expanded=True):
                st.dataframe(pd.DataFrame(problems), hide_index=True, width="stretch")

    with st.container(key="panel_gt_export"):
        step("遺伝子型表の出力")
        e1, e2 = st.columns(2)
        with e1:
            st.download_button("遺伝子型表CSV（AA/AB/BB・ワイド形式）", genotype_table_bytes(scored),
                               "caps_maker_genotypes.csv", "text/csv", width="stretch")
            st.caption("Role・Line・Phenotype列＋マーカー列。Chromosome/Position行、親・F1・後代行。Phenotypeは空欄です。")
        with e2:
            st.download_button("R/qtl形式CSV（A/H/B/C/D）", rqtl_csv_bytes(scored),
                               "caps_maker_rqtl.csv", "text/csv", width="stretch")
            st.caption("後代のみ。位置は物理位置（Mb）で、連鎖地図のcMではありません。")
        legend = "・".join(f"{CALL_LABEL[c]}={c}" for c in "AHBCD")
        st.caption(f"コード対応：{legend}、欠測=-")
    return scored


def narrowing_panel(scored) -> None:
    with st.container(key="panel_narrow"):
        ui.panel_head("組換えと絞り込み",
                      "隣り合うマーカー間で型が変わった後代を組換え個体として抽出します。表現型と照らして原因区間を狭め、"
                      "次の回は「ジェノタイピング」で狭めた区間を高密度にし、組換え個体だけを記録します。")
        if scored is None:
            st.info("「ジェノタイピング」で記入済みの記録シートを読み込むと、ここに組換え個体と区間ごとの組換え数を表示します。")
            return
        breaks, intervals = recombination_report(scored)
        recombinants = sorted({b["個体"] for b in breaks})
        ui.stats([("組換え個体", len(recombinants), "key"), ("組換え", len(breaks), ""),
                  ("組換えのないマーカー間", sum(1 for r in intervals if not r["組換え個体"]), "")])
        if intervals:
            ui.subheading("マーカー間ごとの組換え")
            st.dataframe(pd.DataFrame(intervals), hide_index=True, width="stretch")
            st.caption("「組換え個体」が0のマーカー間は、この集団ではそれ以上分割できません。欠測をはさむ組換えは"
                       "複数のマーカー間に数え、「うち区間が確定」には含めません。優性型の記録に隠れた組換えは数えていません。")
        if breaks:
            left, right = st.columns([2, 1], gap="large")
            with left:
                ui.subheading(f"組換え個体 {len(recombinants)} 個体")
                st.dataframe(pd.DataFrame(breaks), hide_index=True, width="stretch")
            with right:
                ui.subheading("次の回の個体一覧")
                st.caption("「ジェノタイピング」の「一覧を貼り付け」に使えます。")
                st.code("\n".join(recombinants), language=None)
