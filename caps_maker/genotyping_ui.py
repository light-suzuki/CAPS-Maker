"""Genotyping tab: bench sheet, scoring sheet, call QC and genotype-table export."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from . import ui
from .core import table_bytes
from .genotyping import (
    CALL_LABEL, POPULATION_FOR_GOAL, POPULATIONS, bench_sheet, genotype_table_bytes,
    marker_name, qc_report, read_scoring_sheet, rqtl_csv_bytes, scoring_sheet, unique_marker_names,
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


def genotyping_panel(markers: list[dict], goal: str, profile) -> None:
    step = ui.StepCounter()

    step("集団と個体", "交配集団の種類で、親・F1コントロールと分離比の検定内容が決まります。")
    options = list(POPULATIONS)
    # Keyed per design goal so that changing the goal resets this default.
    population = st.selectbox("集団の種類", options, key="gt_population_" + str(options.index(
                                  POPULATION_FOR_GOAL.get(goal, options[0]))),
                              index=options.index(POPULATION_FOR_GOAL.get(goal, options[0])))
    samples = _sample_names()
    include_f1 = st.checkbox("F1（ヘテロ）のコントロール行を入れる", value=POPULATIONS[population] != "RIL",
                             key="gt_include_f1")

    step("マーカーと記録シート",
         "マーカー設計で表示された使える候補から選びます。記録シートには遺伝子型ごとの予測バンドが入り、"
         "そのままゲル判読の早見表になります。")
    if not markers:
        st.info("使えるマーカー候補がまだありません。「マーカー設計」で解析すると、ここで選べます。"
                "記録済みのシートの読み込み（下）は候補なしでも使えます。")
    else:
        names = unique_marker_names(markers)
        # Candidates are ranked; alternative primer pairs for one locus share a
        # base name. One assay per locus is the sensible default for a panel.
        first_per_locus = list({marker_name(m, i): i for i, m in reversed(list(enumerate(markers)))}.values())
        chosen = st.multiselect(
            "記録するマーカー", list(range(len(markers))), key="gt_markers",
            default=sorted(first_per_locus),
            format_func=lambda i: (f"{names[i]}（{markers[i].get('A産物(bp)', '')} bp・"
                                   f"{str(markers[i].get('マーカー型', '')).split('（')[0]}）"))
        st.caption(f"{len(markers)} 候補・{len(first_per_locus)} 座位。初期選択は各座位で最上位の1組です。")
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

    step("記録した遺伝子型の確認", "記入済みの記録シートを読み込み、コントロール・判別できない記録・欠測・分離比を確認します。")
    upload = st.file_uploader("記入済みの記録シート（CSV / TSV）", type=["csv", "tsv", "txt"], key="gt_upload")
    if not upload:
        return
    try:
        scored = read_scoring_sheet(upload.getvalue())
    except (UnicodeDecodeError, ValueError) as exc:
        st.error(str(exc))
        return
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
