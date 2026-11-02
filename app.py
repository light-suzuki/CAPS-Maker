"""Streamlit GUI for two-FASTA CAPS marker design."""
from __future__ import annotations

import hashlib

import pandas as pd
import streamlit as st

from primerblast_oss.gel import analyze_digest_patterns, virtual_gel_svg

from caps_maker.core import (
    align_alleles, available_enzymes, parse_sequence, scan_caps, table_bytes,
)
from caps_maker.primers import design_caps

st.set_page_config(page_title="CAPS-Maker", page_icon="🧬", layout="wide")
st.title("CAPS-Maker")
st.caption("2配列からCAPS候補を探索。制限酵素の指定、Primer-BLAST、CSV/TSV出力に対応します。")

with st.sidebar:
    st.header("解析設定")
    mode = st.radio(
        "解析モード",
        ["制限酵素スキャン", "Primer-BLAST付きマーカー設計"],
        help="スキャンは入力されたDNA全長を仮想消化します。"
             "Primer-BLASTモードでは実際のPCR産物を消化します。",
    )
    all_enzymes = available_enzymes()
    enzyme_mode = st.radio(
        "制限酵素の対象", ["利用可能な全酵素", "指定した酵素だけ"], horizontal=False
    )
    if enzyme_mode == "指定した酵素だけ":
        defaults = [name for name in ("MseI", "HinfI", "XbaI", "EcoRI", "BsmFI")
                    if name in all_enzymes]
        enzymes = st.multiselect(
            "使用する制限酵素（複数選択可）", all_enzymes, default=defaults
        )
    else:
        enzymes = all_enzymes
        st.caption(f"PCRに適した切断情報検証済みの{len(enzymes)}酵素を探索")
    gel_min_gap = st.number_input("最低バンド差の目安 (bp)", 1, 500, 25, 5)
    if mode == "Primer-BLAST付きマーカー設計":
        st.divider()
        st.subheader("PCR / BLAST")
        blast_db_text = st.text_area(
            "BLAST DBのローカルパス（1行につき1つ）",
            placeholder="/home/user/blastdb/pisum_v1a",
            help="makeblastdb済みのローカルDBプレフィックスを指定してください。"
                 "ゲノムFASTAファイルそのものとは異なります。",
        )
        product_min = st.number_input("PCR最小産物長 (bp)", 70, 2000, 200, 10)
        product_max = st.number_input("PCR最大産物長 (bp)", 70, 4000, 800, 10)
        num_pairs = st.number_input("設計する候補プライマーペア", 1, 100, 15)

st.subheader("配列入力")
input_type = st.radio("入力方法", ["配列を貼り付け", "FASTAファイルをアップロード"],
                      horizontal=True)
left, right = st.columns(2)
if input_type == "配列を貼り付け":
    with left:
        raw_a = st.text_area("配列A（例: 野生型）", height=180,
                             placeholder=">wildtype\nATGC...")
    with right:
        raw_b = st.text_area("配列B（例: 変異型）", height=180,
                             placeholder=">mutant\nATGA...")
else:
    with left:
        file_a = st.file_uploader("配列AのFASTA", type=["fa", "fasta", "fna", "txt"],
                                  key="fasta_a")
    with right:
        file_b = st.file_uploader("配列BのFASTA", type=["fa", "fasta", "fna", "txt"],
                                  key="fasta_b")
    raw_a = file_a.getvalue().decode("utf-8-sig") if file_a else ""
    raw_b = file_b.getvalue().decode("utf-8-sig") if file_b else ""

a_record = b_record = alignment = selected_variant = None
if raw_a.strip() and raw_b.strip():
    try:
        a_record = parse_sequence(raw_a)
        b_record = parse_sequence(raw_b)
        alignment = align_alleles(a_record.sequence, b_record.sequence)
        st.caption(
            f"A: {len(a_record.sequence):,} bp / B: {len(b_record.sequence):,} bp / "
            f"検出した多型: {len(alignment.variants):,} 件"
        )
        if any(letter not in "ACGT" for letter in a_record.sequence + b_record.sequence):
            st.warning("NやIUPAC曖昧塩基を含みます。制限酵素の検出には不確実性があります。")
        if mode == "Primer-BLAST付きマーカー設計":
            if alignment.variants:
                selected_index = st.selectbox(
                    "標的多型を選択",
                    options=list(range(len(alignment.variants))),
                    format_func=lambda i: alignment.variants[i].label,
                )
                selected_variant = alignment.variants[selected_index]
            else:
                st.info("配列A/Bは同一です。標的多型がありません。")
    except (UnicodeDecodeError, ValueError) as exc:
        st.error(str(exc))
        a_record = b_record = alignment = None

signature_input = (
    raw_a + "\0" + raw_b + "\0" + mode + "\0" + str(enzymes)
    + "\0" + str(gel_min_gap)
    + ("\0" + str(selected_variant) + "\0" + blast_db_text + "\0"
       + str(product_min) + "\0" + str(product_max) + "\0" + str(num_pairs)
       if mode == "Primer-BLAST付きマーカー設計" else "")
)
signature = hashlib.sha256(signature_input.encode("utf-8")).hexdigest()

run = st.button("CAPSマーカーを探索", type="primary",
                disabled=not (a_record and b_record and enzymes))
if run:
    try:
        if mode == "制限酵素スキャン":
            rows = scan_caps(a_record.sequence, b_record.sequence,
                             enzymes, int(gel_min_gap))
            summary = {"候補数": len(rows)}
            st.session_state.caps_result_mode = mode
        else:
            if selected_variant is None:
                raise ValueError("標的多型を選択してください。")
            if int(product_max) < int(product_min):
                raise ValueError("PCR最大産物長は最小産物長以上にしてください。")
            databases = [line.strip() for line in blast_db_text.splitlines()
                         if line.strip()]
            with st.spinner("Primer3で設計し、BLASTと両アレルの切断を確認中"):
                rows, summary = design_caps(
                    a_record.sequence, b_record.sequence, alignment,
                    selected_variant, enzymes, databases,
                    product_min=int(product_min), product_max=int(product_max),
                    num_pairs=int(num_pairs), min_gap=int(gel_min_gap)
                )
        st.session_state.caps_result_mode = mode
        st.session_state.caps_rows = rows
        st.session_state.caps_summary = summary
        st.session_state.caps_signature = signature
    except (ValueError, RuntimeError, OSError) as exc:
        st.error(str(exc))
        st.session_state.pop("caps_rows", None)

if "caps_rows" in st.session_state:
    st.divider()
    st.subheader("解析結果")
    if st.session_state.get("caps_signature") != signature:
        st.warning("入力や設定が変わっています。以下は前回実行時の結果です。")
    if st.session_state.get("caps_result_mode") == "制限酵素スキャン":
        st.info("これは入力配列全長の仮想消化です。Primer3で決めたPCR産物の"
                "切断パターンではありません。候補を実験に用いる前に産物単位で再評価してください。")
    summary = st.session_state.get("caps_summary", {})
    if "primer3_pairs" in summary:
        st.caption(
            f"Primer3候補 {summary['primer3_pairs']} 組 / "
            f"両アレルで結合確認 {summary['both_alleles_bind']} 組 / "
            f"CAPS候補 {summary['caps_candidates']} 件"
        )
        st.caption("BLASTの判定は指定DBに対するもので、未完全な探索は特異性を保証しません。")
        if not summary["both_alleles_bind"]:
            st.warning("両アレルでプライマー結合を確認できたペアがありません。")
    original_rows = st.session_state.caps_rows
    if not original_rows:
        st.warning("該当する酵素・プライマーの組み合わせはありませんでした。")
    else:
        filter_left, filter_right = st.columns(2)
        with filter_left:
            distinguishable_only = st.checkbox("ゲル判別可能な候補のみ", value=True)
        with filter_right:
            site_count_only = st.checkbox("切断サイト数が異なる候補のみ", value=False)
        visible_rows = [
            row for row in original_rows
            if (not distinguishable_only or row["ゲル判別"] == "可能")
            and (not site_count_only or row["認識サイト数の差"])
        ]
        st.caption(f"表示 {len(visible_rows)} / 全 {len(original_rows)} 件")
        st.dataframe(pd.DataFrame(visible_rows, columns=list(original_rows[0])),
                     use_container_width=True, hide_index=True)
        if not visible_rows:
            st.info("現在のフィルターに一致する候補はありません。")
        else:
            dl1, dl2 = st.columns(2)
            with dl1:
                st.download_button(
                    "表をCSVでダウンロード（Excel対応）",
                    data=table_bytes(visible_rows),
                    file_name="caps_maker_results.csv", mime="text/csv",
                )
            with dl2:
                st.download_button(
                    "表をTSVでダウンロード",
                    data=table_bytes(visible_rows, "\t"),
                    file_name="caps_maker_results.tsv",
                    mime="text/tab-separated-values",
                )

            st.markdown("#### AA / AB / BB の予測ゲル")
            selected_row_index = st.selectbox(
                "ゲル画像に表示する候補",
                list(range(len(visible_rows))),
                format_func=lambda index: (
                    f"{index + 1}. {visible_rows[index]['酵素']} "
                    f"({visible_rows[index]['マーカー型']}; "
                    f"評価 {visible_rows[index]['ゲル評価']})"
                ),
            )
            selected_row = visible_rows[selected_row_index]

            def expand_bands(bands):
                lengths = []
                for token in bands.split(" + "):
                    if "×" in token:
                        size, count = token.split("×")
                        lengths.extend([int(size)] * int(count))
                    else:
                        lengths.append(int(token))
                return lengths

            gel = analyze_digest_patterns(
                expand_bands(selected_row["A断片(bp)"]),
                expand_bands(selected_row["B断片(bp)"]),
            )
            svg = virtual_gel_svg(gel)
            st.image(svg, width=520,
                     caption="M（サイズマーカー）、AA、AB、BB。予測された標的断片のみ。")
            st.download_button("予測ゲルをSVGで保存", data=svg,
                               file_name="caps_maker_gel.svg",
                               mime="image/svg+xml")
            st.caption(
                "このゲルは完全消化・両アレルの増幅を仮定します。"
                "オフターゲットPCR産物の切断断片は重ねていません。"
            )
        with st.expander("解析の注意点"):
            st.markdown(
                "結果はin silico予測で、実際のPCR・制限酵素消化の成功を保証しません。"
                "特に小断片の検出、部分消化、ゲル条件、メチル化、"
                "BLAST DBに含まれないアレルは実験前に確認してください。"
                "ゲル判別可能という評価はゲルモデルの予測です。"
            )
