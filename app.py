"""CAPS-Maker: reference-aware, gel-aware CAPS/dCAPS marker designer."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import streamlit as st

from primerblast_oss.gel import virtual_gel_svg
from caps_maker.batch import choose_sites, discover_sites
from caps_maker.core import (
    align_alleles, available_enzymes, parse_sequence, scan_caps, table_bytes,
)
from caps_maker.existing import analyze_primers, locate_reference_product
from caps_maker.experiments import (
    GelProfile, PROFILE_PRESETS, assess_gel, parse_custom_ladder,
)
from caps_maker.primers import design_caps
from caps_maker.sources import (
    read_sample_variants, sample_names, source_window,
)
from caps_maker.validation import parse_fragment_column

st.set_page_config(page_title="CAPS-Maker", page_icon="🧬", layout="wide")
st.title("CAPS-Maker")
st.caption(
    "2系統の多型からCAPS/dCAPSを設計し、ゲル実験条件と"
    "ゲノム特異性を評価するローカル解析ツール"
)

@st.cache_data(show_spinner=False)
def cached_sample_names(path: str) -> list[str]:
    return sample_names(path)


def parsing_sig(*args) -> str:
    return hashlib.sha256(json.dumps(args, sort_keys=True, default=str).encode()).hexdigest()


def show_result_table(rows: list[dict], mode: str) -> None:
    st.divider()
    st.subheader("マーカー候補一覧")
    if not rows:
        st.warning("条件を満たす候補はありません。")
        return
    counts = {
        "総候補": len(rows),
        "ゲル判別可能": sum(r.get("ゲル判別") == "可能" for r in rows),
        "参照特異性合格": sum(r.get("参照検証採用可") is True for r in rows),
    }
    st.caption(" / ".join(f"{k} {v}件" for k, v in counts.items()))
    col1, col2, col3 = st.columns(3)
    with col1:
        only_scored = st.checkbox("ゲルで判別可能", value=False)
    with col2:
        only_codominant = st.checkbox("共優性のみ", value=False)
    with col3:
        only_validated = st.checkbox("参照特異性合格のみ", value=False)
    visible = [
        row for row in rows
        if (not only_scored or row.get("ゲル判別") == "可能")
        and (not only_codominant or "共優性" in row.get("マーカー型", ""))
        and (not only_validated or row.get("参照検証採用可") is True)
    ]
    st.dataframe(pd.DataFrame(visible, columns=list(rows[0])),
                 use_container_width=True, hide_index=True)
    st.caption(f"{len(visible)} / {len(rows)}件を表示")
    if not visible:
        return
    download1, download2, download3 = st.columns(3)
    with download1:
        st.download_button("結果CSV（Excel対応）", table_bytes(visible),
                           "caps_maker_results.csv", "text/csv")
    with download2:
        st.download_button("結果TSV", table_bytes(visible, "\t"),
                           "caps_maker_results.tsv", "text/tab-separated-values")
    with download3:
        if all("Forward(5'-3')" in row for row in visible):
            ordered = set()
            order_rows = []
            for row in visible:
                key = (row["Forward(5'-3')"], row["Reverse(5'-3')"])
                if key in ordered:
                    continue
                ordered.add(key)
                index = len(ordered)
                order_rows.extend([
                    {"Primer ID": f"CAPS_{index:03d}_F",
                     "Sequence 5'-3'": key[0], "Tm(°C)": row.get("Tm F(°C)", "")},
                    {"Primer ID": f"CAPS_{index:03d}_R",
                     "Sequence 5'-3'": key[1], "Tm(°C)": row.get("Tm R(°C)", "")},
                ])
            st.download_button("プライマー発注CSV", table_bytes(order_rows),
                               "caps_maker_primer_order.csv", "text/csv")

    st.markdown("#### AA / AB / BB 予測ゲル")
    choice = st.selectbox(
        "表示するマーカー",
        list(range(len(visible))),
        format_func=lambda i: f"{i+1}. {visible[i]['酵素']}"
                              f" / {visible[i].get('設計方式', mode)}"
                              f" / {visible[i]['マーカー型']}",
    )
    row = visible[choice]
    try:
        gel = assess_gel(
            parse_fragment_column(row["A断片(bp)"]),
            parse_fragment_column(row["B断片(bp)"]),
            st.session_state.get("active_gel_profile") or GelProfile(),
        )
        svg = virtual_gel_svg(gel)
        st.image(svg, width=520)
        st.download_button("仮想ゲルSVGを保存", svg, "caps_maker_gel.svg",
                           "image/svg+xml")
        st.caption(f"ラダー {gel['ladder']} / ゲル {gel['gel_percent']}% / "
                   f"可視 {gel['visible_range_bp'][0]}～{gel['visible_range_bp'][1]} bp。"
                   "完全消化と両アレルの増幅を仮定。"
                   "副増幅のバンドは仮想画像には重ねていません。")
    except (ValueError, KeyError) as exc:
        st.warning("ゲル表示に必要な情報がありません: " + str(exc))


with st.sidebar:
    st.header("入力・設計")
    input_mode = st.radio("配列の指定方法", [
        "2つのFASTA配列", "参照区間＋2系統VCF",
    ])
    run_mode = st.radio("解析方式", [
        "制限酵素の探索", "新規プライマー設計", "既存プライマーから解析",
        "区間マーカーを一括設計",
    ])
    st.divider()
    st.subheader("制限酵素")
    enzyme_pool = available_enzymes()
    enzyme_mode = st.radio("使用する酵素", ["指定した酵素", "全酵素"],
                          horizontal=True)
    if enzyme_mode == "指定した酵素":
        initial = [name for name in ("MseI", "HinfI", "XbaI", "EcoRI", "BsmFI")
                   if name in enzyme_pool]
        enzymes = st.multiselect("候補酵素（複数選択）",
                                 enzyme_pool, default=initial)
    else:
        enzymes = enzyme_pool
    st.caption(f"対象 {len(enzymes)} 酵素")

    st.subheader("ゲル・ラダー条件")
    preset = st.selectbox("実験プリセット",
                          list(PROFILE_PRESETS) + ["詳細設定"])
    base_profile = PROFILE_PRESETS.get(preset, GelProfile())
    with st.expander("ゲルの詳細条件", expanded=preset == "詳細設定"):
        ladder_choice = st.selectbox(
            "ラダーマーカー",
            ["auto", "100bp", "1kb", "custom"],
            index=["auto", "100bp", "1kb", "custom"].index(base_profile.ladder),
        )
        ladder_text = st.text_input(
            "カスタムラダー（サイズbpをカンマ区切り）",
            value="100, 200, 300, 400, 500, 600, 700, 800, 900, 1000"
        )
        concentrations = [None, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2,
                          1.5, 2.0, 3.0, 4.0, 5.0]
        gel_percentage = st.selectbox(
            "ゲル濃度（%）",
            concentrations, index=concentrations.index(base_profile.gel_percent),
            format_func=lambda value: "自動推奨" if value is None else str(value) + "%",
        )
        min_visible = st.number_input("検出可能な最小断片(bp)",
                                      1, 5000, base_profile.min_visible_bp)
        max_visible = st.number_input("可視上限(bp、0は自動)",
                                      0, 50000, 0)
        min_log_gap = st.number_input(
            "ゲル上の分離閾値(log10)", 0.001, 0.500,
            base_profile.min_log10_gap, 0.005, format="%.3f")
        gel_min_gap = st.number_input("参考：最低バンド差(bp)", 1, 500, 25)
    try:
        gel_profile = GelProfile(
            ladder=ladder_choice,
            custom_bands=parse_custom_ladder(ladder_text) if ladder_choice == "custom" else (),
            gel_percent=gel_percentage,
            min_visible_bp=int(min_visible),
            min_log10_gap=float(min_log_gap),
            max_visible_bp=int(max_visible) if max_visible else None,
        )
        gel_profile.validate()
        st.session_state.active_gel_profile = gel_profile
    except ValueError as exc:
        st.error(str(exc))
        gel_profile = GelProfile()

    with st.expander("PCR / 特異性の詳細設定", expanded=run_mode != "制限酵素の探索"):
        blast_db_text = st.text_area(
            "ローカルBLAST DB（1行に1件）", placeholder="/data/blastdb/cameor_v2",
            help="座標指定のリファレンスと同じアセンブリのDBを最初に指定してください。"
        )
        product_min = st.number_input("増幅産物長 下限(bp)", 70, 3000, 200, 10)
        product_max = st.number_input("増幅産物長 上限(bp)", 70, 5000, 800, 10)
        num_pairs = st.number_input("Primer3の候補ペア数", 1, 100, 12)
        specificity_choice = st.radio(
            "特異性判定方式",
            ["厳密（標的以外の産物を認めない）",
             "ゲル判別型（副産物も消化して評価）"],
        )
        specificity_policy = (
            "strict" if specificity_choice.startswith("厳密") else "gel-tolerant"
        )
        dcaps_choice = st.radio(
            "dCAPS探索", ["CAPSのみ", "CAPSが使えなければdCAPS", "CAPS/dCAPS両方"]
        )
        dcaps_mode = {
            "CAPSのみ": "none",
            "CAPSが使えなければdCAPS": "fallback",
            "CAPS/dCAPS両方": "always",
        }[dcaps_choice]
        max_dcaps_mismatches = st.selectbox("dCAPS人工ミスマッチ数の上限", [1, 2])
        st.caption("ゲル判別型は全副産物を実配列で再構成できた場合のみ合格。"
                   "BLASTの探索未完了・参照座標不明なら採用しません。")
    databases = [line.strip() for line in blast_db_text.splitlines() if line.strip()]

raw_a = raw_b = ""
locus = alignment = selected_variant = None
source_args = None
vcf_data = None

if input_mode == "2つのFASTA配列":
    st.subheader("2系統の配列を入力")
    input_format = st.radio("配列入力", ["貼り付け", "FASTAアップロード"], horizontal=True)
    ca, cb = st.columns(2)
    if input_format == "貼り付け":
        with ca:
            raw_a = st.text_area("配列A（野生型など）", height=160,
                                 placeholder=">allele_A\nATG...")
        with cb:
            raw_b = st.text_area("配列B（変異型など）", height=160,
                                 placeholder=">allele_B\nATG...")
    else:
        with ca:
            fa = st.file_uploader("配列A FASTA", type=["fa", "fasta", "fna", "txt"])
        with cb:
            fb = st.file_uploader("配列B FASTA", type=["fa", "fasta", "fna", "txt"])
        raw_a = fa.getvalue().decode("utf-8-sig") if fa else ""
        raw_b = fb.getvalue().decode("utf-8-sig") if fb else ""
    if raw_a.strip() and raw_b.strip():
        try:
            record_a = parse_sequence(raw_a)
            record_b = parse_sequence(raw_b)
            locus = (record_a.sequence, record_b.sequence)
            alignment = align_alleles(*locus)
            st.caption(f"A={len(locus[0])} bp / B={len(locus[1])} bp / "
                       f"多型 {len(alignment.variants)} 件")
        except ValueError as exc:
            st.error(str(exc))
    if run_mode == "区間マーカーを一括設計":
        st.warning("大量マーカー設計には参照区間＋2系統VCFを選んでください。")
else:
    st.subheader("参照ゲノムの区間と2系統")
    rc1, rc2 = st.columns(2)
    with rc1:
        ref_kind = st.radio("参照配列ソース", ["FASTA", "BLAST DB"], horizontal=True)
        ref_locator = st.text_input("FASTAファイル / BLAST DBのローカルパス",
                                    placeholder="/data/genome.fa")
        chrom = st.text_input("染色体・コンティグID", placeholder="Chr4")
    with rc2:
        vcf_path = st.text_input("2系統を含むVCF（.vcf.gz＋.tbi/.csi）",
                                 placeholder="/data/cross.vcf.gz")
        region_start = st.number_input("参照開始(1-based)", 1, 10**10, 1, step=1)
        region_end = st.number_input("参照終了(1-based, inclusive)",
                                     1, 10**10, 8000, step=1)
    samples = []
    if vcf_path and Path(vcf_path).exists():
        try:
            samples = cached_sample_names(vcf_path)
        except (OSError, ValueError) as exc:
            st.error("VCFのヘッダーを読めません: " + str(exc))
    a_sel, b_sel = st.columns(2)
    with a_sel:
        sample_a = st.selectbox("系統A（VCFサンプル）", samples, index=0 if samples else None)
    with b_sel:
        sample_b = st.selectbox("系統B（VCFサンプル）", samples,
                                index=1 if len(samples) > 1 else 0 if samples else None)
    if (ref_locator and chrom and sample_a and sample_b
            and region_end >= region_start and sample_a != sample_b):
        source_args = (
            ref_kind, ref_locator, chrom, int(region_start), int(region_end),
            vcf_path, sample_a, sample_b,
        )
        st.caption(f"{chrom}:{region_start:,}–{region_end:,} / {sample_a} vs {sample_b}")
        if run_mode not in ("区間マーカーを一括設計", "既存プライマーから解析"):
            if int(region_end) - int(region_start) + 1 <= 8000:
                try:
                    vcf_data = read_sample_variants(vcf_path, chrom,
                                                    int(region_start), int(region_end),
                                                    sample_a, sample_b)
                    win = source_window(ref_kind, ref_locator, chrom,
                                        int(region_start), int(region_end), vcf_data)
                    locus = (win.a, win.b)
                    alignment = align_alleles(*locus)
                    st.session_state.single_locus_window = win
                    st.caption(f"区間内の異なるVCF変異 {sum(v.different for v in vcf_data)}件")
                except (OSError, ValueError, RuntimeError) as exc:
                    st.error(str(exc))
            else:
                st.caption("長い区間です。新規設計では下の標的ゲノム座標を指定するか、"
                           "区間一括設計モードを使用してください。")
    elif source_args is None:
        st.info("参照、染色体、異なる2サンプル、正しい開始/終了座標を指定してください。")

if run_mode == "既存プライマーから解析":
    st.subheader("既存PCRプライマー")
    ec1, ec2 = st.columns(2)
    with ec1:
        existing_forward = st.text_input("Forward（5'→3'）")
    with ec2:
        existing_reverse = st.text_input("Reverse（5'→3'）")

if run_mode == "新規プライマー設計":
    st.subheader("標的多型")
    if source_args and not locus:
        target_ref_pos = st.number_input(
            "標的の参照座標(1-based)", min_value=int(region_start),
            max_value=int(region_end), value=int(region_start))
    elif alignment and alignment.variants:
        index = st.selectbox(
            "標的多型を選択", list(range(len(alignment.variants))),
            format_func=lambda i: alignment.variants[i].label,
        )
        selected_variant = alignment.variants[index]
    else:
        st.caption("配列A/B間に多型があることを確認してください。")

signature = parsing_sig(
    raw_a, raw_b, source_args, run_mode, enzymes,
    gel_profile, gel_min_gap, databases, product_min, product_max,
    num_pairs, dcaps_mode, max_dcaps_mismatches, specificity_policy,
    selected_variant,
)

if run_mode in ("制限酵素の探索", "新規プライマー設計", "既存プライマーから解析"):
    available = bool(enzymes) and (bool(locus) or bool(source_args))
    if st.button("CAPS候補を解析", type="primary", disabled=not available):
        try:
            window = st.session_state.get("single_locus_window") if source_args and locus else None
            if run_mode == "既存プライマーから解析":
                if not existing_forward or not existing_reverse:
                    raise ValueError("Forward/Reverseプライマーを入力してください。")
                evidence = None
                if source_args:
                    exact, evidence = locate_reference_product(
                        existing_forward, existing_reverse, databases, chrom,
                        int(region_start), int(region_end),
                    )
                    vcf_data = read_sample_variants(
                        vcf_path, chrom, exact[0], exact[1], sample_a, sample_b)
                    window = source_window(ref_kind, ref_locator, chrom,
                                           exact[0], exact[1], vcf_data)
                    locus = window.a, window.b
                    alignment = align_alleles(*locus)
                if not locus:
                    raise ValueError("解析する2配列が必要です。")
                rows, summary = analyze_primers(
                    *locus, alignment, existing_forward, existing_reverse,
                    enzymes, databases, min_bp=int(product_min),
                    max_bp=int(product_max), min_gap=int(gel_min_gap),
                    gel_profile=gel_profile,
                    specificity_policy=specificity_policy,
                    window=window, check_evidence=evidence,
                )
            elif run_mode == "制限酵素の探索":
                if not locus:
                    raise ValueError("8,000 bp以下の2配列または参照区間が必要です。")
                rows = scan_caps(*locus, enzymes, int(gel_min_gap),
                                 profile=gel_profile)
                summary = {"全長消化候補": len(rows)}
            else:
                if source_args and not locus:
                    low = max(int(region_start), int(target_ref_pos) - int(product_max))
                    high = min(int(region_end), int(target_ref_pos) + int(product_max))
                    vcf_data = read_sample_variants(
                        vcf_path, chrom, low, high, sample_a, sample_b)
                    window = source_window(
                        ref_kind, ref_locator, chrom, low, high, vcf_data)
                    locus = window.a, window.b
                    alignment = align_alleles(*locus)
                    loc_index = window.ref_to_a.get(int(target_ref_pos) - window.start)
                    if loc_index is None:
                        raise ValueError("標的座標の塩基は系統Aで欠失しています。")
                    candidates = [v for v in alignment.variants
                                  if v.position_a == loc_index]
                    if not candidates:
                        raise ValueError("指定座標に2系統間の差がありません。")
                    selected_variant = candidates[0]
                if not locus or selected_variant is None:
                    raise ValueError("標的多型を指定してください。")
                with st.spinner("Primer3とゲノムBLAST・制限酵素消化を実行中"):
                    rows, summary = design_caps(
                        *locus, alignment, selected_variant, enzymes, databases,
                        product_min=int(product_min), product_max=int(product_max),
                        num_pairs=int(num_pairs), min_gap=int(gel_min_gap),
                        gel_profile=gel_profile, dcaps_mode=dcaps_mode,
                        dcaps_max_mismatches=int(max_dcaps_mismatches),
                        specificity_policy=specificity_policy, window=window,
                    )
            st.session_state["caps_results"] = rows
            st.session_state["caps_summary"] = summary
            st.session_state["caps_signature"] = signature
            st.session_state["caps_result_mode"] = run_mode
        except (OSError, RuntimeError, ValueError, KeyError) as exc:
            st.error(str(exc))
            st.session_state.pop("caps_results", None)

if run_mode == "区間マーカーを一括設計" and source_args:
    st.subheader("区間内のマーカー密度")
    plan1, plan2, plan3 = st.columns(3)
    with plan1:
        panel_mode = st.selectbox("設計密度",
                                  ["全酵素サイト", "一定間隔", "指定本数"])
    with plan2:
        spacing = st.number_input("目標間隔(bp)", 1, 10_000_000, 10000, 1000)
        target_count = st.number_input("目標本数", 1, 1000, 20)
    with plan3:
        max_discovery = st.number_input("探索変異上限", 1, 50000, 1000, 100)
        max_targets = st.number_input("候補サイト数の上限", 1, 1000, 100, 10)
    batch_sig = parsing_sig(source_args, enzymes, panel_mode, spacing,
                            target_count, max_discovery, max_targets, dcaps_mode)
    if st.button("ステップ1：CAPS/dCAPSサイトを先に探索", type="primary"):
        try:
            with st.spinner("区間内の制限酵素サイトを探索中"):
                variants = read_sample_variants(
                    vcf_path, chrom, int(region_start), int(region_end),
                    sample_a, sample_b)
                discovered, stats = discover_sites(
                    ref_kind, ref_locator, chrom,
                    int(region_start), int(region_end), variants, enzymes,
                    include_dcaps=dcaps_mode != "none",
                    max_snps=int(max_discovery),
                )
                planned, panel_stats = choose_sites(
                    discovered, panel_mode, int(region_start), int(region_end),
                    spacing_bp=int(spacing), target_count=int(target_count),
                    limit=int(max_targets),
                )
            st.session_state["batch_variants"] = variants
            st.session_state["batch_sites"] = planned
            st.session_state["batch_sig"] = batch_sig
            st.session_state["batch_stats"] = {**stats, **panel_stats}
        except (OSError, RuntimeError, ValueError, KeyError) as exc:
            st.error(str(exc))
            st.session_state.pop("batch_sites", None)

    if "batch_sites" in st.session_state:
        stats = st.session_state["batch_stats"]
        st.caption(
            f"差のある変異 {stats['n_differential_variants']} 件 / "
            f"認識サイト {stats['n_sites']} 件 / "
            f"設計対象 {stats['selected']} 件"
        )
        if stats.get("truncated"):
            st.warning("設計対象が上限を超えています。上限件数までの結果です。")
        if batch_sig != st.session_state.get("batch_sig"):
            st.warning("入力・解析設定が変更されました。ステップ1を再実行してください。")
        else:
            planned = st.session_state["batch_sites"]
            data = [
                {"標的": x.identifier, "染色体": x.chrom,
                 "変異位置": x.position, "制限酵素サイト": x.site_position,
                 "酵素": x.enzyme, "種類": x.kind}
                for x in planned
            ]
            st.dataframe(pd.DataFrame(data), hide_index=True, use_container_width=True)
            if planned:
                indices = st.multiselect(
                    "プライマーを設計するサイト",
                    list(range(len(planned))),
                    default=list(range(min(10, len(planned)))),
                    format_func=lambda i: (
                        f"{planned[i].chrom}:{planned[i].site_position} "
                        f"{planned[i].enzyme} {planned[i].kind}"
                    ),
                )
                if st.button("ステップ2：選択サイトのPrimer-BLASTを実行",
                             disabled=not indices):
                    all_rows = []
                    failed = []
                    progress = st.progress(0)
                    for ordinal, site_i in enumerate(indices, 1):
                        site = planned[site_i]
                        try:
                            low = max(int(region_start), site.position - int(product_max))
                            high = min(int(region_end), site.position + int(product_max))
                            window = source_window(
                                ref_kind, ref_locator, chrom, low, high,
                                st.session_state["batch_variants"])
                            alignment = align_alleles(window.a, window.b)
                            in_a = window.ref_to_a.get(site.position - low)
                            possible = [v for v in alignment.variants
                                        if in_a is not None and v.position_a == in_a]
                            if not possible:
                                raise ValueError("標的変異を局所配列に配置できません。")
                            variant = possible[0]
                            rows, _ = design_caps(
                                window.a, window.b, alignment, variant,
                                [site.enzyme], databases,
                                product_min=int(product_min), product_max=int(product_max),
                                num_pairs=int(num_pairs), min_gap=int(gel_min_gap),
                                gel_profile=gel_profile,
                                dcaps_mode=("always" if site.kind == "dCAPS" else "none"),
                                dcaps_max_mismatches=int(max_dcaps_mismatches),
                                specificity_policy=specificity_policy, window=window,
                            )
                            rows = [r for r in rows if r.get("設計方式") == site.kind]
                            if not rows:
                                raise ValueError("条件を満たすプライマー・酵素ペアがありません。")
                            for row in rows[:3]:
                                row["染色体"] = chrom
                                row["標的座標(1-based)"] = site.position
                                row["制限酵素サイト(1-based)"] = site.site_position
                                all_rows.append(row)
                        except (OSError, ValueError, RuntimeError) as exc:
                            failed.append({"標的": site.identifier, "理由": str(exc)})
                        progress.progress(ordinal / len(indices))
                    st.session_state["caps_results"] = sorted(
                        all_rows,
                        key=lambda r: (not r.get("参照検証採用可"),
                                       -float(r.get("総合設計スコア", 0))),
                    )
                    st.session_state["caps_summary"] = {
                        "attempted": len(indices), "failed": failed,
                        "result_rows": len(all_rows),
                    }
                    st.session_state["caps_signature"] = batch_sig
                    st.session_state["caps_result_mode"] = run_mode

if "caps_results" in st.session_state:
    if run_mode == st.session_state.get("caps_result_mode"):
        summary = st.session_state.get("caps_summary") or {}
        with st.expander("実行サマリー・未成功サイト", expanded=False):
            st.json(summary)
            if summary.get("failed"):
                st.download_button(
                    "未成功サイトCSV", table_bytes(summary["failed"]),
                    "caps_maker_unresolved.csv", "text/csv",
                )
        show_result_table(st.session_state["caps_results"], run_mode)
    else:
        st.caption("別の解析モードで実行された結果があります。")

st.info(
    "スコアは実験成功確率ではなく設計上の優先順位です。"
    "CAPS探索単独は入力DNAの全長消化でありPCR産物の結果ではありません。"
    "dCAPSはプライマー改変によって実際のPCR産物配列が変わるため、"
    "提案されたプライマーと酵素で必ず実験検証してください。"
)
