"""Named local resources, enzyme sets and reusable design settings."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from .catalog import catalog_path, load_catalog, local_path, upsert_entry, validate_resource
from .core import ENZYME_METADATA, available_enzymes
from .design_settings import BUILTINS, apply_settings, snapshot_settings


def template_controls() -> None:
    saved = load_catalog()["design_templates"]
    options = {"builtin:" + k: v for k, v in BUILTINS.items()}
    options.update({"saved:" + entry["id"]: entry for entry in saved})
    # Applied in a callback, so the controls may sit anywhere on the page.
    def apply() -> None:
        try:
            apply_settings(options[st.session_state["template_choice"]]["settings"], st.session_state)
            st.session_state["_template_msg"] = ("ok", "適用しました。場所・系統・DBを選んで解析してください。")
        except (TypeError, ValueError) as exc:
            st.session_state["_template_msg"] = ("error", str(exc))

    st.selectbox("テンプレートを選ぶ", list(options), key="template_choice",
                 format_func=lambda v: options[v]["name"])
    st.button("このテンプレートを適用", width="stretch", type="primary", on_click=apply)
    tone, message = st.session_state.pop("_template_msg", (None, None))
    if tone == "ok":
        st.success(message)
    elif tone:
        st.error(message)
    st.divider()
    with st.form("save_design_template", border=False):
        name = st.text_input("設定を保存する名前", placeholder="例：F2用・50 kb間隔")
        submit = st.form_submit_button("現在の設定をテンプレート保存", width="stretch")
    st.caption("条件だけを保存します。DB・座標・試料名・プライマー配列は含めません。")
    if submit:
        try:
            old = next((e for e in saved if e["name"] == name.strip()), None)
            entry = {"name": name, "settings": snapshot_settings(st.session_state)}
            if old:
                entry["id"] = old["id"]
            upsert_entry("design_templates", entry)
            st.success("保存しました。次回もこの端末で選べます。")
        except (OSError, ValueError) as exc:
            st.error(str(exc))


def resource_picker(kind: str, label: str, key: str) -> str:
    entries = [e for e in load_catalog()["databases"] if e["kind"] == kind]
    by_id = {e["id"]: e for e in entries}
    chosen = st.selectbox(label, ["manual", *by_id], key=key, index=1 if entries else 0,
                          format_func=lambda v: "パスを直接指定" if v == "manual" else by_id[v]["name"])
    entry = by_id.get(chosen)
    if entry:
        st.caption(f"{entry['path']}" + (f" / {entry['assembly']}" if entry.get("assembly") else ""))
        return entry["path"]
    value = st.text_input(kind + "のローカルパス", key=key + "_manual",
                          help="DB管理タブで登録すると、名前で選べます。")
    return local_path(value) if value.strip() else ""


def database_manager() -> None:
    st.subheader("名前で使うローカルデータ")
    st.caption("BLAST DB・参照FASTA・VCFを登録します。元のファイルは移動しません。")
    entries = load_catalog()["databases"]
    by_id = {e["id"]: e for e in entries}
    left, right = st.columns([1, 1.2], gap="large")
    with left:
        editing = st.selectbox("登録・更新するデータ", ["new", *by_id],
                               format_func=lambda v: "新しく登録" if v == "new" else by_id[v]["name"])
        old = by_id.get(editing, {})
        with st.form("resource_" + editing):
            name = st.text_input("登録名", value=old.get("name", ""), placeholder="例：参照ゲノム v1 特異性検索")
            kinds = ["BLAST DB", "参照FASTA", "VCF"]
            kind = st.selectbox("データ種別", kinds, index=kinds.index(old.get("kind", "BLAST DB")))
            path = st.text_input("ファイルパス / DB接頭辞", value=old.get("path", ""),
                                 placeholder="C:\\data\\blastdb\\reference または /data/reference")
            assembly = st.text_input("アセンブリ・版（任意）", value=old.get("assembly", ""))
            st.caption("BLASTはフォルダーではなくDB接頭辞。FASTAは.fai、VCFは.tbi/.csiが必要です。")
            submit = st.form_submit_button("確認して保存", type="primary")
        if submit:
            try:
                if not name.strip():
                    raise ValueError("名前を入力してください。")
                checked, detail = validate_resource(kind, path)
                entry = {"name": name, "kind": kind, "path": checked,
                         "assembly": assembly.strip(), "check": detail}
                if old:
                    entry["id"] = old["id"]
                upsert_entry("databases", entry)
                st.success("登録しました。マーカー設計タブで名前を選べます。")
            except (OSError, RuntimeError, ValueError) as exc:
                st.error(str(exc))
    with right:
        st.markdown("#### 登録済みデータ")
        current = load_catalog()["databases"]
        if current:
            st.dataframe(pd.DataFrame([
                {"名前": e["name"], "種類": e["kind"], "アセンブリ・版": e.get("assembly", ""),
                 "パス": e["path"], "確認内容": e.get("check", "未確認")} for e in current
            ]), hide_index=True, width="stretch")
        else:
            st.info("まだ登録がありません。左のフォームから名前とパスを登録してください。")
        st.caption("登録時の確認は読み込みと索引までです。参照とDBの同一性・特異性は解析時に確認します。")
    with st.expander("ローカル保存先"):
        st.code(str(catalog_path()), language=None)
        st.caption("登録名・パスと酵素セットをこの端末に保存します。")


def enzyme_manager() -> list[str]:
    st.subheader("使える酵素から設計する")
    st.caption("研究室の酵素セットを保存し、解析に使う酵素を選びます。価格・反応温度は自動判定しません。")
    pool = available_enzymes()
    saved = load_catalog()["enzyme_sets"]
    by_id = {e["id"]: e for e in saved}
    with st.container(key="seg_enzyme_mode"):
        mode = st.radio("使用する酵素", ["指定した酵素", "保存した酵素セット", "全酵素"],
                        horizontal=True, key="enzyme_mode")
    left, right = st.columns([1.1, 1], gap="large")
    with left:
        preset = None
        if mode == "保存した酵素セット":
            if saved:
                preset_id = st.selectbox("酵素セット", list(by_id), format_func=lambda v: by_id[v]["name"])
                preset = by_id[preset_id]
            else:
                st.info("保存済みセットがありません。酵素を選んで下のフォームで保存できます。")
        if mode == "全酵素":
            excluded = st.multiselect("除外する酵素", pool, key="enzyme_excluded")
            enzymes = [name for name in pool if name not in excluded]
        else:
            initial = preset["enzymes"] if preset else ["MseI", "HinfI", "XbaI", "EcoRI", "BsmFI"]
            enzymes = st.multiselect("候補酵素（複数選択）", pool,
                                      default=[e for e in initial if e in pool],
                                      key="enzyme_selection_" + (preset["id"] if preset else "manual"))
        st.metric("解析に使う酵素", f"{len(enzymes)} 種")
        with st.form("save_enzyme_set"):
            set_name = st.text_input("酵素セット名", placeholder="例：研究室の保有酵素")
            save = st.form_submit_button("この酵素セットを保存", disabled=not enzymes)
        if save:
            try:
                old = next((e for e in saved if e["name"] == set_name.strip()), None)
                entry = {"name": set_name, "enzymes": sorted(enzymes)}
                if old:
                    entry["id"] = old["id"]
                upsert_entry("enzyme_sets", entry)
                st.success("酵素セットを保存しました。")
            except (OSError, ValueError) as exc:
                st.error(str(exc))
    with right:
        st.markdown("#### 選択中の酵素")
        st.dataframe(pd.DataFrame([
            {"酵素": e, "認識配列": ENZYME_METADATA[e].recognition,
             "切断位置 上鎖": ENZYME_METADATA[e].top_cut,
             "切断位置 下鎖": ENZYME_METADATA[e].bottom_cut} for e in enzymes
        ]), hide_index=True, width="stretch", height=300)
        st.caption("切断位置は認識配列の開始を0とする鎖ごとの位置です。")
    st.session_state["active_enzymes"] = list(enzymes)
    return enzymes


def primary_databases() -> list[str]:
    primary = resource_picker("BLAST DB", "設計参照に対応するBLAST DB（最初に検索）", "primary_db")
    extras = [e for e in load_catalog()["databases"] if e["kind"] == "BLAST DB" and e["path"] != primary]
    by_id = {e["id"]: e for e in extras}
    with st.expander("追加のゲノムも検索する（任意）"):
        selected = st.multiselect("追加で確認する登録DB", list(by_id),
                                  format_func=lambda v: by_id[v]["name"], key="extra_dbs")
        manual = st.text_area("追加DBのパス（任意、1行に1件）", key="extra_db_paths")
    values = [primary, *(by_id[v]["path"] for v in selected), *manual.splitlines()]
    return list(dict.fromkeys(local_path(v) for v in values if v.strip()))
