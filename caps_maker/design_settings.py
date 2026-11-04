"""Reusable design recipes. Coordinates, sample IDs and DB paths are not recipes."""
from __future__ import annotations

from copy import deepcopy

from .genotyping import GOALS


DEFAULTS = {
    "input_mode": "参照区間＋2系統VCF", "run_mode": "新規プライマー設計",
    "product_min": 200, "product_max": 800, "num_pairs": 12,
    "specificity_choice": "厳密（標的以外の産物を認めない）",
    "dcaps_choice": "CAPSが使えなければdCAPS", "max_dcaps_mismatches": 1,
    "goal": "AA・AB・BBを区別（ヘテロも判別）", "preset": "標準（2%アガロース）",
    "panel_mode": "一定間隔", "spacing": 10000, "target_count": 20,
    "max_discovery": 1000, "max_targets": 100, "min_spacing": 0, "radius": 0,
    "enzyme_mode": "指定した酵素", "enzyme_selection_manual": ["MseI", "HinfI", "XbaI", "EcoRI", "BsmFI"],
    "enzyme_excluded": [], "primer_limit": False,
}
BUILTINS = {
    "single": {"name": "単一座位・標準", "settings": {}},
    "spaced": {"name": "区間に配置・目標100 kb", "settings": {
        "run_mode": "区間マーカーを一括設計", "spacing": 100000, "min_spacing": 50000}},
    "count": {"name": "区間から最大10座位", "settings": {
        "run_mode": "区間マーカーを一括設計", "panel_mode": "指定本数", "target_count": 10}},
    "existing": {"name": "既存プライマーを確認", "settings": {"run_mode": "既存プライマーから解析"}},
    "backcross_a": {"name": "戻し交雑（反復親A）用", "settings": {"goal": "AAとABを区別（Aへの戻し交雑）"}},
    "ril": {"name": "RIL・DH用", "settings": {"goal": "AAとBBを区別（ヘテロ判別は不要）"}},
}
RANGES = {"product_min": (70, 3000), "product_max": (70, 5000), "num_pairs": (1, 100),
          "spacing": (1, 10000000), "target_count": (1, 1000), "max_discovery": (1, 50000),
          "max_targets": (1, 1000), "min_spacing": (0, 10000000), "radius": (0, 10000000),
          "max_dcaps_mismatches": (1, 2)}
GEL_KEYS = ("ladder_choice", "ladder_text", "gel_percentage", "min_visible", "max_visible",
            "min_log_gap", "gel_min_gap")


def reset_gel_settings(state) -> None:
    for key in GEL_KEYS:
        if key in state:
            del state[key]


def apply_settings(settings: dict, state) -> None:
    merged = {**deepcopy(DEFAULTS), **{k: deepcopy(v) for k, v in settings.items() if k in DEFAULTS}}
    for key, (low, high) in RANGES.items():
        value = merged[key]
        if not isinstance(value, int) or isinstance(value, bool) or not low <= value <= high:
            raise ValueError(f"テンプレートの{key}が範囲外です。")
    if merged["product_max"] < merged["product_min"]:
        raise ValueError("PCR産物長の上限が下限より小さいテンプレートです。")
    enums = {
        "input_mode": ("参照区間＋2系統VCF", "2つのFASTA配列"),
        "run_mode": ("新規プライマー設計", "既存プライマーから解析", "区間マーカーを一括設計", "制限酵素の探索"),
        "panel_mode": ("一定間隔", "指定本数", "全酵素サイト"),
        "specificity_choice": ("厳密（標的以外の産物を認めない）", "ゲル判別型（副産物も消化して評価）"),
        "dcaps_choice": ("CAPSのみ", "CAPSが使えなければdCAPS", "CAPS/dCAPS両方"),
        "goal": tuple(GOALS),
        "enzyme_mode": ("指定した酵素", "保存した酵素セット", "全酵素"),
    }
    from .experiments import PROFILE_PRESETS
    enums["preset"] = (*PROFILE_PRESETS, "詳細設定")
    if any(merged[k] not in options for k, options in enums.items()):
        raise ValueError("テンプレートに未対応の選択項目があります。")
    gel_values = {k: v for k, v in settings.get("gel_options", {}).items() if k in GEL_KEYS}
    from .core import available_enzymes
    for key in ("enzyme_selection_manual", "enzyme_excluded"):
        if not isinstance(merged[key], list) or set(merged[key]) - set(available_enzymes()):
            raise ValueError("テンプレートに未対応の酵素があります。")
    if gel_values:
        from .experiments import GelProfile, parse_custom_ladder
        base = PROFILE_PRESETS.get(merged["preset"], GelProfile())
        ladder = gel_values.get("ladder_choice", base.ladder)
        profile = GelProfile(
            ladder=ladder,
            custom_bands=parse_custom_ladder(gel_values.get("ladder_text", "")) if ladder == "custom" else (),
            gel_percent=gel_values.get("gel_percentage", base.gel_percent),
            min_visible_bp=gel_values.get("min_visible", base.min_visible_bp),
            max_visible_bp=gel_values.get("max_visible") or None,
            min_log10_gap=gel_values.get("min_log_gap", base.min_log10_gap))
        profile.validate()
        for key, low, high in (("min_visible", 1, 5000), ("max_visible", 0, 50000), ("gel_min_gap", 1, 500)):
            if key in gel_values and not low <= gel_values[key] <= high:
                raise ValueError("テンプレートのゲル条件が範囲外です。")
    for key, value in merged.items():
        state[key] = value
    # Rebuild the gel widgets from the selected preset rather than retaining old values.
    reset_gel_settings(state)
    for key, value in gel_values.items():
        state[key] = value


def snapshot_settings(state) -> dict:
    saved = {key: deepcopy(state.get(key, default)) for key, default in DEFAULTS.items()}
    saved["enzyme_mode"] = "指定した酵素"
    saved["enzyme_selection_manual"] = list(state.get("active_enzymes", saved["enzyme_selection_manual"]))
    saved["gel_options"] = {key: deepcopy(state[key]) for key in GEL_KEYS if key in state}
    return saved
