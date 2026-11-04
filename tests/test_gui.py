"""Headless Streamlit smoke test: inputs, enzyme scan, and downloads."""
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest


def find_by_label(items, label):
    return next(item for item in items if item.label == label)


def test_gui_accepts_two_fasta_but_does_not_endorse_unchecked_scan(tmp_path, monkeypatch):
    monkeypatch.setenv("CAPS_MAKER_HOME", str(tmp_path / "settings"))
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=25).run()
    assert not at.exception
    find_by_label(at.radio, "配列の指定方法").set_value("2つのFASTA配列")
    find_by_label(at.radio, "解析方式").set_value("制限酵素の探索")
    at.run()

    a = "C" * 60 + "GAATTC" + "C" * 60
    b = "C" * 60 + "GACTTC" + "C" * 60
    find_by_label(at.text_area, "配列A（野生型など）").input(">allele_a\n" + a)
    find_by_label(at.text_area, "配列B（変異型など）").input(">allele_b\n" + b)
    at.run()
    assert not at.exception
    find_by_label(at.button, "CAPS候補を解析").click().run()
    assert not at.exception
    assert at.session_state["caps_results"]
    assert any(row["酵素"] == "EcoRI" for row in at.session_state["caps_results"])
    assert any("増幅・特異性が未評価" in item.value for item in at.info)
    assert not any(item.label == "結果CSV（Excel対応）" for item in at.get("download_button"))


def test_register_reference_is_persistent_and_selectable(tmp_path, monkeypatch):
    import pysam
    from caps_maker.catalog import load_catalog

    monkeypatch.setenv("CAPS_MAKER_HOME", str(tmp_path / "settings"))
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">Chr1\nACGTACGT\n")
    pysam.faidx(str(fasta))
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=25).run()
    assert [tab.label for tab in at.tabs] == ["マーカー設計", "ジェノタイピング", "組換え・絞り込み", "設定",
                                            "参照データ・DB", "制限酵素", "ゲル・判定条件"]
    find_by_label(at.text_input, "登録名").input("テスト参照")
    find_by_label(at.selectbox, "データ種別").set_value("参照FASTA")
    find_by_label(at.text_input, "ファイルパス / DB接頭辞").input(str(fasta))
    find_by_label(at.button, "確認して保存").click().run()
    assert not at.exception
    assert load_catalog()["databases"][0]["name"] == "テスト参照"
    reopened = AppTest.from_file(str(app), default_timeout=25).run()
    assert not reopened.exception
    selector = find_by_label(reopened.selectbox, "登録した参照データ")
    assert "テスト参照" in selector.options
    assert selector.value != "manual"


def test_enzyme_set_survives_new_session(tmp_path, monkeypatch):
    monkeypatch.setenv("CAPS_MAKER_HOME", str(tmp_path / "settings"))
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=25).run()
    find_by_label(at.multiselect, "候補酵素（複数選択）").set_value(["EcoRI", "HinfI"])
    find_by_label(at.text_input, "酵素セット名").input("手持ちテスト")
    find_by_label(at.button, "この酵素セットを保存").click().run()
    assert not at.exception
    reopened = AppTest.from_file(str(app), default_timeout=25).run()
    find_by_label(reopened.radio, "使用する酵素").set_value("保存した酵素セット").run()
    assert not reopened.exception
    assert find_by_label(reopened.multiselect, "候補酵素（複数選択）").value == ["EcoRI", "HinfI"]


def test_design_template_changes_route_spacing_and_survives_session(tmp_path, monkeypatch):
    from caps_maker.catalog import load_catalog

    monkeypatch.setenv("CAPS_MAKER_HOME", str(tmp_path / "settings"))
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=25).run()
    find_by_label(at.selectbox, "テンプレートを選ぶ").set_value("builtin:spaced")
    find_by_label(at.button, "このテンプレートを適用").click().run()
    assert not at.exception
    assert find_by_label(at.radio, "解析方式").value == "区間マーカーを一括設計"
    assert find_by_label(at.number_input, "目標間隔(bp)").value == 100000
    find_by_label(at.number_input, "目標間隔(bp)").set_value(50000)
    find_by_label(at.number_input, "増幅産物長 上限(bp)").set_value(600)
    find_by_label(at.text_input, "設定を保存する名前").input("50 kb test")
    find_by_label(at.button, "現在の設定をテンプレート保存").click().run()
    assert not at.exception
    saved = load_catalog()["design_templates"][0]
    reopened = AppTest.from_file(str(app), default_timeout=25).run()
    find_by_label(reopened.selectbox, "テンプレートを選ぶ").set_value("saved:" + saved["id"])
    find_by_label(reopened.button, "このテンプレートを適用").click().run()
    assert not reopened.exception
    assert find_by_label(reopened.number_input, "目標間隔(bp)").value == 50000
    assert find_by_label(reopened.number_input, "増幅産物長 上限(bp)").value == 600
    find_by_label(reopened.selectbox, "テンプレートを選ぶ").set_value("builtin:existing")
    find_by_label(reopened.button, "このテンプレートを適用").click().run()
    assert not reopened.exception
    assert not any(item.label == "染色体・コンティグID" for item in reopened.text_input)
    assert any(item.label == "Forward（5'→3'）" for item in reopened.text_input)
    find_by_label(reopened.selectbox, "実験プリセット").set_value("小断片重視（3%アガロース）").run()
    assert not reopened.exception
    assert find_by_label(reopened.selectbox, "ゲル濃度（%）").value == 3.0
    assert find_by_label(reopened.number_input, "検出可能な最小断片(bp)").value == 25


def test_changing_long_interval_target_invalidates_previous_results(tmp_path, monkeypatch):
    import pysam
    from caps_maker.catalog import upsert_entry
    import caps_maker.primers

    monkeypatch.setenv("CAPS_MAKER_HOME", str(tmp_path / "settings"))
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 3000 + "\n")
    pysam.faidx(str(fasta))
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=12000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        "chr1\t4001\t.\tA\tC\t60\tPASS\t.\tGT\t0/0\t1/1\n"
        "chr1\t7001\t.\tA\tC\t60\tPASS\t.\tGT\t0/0\t1/1\n"
    )
    vcfgz = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(vcf), str(vcfgz))
    pysam.tabix_index(str(vcfgz), preset="vcf")
    for name, kind, path in (("Reference", "参照FASTA", fasta),
                             ("Variants", "VCF", vcfgz), ("Search DB", "BLAST DB", tmp_path / "db")):
        upsert_entry("databases", {"name": name, "kind": kind, "path": str(path)})
    # Only the GUI's result invalidation is under test; no PCR claim is made.
    monkeypatch.setattr(caps_maker.primers, "design_caps", lambda *a, **k: ([], {}))
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=25).run()
    find_by_label(at.text_input, "染色体・コンティグID").input("chr1")
    find_by_label(at.number_input, "参照終了(1-based, inclusive)").set_value(10000)
    at.run()
    find_by_label(at.number_input, "標的の参照座標(1-based)").set_value(4001).run()
    find_by_label(at.button, "CAPS候補を解析").click().run()
    assert not at.exception
    assert at.session_state["caps_config"]["target_reference_position"] == 4001
    assert at.session_state["caps_config"]["target_variant"]
    find_by_label(at.number_input, "標的の参照座標(1-based)").set_value(7001).run()
    assert not at.exception
    assert any("入力または解析条件が変わりました" in item.value for item in at.info)
    assert not any(item.value == "使えるマーカー候補" for item in at.subheader)


@pytest.mark.integration
def test_named_resources_to_real_usable_pcr_results(tmp_path, monkeypatch):
    import random
    import shutil
    import subprocess
    import pysam
    from caps_maker.catalog import upsert_entry
    from caps_maker.genotyping import usable_marker

    if not all(shutil.which(p) for p in ("makeblastdb", "primer3_core", "blastn", "blastdbcmd")):
        pytest.skip("requires Primer3 and BLAST+")
    monkeypatch.setenv("CAPS_MAKER_HOME", str(tmp_path / "settings"))
    rng = random.Random(1841)
    bases = [rng.choice("ACGT") for _ in range(850)]
    bases[418:424] = list("GAATTC")
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">Chr4\n" + "".join(bases) + "\n")
    pysam.faidx(str(fasta))
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=Chr4,length=850>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        "Chr4\t421\t.\tA\tC\t60\tPASS\t.\tGT\t0/0\t1/1\n"
    )
    vcfgz = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(vcf), str(vcfgz))
    pysam.tabix_index(str(vcfgz), preset="vcf")
    db = tmp_path / "reference"
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-parse_seqids", "-out", str(db)],
                   check=True, capture_output=True)
    for name, kind, path in (("Test genome", "参照FASTA", fasta),
                             ("Test variants", "VCF", vcfgz), ("Test search DB", "BLAST DB", db)):
        upsert_entry("databases", {"name": name, "kind": kind, "path": str(path)})
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=30).run()
    find_by_label(at.text_input, "染色体・コンティグID").input("Chr4")
    find_by_label(at.number_input, "参照終了(1-based, inclusive)").set_value(850)
    find_by_label(at.radio, "dCAPS探索").set_value("CAPSのみ")
    at.run()
    find_by_label(at.button, "CAPS候補を解析").click().run()
    assert not at.exception
    assert any(usable_marker(r) for r in at.session_state["caps_results"])
    assert any("条件を満たす" in c.value for c in at.caption)
    tables = [frame.value for frame in at.dataframe if "特異性区分" in frame.value.columns]
    assert tables and all(str(v).startswith("A") for table in tables for v in table["特異性区分"])
    find_by_label(at.number_input, "増幅産物長 上限(bp)").set_value(600).run()
    assert not at.exception
    assert not any("特異性区分" in frame.value.columns for frame in at.dataframe)
    passing = next(r for r in at.session_state["caps_results"] if usable_marker(r))
    find_by_label(at.radio, "解析方式").set_value("既存プライマーから解析").run()
    assert not any(item.label == "染色体・コンティグID" for item in at.text_input)
    for swapped in (False, True):
        f, r = passing["Forward(5'-3')"], passing["Reverse(5'-3')"]
        if swapped:
            f, r = r, f
        find_by_label(at.text_input, "Forward（5'→3'）").input(f)
        find_by_label(at.text_input, "Reverse（5'→3'）").input(r)
        at.run()
        find_by_label(at.button, "CAPS候補を解析").click().run()
        assert not at.exception
        assert not at.error
        assert any(usable_marker(r) for r in at.session_state["caps_results"])
        assert at.session_state["caps_config"]["resolved_target"]["chrom"] == "Chr4"
        assert at.session_state["caps_summary"]["input_primers_swapped"] is swapped


def test_population_goal_and_genotyping_tab_render(tmp_path, monkeypatch):
    monkeypatch.setenv("CAPS_MAKER_HOME", str(tmp_path / "settings"))
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=25).run()
    goal = find_by_label(at.radio, "識別したい遺伝子型")
    assert len(goal.options) == 4
    find_by_label(at.selectbox, "テンプレートを選ぶ").set_value("builtin:backcross_a")
    find_by_label(at.button, "このテンプレートを適用").click().run()
    assert not at.exception
    assert find_by_label(at.radio, "識別したい遺伝子型").value == "AAとABを区別（Aへの戻し交雑）"
    assert find_by_label(at.selectbox, "集団の種類").value == "戻し交雑（反復親A）"
    assert any("使えるマーカー候補がまだありません" in item.value for item in at.info)
