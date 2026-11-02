"""Headless Streamlit smoke test: inputs, enzyme scan, and downloads."""
from pathlib import Path

from streamlit.testing.v1 import AppTest


def find_by_label(items, label):
    return next(item for item in items if item.label == label)


def test_gui_accepts_two_fasta_and_shows_caps_table():
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=25).run()
    assert not at.exception

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
