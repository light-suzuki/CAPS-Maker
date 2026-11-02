"""Headless GUI smoke test: no browser is required."""
from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_gui_accepts_two_fasta_and_shows_caps_table():
    app = Path(__file__).resolve().parents[1] / "app.py"
    at = AppTest.from_file(str(app), default_timeout=20).run()
    assert not at.exception

    a = "C" * 60 + "GAATTC" + "C" * 60
    b = "C" * 60 + "GACTTC" + "C" * 60
    at.text_area[0].input(">allele_a\n" + a)
    at.text_area[1].input(">allele_b\n" + b)
    at.run()
    assert not at.exception
    at.button[0].click().run()
    assert not at.exception
    assert at.session_state.caps_rows
    assert any(row["酵素"] == "EcoRI" for row in at.session_state.caps_rows)
