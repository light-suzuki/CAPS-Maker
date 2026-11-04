"""Presentation helpers for the Streamlit GUI.

The visual language follows the Digital Agency Design System (デジタル庁デザインシステム).
Colour, type, radius and focus values are taken from ``@digital-go-jp/design-tokens``
v2.1.0 (MIT License, Copyright (c) 2023 デジタル庁). Only presentation lives here.

Styling hooks are Streamlit keyed containers: ``st.container(key="seg_…")`` renders a
short radio as a segmented control, ``st.container(key="cards_…")`` renders a radio
as selectable cards and ``st.container(key="panel_…")`` is a white panel on the grey page
(``panel_run`` also sticks to the top while the page scrolls).
"""
from __future__ import annotations

from html import escape
from typing import Sequence

import streamlit as st

_TOKENS = """
:root {
  --dads-key-50: #e8f1fe; --dads-key-100: #d9e6ff; --dads-key-200: #c5d7fb;
  --dads-key-800: #0031d8; --dads-key-900: #0017c1; --dads-key-1000: #00118f;
  --dads-key-1200: #000060;
  --dads-gray-50: #f2f2f2; --dads-gray-100: #e6e6e6; --dads-gray-200: #cccccc;
  --dads-gray-420: #949494; --dads-gray-536: #767676; --dads-gray-600: #666666;
  --dads-gray-700: #4d4d4d; --dads-gray-800: #333333; --dads-gray-900: #1a1a1a;
  --dads-success: #197a4b; --dads-success-bg: #e6f5ec;
  --dads-error: #ce0000; --dads-error-bg: #fdeeee;
  --dads-warning: #927200; --dads-warning-bg: #fbf5e0;
  --dads-focus: #ffd43d; --dads-radius: 8px;
}
"""

_SEG = '[class*="st-key-seg_"]'  # must stay a single selector: it prefixes descendant rules
_CARDS = '[class*="st-key-cards_"]'

_CSS = f"""
html, body {{ letter-spacing: .02em; }}
.stApp {{ background: #f3f5f8; }}
[data-testid="stMainBlockContainer"] {{ max-width: 1280px; padding: 1rem 2rem 4rem; }}
[data-testid="stHeader"] {{ background: transparent; height: 0; }}
[data-testid="stSidebarCollapsedControl"], [data-testid="stSidebar"] {{ display: none; }}
p, li {{ line-height: 1.7; }}
h3 {{ font-size: 1.25rem !important; font-weight: 700 !important; padding-top: .75rem !important; }}
h4 {{ font-size: 1.0625rem !important; font-weight: 700 !important; }}

/* app bar */
.caps-bar {{ display: flex; align-items: center; justify-content: space-between; gap: 1rem 2rem;
  flex-wrap: wrap; padding: 1rem 1.5rem; border-radius: 14px; color: #fff;
  background: linear-gradient(120deg, var(--dads-key-1200) 0%, var(--dads-key-1000) 45%, var(--dads-key-900) 100%); }}
.caps-bar__brand {{ display: flex; align-items: center; gap: .875rem; }}
.caps-bar__logo {{ flex: none; width: 40px; height: 40px; }}
.caps-bar__name {{ font-size: 1.375rem; font-weight: 700; line-height: 1.2; }}
.caps-bar__sub {{ font-size: .8125rem; opacity: .82; line-height: 1.5; }}
.caps-ctx {{ display: flex; flex-wrap: wrap; gap: .5rem; }}
.caps-ctx__item {{ display: flex; flex-direction: column; gap: .1rem; padding: .4rem .8rem;
  border-radius: 8px; background: rgba(255,255,255,.1); border: 1px solid rgba(255,255,255,.18); min-width: 6.5rem; }}
.caps-ctx__k {{ font-size: .6875rem; letter-spacing: .08em; opacity: .75; line-height: 1.3; }}
.caps-ctx__v {{ font-size: .875rem; font-weight: 700; line-height: 1.35; white-space: nowrap; }}
.caps-ctx__item--off .caps-ctx__v {{ color: var(--dads-focus); }}
@media (max-width: 640px) {{
  [data-testid="stMainBlockContainer"] {{ padding: .5rem 1rem 3rem; }}
  .caps-bar {{ padding: .875rem 1rem; }}
  .caps-bar__logo, .caps-bar__sub {{ display: none; }}
  .caps-ctx__item {{ min-width: 0; }}
}}

/* panels: st.container(key="panel_…") */
[class*="st-key-panel_"] {{ background: #fff; border: 1px solid #e3e6eb; border-radius: 14px;
  padding: 1.25rem 1.5rem 1.5rem; box-shadow: 0 1px 2px rgba(16,24,40,.04); }}
[data-testid="stColumn"]:has(.st-key-panel_run) {{ position: sticky; top: 1rem; align-self: flex-start; }}
@media (max-width: 640px) {{
  [class*="st-key-panel_"] {{ padding: 1rem 1rem 1.25rem; border-radius: 12px; }}
  [data-testid="stColumn"]:has(.st-key-panel_run) {{ position: static; }}
}}
.caps-ph {{ display: flex; align-items: flex-start; gap: .75rem; margin-bottom: .25rem; }}
.caps-ph__num {{ flex: none; width: 1.75rem; height: 1.75rem; border-radius: 999px; margin-top: .05rem;
  display: inline-flex; align-items: center; justify-content: center;
  background: var(--dads-key-900); color: #fff; font-weight: 700; font-size: .875rem; }}
.caps-ph__title {{ font-size: 1.125rem; font-weight: 700; line-height: 1.5; color: var(--dads-gray-900); }}
.caps-ph__lead {{ font-size: .875rem; color: var(--dads-gray-700); line-height: 1.6; }}
.caps-lead {{ font-size: .9375rem; color: var(--dads-gray-700); line-height: 1.6; margin: -.25rem 0 .25rem; }}

.caps-sub {{ font-size: .9375rem; font-weight: 700; margin: .75rem 0 0; color: var(--dads-gray-800); }}

/* radio as selectable cards  (DOM: radiogroup > div[data-selected] > label + caption) */
{_CARDS} [data-testid="stElementContainer"], {_CARDS} [data-testid="stRadio"] {{ width: 100%; }}
{_CARDS} [role="radiogroup"] {{ display: grid; gap: .75rem; width: 100%;
  grid-template-columns: repeat(2, minmax(0, 1fr)); }}
@media (max-width: 640px) {{ {_CARDS} [role="radiogroup"] {{ grid-template-columns: 1fr; }} }}
{_CARDS} [role="radiogroup"] > div, {_CARDS} [data-testid="stRadioOption"] {{ width: 100%; }}
{_CARDS} [data-testid="stRadioOption"] p {{ text-align: left; }}
{_CARDS} [role="radiogroup"] > div {{ margin: 0; padding: .875rem 1rem 1rem; position: relative;
  border: 1px solid var(--dads-gray-420); border-radius: var(--dads-radius); background: #fff;
  transition: border-color .12s, background-color .12s; }}
{_CARDS} [role="radiogroup"] > div:hover {{ border-color: var(--dads-key-900); background: var(--dads-key-50); }}
{_CARDS} [role="radiogroup"] > div[data-selected="true"] {{ background: var(--dads-key-50);
  border-color: var(--dads-key-900); box-shadow: inset 0 0 0 1px var(--dads-key-900); }}
/* the whole card is the click target */
{_CARDS} [data-testid="stRadioOption"], {_CARDS} [data-testid="stRadioOption"] > div {{ position: static !important; }}
{_CARDS} [data-testid="stRadioOption"]::after {{ content: ""; position: absolute; inset: 0; z-index: 1;
  cursor: pointer; }}
{_CARDS} [data-testid="stRadioOption"] p {{ font-weight: 700; font-size: 1rem; color: var(--dads-gray-900); }}
{_CARDS} [data-selected="true"] [data-testid="stRadioOption"] p {{ color: var(--dads-key-1000); }}
{_CARDS} [data-testid="stRadioCaption"] {{ padding-left: 1.75rem; }}
{_CARDS} [data-testid="stCaptionContainer"] p {{ font-weight: 400 !important; font-size: .875rem !important;
  color: var(--dads-gray-700) !important; line-height: 1.6; }}

/* short radio as segmented control */
{_SEG} [role="radiogroup"] {{ gap: 0; display: inline-flex; flex-wrap: wrap; width: fit-content;
  border: 1px solid var(--dads-key-900); border-radius: var(--dads-radius); overflow: hidden; }}
{_SEG} [role="radiogroup"] > div {{ margin: 0; background: #fff; border-right: 1px solid var(--dads-key-900); }}
{_SEG} [role="radiogroup"] > div:last-child {{ border-right: 0; }}
{_SEG} [data-testid="stRadioOption"] {{ margin: 0; padding: .5rem 1.125rem; cursor: pointer; }}
{_SEG} [data-testid="stRadioOption"] > div > div:first-child {{ display: none; }}
{_SEG} [data-testid="stRadioOption"] p {{ color: var(--dads-key-900); }}
{_SEG} [role="radiogroup"] > div:hover {{ background: var(--dads-key-50); }}
{_SEG} [role="radiogroup"] > div[data-selected="true"] {{ background: var(--dads-key-900); }}
{_SEG} [data-selected="true"] [data-testid="stRadioOption"] p {{ color: #fff; font-weight: 700; }}

/* form controls */
[data-testid="stTextAreaRootElement"], [data-testid="stTextInputRootElement"],
[data-testid="stNumberInputContainer"], [data-testid="stSelectbox"] div:has(> input),
[data-testid="stMultiSelect"] div:has(> input) {{
  border-color: var(--dads-gray-600) !important; background-color: #fff !important; }}
[data-testid="stNumberInputContainer"] button {{ background: #fff; }}
[data-testid="stTextAreaRootElement"]:hover, [data-testid="stTextInputRootElement"]:hover,
[data-testid="stNumberInputContainer"]:hover, [data-testid="stSelectbox"] div:has(> input):hover {{
  border-color: var(--dads-gray-900) !important; }}
[data-testid="stTextInputRootElement"]:has(input:disabled),
[data-testid="stNumberInputContainer"]:has(input:disabled),
[data-testid="stTextAreaRootElement"]:has(textarea:disabled) {{
  background-color: var(--dads-gray-50) !important; border-color: var(--dads-gray-200) !important; }}
[data-testid="stWidgetLabel"] p {{ font-weight: 700; color: var(--dads-gray-800); }}
[data-baseweb="tag"] {{ background: var(--dads-key-100) !important; border-radius: 6px !important; }}
[data-baseweb="tag"] span {{ color: var(--dads-key-1000) !important; font-weight: 700; }}
textarea {{ font-family: 'Noto Sans Mono', ui-monospace, monospace !important; font-size: .875rem !important; }}

/* buttons */
[data-testid="stBaseButton-primary"], [data-testid="stBaseButton-primaryFormSubmit"] {{
  min-height: 3rem; padding: .5rem 1.75rem; font-weight: 700; font-size: 1rem;
  background: var(--dads-key-900); border: 1px solid var(--dads-key-900); color: #fff; }}
[data-testid="stBaseButton-primary"]:hover, [data-testid="stBaseButton-primaryFormSubmit"]:hover {{
  background: var(--dads-key-1000); border-color: var(--dads-key-1000); color: #fff;
  text-decoration: underline; text-underline-offset: 3px; }}
[data-testid="stBaseButton-primary"]:active {{ background: var(--dads-key-1200); }}
[data-testid="stBaseButton-primary"]:disabled {{ background: var(--dads-gray-200);
  border-color: var(--dads-gray-200); color: var(--dads-gray-600); text-decoration: none; }}
[data-testid="stBaseButton-secondary"], [data-testid="stBaseButton-secondaryFormSubmit"] {{
  min-height: 2.75rem; font-weight: 700; background: #fff;
  border: 1px solid var(--dads-key-900); color: var(--dads-key-900); }}
[data-testid="stBaseButton-secondary"]:hover, [data-testid="stBaseButton-secondaryFormSubmit"]:hover {{
  background: var(--dads-key-50); border-color: var(--dads-key-1000); color: var(--dads-key-1000);
  text-decoration: underline; text-underline-offset: 3px; }}
[data-testid="stBaseButton-secondary"]:disabled {{ border-color: var(--dads-gray-200);
  color: var(--dads-gray-420); background: #fff; text-decoration: none; }}
[data-testid="stDownloadButton"] button {{ width: 100%; justify-content: flex-start; }}

/* keyboard focus: black outline with yellow ring */
:focus-visible, label:has(input:focus-visible) {{ outline: 4px solid #000 !important;
  outline-offset: 2px !important; box-shadow: 0 0 0 2px var(--dads-focus) !important; }}

/* tabs (react-aria DOM): the top level inside key "mainnav" is the main navigation */
[role="tablist"] {{ gap: .25rem; border-bottom: 1px solid var(--dads-gray-200); }}
[data-testid="stTab"] {{ padding: .625rem 1rem; height: auto; }}
[data-testid="stTab"] p {{ font-size: .9375rem; color: var(--dads-gray-700); }}
[data-testid="stTab"][aria-selected="true"] p {{ color: var(--dads-key-900); font-weight: 700; }}
[data-testid="stTab"]:hover p {{ color: var(--dads-key-1000); text-decoration: underline; }}
[data-testid="stTab"] .react-aria-SelectionIndicator {{ height: 3px; background: var(--dads-key-900); }}
.st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] {{ margin: .75rem 0 1.25rem; padding: .3rem; gap: .25rem; border: 1px solid #e3e6eb;
  border-radius: 12px; background: #fff; width: fit-content; max-width: 100%; overflow-x: auto;
  box-shadow: 0 1px 2px rgba(16,24,40,.04); }}
.st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"] {{ padding: .55rem 1.25rem; border-radius: 9px; margin: 0; }}
.st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"] p {{ font-size: 1rem; font-weight: 700; color: var(--dads-gray-700); }}
.st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"]:hover {{ background: var(--dads-key-50); }}
.st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"][aria-selected="true"] {{ background: var(--dads-key-900); }}
.st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"][aria-selected="true"] p {{ color: #fff; text-decoration: none; }}
.st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"] .react-aria-SelectionIndicator {{ display: none; }}
@media (max-width: 640px) {{
  .st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"] {{ padding: .45rem .8rem; }}
  .st-key-mainnav > [data-testid="stTabs"] > div > [role="tablist"] > [data-testid="stTab"] p {{ font-size: .875rem; }}
}}

/* banners and panels */
[data-testid="stAlertContainer"] {{ border-radius: var(--dads-radius); border: 1px solid currentColor; }}
[data-testid="stExpander"] details {{ border-color: var(--dads-gray-200); border-radius: var(--dads-radius); }}
[data-testid="stExpander"] summary p {{ font-weight: 700; }}
[data-testid="stForm"] {{ border-color: var(--dads-gray-200); border-radius: var(--dads-radius); }}

/* stat cards */
.caps-stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: .75rem; margin: .25rem 0 1rem; }}
.caps-stat {{ border: 1px solid var(--dads-gray-200); border-radius: var(--dads-radius);
  padding: .75rem 1rem; background: #fff; }}
.caps-stat__label {{ font-size: .8125rem; color: var(--dads-gray-700); line-height: 1.5; }}
.caps-stat__value {{ font-size: 1.75rem; font-weight: 700; line-height: 1.3;
  font-variant-numeric: tabular-nums; color: var(--dads-gray-900); }}
.caps-stat__unit {{ font-size: .875rem; font-weight: 400; margin-left: .2rem; color: var(--dads-gray-700); }}
.caps-stat--key {{ border-color: var(--dads-key-900); background: var(--dads-key-50); }}
.caps-stat--key .caps-stat__value {{ color: var(--dads-key-900); }}
.caps-stat--ok .caps-stat__value {{ color: var(--dads-success); }}

/* status tags and readiness checklist */
.caps-chips {{ display: flex; flex-wrap: wrap; gap: .5rem; margin: .25rem 0; }}
.caps-tag {{ display: inline-flex; align-items: center; gap: .375rem; font-size: .875rem; line-height: 1.4;
  padding: .3rem .7rem; border-radius: 6px; background: var(--dads-gray-50); color: var(--dads-gray-800); }}
.caps-tag--ok {{ background: var(--dads-success-bg); color: var(--dads-success); }}
.caps-tag--key {{ background: var(--dads-key-50); color: var(--dads-key-1000); }}
.caps-check {{ list-style: none; margin: .25rem 0 1rem; padding: 0; display: grid; gap: .375rem; }}
.caps-check li {{ display: flex; gap: .625rem; align-items: baseline; font-size: .9375rem; line-height: 1.6; }}
.caps-check__icon {{ flex: none; width: 1.25rem; height: 1.25rem; border-radius: 999px; font-size: .75rem;
  font-weight: 700; display: inline-flex; align-items: center; justify-content: center;
  transform: translateY(.2rem); }}
.caps-check__ok .caps-check__icon {{ background: var(--dads-success); color: #fff; }}
.caps-check__ng .caps-check__icon {{ background: #fff; color: var(--dads-error);
  box-shadow: inset 0 0 0 2px var(--dads-error); }}
.caps-check__ng {{ color: var(--dads-gray-900); }}
.caps-check small {{ color: var(--dads-gray-700); font-size: .8125rem; }}

/* run summary */
.caps-dl {{ margin: .25rem 0 .75rem; display: grid; grid-template-columns: auto 1fr; column-gap: 1rem; }}
.caps-dl dt {{ font-size: .8125rem; color: var(--dads-gray-700); line-height: 1.5; padding: .4rem 0;
  border-bottom: 1px solid var(--dads-gray-100); }}
.caps-dl dd {{ margin: 0; font-size: .875rem; font-weight: 700; color: var(--dads-gray-900);
  line-height: 1.5; padding: .4rem 0; border-bottom: 1px solid var(--dads-gray-100);
  text-align: right; word-break: break-all; }}
.caps-note {{ font-size: .8125rem; line-height: 1.7; color: var(--dads-gray-700);
  border-left: 4px solid var(--dads-gray-420); padding: .25rem 0 .25rem .75rem; }}
.caps-footer {{ margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--dads-gray-200);
  font-size: .8125rem; color: var(--dads-gray-700); line-height: 1.7; }}
[data-testid="stImage"] img {{ border-radius: 6px; }}
"""

_LOGO = """<svg class="caps-bar__logo" viewBox="0 0 44 44" aria-hidden="true">
<rect width="44" height="44" rx="10" fill="#fff" fill-opacity=".14"/>
<rect x="8" y="8" width="8" height="28" rx="2" fill="#fff" opacity=".16"/>
<rect x="18" y="8" width="8" height="28" rx="2" fill="#fff" opacity=".16"/>
<rect x="28" y="8" width="8" height="28" rx="2" fill="#fff" opacity=".16"/>
<rect x="8" y="14" width="8" height="2.5" rx="1" fill="#fff"/>
<rect x="18" y="14" width="8" height="2.5" rx="1" fill="#fff"/>
<rect x="18" y="22" width="8" height="2.5" rx="1" fill="#ffd43d"/>
<rect x="18" y="29" width="8" height="2.5" rx="1" fill="#ffd43d"/>
<rect x="28" y="22" width="8" height="2.5" rx="1" fill="#fff"/>
<rect x="28" y="29" width="8" height="2.5" rx="1" fill="#fff"/>
</svg>"""


def inject_styles() -> None:
    st.html(f"<style>{_TOKENS}{_CSS}</style>")


def header(context: Sequence[tuple[str, str, bool]] = ()) -> str:
    """App bar HTML; context items are (label, value, ok)."""
    items = "".join(
        f'<div class="caps-ctx__item{"" if ok else " caps-ctx__item--off"}">'
        f'<span class="caps-ctx__k">{escape(k)}</span><span class="caps-ctx__v">{escape(v)}</span></div>'
        for k, v, ok in context)
    return ('<div class="caps-bar"><div class="caps-bar__brand">' + _LOGO +
            '<div><div class="caps-bar__name">CAPS-Maker</div>'
            '<div class="caps-bar__sub">ジェノタイピング用 CAPS / dCAPS マーカー・ローカル in silico 予測</div>'
            f'</div></div><div class="caps-ctx">{items}</div></div>')


def panel_head(title: str, lead: str | None = None, num: int | None = None) -> None:
    num_html = f'<span class="caps-ph__num">{num}</span>' if num is not None else ""
    lead_html = f'<div class="caps-ph__lead">{escape(lead)}</div>' if lead else ""
    st.html(f'<div class="caps-ph">{num_html}<div><div class="caps-ph__title">{escape(title)}</div>'
            f'{lead_html}</div></div>')


def lead(text: str) -> None:
    st.html(f'<div class="caps-lead">{escape(text)}</div>')


class StepCounter:
    """Numbered panel headings, so a hidden step never leaves a gap."""

    def __init__(self) -> None:
        self.n = 0

    def __call__(self, title: str, lead: str | None = None) -> None:
        self.n += 1
        panel_head(title, lead, self.n)


def subheading(text: str) -> None:
    st.html(f'<div class="caps-sub">{escape(text)}</div>')


def chips(items: Sequence[tuple[str, str]]) -> None:
    """Status tags; each item is (plain text, tone) with tone in '', 'ok', 'key'."""
    parts = [f'<span class="caps-tag{" caps-tag--" + tone if tone else ""}">{escape(text)}</span>'
             for text, tone in items]
    if parts:
        st.html('<div class="caps-chips">' + "".join(parts) + "</div>")


def checklist(items: Sequence[tuple[bool, str, str]]) -> None:
    """Prerequisites before running; each item is (ok, label, detail)."""
    rows = []
    for ok, label, detail in items:
        cls, icon = ("caps-check__ok", "✓") if ok else ("caps-check__ng", "!")
        extra = f" <small>{escape(detail)}</small>" if detail else ""
        rows.append(f'<li class="{cls}"><span class="caps-check__icon" aria-hidden="true">{icon}</span>'
                    f'<span><span class="sr-only">{"完了" if ok else "未完了"}：</span>'
                    f'{escape(label)}{extra}</span></li>')
    st.html('<ul class="caps-check">' + "".join(rows) + "</ul>"
            '<style>.sr-only{position:absolute;width:1px;height:1px;overflow:hidden;'
            'clip:rect(0 0 0 0);white-space:nowrap}</style>')


def stats(items: Sequence[tuple[str, int | str, str]]) -> None:
    """Stat cards; each item is (label, value, tone) with tone in '', 'key', 'ok'."""
    cards = []
    for label, value, tone in items:
        cls = "caps-stat" + (f" caps-stat--{tone}" if tone else "")
        cards.append(f'<div class="{cls}"><div class="caps-stat__label">{escape(label)}</div>'
                     f'<div class="caps-stat__value">{escape(str(value))}'
                     '<span class="caps-stat__unit">件</span></div></div>')
    st.html('<div class="caps-stats">' + "".join(cards) + "</div>")


def summary(rows: Sequence[tuple[str, str]]) -> None:
    body = "".join(f"<dt>{escape(k)}</dt><dd>{escape(v)}</dd>" for k, v in rows)
    st.html(f'<dl class="caps-dl">{body}</dl>')


def footer(text: str) -> None:
    st.html(f'<div class="caps-footer">{escape(text)}</div>')
