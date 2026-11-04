"""Local resource, enzyme-set and design-template catalog."""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
from uuid import uuid4


def catalog_path() -> Path:
    root = os.environ.get("CAPS_MAKER_HOME")
    return (Path(root).expanduser() if root else Path.home() / ".codex" / "caps-maker") / "catalog.json"


def load_catalog(path: Path | None = None) -> dict:
    path = path or catalog_path()
    if not path.exists():
        return {"schema": 1, "databases": [], "enzyme_sets": [], "design_templates": []}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError("保存データを読めません。既存ファイルは上書きしていません。") from exc
    if (not isinstance(data, dict) or data.get("schema") != 1
            or any(not isinstance(data.get(k), list) for k in ("databases", "enzyme_sets"))):
        raise ValueError("保存データの形式が未対応です。既存ファイルは上書きしていません。")
    data.setdefault("design_templates", [])
    if not isinstance(data["design_templates"], list):
        raise ValueError("設計テンプレートの保存形式が不正です。")
    return data


def save_catalog(data: dict, path: Path | None = None) -> None:
    path = path or catalog_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=".catalog-", suffix=".json", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(data, handle, ensure_ascii=False, indent=2)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def upsert_entry(section: str, entry: dict, path: Path | None = None) -> dict:
    if section not in ("databases", "enzyme_sets", "design_templates"):
        raise ValueError("未対応の登録先です。")
    entry = dict(entry)
    entry["name"] = str(entry.get("name", "")).strip()
    if not entry["name"]:
        raise ValueError("名前を入力してください。")
    data = load_catalog(path)
    entry.setdefault("id", uuid4().hex)
    for old in data[section]:
        if (old["id"] != entry["id"] and old["name"].casefold() == entry["name"].casefold()):
            raise ValueError("同じ名前が登録済みです。更新対象を選ぶか、別の名前にしてください。")
    data[section] = [old for old in data[section] if old["id"] != entry["id"]] + [entry]
    save_catalog(data, path)
    return entry


def local_path(value: str) -> str:
    value = value.strip().strip('"').strip("'")
    if not value:
        raise ValueError("パスを入力してください。")
    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", value):
        value = "/mnt/" + value[0].lower() + "/" + value[3:].replace("\\", "/")
    return str(Path(value).expanduser())


def validate_resource(kind: str, value: str) -> tuple[str, str]:
    """Check actual readability/indexes, without asserting assembly identity."""
    path = local_path(value)
    if kind == "BLAST DB":
        try:
            result = subprocess.run(["blastdbcmd", "-db", path, "-dbtype", "nucl", "-info"],
                                    capture_output=True, text=True, timeout=20)
        except FileNotFoundError as exc:
            raise ValueError("blastdbcmdが見つかりません。BLAST+を確認してください。") from exc
        except subprocess.TimeoutExpired as exc:
            raise ValueError("DBの確認が20秒以内に完了しませんでした。") from exc
        if result.returncode:
            raise ValueError("BLAST DBを開けません。フォルダーではなくDBの接頭辞を指定してください。\n"
                             + result.stderr.strip()[-500:])
        return path, "nucleotide DBを確認"
    if kind == "参照FASTA":
        import pysam
        if not Path(path + ".fai").exists():
            raise ValueError("参照FASTAの.fai索引がありません。")
        with pysam.FastaFile(path) as handle:
            if not handle.references:
                raise ValueError("参照FASTAに配列がありません。")
            handle.fetch(handle.references[0], 0, 1)
            return path, f"FASTA・索引を確認 / {len(handle.references)}配列"
    if kind == "VCF":
        import pysam
        if not any(Path(path + suffix).exists() for suffix in (".tbi", ".csi")):
            raise ValueError("VCFの.tbiまたは.csi索引がありません。")
        with pysam.VariantFile(path) as handle:
            if not handle.header.samples or not handle.header.contigs:
                raise ValueError("VCFにサンプルまたは染色体ヘッダーがありません。")
            next(handle.fetch(next(iter(handle.header.contigs)), 0, 1), None)
            return path, f"VCF・索引を確認 / {len(handle.header.samples)}サンプル"
    raise ValueError("未対応のデータ種別です。")
