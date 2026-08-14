"""Release integrity checks for the RESCUE repository."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

FORBIDDEN_PATTERNS = [
    r"(^|[\s\"'=])/(home|data)/[A-Za-z0-9._-]+",
    r"(^|[\s\"'=])/add[0-9]+/users/",
    r"10\.[0-9]+\.[0-9]+\.[0-9]+",
    r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}",
    r"\." + "codex",
    r"\." + "claude",
    "paper" + r"-backup",
    "CUDA" + r"_VISIBLE_DEVICES",
    r"\b" + ("S" + "LURM") + r"\b",
    r"\b" + ("s" + "batch") + r"\b",
    r"conda run -n [A-Za-z0-9._-]+",
    "%" + "F %" + "T",
    r"date \+%s",
]

REQUIRED_FILES = [
    "README.md",
    "pyproject.toml",
    "requirements.txt",
    "docs/DATA.md",
    "scripts/prepare_data/README.md",
    "data/README.md",
    "data/v3/mouse_paralog/proteins.fa",
    "configs/exp1_tools.yaml",
    "scripts/rescue_paths.py",
    "scripts/stage0/normalize.py",
    "scripts/exp1/common.py",
    "scripts/exp5/common.py",
    "scripts/tprs/phase5_run.py",
    "results/exp1/summary.tsv",
    "results/exp5/cell/per_cell.tsv",
    "results/tprs/phase5/ablation_summary_all_species.json",
    "results/tprs/phase5/heuristic_baselines.json",
    "results/tprs/phase5/feature_importance_phase5.json",
]

EXCLUDED_SUBSTRINGS = [
    "genome" + "_" + "first",
    "label" + "_" + "free",
    "exp" + "9",
    "deploy" + "able",
    "." + "codex",
    "." + "claude",
]


def iter_text_files():
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix()
        if ".git/" in rel:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        yield rel, text


def test_required_files_exist():
    missing = [p for p in REQUIRED_FILES if not (ROOT / p).exists()]
    assert not missing, f"missing required release files: {missing}"


def test_no_excluded_release_paths():
    all_paths = [p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")]
    offenders = [
        p for p in all_paths
        if any(fragment in p for fragment in EXCLUDED_SUBSTRINGS)
    ]
    assert not offenders, f"excluded paths present: {offenders[:20]}"


def test_no_private_patterns_in_text_files():
    hits = []
    compiled = [re.compile(p) for p in FORBIDDEN_PATTERNS]
    for rel, text in iter_text_files():
        for lineno, line in enumerate(text.splitlines(), 1):
            for pattern in compiled:
                if pattern.search(line):
                    hits.append(f"{rel}:{lineno}: {line[:160]}")
    assert not hits, "private or machine-specific patterns found:\n" + "\n".join(hits[:50])


def test_summary_artifacts_parse():
    exp1 = ROOT / "results/exp1/summary.tsv"
    exp5 = ROOT / "results/exp5/cell/per_cell.tsv"

    with exp1.open() as fh:
        exp1_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(exp1_rows) >= 40
    assert {"tool", "species", "regime", "exon_f1_mean"} <= set(exp1_rows[0])

    with exp5.open() as fh:
        exp5_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(exp5_rows) >= 100
    assert {"tool", "species", "regime", "padding", "mean_exon_f1"} <= set(exp5_rows[0])

    for rel in [
        "results/tprs/phase5/ablation_summary_all_species.json",
        "results/tprs/phase5/heuristic_baselines.json",
        "results/tprs/phase5/feature_importance_phase5.json",
    ]:
        data = json.loads((ROOT / rel).read_text())
        assert isinstance(data, dict) and data


def main() -> int:
    test_required_files_exist()
    test_no_excluded_release_paths()
    test_no_private_patterns_in_text_files()
    test_summary_artifacts_parse()
    print("release integrity checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
