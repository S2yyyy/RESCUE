"""Repository-local paths and external-tool lookup for RESCUE scripts.

The release defaults to paths inside the checkout. Large reference files and
tool installations can be supplied with environment variables documented in
the README.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path


def _path_from_env(name: str, default: Path) -> Path:
    return Path(os.environ.get(name, default)).expanduser().resolve()


PROJECT_ROOT = _path_from_env(
    "RESCUE_REPO", Path(__file__).resolve().parents[1]
)
DATA_ROOT = _path_from_env("RESCUE_DATA_ROOT", PROJECT_ROOT / "data")
RESULTS_ROOT = _path_from_env("RESCUE_RESULTS_ROOT", PROJECT_ROOT / "results")
TOOLS_ROOT = _path_from_env("RESCUE_TOOLS_ROOT", PROJECT_ROOT / "tools")

PROCESSED_DATA = _path_from_env("RESCUE_PROCESSED_DATA", DATA_ROOT / "processed")
V3_DATA = _path_from_env("RESCUE_V3_DATA", DATA_ROOT / "v3")
V3_RESULTS = _path_from_env("RESCUE_V3_RESULTS", RESULTS_ROOT / "v3")

GENOMES = _path_from_env("RESCUE_GENOME_DIR", PROCESSED_DATA / "genomes")
GENE_STRUCTURES = _path_from_env(
    "RESCUE_GENE_STRUCTURE_DIR", PROCESSED_DATA / "gene_structures"
)
PROTEOMES = _path_from_env("RESCUE_PROTEOME_DIR", PROCESSED_DATA / "proteins")
ORTHOLOGS = _path_from_env("RESCUE_ORTHOLOG_DIR", PROCESSED_DATA / "orthologs")

CONDA = os.environ.get("CONDA_EXE") or shutil.which("conda") or "conda"


def binary(env_var: str, default_name: str) -> str:
    """Return an external executable from env, PATH, or its bare command name."""
    return os.environ.get(env_var) or shutil.which(default_name) or default_name


MINIPROT = binary("MINIPROT_BIN", "miniprot")
SAMTOOLS = binary("SAMTOOLS_BIN", "samtools")
TBLASTN = binary("TBLASTN_BIN", "tblastn")
PERL = binary("PERL_BIN", "perl")

