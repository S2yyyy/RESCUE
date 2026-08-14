#!/usr/bin/env bash
# Verify a regenerated data tree against the checksums of the data used in
# the paper. Run from anywhere; checks files under $RESCUE_DATA_ROOT
# (default: ./data).
set -euo pipefail
DATA_ROOT="${RESCUE_DATA_ROOT:-$PWD/data}"
LIST="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/checksums.md5"
cd "$DATA_ROOT"
grep -v '^#' "$LIST" | md5sum -c -
