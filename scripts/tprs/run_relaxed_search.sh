#!/usr/bin/env bash
# TPRS Phase 1 — Broad-recall A-stage with relaxed tool parameters.
# Runs only paralog regime (1:1 is saturated by miniprot family).
#
# Tools:
#   miniprot-ultra  — zero threshold, 500 max hits, 500kb max intron
#   diamond-relaxed — max-target-seqs 5000 (up from 1000)
#   mmseqs-relaxed  — max-seqs 5000, sensitivity -s 1.0 (down from 5.7)
#   last-relaxed    — E-value 1.0 (up from 0.05), frameshift cost 20 (up from 15)
#
# Output: results/tprs/raw/<tool>.<species>.paralog.<ext>
set -euo pipefail

REPO="${RESCUE_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
DATA_ROOT="${RESCUE_DATA_ROOT:-"$REPO/data"}"
GENOMES="${RESCUE_GENOME_DIR:-"$DATA_ROOT/processed/genomes"}"
DBROOT=$REPO/results/exp3/db
RAW=$REPO/results/tprs/raw
mkdir -p "$RAW"

THREADS="${THREADS:-16}"
DIAMOND_BIN="${DIAMOND_BIN:-diamond}"
MMSEQS_BIN="${MMSEQS_BIN:-mmseqs}"
LASTAL_BIN="${LASTAL_BIN:-lastal}"
SPECIES=(mouse dog chicken zebrafish)

log() {
  local tool="$1"; shift
  printf '[RESCUE] %s %s\n' "$tool" "$*" | tee -a "$RAW/${tool}.log"
}

# --------------- miniprot ultra-relaxed ---------------
run_miniprot_ultra_cell() {
  local sp=$1
  local mpi="$GENOMES/${sp}.mpi"
  local proteins="$REPO/data/v3/${sp}_paralog/proteins.fa"
  local out="$RAW/miniprot_ultra.${sp}.paralog.gff"
  if [ -s "$out" ]; then
    log miniprot_ultra "skip ${sp} (exists, $(wc -l <"$out") lines)"
    return
  fi
  if [ ! -f "$mpi" ]; then
    log miniprot_ultra "skip ${sp} (no index at $mpi)"
    return
  fi
  local t=$SECONDS
  log miniprot_ultra "begin ${sp}"
  miniprot --gff -t "$THREADS" -N 500 -p 0.0 --outs=0.0 --outc=0.0 -G 500000 \
    "$mpi" "$proteins" > "$out" 2>>"$RAW/miniprot_ultra.log"
  log miniprot_ultra "done ${sp} in $(( SECONDS - t )) s; lines=$(wc -l <"$out")"
}

# --------------- miniprot relaxed (practical sweet-spot) ---------------
run_miniprot_relaxed_cell() {
  local sp=$1
  local mpi="$GENOMES/${sp}.mpi"
  local proteins="$REPO/data/v3/${sp}_paralog/proteins.fa"
  local out="$RAW/miniprot_relaxed.${sp}.paralog.gff"
  if [ -s "$out" ]; then
    log miniprot_relaxed "skip ${sp} (exists, $(wc -l <"$out") lines)"
    return
  fi
  if [ ! -f "$mpi" ]; then
    log miniprot_relaxed "skip ${sp} (no index at $mpi)"
    return
  fi
  local t=$SECONDS
  log miniprot_relaxed "begin ${sp}"
  miniprot --gff -t "$THREADS" -N 50 --outn=50 --outs=0.2 --outc=0.05 \
    "$mpi" "$proteins" > "$out" 2>>"$RAW/miniprot_relaxed.log"
  log miniprot_relaxed "done ${sp} in $(( SECONDS - t )) s; lines=$(wc -l <"$out")"
}

# --------------- DIAMOND relaxed ---------------
run_diamond_relaxed_cell() {
  local sp=$1
  local db="$DBROOT/diamond_protdb/${sp}.paralog"
  local genome="$GENOMES/${sp}.fa"
  local out="$RAW/diamond_relaxed.${sp}.paralog.tsv"
  if [ -s "$out" ]; then
    log diamond_relaxed "skip ${sp} (exists, $(wc -l <"$out") lines)"
    return
  fi
  if [ ! -f "${db}.dmnd" ]; then
    log diamond_relaxed "skip ${sp} (no DB at $db)"
    return
  fi
  local t=$SECONDS
  log diamond_relaxed "begin ${sp}"
  "$DIAMOND_BIN" blastx \
    --db "$db" \
    --query "$genome" \
    --frameshift 15 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    -o "$out" \
    --threads "$THREADS" --max-target-seqs 5000 --evalue 1e-3 --quiet \
    >>"$RAW/diamond_relaxed.log" 2>&1
  log diamond_relaxed "done ${sp} in $(( SECONDS - t )) s; lines=$(wc -l <"$out")"
}

# --------------- MMseqs2 relaxed ---------------
run_mmseqs_relaxed_cell() {
  local sp=$1
  local db="$DBROOT/mmseqs2/${sp}/genomedb"
  local proteins="$REPO/data/v3/${sp}_paralog/proteins.fa"
  local out="$RAW/mmseqs_relaxed.${sp}.paralog.tsv"
  local tmp="$RAW/mmseqs_relaxed.${sp}.paralog.tmp"
  if [ -s "$out" ]; then
    log mmseqs_relaxed "skip ${sp} (exists, $(wc -l <"$out") lines)"
    return
  fi
  if [ ! -f "$db" ]; then
    log mmseqs_relaxed "skip ${sp} (no DB at $db)"
    return
  fi
  rm -rf "$tmp"; mkdir -p "$tmp"
  local t=$SECONDS
  log mmseqs_relaxed "begin ${sp}"
  "$MMSEQS_BIN" easy-search \
    "$proteins" "$db" "$out" "$tmp" \
    -s 1.0 --search-type 2 --threads "$THREADS" \
    --max-seqs 5000 --format-mode 0 \
    >>"$RAW/mmseqs_relaxed.log" 2>&1
  rm -rf "$tmp"
  log mmseqs_relaxed "done ${sp} in $(( SECONDS - t )) s; lines=$(wc -l <"$out")"
}

# --------------- LAST relaxed ---------------
run_last_relaxed_cell() {
  local sp=$1
  local db="$DBROOT/last_protdb/${sp}.paralog"
  local genome="$GENOMES/${sp}.fa"
  local out="$RAW/last_relaxed.${sp}.paralog.tab"
  if [ -s "$out" ]; then
    log last_relaxed "skip ${sp} (exists, $(wc -l <"$out") lines)"
    return
  fi
  if [ ! -f "${db}.prj" ]; then
    log last_relaxed "skip ${sp} (no DB at $db)"
    return
  fi
  local t=$SECONDS
  log last_relaxed "begin ${sp}"
  "$LASTAL_BIN" -F20 -P "$THREADS" -E1.0 -fBlastTab \
    "$db" "$genome" \
    >"$out" 2>>"$RAW/last_relaxed.log"
  log last_relaxed "done ${sp} in $(( SECONDS - t )) s; lines=$(wc -l <"$out")"
}

# --------------- Driver ---------------
run_tool_all_paralog() {
  local tool_func=$1
  for sp in "${SPECIES[@]}"; do
    "$tool_func" "$sp"
  done
  log "${tool_func#run_}" "ALL DONE"
}

case "${1:-help}" in
  miniprot-ultra)    run_tool_all_paralog run_miniprot_ultra_cell ;;
  miniprot-relaxed)  run_tool_all_paralog run_miniprot_relaxed_cell ;;
  diamond)           run_tool_all_paralog run_diamond_relaxed_cell ;;
  mmseqs)            run_tool_all_paralog run_mmseqs_relaxed_cell ;;
  last)              run_tool_all_paralog run_last_relaxed_cell ;;
  all)
    run_tool_all_paralog run_miniprot_ultra_cell
    run_tool_all_paralog run_diamond_relaxed_cell
    run_tool_all_paralog run_mmseqs_relaxed_cell
    run_tool_all_paralog run_last_relaxed_cell
    ;;
  help|*) echo "usage: $0 {miniprot-ultra|miniprot-relaxed|diamond|mmseqs|last|all}"; exit 1 ;;
esac
