"""Per-tool parsers.

Every parser is a pure reader: it takes the tool's native output file and
yields either raw HSPs (HSP-only tools) or locus-level Candidates (tools
with native locus output). No tool code is touched.

Parser contract
---------------
- HSP parsers:       iter_hsps_*(path) -> Iterator[(query_id, HSP)]
- Native parsers:    iter_candidates_*(path, query_lengths) -> Iterator[(query_id, Candidate)]

The driver in `normalize.py` groups HSPs by query, calls aggregate_hsps,
then feeds everything through finalize_candidates.

All genomic coordinates are 0-based half-open in the emitted records.
Tool-specific coordinate conventions are converted in each parser.
"""

from __future__ import annotations

import re
from typing import Iterable, Iterator

from .common import HSP, Candidate, query_coverage


# ===========================================================================
# HSP-only tools (fed into aggregate_hsps downstream)
# ===========================================================================

def iter_hsps_blast_tab(path: str) -> Iterator[tuple[str, HSP]]:
    """NCBI tblastn / DIAMOND / any `-outfmt 6` producer.

    Expected columns (standard 12):
      qseqid sseqid pident length mismatch gapopen
      qstart qend sstart send evalue bitscore

    qstart/qend are 1-based inclusive on the protein; sstart/send may appear
    reversed when the hit is on the minus strand of the subject DNA.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qid, sid = f[0], f[1]
            qstart = int(f[6]) - 1   # → 0-based
            qend = int(f[7])         # inclusive → half-open
            sstart = int(f[8])
            send = int(f[9])
            bitscore = float(f[11])
            if sstart <= send:
                strand = "+"
                s_lo, s_hi = sstart - 1, send
            else:
                strand = "-"
                s_lo, s_hi = send - 1, sstart
            yield qid, HSP(
                chrom=sid, start=s_lo, end=s_hi, strand=strand,
                bitscore=bitscore, q_start=qstart, q_end=qend,
            )


def iter_hsps_last_tab(path: str) -> Iterator[tuple[str, HSP]]:
    """LAST `-f BlastTab` (protein-to-6-frame) output.

    Columns identical to BLAST tab12. Protein is in the query column; the
    translated DNA is in the subject column with 1-based coordinates that
    may be reversed on minus strand.
    """
    yield from iter_hsps_blast_tab(path)


def iter_hsps_blast_tab_swapped(path: str) -> Iterator[tuple[str, HSP]]:
    """BLAST tab12 where DB is the protein side (i.e., genome was the query).

    Used for DIAMOND blastx (genome.fa as query, proteins.dmnd as db) and for
    LAST runs where the protein database is queried with the genome.

    Field semantics with this orientation:
      f[0] qseqid = chromosome name      ← becomes HSP.chrom (and BED `score` is per protein)
      f[1] sseqid = protein id           ← becomes the yielded query_id
      f[6] qstart, f[7] qend = DNA coords (genome)
      f[8] sstart, f[9] send = protein coords
      f[11] bitscore
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            chrom = f[0]
            protein_id = f[1]
            genome_a = int(f[6])
            genome_b = int(f[7])
            prot_a = int(f[8])
            prot_b = int(f[9])
            bitscore = float(f[11])
            if genome_a <= genome_b:
                strand = "+"
                g_lo, g_hi = genome_a - 1, genome_b
            else:
                strand = "-"
                g_lo, g_hi = genome_b - 1, genome_a
            qstart = min(prot_a, prot_b) - 1
            qend = max(prot_a, prot_b)
            if qstart >= qend:
                continue
            yield protein_id, HSP(
                chrom=chrom, start=g_lo, end=g_hi, strand=strand,
                bitscore=bitscore, q_start=qstart, q_end=qend,
            )


def iter_hsps_blast_report(path: str) -> Iterator[tuple[str, HSP]]:
    """genBlastA-emitted `.blast.report` (one tabular row per tblastn HSP).

    Format (one header line, then per-query blocks separated by `>QID` lines):
      HSP_ID  Protein_length  HSP_Score  E-value  P-value  P_Idn  Chr
      Chr_start  Chr_end  Chr_length  chr_strand  query_start  query_end  query_length

    chr_strand is 1 or -1 (genBlastA convention). All coordinates 1-based
    inclusive. Empty `>QID` blocks (no HSPs) are tolerated.
    """
    cur_qid: str | None = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("HSP_ID"):
                continue
            if line.startswith(">"):
                cur_qid = line[1:].split()[0]
                continue
            f = line.split("\t")
            if len(f) < 14 or cur_qid is None:
                continue
            try:
                bitscore = float(f[2])
                chrom = f[6]
                cs = int(f[7])
                ce = int(f[8])
                strand_int = int(f[10])
                qs = int(f[11])
                qe = int(f[12])
            except ValueError:
                continue
            strand = "+" if strand_int >= 0 else "-"
            g_lo = min(cs, ce) - 1
            g_hi = max(cs, ce)
            qstart = min(qs, qe) - 1
            qend = max(qs, qe)
            if g_lo >= g_hi or qstart >= qend:
                continue
            yield cur_qid, HSP(
                chrom=chrom, start=g_lo, end=g_hi, strand=strand,
                bitscore=bitscore, q_start=qstart, q_end=qend,
            )


def iter_hsps_blat_psl(path: str) -> Iterator[tuple[str, HSP]]:
    """BLAT PSL output (protein-mode `blat -t=dnax -q=prot`).

    PSL columns (0-indexed):
      0  matches
      1  misMatches
      ...
      8  strand   — 2-char for prot→dnax, e.g. '++' or '+-'; subject strand is 2nd
      9  qName
      10 qSize
      11 qStart (0-based)
      12 qEnd   (half-open)
      13 tName
      14 tSize
      15 tStart (0-based)
      16 tEnd
      17 blockCount
      18 blockSizes
      19 qStarts
      20 tStarts

    BLAT has no bitscore; we use `matches − misMatches` as the rank signal,
    which aligns with BLAT's own `pslReps -nohead -minCover=0.5` scoring.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 21:
                continue
            # Skip header rows: f[0] should be the integer 'matches' count.
            if not f[0].lstrip().isdigit():
                continue
            matches = int(f[0])
            mismatches = int(f[1])
            strand_field = f[8]
            strand = strand_field[-1] if len(strand_field) >= 2 else strand_field
            qid = f[9]
            qstart = int(f[11])
            qend = int(f[12])
            tname = f[13]
            tstart = int(f[15])
            tend = int(f[16])
            yield qid, HSP(
                chrom=tname, start=tstart, end=tend, strand=strand,
                bitscore=float(matches - mismatches),
                q_start=qstart, q_end=qend,
            )


# ===========================================================================
# Tools with native locus-level output (pass-through)
# ===========================================================================

def iter_candidates_miniprot_paf(path: str,
                                 query_lengths: dict[str, int]
                                 ) -> Iterator[tuple[str, Candidate]]:
    """miniprot `--outn N` PAF output (one line = one alignment).

    Accepts both raw PAF and miniprot's GFF3 with embedded ``##PAF\\t...``
    comment lines (the format produced by ``miniprot --gff``). The first
    field after stripping ``##PAF\\t`` is the standard 12-column PAF.

    PAF columns (0-indexed):
       0 qname       (protein id)
       1 qlen        (residues)
       2 qstart      0-based
       3 qend        half-open
       4 strand      + / -
       5 tname       (chrom)
       6 tlen
       7 tstart      0-based
       8 tend
       9 mlen        (matching residues)
      10 blen        (block length)
      11 mapq
      ... tag fields: ms:i:<chain_score>, AS:i:<align_score>, ...

    Score = `ms:i:` (chaining score) if present, else `AS:i:`, else `mlen`.
    We do not re-rank; we preserve whatever miniprot itself used to order
    the candidates via a higher-is-better numeric key.
    """
    tag_re = re.compile(r"(ms|AS):i:(-?\d+)")
    for line in open(path):
        if not line.strip():
            continue
        if line.startswith("##PAF\t"):
            line = line[len("##PAF\t"):]
        elif line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        qid = f[0]
        qlen = int(f[1]) if f[1].isdigit() else query_lengths.get(qid, 0)
        qstart = int(f[2])
        qend = int(f[3])
        strand = f[4]
        chrom = f[5]
        tstart = int(f[7])
        tend = int(f[8])
        mlen = int(f[9])

        score = float(mlen)
        for tag, val in tag_re.findall("\t".join(f[12:])):
            if tag == "ms":
                score = float(val)
                break
            if tag == "AS":
                score = float(val)

        qlen_eff = qlen if qlen > 0 else max(1, qend - qstart)
        qcov = (qend - qstart) / qlen_eff
        yield qid, Candidate(
            chrom=chrom, start=tstart, end=tend, strand=strand,
            score=score, q_start=qstart, q_end=qend, q_cov=min(1.0, qcov),
        )


def iter_candidates_spaln_gff(path: str,
                              query_lengths: dict[str, int]
                              ) -> Iterator[tuple[str, Candidate]]:
    """Spaln3 GFF3 output; one ``mRNA`` feature per candidate.

    The attribute field carries either ``Target=<qid> <qstart> <qend>`` or
    just ``Target=<qid>``. When the per-query coords are absent on the
    mRNA record we recover them from the union of child ``cds`` records
    that share the mRNA's ``ID=`` via ``Parent=``. Top-level ``gene``
    records are skipped because Spaln3 always emits a redundant mRNA at
    the same coordinates and that would double-count candidates.
    """
    attr_target_full = re.compile(r"Target=([^ ;]+)\s+(\d+)\s+(\d+)")
    attr_target_id   = re.compile(r"Target=([^ ;\t]+)")
    attr_id          = re.compile(r"\bID=([^;]+)")
    attr_parent      = re.compile(r"\bParent=([^;]+)")

    pending: dict[str, dict] = {}
    cds_qbounds: dict[str, tuple[int, int]] = {}
    pending_order: list[str] = []

    out: list[tuple[str, Candidate]] = []
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        ftype = f[2]
        chrom = f[0]
        start = int(f[3]) - 1
        end = int(f[4])
        score_str = f[5]
        strand = f[6]
        attr = f[8]

        if ftype == "mRNA":
            m_full = attr_target_full.search(attr)
            if m_full:
                qid = m_full.group(1)
                qstart = int(m_full.group(2)) - 1
                qend = int(m_full.group(3))
                qlen = query_lengths.get(qid, max(1, qend - qstart))
                try:
                    score = float(score_str)
                except ValueError:
                    score = float(end - start)
                out.append((qid, Candidate(
                    chrom=chrom, start=start, end=end, strand=strand,
                    score=score, q_start=qstart, q_end=qend,
                    q_cov=min(1.0, (qend - qstart) / max(1, qlen)),
                )))
                continue
            m_id_only = attr_target_id.search(attr)
            if not m_id_only:
                continue
            qid = m_id_only.group(1)
            try:
                score = float(score_str)
            except ValueError:
                score = float(end - start)
            m_recid = attr_id.search(attr)
            if not m_recid:
                continue
            rec_id = m_recid.group(1)
            pending[rec_id] = {
                "qid": qid, "chrom": chrom, "start": start, "end": end,
                "strand": strand, "score": score,
            }
            pending_order.append(rec_id)
        elif ftype.lower() == "cds":
            m_p = attr_parent.search(attr)
            if not m_p:
                continue
            parents = m_p.group(1).split(",")
            m_full = attr_target_full.search(attr)
            if not m_full:
                continue
            qstart = int(m_full.group(2)) - 1
            qend = int(m_full.group(3))
            for parent_id in parents:
                lo_hi = cds_qbounds.get(parent_id)
                if lo_hi is None:
                    cds_qbounds[parent_id] = (qstart, qend)
                else:
                    cds_qbounds[parent_id] = (min(lo_hi[0], qstart),
                                              max(lo_hi[1], qend))

    for rec_id in pending_order:
        rec = pending[rec_id]
        qstart_lo, qend_hi = cds_qbounds.get(rec_id, (None, None))
        qid = rec["qid"]
        qlen = query_lengths.get(qid, 0)
        if qstart_lo is not None and qend_hi is not None:
            qstart, qend = qstart_lo, qend_hi
        elif qlen > 0:
            qstart, qend = 0, qlen
        else:
            qstart, qend = 0, 1
        out.append((qid, Candidate(
            chrom=rec["chrom"], start=rec["start"], end=rec["end"],
            strand=rec["strand"], score=rec["score"],
            q_start=qstart, q_end=qend,
            q_cov=min(1.0, (qend - qstart) / max(1, qlen if qlen > 0 else qend - qstart)),
        )))
    yield from out


def iter_candidates_gth_gff(path: str,
                            query_lengths: dict[str, int]
                            ) -> Iterator[tuple[str, Candidate]]:
    """GenomeThreader GFF3; same convention as Spaln for our purposes."""
    yield from iter_candidates_spaln_gff(path, query_lengths)


def iter_candidates_genblasta(path: str,
                              query_lengths: dict[str, int]
                              ) -> Iterator[tuple[str, Candidate]]:
    """genBlastA cluster output as embedded in genblastg work-dir files.

    The genblastg pipeline writes its phase-1 (genBlastA HSP-grouping) results
    inside a header block:
        //*****************START******************//
        //   candidate genes: (in ranked order)   //
        //for query: <QID>                        //
        //****************************************//
        <QID>|<chrom>:<start>..<end>|<strand>|gene cover:<N>(<P>%)|score:<N>|rank:<N>
        ...
        <blank line>
        //******************END*******************//

    We parse only the phase-1 block (between START header and END marker) so
    that genBlastG phase-2 refined coords don't leak in.

    Coordinates are 1-based inclusive; strand is `+` or `-`.
    """
    row_re = re.compile(
        r"^(?P<qid>\S+?)\|(?P<chrom>\S+?):"
        r"(?P<start>\d+)\.\.(?P<end>\d+)\|"
        r"(?P<strand>[+-])\|"
        r"gene\s*cover:(?P<cov>\d+)\((?P<cov_pct>\d+(?:\.\d+)?)%\)\|"
        r"score:(?P<score>-?\d+(?:\.\d+)?)\|"
        r"rank:(?P<rank>\d+)"
    )
    in_phase1 = False
    for line in open(path):
        if "candidate genes" in line:
            in_phase1 = True
            continue
        if "END***" in line or "phase 2" in line:
            in_phase1 = False
            continue
        if not in_phase1:
            continue
        m = row_re.search(line)
        if not m:
            continue
        qid = m.group("qid")
        chrom = m.group("chrom")
        strand = m.group("strand")
        start = int(m.group("start")) - 1
        end = int(m.group("end"))
        score = float(m.group("score"))
        cov_pct = float(m.group("cov_pct")) / 100.0
        qlen = query_lengths.get(qid, 0)
        qend = int(round(cov_pct * qlen)) if qlen else 0
        yield qid, Candidate(
            chrom=chrom, start=start, end=end, strand=strand,
            score=score, q_start=0, q_end=qend, q_cov=min(1.0, cov_pct),
        )


def iter_candidates_scipio_yaml(path: str,
                                query_lengths: dict[str, int]
                                ) -> Iterator[tuple[str, Candidate]]:
    """Scipio v1.4 YAML output (concatenated or single-protein).

    Scipio emits a top-level mapping ``ENSPxxxxxxxxxxxxx:`` whose value is a
    list. Each list item is a separate candidate locus and contributes one
    BED row. The chrom field comes from the ``target:`` line and is the
    first whitespace-separated token (strip any FASTA description).

    The parser is line-based to avoid a PyYAML dependency. State is flushed
    when:
      - a new ``- <key>: ...`` list item is encountered (sibling candidate),
      - a new ``ENSPxxxxxxxxxxxxx:`` top-level key starts (new query), or
      - a YAML document separator ``---`` is seen.
    """
    state: dict[str, object] = {
        "qid": None,
        "chrom": None,
        "strand": None,
        "dna_start": None,
        "dna_end": None,
        "prot_start": None,
        "prot_end": None,
        "matches": None,
        "plen": None,
        "in_candidate": False,
    }

    def consume_kv(key: str, val: str) -> None:
        if key == "target":
            tok = val.split()
            if tok:
                state["chrom"] = tok[0]
        elif key == "strand":
            state["strand"] = "+" if val in ("+", "forward", "1") else "-"
        elif key in ("dna_start", "dna_end", "prot_start", "prot_end", "prot_len"):
            try:
                state[key if key != "prot_len" else "plen"] = int(val)
            except ValueError:
                pass
        elif key == "matches":
            try:
                state["matches"] = float(val)
            except ValueError:
                pass

    def flush(out: list) -> None:
        if (state["in_candidate"]
                and state["qid"] and state["chrom"]
                and state["dna_start"] is not None
                and state["dna_end"] is not None):
            qid = state["qid"]
            chrom = state["chrom"]
            dna_start = state["dna_start"]
            dna_end = state["dna_end"]
            prot_start = state["prot_start"] or 0
            prot_end = state["prot_end"] or 0
            matches = state["matches"] or 0.0
            plen = state["plen"] or query_lengths.get(qid, max(1, prot_end - prot_start))
            out.append((qid, Candidate(
                chrom=chrom,
                start=min(dna_start, dna_end) - 1,
                end=max(dna_start, dna_end),
                strand=state["strand"] or "+",
                score=float(matches),
                q_start=prot_start,
                q_end=prot_end if prot_end > 0 else plen,
                q_cov=min(1.0, ((prot_end or plen) - prot_start) / max(1, plen)),
            )))
        state["chrom"] = None
        state["strand"] = None
        state["dna_start"] = None
        state["dna_end"] = None
        state["prot_start"] = None
        state["prot_end"] = None
        state["matches"] = None
        state["plen"] = None
        state["in_candidate"] = False

    out: list[tuple[str, Candidate]] = []

    for line in open(path):
        stripped = line.rstrip("\n")
        if not stripped or stripped.lstrip().startswith("#"):
            continue
        if stripped.startswith("---"):
            flush(out)
            state["qid"] = None
            continue
        if (not stripped.startswith(" ") and not stripped.startswith("\t")
                and stripped.endswith(":") and ":" in stripped):
            flush(out)
            raw_qid = stripped.rstrip(":").strip() or None
            # Scipio emits alternative loci as `<qid>_(<n>):`; strip the
            # parenthetical suffix so candidates feed back into the original
            # query bucket.
            if raw_qid is not None:
                m_alt = re.match(r"^(?P<qid>.+?)_\(\d+\)$", raw_qid)
                if m_alt:
                    raw_qid = m_alt.group("qid")
            state["qid"] = raw_qid
            continue
        s_lstrip = stripped.lstrip()
        if s_lstrip.startswith("- "):
            flush(out)
            state["in_candidate"] = True
            after = s_lstrip[2:]
            if ":" in after:
                key, _, val = after.partition(":")
                consume_kv(key.strip(), val.strip().strip('"\''))
            continue
        if not state["in_candidate"]:
            continue
        if ":" not in s_lstrip:
            continue
        key, _, val = s_lstrip.partition(":")
        consume_kv(key.strip(), val.strip().strip('"\''))
    flush(out)
    yield from out


def iter_candidates_genblastg(path: str,
                              query_lengths: dict[str, int]
                              ) -> Iterator[tuple[str, Candidate]]:
    """genBlastG concatenated stdout. The candidate-list block is emitted
    once per query and looks like:

        //*****************START******************//
        //   candidate genes: (in ranked order)   //
        //for query: ENSP00000378787              //
        //****************************************//
        ENSP00000378787|GL456221.1:78825..79442|-|gene cover:206(51.5%)|score:-4.62625|rank:1
        HSP_ID[6]:(78825-79442);query:(46-251); pid: 20.3883
        ...

    The pipe-delimited candidate line carries the chrom, span, strand,
    cover percent, score, and rank for each genBlastG locus. We parse
    those lines only — the HSP detail rows are ignored.

    Coordinates in the pipe line are 1-based inclusive (genBlastG's BLAST
    convention); converted to 0-based half-open here.
    """
    line_re = re.compile(
        r"^(?P<qid>\S+?)\|"
        r"(?P<chrom>[^|:]+?):(?P<start>\d+)\.\.(?P<end>\d+)\|"
        r"(?P<strand>[+-])\|"
        r"gene cover:(?P<cov_n>\d+)\((?P<cov_pct>[\d.]+)%\)\|"
        r"score:(?P<score>-?[\d.eE+-]+)\|"
        r"rank:(?P<rank>\d+)"
    )
    for line in open(path):
        m = line_re.search(line)
        if not m:
            continue
        qid = m.group("qid")
        chrom = m.group("chrom")
        start = int(m.group("start")) - 1   # 1-based inclusive -> 0-based half-open
        end = int(m.group("end"))
        strand = m.group("strand")
        cov_pct = float(m.group("cov_pct")) / 100.0
        try:
            score = float(m.group("score"))
        except ValueError:
            score = 0.0
        # genBlastG scores are negative-better; flip sign so higher = better.
        score_higher_better = -score
        qlen = query_lengths.get(qid, 0)
        qend = int(round(cov_pct * qlen)) if qlen > 0 else 0
        yield qid, Candidate(
            chrom=chrom, start=start, end=end, strand=strand,
            score=score_higher_better, q_start=0, q_end=qend,
            q_cov=min(1.0, cov_pct),
        )


def iter_candidates_procompart(path: str,
                               query_lengths: dict[str, int]
                               ) -> Iterator[tuple[str, Candidate]]:
    """ProCompart output.

    Each row after ProCompart is essentially a compartment record:
      query subject strand compartment_id q_start q_end s_start s_end score

    The exact column layout is taken from pegi3s/prosplign-procompart
    pipeline logs: tab-separated with the columns above.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            qid = f[0]
            chrom = f[1]
            strand = f[2] if f[2] in ("+", "-") else "+"
            qstart = int(f[4]) - 1
            qend = int(f[5])
            sstart = int(f[6])
            send = int(f[7])
            lo = min(sstart, send) - 1
            hi = max(sstart, send)
            score = float(f[8])
            qlen = query_lengths.get(qid, max(1, qend - qstart))
            yield qid, Candidate(
                chrom=chrom, start=lo, end=hi, strand=strand,
                score=score, q_start=qstart, q_end=qend,
                q_cov=min(1.0, (qend - qstart) / qlen),
            )
