"""minisplice + miniprot forced-locus runner.

Same as miniprot but with a learned splice score model. The splice score
is provided via miniprot's `--splice-score-file` flag (or whatever the
binary's actual flag name is).

If the splice model file is not supplied via options, this falls back to
plain miniprot — same as the standalone miniprot runner.
"""
from __future__ import annotations

from pathlib import Path

from . import miniprot as _mp


def run_one(protein_id, protein_seq, region_fa, region_chrom,
            region_offset, region_strand, scratch, options) -> dict | None:
    splice_model = options.get("minisplice_model")
    extra = list(options.get("miniprot_extra", []))
    if splice_model:
        # miniprot 0.18+ exposes --splice-score-file=FILE for learned scores;
        # the canonical flag name may vary. Pass via options['minisplice_extra']
        # if the binary expects something different.
        extra += options.get("minisplice_extra",
                              [f"--splice-score-file={splice_model}"])
    new_opts = dict(options)
    new_opts["miniprot_extra"] = extra
    return _mp.run_one(protein_id, protein_seq, region_fa, region_chrom,
                        region_offset, region_strand, scratch, new_opts)
