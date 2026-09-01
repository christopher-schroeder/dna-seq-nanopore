#!/usr/bin/env python
"""Check whether the input carries base modification tags."""

import os
import sys

import pysam

MOD_TAGS = {"mm", "ml"}
MAX_READS = 10000

valid_reads = 0
with pysam.AlignmentFile(
    snakemake.input.xam, reference_filename=snakemake.input.ref
) as xam:
    for i, alignment in enumerate(xam):
        tags = {tag.lower() for tag, _ in alignment.get_tags()}
        if MOD_TAGS <= tags:
            valid_reads += 1
            break
        if i >= MAX_READS - 1:
            break

if valid_reads == 0:
    sys.stderr.write(
        f"{snakemake.input.xam} carries no MM/ML tags in its first "
        f"{MAX_READS} reads: it is not a modified-base alignment.\n"
    )
    sys.exit(os.EX_DATAERR)

open(snakemake.output.check, mode="w").close()
