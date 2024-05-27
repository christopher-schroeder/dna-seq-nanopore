#!/usr/bin/env python
"""Check whether the input is a modbam."""

import pysam
import os
import sys

valid_reads = 0
fields = ['mm', 'ml']
for i, alignment in enumerate(pysam.AlignmentFile(snakemake.input.xam)):
    print(alignment)
    n_tags = len([
        tag for (tag, val) in alignment.get_tags() if tag.lower() in fields])
    if n_tags == 2:
        valid_reads += 1
        break
    if i >= 9999:
        break

if valid_reads == 0:
    sys.exit(os.EX_DATAERR)

open(snakemake.output.check, mode='w').close()