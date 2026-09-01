"""Pair up the two strand rows modkit emits per CpG in a bedMethyl file."""

import sys

filename = sys.argv[1] if len(sys.argv) > 1 else snakemake.input.bed

last_split = None
for i, curr in enumerate(open(filename)):
    curr_split = curr.strip().split("\t")
    if i % 2 == 0:
        last_split = curr_split
    else:
        # the two rows of a pair must describe the same interval
        assert curr_split[0] == last_split[0], (last_split, curr_split)
        assert curr_split[1] == last_split[1], (last_split, curr_split)
        assert curr_split[2] == last_split[2], (last_split, curr_split)
        last_values = last_split[-1].split(" ")
        curr_values = curr_split[-1].split(" ")
        print("last", last_values)
        print("curr", curr_values)
