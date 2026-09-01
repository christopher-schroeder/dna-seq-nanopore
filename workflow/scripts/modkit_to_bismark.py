import gzip

with open(snakemake.input.bed) as f, gzip.open(snakemake.output.cov, "wt") as out:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, mod = fields[0], fields[1], fields[2], fields[3]
        if mod != "m":  # only 5mC; skip 5hmC etc.
            continue
        parts = fields[9].split()
        nvalid = int(parts[0])
        nmod   = int(parts[2])
        ncan   = int(parts[3])
        if nvalid == 0:
            continue
        pct = nmod / nvalid * 100
        # bismark cov: chrom  start(1-based)  end  pct_meth  count_meth  count_unmeth
        out.write(f"{chrom}\t{int(start)+1}\t{end}\t{pct:.6f}\t{nmod}\t{ncan}\n")