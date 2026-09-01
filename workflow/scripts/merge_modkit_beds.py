from collections import defaultdict

inputs = [
    snakemake.input.bed_1,
    snakemake.input.bed_2,
    snakemake.input.bed_ungrouped,
]

agg = {}

for path in inputs:
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end, mod, _score, strand, s2, e2, color = fields[:9]
            # fields[9] is space-separated: nvalid pct nmod ncan noth ndel nfail ndiff nnocall
            parts = fields[9].split()
            nvalid  = int(parts[0])
            nmod    = int(parts[2])
            ncan    = int(parts[3])
            noth    = int(parts[4])
            ndel    = int(parts[5])
            nfail   = int(parts[6])
            ndiff   = int(parts[7])
            nnocall = int(parts[8])

            key = (chrom, int(start), int(end), mod)
            if key not in agg:
                agg[key] = {
                    "strand": strand, "s2": s2, "e2": e2, "color": color,
                    "nvalid": 0, "nmod": 0, "ncan": 0, "noth": 0,
                    "ndel": 0, "nfail": 0, "ndiff": 0, "nnocall": 0,
                }
            a = agg[key]
            a["nvalid"]  += nvalid
            a["nmod"]    += nmod
            a["ncan"]    += ncan
            a["noth"]    += noth
            a["ndel"]    += ndel
            a["nfail"]   += nfail
            a["ndiff"]   += ndiff
            a["nnocall"] += nnocall

with open(snakemake.output.bed, "w") as out:
    for (chrom, start, end, mod) in sorted(agg.keys()):
        a = agg[(chrom, start, end, mod)]
        pct = (a["nmod"] / a["nvalid"] * 100) if a["nvalid"] > 0 else 0.0
        col10plus = f'{a["nvalid"]} {pct:.2f} {a["nmod"]} {a["ncan"]} {a["noth"]} {a["ndel"]} {a["nfail"]} {a["ndiff"]} {a["nnocall"]}'
        out.write(
            f'{chrom}\t{start}\t{end}\t{mod}\t{a["nvalid"]}\t{a["strand"]}\t'
            f'{a["s2"]}\t{a["e2"]}\t{a["color"]}\t{col10plus}\n'
        )