import pandas as pd
from collections import defaultdict

filename = "/projects/humgen/science/depienne/project818/dna-seq-nanopore/test.bed"

df = pd.read_csv(filename, delimiter="\t", na_values=["-"], low_memory=False)

allele1 = defaultdict(list)
allele2 = defaultdict(list)

for _, row in df.iterrows():
    key =  tuple(row[["#chrom", "start", "end", "repeat_unit"]])
    allele1[key].append(row["allele1:size"])
    allele2[key].append(row["allele2:size"])

print(allele1)
print(allele2)