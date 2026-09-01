import sys
from collections import defaultdict

import pandas as pd

filename = sys.argv[1]

df = pd.read_csv(filename, delimiter="\t", na_values=["-"], low_memory=False)

allele1 = defaultdict(list)
allele2 = defaultdict(list)

for _, row in df.iterrows():
    key =  tuple(row[["#chrom", "start", "end", "repeat_unit"]])
    allele1[key].append(row["allele1:size"])
    allele2[key].append(row["allele2:size"])

print(allele1)
print(allele2)