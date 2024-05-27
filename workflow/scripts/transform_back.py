import argparse
import pysam


parser = argparse.ArgumentParser()
parser.add_argument("vcf")
args = parser.parse_args()

with pysam.VariantFile(args.vcf, "r") as f:
    with pysam.VariantFile("-", "w", header=f.header) as o:
        for record in f:
            if "OSR" in record.info:
                record.ref = record.info["OSR"]
                del record.info["OSR"]
            if "OSA" in record.info:
                record.alts = (record.info["OSA"], )
                del record.info["OSA"]
            o.write(record)