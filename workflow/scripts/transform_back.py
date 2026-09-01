import pysam

with pysam.VariantFile(snakemake.input.calls, "r") as f:
    with pysam.VariantFile(snakemake.output.calls, "w", header=f.header) as o:
        for record in f:
            if "SEQ" in record.info:
                record.alts = (record.info["SEQ"],)
                del record.info["SEQ"]
            o.write(record)