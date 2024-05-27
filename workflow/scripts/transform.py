import pysam

#reference = pysam.FastaFile(snakemake.input.reference)

with pysam.VariantFile(snakemake.input.calls, "r") as f:
    f.header.add_meta('INFO', items=[('ID',"SEQ"), ('Number',1), ('Type','String'), ('Description','Original Sequence')])
    with pysam.VariantFile(snakemake.output.calls, "w", header=f.header) as o:
        for record in f:
            svtype = record.info["SVTYPE"]
            if svtype == "INS":
                record.info["SEQ"] = record.alts[0]
                record.alts = ("<INS>",)
            elif svtype == "DEL":
                record.info["SEQ"] = record.ref
                record.alts = ("<DEL>",)
            # elif svtype == "BND":
            #     record.ref = base
            #     record.alts = (alt.replace("N", base),)
            #     record.pos = record.pos - 1
            elif svtype == "INV":
                pass
            elif svtype == "DUP":
                pass

            o.write(record)