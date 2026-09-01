import argparse
import pysam
from collections import defaultdict

# controls = pd.read_csv("config/vcfs.tsv", sep="\t")
controls = snakemake.params.controls

with pysam.VariantFile(snakemake.input.variants) as f:
    header = f.header

    header.add_meta("INFO", items=[('ID',"CTRL_AF"), ('Number',1), ('Type','Float'), ('Description', 'The allele frequency from control samples.')])
    header.add_meta("INFO", items=[('ID',"CTRL_GT_HOM_WT"), ('Number',1), ('Type','Float'), ('Description', 'Number of controls with unknown variant genotype.')])
    header.add_meta("INFO", items=[('ID',"CTRL_GT_HET"), ('Number',1), ('Type','Float'), ('Description', 'Number of controls with heterozygous genotype.')])
    header.add_meta("INFO", items=[('ID',"CTRL_GT_HOM"), ('Number',1), ('Type','Float'), ('Description', 'Number of controls with homozygous genotype.')])

    # header.add_meta("INFO", items=[('ID',"CTRL_AC_AFR"), ('Number',1), ('Type','Float'), ('Description', 'The allele .')])
    with pysam.VariantFile(snakemake.output.variants, "w", header=f.header) as o:
        samples = list(f.header.samples)
        # # this seems to be too much work, but it ensures the correct order of the list
        # x = [(s, ethnicity[suffix]) for s in samples if (suffix:=s.split("_", 1)[1]) in ethnicity]
        # control_samples, ethnicity = zip(*x)

        control_samples = [s for s in samples if s.endswith("_SAMPLE")]

        for record in f:
            population = 0
            allele = 0
            gt_hom_wt = 0
            gt_het = 0
            gt_hom = 0

            for s in control_samples:
                gt = record.samples[s]["GT"]
                
                if gt == (None, None):
                    gt_hom_wt += 1
                    continue

                population += 1

                if gt[0] == 1:
                    allele += 1
                if gt[1] == 1:
                    allele += 1

                if gt == (1,0) or gt == (0, 1):
                    gt_het += 1
                    continue

                if gt == (1,1):
                    gt_hom += 1

            ref_freq = ((2 * gt_hom_wt) + gt_het) / (2 * (gt_hom_wt + gt_het + gt_hom))
            alt_freq = 1 - ref_freq
            record.info["CTRL_AF"] = alt_freq
            record.info["CTRL_GT_HOM_WT"] = gt_hom_wt
            record.info["CTRL_GT_HET"] = gt_het
            record.info["CTRL_GT_HOM"] = gt_hom
            o.write(record)
            # print(record.chrom, record.pos, alt_freq)
# for record in f:
    # print(record)