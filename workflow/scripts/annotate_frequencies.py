"""Annotate merged SV calls with genotype counts from the control cohort.

Jasmine prefixes every input file's sample column with its index in the file
list, so the control columns are the ones named "<n>_SAMPLE"; the first file is
the cohort of interest and is left out of the counts.
"""

import pysam

controls = snakemake.params.controls

with pysam.VariantFile(snakemake.input.variants) as f:
    header = f.header

    header.add_meta("INFO", items=[('ID', "CTRL_AF"), ('Number', 1), ('Type', 'Float'), ('Description', 'Alternate allele frequency across the called control genotypes.')])
    header.add_meta("INFO", items=[('ID', "CTRL_GT_HOM_WT"), ('Number', 1), ('Type', 'Integer'), ('Description', 'Number of controls genotyped homozygous reference.')])
    header.add_meta("INFO", items=[('ID', "CTRL_GT_HET"), ('Number', 1), ('Type', 'Integer'), ('Description', 'Number of controls genotyped heterozygous.')])
    header.add_meta("INFO", items=[('ID', "CTRL_GT_HOM"), ('Number', 1), ('Type', 'Integer'), ('Description', 'Number of controls genotyped homozygous alternate.')])
    header.add_meta("INFO", items=[('ID', "CTRL_GT_MISSING"), ('Number', 1), ('Type', 'Integer'), ('Description', 'Number of controls without a genotype call at this site.')])
    header.add_meta("INFO", items=[('ID', "CTRL_AN"), ('Number', 1), ('Type', 'Integer'), ('Description', 'Number of called alleles across the controls (the CTRL_AF denominator).')])

    with pysam.VariantFile(snakemake.output.variants, "w", header=f.header) as o:
        samples = list(f.header.samples)

        control_samples = [s for s in samples if s.endswith("_SAMPLE")]
        if not control_samples:
            raise ValueError(
                "no control columns (<n>_SAMPLE) found in "
                f"{snakemake.input.variants}; got {samples}"
            )

        for record in f:
            alt_alleles = 0
            called_alleles = 0
            gt_hom_wt = 0
            gt_het = 0
            gt_hom = 0
            gt_missing = 0

            for s in control_samples:
                gt = record.samples[s]["GT"]

                # a genotype is missing if it is absent or every allele is "."
                if not gt or all(allele is None for allele in gt):
                    gt_missing += 1
                    continue

                # only alleles that were actually called carry information
                called = [allele for allele in gt if allele is not None]
                n_alt = sum(1 for allele in called if allele > 0)
                called_alleles += len(called)
                alt_alleles += n_alt

                if n_alt == 0:
                    gt_hom_wt += 1
                elif n_alt == len(called):
                    gt_hom += 1
                else:
                    gt_het += 1

            record.info["CTRL_AF"] = (alt_alleles / called_alleles) if called_alleles else 0.0
            record.info["CTRL_AN"] = called_alleles
            record.info["CTRL_GT_HOM_WT"] = gt_hom_wt
            record.info["CTRL_GT_HET"] = gt_het
            record.info["CTRL_GT_HOM"] = gt_hom
            record.info["CTRL_GT_MISSING"] = gt_missing
            o.write(record)
