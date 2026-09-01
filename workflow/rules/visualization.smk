rule snps_table:
    input:
        "results/snps/{group}.annotated.inhouse.bcf"
    output:
        "results/tables/{group}.snps.tsv"
    params:
        expression=lambda wc: "CHROM, POS, REF, ALT, QUAL, CSQ['SYMBOL'], CSQ['Consequence'], CSQ['IMPACT'], CSQ['Feature'], INFO['inhouse_AF'], INFO['inhouse_AC'], INFO['inhouse_nhomalt'], INFO['gnomad_AF'], INFO['gnomad_AC'], INFO['gnomad_nhomalt'], INFO['gnomad_AF_XX'], INFO['gnomad_AC_XX'], INFO['gnomad_nhomalt_XX']," + ", ".join(f"FORMAT['GT']['{s}'], FORMAT['AD']['{s}'][0], FORMAT['DP']['{s}']" for s in get_group_samples(wc.group))
    conda:
        "../envs/vembrane.yaml"
    benchmark:
        "benchmarks/table/{group}.snps.txt"
    resources:
        mem_mb=256
    group:
        "group"
    shell:
        """vembrane table --overwrite-number-format GT=2 --overwrite-number-format AD=2 --annotation-key CSQ "{params.expression}" {input} > {output}"""


rule snps_table_maf:
    input:
        "results/snps/{group}.annotated.inhouse.maf.{maf}.bcf"
    output:
        "results/tables/{group}.snps.maf.{maf}.tsv"
    params:
        expression=lambda wc: "CHROM, POS, REF, ALT, QUAL, CSQ['SYMBOL'], CSQ['Consequence'], CSQ['IMPACT'], CSQ['Feature'], INFO['inhouse_AF'], INFO['inhouse_AC'], INFO['inhouse_nhomalt'], INFO['gnomad_AF'], INFO['gnomad_AC'], INFO['gnomad_nhomalt'], INFO['gnomad_AF_XX'], INFO['gnomad_AC_XX'], INFO['gnomad_nhomalt_XX']," + ", ".join(f"FORMAT['GT']['{s}'], FORMAT['AD']['{s}'][0], FORMAT['DP']['{s}']" for s in get_group_samples(wc.group))
    conda:
        "../envs/vembrane.yaml"
    benchmark:
        "benchmarks/table_snps_maf/{group}.{maf}.txt"
    resources:
        mem_mb=256
    group:
        "group"
    shell:
        """vembrane table --overwrite-number-format GT=2 --overwrite-number-format AD=2 --annotation-key CSQ "{params.expression}" {input} > {output}"""


rule sv_table:
    input:
        "results/sv/{group}.annotated.repeats.bcf"
    output:
        "results/tables/{group}.sv.tsv"
    params:
        expression=lambda wc: f"CHROM, POS, INFO['SVTYPE'], INFO['SVLEN'], QUAL, CSQ['SYMBOL'], CSQ['Consequence'], CSQ['IMPACT'], CSQ['Feature'], INFO.get('gnomad_AF', ''), INFO.get('gnomad_AC', ''), INFO.get('gnomad_nhomalt', ''), INFO.get('CTRL_AF', ''), INFO.get('CTRL_GT_HOM_WT', ''), INFO.get('CTRL_GT_HET', ''), INFO.get('CTRL_GT_HOM', ''), INFO.get('CTRL_GT_MISSING', ''), INFO.get('SR_LOCATION', ''), INFO.get('SR_PERIOD'), INFO.get('SR_COPYNUMBER'), INFO.get('SR_CONSENSUS_SIZE'), INFO.get('SR_PER_MATCH'), INFO.get('SR_SEQUENCE'), " + ", ".join(f"FORMAT['GT']['{sample}'], FORMAT['DV']['{sample}'], FORMAT['DR']['{sample}'] + FORMAT['DV']['{sample}']" for sample in get_group_samples(wc.group))
    conda:
        "../envs/vembrane.yaml"
    resources:
        mem_mb=256
    shell:
        """vembrane table --overwrite-number-format GT=2 --annotation-key CSQ "{params.expression}" {input} > {output}"""



rule str_table:
    input:
        "results/str/{sample}.annotated.vep.bcf"
    output:
        "results/tables/{sample}.str.tsv"
    params:
        expression="CHROM, POS, ALT, CSQ['SYMBOL'], CSQ['Consequence'], CSQ['IMPACT'], CSQ['Feature'], FORMAT['AS1'][SAMPLES[0]], FORMAT['ACN1'][SAMPLES[0]], FORMAT['ASP1'][SAMPLES[0]], FORMAT['AS2'][SAMPLES[0]], FORMAT['ACN2'][SAMPLES[0]], FORMAT['ASP2'][SAMPLES[0]]"
    conda:
        "../envs/vembrane.yaml"
    resources:
        mem_mb=256
    shell:
        """vembrane table --overwrite-number-format GT=2 --annotation-key CSQ "{params.expression}" {input} > {output}"""
