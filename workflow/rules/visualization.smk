rule table_sv:
    input:
        "results/sv/{group}.annotated.repeats.bcf"
    output:
        "results/tables/{group}.sv.tsv"
    params:
        expression=lambda wc: f"CHROM, POS, INFO['SVTYPE'], INFO['SVLEN'], QUAL, CSQ['SYMBOL'], CSQ['Consequence'], CSQ['IMPACT'], CSQ['Feature'], INFO.get('gnomadAF', ''), INFO.get('gnomadAC', ''), INFO.get('gnomadnhomalt', ''), INFO.get('gnomadAF_XX', ''), INFO.get('gnomadAC_XX', ''), INFO.get('gnomadnhomalt_XX', ''), INFO.get('gnomadAF_XY', ''), INFO.get('gnomadAC_XY', ''), INFO.get('gnomadnhomalt_XY', ''), INFO.get('SR_LOCATION', ''), INFO.get('SR_PERIOD'), INFO.get('SR_COPYNUMBER'), INFO.get('SR_CONSENSUS_SIZE'), INFO.get('SR_PER_MATCH'), INFO.get('SR_SEQUENCE'), " + ", ".join(f"FORMAT['GT']['{sample}'], FORMAT['DV']['{sample}'], FORMAT['DR']['{sample}'] + FORMAT['DV']['{sample}']" for sample in get_group_samples(wc.group))
    conda:
        "../envs/vembrane.yaml"
    resources:
        mem_mb=256
    shell:
        """vembrane table --overwrite-number-format GT=2 --annotation-key CSQ "{params.expression}" {input} > {output}"""


rule table_str:
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