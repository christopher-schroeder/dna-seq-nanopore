rule align_and_qsFilter:
    threads: 64
    resources:
        mem_mb =160000
    input:
        reads="results/basecalls_sample/{sample}.ubam",
        reference=REFERENCE,
        mmi_reference=f"{REFERENCE}.mmi"
    output:
        passed="results/alignment/{sample}.cram",
        passed_index="results/alignment/{sample}.cram.crai",
        failed="results/alignment_failed/{sample}.cram",
    params:
       qscore_filter = 10
    conda:
        "../envs/minimap2.yaml"
    group:
        lambda wc: wc.sample
    shell:
        """
        samtools bam2fq -@ {threads} -T 1 {input.reads} \
        | minimap2 -y -t {threads} -ax map-ont {input.mmi_reference} - \
        | samtools sort -@ {threads} -m 2G \
        | samtools view -@ {threads} -e '[qs] >= {params.qscore_filter}' --output {output.passed} --unoutput {output.failed} -O cram,embed_ref --reference {input.reference} --write-index -
        """