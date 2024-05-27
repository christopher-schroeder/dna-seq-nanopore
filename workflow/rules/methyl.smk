
# Check that the bam has modifications
rule validate_modbam:
    input:
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
        ref=REFERENCE,
    output:
        check="results/checks/validate_modbam/{sample}.txt",
    params:
        basedir=BASEDIR
    conda:
        "../envs/samtools.yaml"
    group:
        lambda wc: f"{wc.sample}"
    script:
        "{params.basedir}/scripts/check_valid_modbam.py"

rule modkit_phase:
    threads: 64
    input:
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
        ref=REFERENCE,
        check="results/checks/validate_modbam/{sample}.txt",
    output:
        bed_1="results/methylation/{sample}_1.bed",
        bed_2="results/methylation/{sample}_2.bed",
        bed_ungrouped="results/methylation/{sample}_ungrouped.bed",
    params:
        outdir="results/methylation/",
        modkit_args="--combine-strands --cpg",
        modkit=f"{BASEDIR}/tools/modkit/modkit",
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        {params.modkit} pileup \\
            {input.xam} \\
            {params.outdir} \\
            --ref {input.ref} \\
            --partition-tag HP \\
            --prefix {wildcards.sample} \\
            --threads {threads} {params.modkit_args}
        """

        # for i in `ls {params.outdir}/`; do
        #     bgzip {params.outdir}/*.bed
        # done