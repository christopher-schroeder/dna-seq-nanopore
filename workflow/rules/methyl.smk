
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

rule modkit_merge:
    threads: 1
    input:
        bed_1="results/methylation/{sample}_1.bed",
        bed_2="results/methylation/{sample}_2.bed",
        bed_ungrouped="results/methylation/{sample}_ungrouped.bed",
    output:
        bed="results/methylation/{sample}_merged.bed",
    group:
        lambda wc: f"{wc.sample}"
    resources:
        mem_mb=100000,
    script:
        "../scripts/merge_modkit_beds.py"


rule methbat_pileup:
    threads: 64
    input:
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
    output:
        bed="results/methbat/pileup/{sample}.5mC.bed.gz",
        tbi="results/methbat/pileup/{sample}.5mC.bed.gz.tbi",
    params:
        prefix="results/methbat/pileup/{sample}",
    conda:
        "../envs/methbat.yaml"
    shell:
        """
        methbat pileup \
            --input-bam {input.xam} \
            --output-prefix {params.prefix} \
            --threads {threads}
        """


rule methbat_profile:
    threads: 1
    input:
        bed="results/methbat/pileup/{sample}.5mC.bed.gz",
        tbi="results/methbat/pileup/{sample}.5mC.bed.gz.tbi",
        regions=config["methbat_regions"],
    output:
        profile="results/methbat/profiles/{sample}.profile.tsv",
    conda:
        "../envs/methbat.yaml"
    shell:
        """
        methbat profile \
            --input-pileup {input.bed} \
            --input-regions {input.regions} \
            --output-region-profile {output.profile}
        """


rule methbat_collection:
    input:
        profiles=expand("results/methbat/profiles/{sample}.profile.tsv", sample=units["sample_name"]),
    output:
        collection="results/methbat/collection.tsv",
    run:
        with open(output.collection, "w") as f:
            f.write("id\tfilename\n")
            for profile in input.profiles:
                sample = os.path.basename(profile).replace(".profile.tsv", "")
                f.write(f"{sample}\t{profile}\n")


rule methbat_build:
    threads: 1
    input:
        collection="results/methbat/collection.tsv",
    output:
        cohort="results/methbat/cohort.profile.tsv",
    conda:
        "../envs/methbat.yaml"
    shell:
        """
        methbat build \
            --input-collection {input.collection} \
            --output-profile {output.cohort}
        """


rule methbat_outliers:
    threads: 1
    input:
        bed="results/methbat/pileup/{sample}.5mC.bed.gz",
        tbi="results/methbat/pileup/{sample}.5mC.bed.gz.tbi",
        cohort="results/methbat/cohort.profile.tsv",
    output:
        outliers="results/methbat/{sample}.outliers.tsv",
    conda:
        "../envs/methbat.yaml"
    shell:
        """
        methbat profile \
            --input-pileup {input.bed} \
            --input-regions {input.cohort} \
            --output-region-profile {output.outliers}
        """
