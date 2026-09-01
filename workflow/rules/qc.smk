## ---------------------------------------------------------------------------
## Read-level nanopore QC
## ---------------------------------------------------------------------------

rule nanoplot:
    # NanoStats.txt is what MultiQC's "nanostat" module parses; the plots and the
    # HTML report are a by-product that is useful on its own.
    threads:
        4
    input:
        xam="results/alignment/{sample}.cram",
        xai="results/alignment/{sample}.cram.crai",
        reference=REFERENCE,
    output:
        stats="results/qc/nanoplot/{sample}/{sample}.NanoStats.txt",
        report="results/qc/nanoplot/{sample}/{sample}.NanoPlot-report.html",
    params:
        outdir="results/qc/nanoplot/{sample}",
        prefix="{sample}.",
    log:
        "logs/nanoplot/{sample}.log"
    benchmark:
        "benchmarks/nanoplot/{sample}.txt"
    conda:
        "../envs/nanoplot.yaml"
    resources:
        mem_mb=64000,
    shell:
        """
        export REF_PATH={input.reference}
        # --tsv_stats emits the "Metrics<TAB>dataset" flavour of NanoStats.txt,
        # which is the format MultiQC expects; --no_static skips the matplotlib
        # renderings, which are slow and fragile at WGS read counts.
        (NanoPlot \
            --cram {input.xam} \
            --threads {threads} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --tsv_stats \
            --N50 \
            --no_static) > {log} 2>&1
        """


## ---------------------------------------------------------------------------
## Alignment QC
## ---------------------------------------------------------------------------

rule samtools_stats:
    input:
        bam= "results/alignment/{sample}.cram",
        ref= REFERENCE,
    output:
        "results/qc/samtools_stats/{sample}.txt",
    params:
        extra="",  # Optional: extra arguments.
    log:
        "logs/samtools_stats/{sample}.log",
    benchmark:
        "benchmarks/samtools_stats/{sample}.log",
    resources:
        mem_mb=1024,
    # group:
    #     "group"
    wrapper:
        "v2.0.0/bio/samtools/stats"


# flagstat and idxstats read flags and the index only, never sequence, so unlike
# `samtools stats` they neither need nor accept a --reference for CRAM input.
rule samtools_flagstat:
    threads:
        4
    input:
        xam="results/alignment/{sample}.cram",
        xai="results/alignment/{sample}.cram.crai",
    output:
        "results/qc/samtools_flagstat/{sample}.flagstat",
    log:
        "logs/samtools_flagstat/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=1024,
    shell:
        "(samtools flagstat -@ {threads} {input.xam} > {output}) 2> {log}"


rule samtools_idxstats:
    threads:
        1
    input:
        xam="results/alignment/{sample}.cram",
        xai="results/alignment/{sample}.cram.crai",
    output:
        "results/qc/samtools_idxstats/{sample}.idxstats",
    log:
        "logs/samtools_idxstats/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=1024,
    shell:
        "(samtools idxstats {input.xam} > {output}) 2> {log}"


rule qualimap:
    input:
        bam="results/alignment/{sample}.cram",
    output:
        outdir=directory("results/qc/qualimap/{sample}"),
    log:
        "logs/qualimap/bamqc/{sample}.log"
    benchmark:
        "benchmarks/qualimap/bamqc/{sample}.txt"
    conda:
        "../envs/qualimap.yaml"
    threads:
        33
    resources:
        mem_mb=90000
    shell:
        """
        samtools view -b {input.bam} | \
        qualimap bamqc \
        -bam /dev/stdin \
        -outdir {output.outdir} \
        -nt {threads} \
        --java-mem-size=80G
        """


## ---------------------------------------------------------------------------
## Pedigree / sex concordance (peddy)
## ---------------------------------------------------------------------------

# PED sex codes; anything unrecognised (or missing) becomes 0 = unknown, which
# peddy still uses to *infer* sex, it just cannot flag a mismatch.
PED_SEX_CODES = {
    "1": "1", "m": "1", "male": "1",
    "2": "2", "f": "2", "female": "2",
}


# samples.tsv is written with the short column names; the PED spellings are
# accepted as aliases so either sheet layout works.
PED_COLUMN_ALIASES = {
    "paternal_id": ("paternal_id", "father"),
    "maternal_id": ("maternal_id", "mother"),
    "sex": ("sex",),
    "phenotype": ("phenotype",),
}


def get_ped_field(sample, column, default="0"):
    """Read an optional pedigree column out of config/samples.tsv."""
    for candidate in PED_COLUMN_ALIASES.get(column, (column,)):
        if candidate not in samples.columns:
            continue
        value = samples.at[sample, candidate]
        if pd.isna(value) or not str(value).strip():
            continue
        return str(value).strip()
    return default


rule peddy_ped:
    # samples.tsv only has to carry sample_name and group; the optional columns
    # father/paternal_id, mother/maternal_id, sex and phenotype are used when
    # present.
    input:
        "config/samples.tsv",
    output:
        ped="results/qc/peddy/{group}.ped",
    run:
        with open(output.ped, "w") as out:
            print("#family_id", "sample_id", "paternal_id", "maternal_id", "sex", "phenotype", sep="\t", file=out)
            for sample in get_group_samples(wildcards.group):
                sex = get_ped_field(sample, "sex").lower()
                print(
                    wildcards.group,
                    sample,
                    get_ped_field(sample, "paternal_id"),
                    get_ped_field(sample, "maternal_id"),
                    PED_SEX_CODES.get(sex, "0"),
                    get_ped_field(sample, "phenotype", "-9"),
                    sep="\t",
                    file=out,
                )


rule peddy_vcf:
    # Deliberately re-merged from the per-sample calls with -0 instead of reusing
    # results/snps/{group}.bcf: plain `bcftools merge` leaves a sample as ./. at
    # every site the *other* samples carry, and peddy would read those as missing
    # rather than hom-ref, wrecking the relatedness estimates.
    threads:
        4
    input:
        vcfs=lambda wc: expand("results/snps_sample/{sample}.vcf.gz", sample=get_group_samples(wc.group)),
        indices=lambda wc: expand("results/snps_sample/{sample}.vcf.gz.tbi", sample=get_group_samples(wc.group)),
    output:
        vcf="results/qc/peddy/{group}.snps.vcf.gz",
        tbi="results/qc/peddy/{group}.snps.vcf.gz.tbi",
    params:
        force_single=lambda w, input: "--force-single" if len(input.vcfs) == 1 else "",
    log:
        "logs/peddy/{group}.merge.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        (bcftools merge -0 -m none {params.force_single} --threads {threads} -O u {input.vcfs} \
            | bcftools view -v snps -m 2 -M 2 --threads {threads} -O z -o {output.vcf} --write-index=tbi) 2> {log}
        """


rule peddy:
    threads:
        8
    input:
        vcf="results/qc/peddy/{group}.snps.vcf.gz",
        tbi="results/qc/peddy/{group}.snps.vcf.gz.tbi",
        ped="results/qc/peddy/{group}.ped",
    output:
        html="results/qc/peddy/{group}.html",
        ped_check="results/qc/peddy/{group}.ped_check.csv",
        sex_check="results/qc/peddy/{group}.sex_check.csv",
        het_check="results/qc/peddy/{group}.het_check.csv",
        background_pca="results/qc/peddy/{group}.background_pca.json",
        peddy_ped="results/qc/peddy/{group}.peddy.ped",
    params:
        prefix="results/qc/peddy/{group}",
        # peddy ships GRCh37 sites by default; "hg38" selects its GRCH38.sites,
        # which -- like this pipeline's reference -- is not chr-prefixed.
        sites="hg38",
    log:
        "logs/peddy/{group}.log"
    benchmark:
        "benchmarks/peddy/{group}.txt"
    conda:
        "../envs/peddy.yaml"
    resources:
        mem_mb=16000,
    shell:
        """
        (peddy --plot \
            --procs {threads} \
            --sites {params.sites} \
            --prefix {params.prefix} \
            {input.vcf} \
            {input.ped}) > {log} 2>&1
        """


## ---------------------------------------------------------------------------
## Aggregation
## ---------------------------------------------------------------------------

def multiqc_input(wildcards):
    paths = [
        *expand("results/qc/nanoplot/{sample}/{sample}.NanoStats.txt", sample=sample_names),
        *expand("results/qc/samtools_stats/{sample}.txt", sample=sample_names),
        *expand("results/qc/samtools_flagstat/{sample}.flagstat", sample=sample_names),
        *expand("results/qc/samtools_idxstats/{sample}.idxstats", sample=sample_names),
        *expand("results/qc/qualimap/{sample}", sample=sample_names),
        *expand("results/mosdepth/{sample}.mosdepth.global.dist.txt", sample=sample_names),
        *expand("results/mosdepth/{sample}.mosdepth.summary.txt", sample=sample_names),
    ]
    if config.get("peddy", True):
        paths += expand("results/qc/peddy/{group}.peddy.ped", group=groups)
        paths += expand("results/qc/peddy/{group}.sex_check.csv", group=groups)
        paths += expand("results/qc/peddy/{group}.het_check.csv", group=groups)
        paths += expand("results/qc/peddy/{group}.ped_check.csv", group=groups)
        paths += expand("results/qc/peddy/{group}.background_pca.json", group=groups)
    return paths


rule multiqc:
    input:
        analysis=multiqc_input,
        config=f"{BASEDIR}/data/multiqc_config.yaml",
    output:
        report="results/qc/multiqc.html",
        data=directory("results/qc/multiqc_data"),
    log:
        "logs/multiqc.log"
    benchmark:
        "benchmarks/multiqc.txt"
    conda:
        "../envs/qc.yaml"
    resources:
        mem_mb=10000,
    shell:
        """
        (multiqc \
            --force \
            --config {input.config} \
            --outdir results/qc \
            --filename multiqc.html \
            --data-dir \
            {input.analysis}) > {log} 2>&1
        """
