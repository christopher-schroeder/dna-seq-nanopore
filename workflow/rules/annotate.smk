rule get_vep_cache:
    output:
        directory("results/resources/vep/cache"),
    params:
        species="homo_sapiens",
        build="GRCh38",
        release="111",
    log:
        "logs/vep/cache.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v4.3.0/bio/vep/cache"

rule download_vep_plugins:
    output:
        directory("results/resources/vep/plugins")
    params:
        release=100
    wrapper:
        "v3.3.6/bio/vep/plugins"


rule sort_str:
    input:
        bed="{x}/{sample}.bed",
    output:
        bed="{x}/{sample}.sorted.bed",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools sort -i {input} > {output}"


rule normalize:
    threads:
        8
    input:
        variants="results/{calls}/{group}.bcf",
        reference=REFERENCE,
    output:
        "results/{calls}/{group}.norm.bcf",
    log:
        "logs/normalize/{calls}/{group}.log"
    benchmark:
        "benchmarks/normalize/{calls}/{group}.txt"
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=2048
    group:
        lambda wc: wc.group
    shell:
        "(bcftools norm -a -m -any -f {input.reference} --atom-overlaps . --threads {threads} -Ob -c w {input.variants} > {output}) 2> {log}"


rule annotate_snps_vep:
    threads:
        16
    input:
        calls="results/snps/{group}.norm.bcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls="results/snps/{group}.annotated.vep.vcf.gz",
        stats="results/snps/{group}.annotated.vep.stats.html",
    params:
        plugins=[],
        extra="--symbol"
    log:
        "logs/annotate/{group}.vep.log",
    wrapper:
        "v3.3.6/bio/vep/annotate"



rule annotate_snps_gnomad:
    threads: 2
    input:
        calls="results/snps/{group}.annotated.vep.vcf.gz",
        calls_index="results/snps/{group}.annotated.vep.vcf.gz.tbi",
        database="resources/gnomad.genomes.v4.1.sites.vcf.gz",
        database_index="resources/gnomad.genomes.v4.1.sites.vcf.gz.tbi",
    output:
        call="results/snps/{group}.annotated.gnomad.bcf",
        call_index="results/snps/{group}.annotated.gnomad.bcf.csi",
    params:
        info="gnomad_AN:=INFO/AN,gnomad_AF:=INFO/AF,gnomad_AC:=INFO/AC,gnomad_nhomalt:=INFO/nhomalt,gnomad_AN_XX:=INFO/AN_XX,gnomad_AC_XX:=INFO/AC_XX,gnomad_AF_XX:=INFO/AF_XX,gnomad_nhomalt_XX:=INFO/nhomalt_XX"
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/annotate/{group}.gnomad.log"
    shell:
        "bcftools annotate -a {input.database} {input.calls} -c CHROM,POS,REF,ALT,{params.info} --threads {threads} -O b -o {output.call} --write-index"



rule sv_sort_vcf:
    threads:
        1
    input:
        vcf="results/sv/{group}.control_annotated.vcf",
    output:
        vcf="results/sv/{group}.sorted.vcf",
    conda:
        "../envs/filtercalls.yaml"
    group:
        lambda wc: wc.group
    shell:
        "bcftools sort {input.vcf} > {output.vcf}"


def sv_transform_input(wc):
    if config.get("jasmine", True):
        return f"results/sv/{wc.group}.jasmine_fix.vcf"
    else:
        return f"results/sv/{wc.group}.filtered.vcf"


rule sv_transform:
    input:
        calls=sv_transform_input,
        reference=REFERENCE,
    output:
        calls="results/sv/{group}.transformed.bcf",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/transform.py"


def get_vep_gff():
    if (vep:=config.get("vep", None)):
        return {"gff": vep}
    return {}


rule sv_annotate_vep:
    threads:
        4
    input:
        **get_vep_gff(),
        calls="results/sv/{group}.transformed.bcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
        fasta=REFERENCE,
    output:
        calls="results/sv/{group}.annotated.vep.bcf",
        stats="results/sv/{group}.annotated.vep.stats.html",
    params:
        plugins=[],
        extra="--symbol"
    log:
        "logs/annotate/sv/{group}.vep.log",
    wrapper:
        "v3.5.2/bio/vep/annotate"


rule sv_transform_back:
    input:
        calls="results/sv/{group}.annotated.vep.bcf",
    output:
        calls="results/sv/{group}.annotated.vep.back.bcf",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/transform_back.py"


rule sv_annotate_gnomad:
    threads: 2
    input:
        calls="results/sv/{group}.annotated.vep.back.bcf",
        calls_index="results/sv/{group}.annotated.vep.back.bcf.csi",
        database="resources/gnomad.v4.1.sv.sites.no_chr.vcf.gz",
        database_index="resources/gnomad.v4.1.sv.sites.no_chr.vcf.gz.tbi",
    output:
        call="results/sv/{group}.annotated.gnomad.bcf",
        call_index="results/sv/{group}.annotated.gnomad.bcf.csi",
    params:
        info="gnomad_AN:=INFO/AN,gnomad_AF:=INFO/AF,gnomad_AC:=INFO/AC,gnomad_nhomalt:=INFO/N_HOMALT,gnomad_AN_XX:=INFO/AN_XX,gnomad_AC_XX:=INFO/AC_XX,gnomad_AF_XX:=INFO/AF_XX,gnomad_nhomalt_XX:=INFO/N_HOMALT_XX"
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/annotate/sv/{group}.gnomad.log"
    shell:
        "bcftools annotate -a {input.database} {input.calls} -c CHROM,POS,REF,ALT,{params.info} --threads {threads} -O b -o {output.call} --write-index"


rule sv_annotate_repeats:
    threads:
        8
    input:
        calls="results/sv/{group}.annotated.gnomad.bcf",
        repeats="/projects/humgen/pipelines/dna-seq-nanopore/workflow/data/simple_repeats.tsv",
        header="/projects/humgen/pipelines/dna-seq-nanopore/workflow/data/simple_repeats.hdr.txt"
    output:
        calls="results/sv/{group}.annotated.repeats.bcf",
        calls_index="results/sv/{group}.annotated.repeats.bcf.csi",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/annotate/sv/{group}.repeats.log"
    shell:
        "bcftools annotate -a {input.repeats} {input.calls} -c CHROM,FROM,TO,SR_LOCATION,SR_PERIOD,SR_COPYNUMBER,SR_CONSENSUS_SIZE,SR_PER_MATCH,SR_SEQUENCE -h {input.header} -o {output.calls} --write-index --threads {threads}"


# In-house allele counts from the database built by the dna-seq-inhouse-db
# pipeline (INHOUSE_N, INHOUSE_CARRIERS, INHOUSE_HET, INHOUSE_HOM, INHOUSE_LOWQUAL
# and their _FRAC fractions). Its records are biallelic and left-aligned against
# the same GRCh38 reference, matching the output of `normalize`. Variants absent
# from the database get no INHOUSE_* tags, i.e. zero carriers.
INHOUSE_DB = config.get(
    "inhouse_db",
    "/projects/humgen/core/database_inhouse_nanopore/results/inhouse_db/inhouse.bcf",
)


rule annotate_inhouse:
    threads:
        8
    input:
        calls="results/snps/{group}.annotated.cadd.bcf",
        calls_index="results/snps/{group}.annotated.cadd.bcf.csi",
        database=INHOUSE_DB,
        database_index=f"{INHOUSE_DB}.csi",
    output:
        call="results/snps/{group}.annotated.inhouse.bcf",
        call_index="results/snps/{group}.annotated.inhouse.bcf.csi",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/annotate/{group}.inhouse.log"
    resources:
        mem_mb=4000
    params:
        info=",".join(
            f"INFO/INHOUSE_{tag}{frac}"
            for tag in ["N", "CARRIERS", "HET", "HOM", "LOWQUAL"]
            for frac in (["", "_FRAC"] if tag != "N" else [""])
        ),
    shell:
        "bcftools annotate -a {input.database} -c {params.info} --pair-logic exact"
        " --threads {threads} -O b -o {output.call} --write-index {input.calls} 2> {log}"
