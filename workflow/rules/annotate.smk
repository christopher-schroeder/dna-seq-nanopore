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
        "v3.3.5/bio/vep/cache"

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


rule annotate_str_vep:
    threads:
        4
    input:
        calls="results/str/{sample}.vcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls="results/str/{sample}.annotated.vep.bcf",
        stats="results/str/{sample}.annotated.vep.stats.html",
    params:
        plugins=[],
        extra="--everything"
    log:
        "logs/annotate_str/{sample}.vep.log",
    wrapper:
        "v3.5.2/bio/vep/annotate"


rule annotate_snps:
    threads:
        4
    input:
        calls="results/snps/{sample}.vcf.gz",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls="results/snps/{sample}.annotated.vcf.gz",
        stats="results/snps/{sample}.annotated.stats.html",
    params:
        plugins=[],
        extra="--symbol"
    log:
        "logs/annotate_snps/{sample}.log",
    wrapper:
        "v3.3.6/bio/vep/annotate"


rule transform_sv:
    input:
        calls="results/sv/{group}.sorted.vcf",
        reference=REFERENCE,
    output:
        calls="results/sv/{group}.transformed.bcf",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/transform.py"


rule annotate_sv_vep:
    threads:
        4
    input:
        calls="results/sv/{group}.transformed.bcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls="results/sv/{group}.annotated.vep.bcf",
        stats="results/sv/{group}.annotated.vep.stats.html",
    params:
        plugins=[],
        extra="--everything"
    log:
        "logs/annotate_sv/{group}.vep.log",
    wrapper:
        "v3.5.2/bio/vep/annotate"


rule annotate_sv_snpsift:
    input:
        call="results/sv/{group}.annotated.vep.bcf",
        database="/projects/humgen/science/resources/gnomad4.vcf.gz",
        database_index="/projects/humgen/science/resources/gnomad4.vcf.gz.tbi",
    output:
        call="results/sv/{group}.annotated.snpsift.bcf",
    params:
        extra="-name gnomad"
    log:
        "logs/annotate_gnomad_sv/{group}.snpsift.log",
    resources:
        mem_mb=30000
    threads: 8
    wrapper:
        "v3.5.0/bio/snpsift/annotate"


rule annotate_sv_repeats:
    input:
        call="results/sv/{group}.annotated.snpsift.bcf",
        repeats="/projects/humgen/pipelines/dna-seq-nanopore/workflow/data/simple_repeats.tsv",
        header="/projects/humgen/pipelines/dna-seq-nanopore/workflow/data/simple_repeats.hdr.txt"
    output:
        call="results/sv/{group}.annotated.repeats.bcf",
    log:
        "logs/annotate_gnomad_sv/{group}.repeats.log",
    params:
        fields="CHROM,FROM,TO,SR_LOCATION,SR_PERIOD,SR_COPYNUMBER,SR_CONSENSUS_SIZE,SR_PER_MATCH,SR_SEQUENCE"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools annotate -a {input.repeats} -c {params.fields} {input.call} -h {input.header} -o {output} -O b"