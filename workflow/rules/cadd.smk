# CADD v1.7 scoring of the small variant calls, mirroring the dna-seq-strelka
# pipeline. CADD.sh drives its own snakemake pipeline (with its own conda/
# apptainer environments), so it only needs the plain vcf.gz of the calls.
CADD_SCRIPT = config.get(
    "cadd_script", "/projects/humgen/tools/CADD-scripts-1.7.3/CADD.sh"
)
CADD_GENOME_BUILD = config.get("cadd_genome_build", "GRCh38")
# The CADD installation holds the scripts, models and annotation data that its
# own containers read via absolute paths. Apptainer is configured with a minimal
# bind list here, so the install tree has to be bound explicitly or the inner
# jobs fail with "No such file or directory" on $CADD/src/scripts/*.
CADD_DIR = os.path.dirname(CADD_SCRIPT)
# CADD.sh stages the calls and all of its intermediates in $TMPDIR; point it at
# node-local scratch rather than the (small) default /tmp.
CADD_TMPDIR = config.get("cadd_tmpdir", "/local/tmp")


rule bcf_to_vcf:
    threads:
        4
    input:
        "results/{type}/{group}.annotated.gnomad.bcf",
    output:
        temp("results/{type}/{group}.annotated.gnomad.vcf.gz"),
    wildcard_constraints:
        type="snps|sv",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/cadd/{type}/{group}.bcf_to_vcf.log"
    resources:
        mem_mb=16000
    shell:
        "bcftools view -O z -o {output} --threads {threads} {input} 2> {log}"


rule filter_primary_contigs:
    # CADD only ships annotations for the primary assembly, so scaffolds (and
    # the mitochondrion) have to be dropped before scoring. They are kept in
    # the annotated output, just without a CADD score.
    threads:
        2
    input:
        "results/{type}/{group}.annotated.gnomad.vcf.gz",
    output:
        temp("results/{type}/{group}.annotated.gnomad.primary.vcf.gz"),
    wildcard_constraints:
        type="snps|sv",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/cadd/{type}/{group}.primary_contigs.log"
    resources:
        mem_mb=4000
    shell:
        "bcftools view -t $(echo {{1..22}} X Y | tr ' ' ',') -O z -o {output} {input} 2> {log}"


rule cadd:
    threads:
        64
    input:
        "results/{type}/{group}.annotated.gnomad.primary.vcf.gz",
    output:
        "results/{type}/{group}.cadd.tsv.gz",
    wildcard_constraints:
        type="snps|sv",
    # CADD.sh starts its own apptainer containers, so this rule must run on the
    # host: container.sif ships no apptainer, and neither the CADD install nor
    # its images are visible from inside it. Ignored when running --sdm conda.
    container:
        None
    conda:
        "../envs/cadd.yaml"
    log:
        "logs/cadd/{type}/{group}.log"
    benchmark:
        "benchmarks/cadd/{type}/{group}.txt"
    resources:
        mem_mb=160000
    # A run takes days; a retry repeats all of it and can overlap a still
    # running attempt if job status is misreported.
    retries: 0
    shell:
        "TMPDIR={CADD_TMPDIR} {CADD_SCRIPT} -c {threads} -g {CADD_GENOME_BUILD}"
        " -r '--bind {CADD_TMPDIR} --bind {CADD_DIR}'"
        " -o {output} {input} &> {log}"


rule annotate_cadd:
    threads:
        2
    input:
        calls="results/{type}/{group}.annotated.gnomad.bcf",
        calls_index="results/{type}/{group}.annotated.gnomad.bcf.csi",
        tsv="results/{type}/{group}.cadd.tsv.gz",
    output:
        tsv_index=temp("results/{type}/{group}.cadd.tsv.gz.tbi"),
        call="results/{type}/{group}.annotated.cadd.bcf",
    wildcard_constraints:
        type="snps|sv",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/cadd/{type}/{group}.annotate.log"
    resources:
        mem_mb=16000
    shell:
        """
        (tabix -s 1 -b 2 -e 2 -f {input.tsv}
        bcftools annotate -a {input.tsv} -c CHROM,POS,REF,ALT,CADD_RAW,CADD_PHRED \
            -h <(printf '##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description="CADD phred score">\\n##INFO=<ID=CADD_RAW,Number=1,Type=Float,Description="CADD raw score">\\n') \
            --threads {threads} -O b -o {output.call} {input.calls}) 2> {log}
        """
