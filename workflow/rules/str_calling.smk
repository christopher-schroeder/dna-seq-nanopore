rule call_str:
    #first subset the repeats BED file, then use this to do straglr genotyping
    threads:
        64
    input:
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
        reference=REFERENCE,
    output:
        # vcf="results/str/{sample}.{chrom}.vcf.gz",
        tsv="results/str_single/{sample}.tsv",
        bed="results/str_single/{sample}.bed",
    params:
        prefix="results/str_single/{sample}"
        # --sex !{params.sex}
    conda:
        "../envs/straglr.yaml"
    shell:
        "straglr.py {input.xam} {input.reference} {params.prefix} --nprocs {threads} --min_support 1 --min_cluster_size 1"
        # '''
        # { grep !{chr} -Fw !{repeat_bed} || true; } > repeats_subset.bed
        #     straglr-genotype --loci repeats_subset.bed \
        #         --sample !{params.sample_name} \
        #         --tsv {output.tsv} \
        #         -v !{chr}_tmp.vcf \
        #         {input.reference} \
        #         --min_support 1 \
        #         --min_cluster_size 1
        #     cat !{chr}_tmp.vcf | vcfstreamsort | bgziptabix !{chr}_straglr.vcf.gz
        # '''


rule merged_str_bed:
    input:
        expand("results/str_single/{sample}.bed", sample=sample_names)
    output:
        "results/str_candidates/merged.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "cat {input} | cut -f 1-4 |  sort -u | bedtools sort > {output}"


rule genotype_str:
    threads:
        64
    input:
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
        reference=REFERENCE,
        bed="results/str_candidates/merged.bed",
    output:
        tsv="results/str/{sample}.tsv",
        bed="results/str/{sample}.bed",
    params:
        prefix="results/str/{sample}"
    conda:
        "../envs/straglr.yaml"
    shell:
        "straglr.py {input.xam} {input.reference} {params.prefix} --nprocs {threads} --min_support 1 --min_cluster_size 1 --loci {input.bed}"


rule str_to_vcf:
    input:
        bed="results/str/{sample}.bed",
    output:
        vcf="results/str/{sample}.vcf",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "python /projects/humgen/pipelines/dna-seq-nanopore/workflow/tools/strling_to_vcf.py {input.bed} {wildcards.sample} | bcftools sort > {output.vcf}"
