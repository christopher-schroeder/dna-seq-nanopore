rule filter_bam:
    threads:
        16
    input:
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
        reference=REFERENCE,
    output:
        xam="results/filtered/{sample}.cram",
        xai="results/filtered/{sample}.cram.crai",
    conda:
        "../envs/samtools.yaml"
    group:
        lambda wc: wc.sample
    shell:
        "samtools view -@ {threads} -F 2308 -o {output.xam}##idx##{output.xai} -O cram,embed_ref --write-index --reference {input.reference} {input.xam}"


rule sniffles2_snf:
    threads:
        64
    input:
        xam="results/filtered/{sample}.cram",
        xai="results/filtered/{sample}.cram.crai",
        reference=REFERENCE,
    output:
        snf="results/sv_snf/{sample}.snf",
    params:
        min_sv_length=30,
        cluster_merge_pos=0,
    conda:
        "../envs/sniffles.yaml"
    group:
        lambda wc: wc.sample
    shell:
        """
        sniffles \
            --threads {threads} \
            --sample-id {wildcards.sample} \
            --output-rnames \
            --minsvlen {params.min_sv_length} \
            --cluster-merge-pos {params.cluster_merge_pos} \
            --input {input.xam} \
            --long-ins-length 10000 \
            --phase \
            --snf {output.snf} \
            --reference {input.reference}
        """


rule sniffles2:
    threads:
        4
    input:
        snfs=lambda wc: expand("results/sv_snf/{sample}.snf", sample=get_group_samples(wc.group)),
        reference=REFERENCE,
    output:
        vcf="results/sv/{group}.vcf",
    params:
        tmp_vcf="results/sv/{group}.tmp.vcf",
        min_sv_length=30,
        cluster_merge_pos=0,
    conda:
        "../envs/sniffles.yaml"
    shell:
        """
        sniffles \
            --threads {threads} \
            --output-rnames \
            --minsvlen {params.min_sv_length} \
            --cluster-merge-pos {params.cluster_merge_pos} \
            --input {input.snfs} \
            --long-ins-length 10000 \
            --phase \
            --vcf {output.vcf} \
            --reference {input.reference}
        sed '/.:0:0:0:NULL/d' {output.vcf} > {params.tmp_vcf}
        mv {params.tmp_vcf} {output.vcf}
        """


rule filterCalls:
    threads: 1
    input:
        vcf="results/sv/{group}.vcf",
        mosdepth_summary=lambda wc: expand("results/mosdepth/{sample}.mosdepth.summary.txt", sample=get_group_samples(wc.group)),
        target_bed="results/resources/all_chromosomes.bed",
    output:
        vcf="results/sv/{group}.filtered.vcf",
        script="results/scripts/{group}.filtercalls.sh",
    params:
        basedir=BASEDIR,
        min_read_support="auto",
        min_read_support_limit=2,
    conda:
        "../envs/filtercalls.yaml"
    group:
        lambda wc: wc.group
    shell:
        """
        python {params.basedir}/scripts/get_filter_calls_command.py \
            --target_bedfile {input.target_bed} \
            --vcf {input.vcf} \
            --depth_summary {input.mosdepth_summary} \
            --min_read_support {params.min_read_support} \
            --min_read_support_limit {params.min_read_support_limit} > {output.script}

        sh {output.script} > {output.vcf}
        """

rule sortVCF:
    threads:
        1
    input:
        vcf="results/sv/{group}.filtered.vcf",
    output:
        vcf="results/sv/{group}.sorted.vcf",
    conda:
        "../envs/filtercalls.yaml"
    group:
        lambda wc: wc.group
    shell:
        "bcftools sort {input.vcf} > {output.vcf}"


# rule indexVCF:
#     threads: 1
#     input:
#         vcf="results/sv/{sample}.sorted.vcf",
#     output:
#         vcf="results/sv/{sample}.vcf.gz",
#         vcf_index="results/sv/{sample}.vcf.gz.tbi"
#     conda:
#         "../envs/vcflib.yaml"
#     group:
#         lambda wc: wc.sample
#     shell:
#         "cat {input.vcf} | bgziptabix {output.vcf}"
