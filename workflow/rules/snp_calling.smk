checkpoint make_chunks:
    threads:
        1
    input:
        xam="results/alignment/{sample}.cram",
        xam_index="results/alignment/{sample}.cram.crai",
        model="results/resources/r1041_e82_400bps_hac_v410.tar.gz", #lambda wc: f"results/resources/{select_model("dna_r10.4.1_e8.2_400bps_hac@v4.1.0")}.tar.gz",
        reference=REFERENCE,
    output:
        contigs="results/clair3/{sample}/tmp/CONTIGS",
        chunks="results/clair3/{sample}/tmp/CHUNK_LIST",
        cmd="results/clair3/{sample}/tmp/CMD",
        # pileup=directory("results/clair3/{sample}/tmp/pileup_output"),
        merge=directory("results/clair3/{sample}/tmp/merge_output"),
        phase=directory("results/clair3/{sample}/tmp/phase_output"),
        gvcf_tmp_path=directory("results/clair3/{sample}/tmp/gvcf_tmp_output"),
        full_alignment=directory("results/clair3/{sample}/tmp/full_alignment_output"),
        phase_vcf=directory("results/clair3/{sample}/tmp/phase_output/phase_vcf"),
        phase_bam=directory("results/clair3/{sample}/tmp/phase_output/phase_bam"),
        candidate_bed=directory("results/clair3/{sample}/tmp/full_alignment_output/candidate_bed"),
    params:
        clair3_prefix="results/clair3/{sample}/",
        sample_name="{sample}",
        vcf_fn="EMPTY",
        ctg_name="EMPTY",
        include_all_ctgs="false",
        min_qual=2,
        var_pct_full=0.7,
        ref_pct_full=0.1,
        snp_min_af=0.08,
        indel_min_af=0.15,
        min_contig_size=0,
        bed_args="",
        chunk_size=5000000,
        model_path="results/resources/r1041_e82_400bps_hac_v410",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        # CW-2456: save command line to add to VCF file (very long command...)
        echo "run_clair3.sh --bam_fn={input.xam} --ref_fn={input.reference} --vcf_fn={params.vcf_fn} --output=clair_output --platform=ont --sample_name={params.sample_name} --model_path={params.model_path} --ctg_name={params.ctg_name} --include_all_ctgs={params.include_all_ctgs} --chunk_num=0 --chunk_size={params.chunk_size} --qual={params.min_qual} --var_pct_full={params.var_pct_full} --ref_pct_full={params.ref_pct_full} --snp_min_af={params.snp_min_af} --indel_min_af={params.indel_min_af} --min_contig_size={params.min_contig_size}" > {output.cmd}
        # CW-2456: prepare other inputs normally
        python $(which clair3.py) CheckEnvs \
            --bam_fn {input.xam} \
            {params.bed_args} \
            --output_fn_prefix {params.clair3_prefix} \
            --ref_fn {input.reference} \
            --vcf_fn {params.vcf_fn} \
            --ctg_name {params.ctg_name} \
            --chunk_num 0 \
            --chunk_size {params.chunk_size} \
            --include_all_ctgs {params.include_all_ctgs} \
            --threads {threads}  \
            --qual {params.min_qual} \
            --sampleName {params.sample_name} \
            --var_pct_full {params.var_pct_full} \
            --ref_pct_full {params.ref_pct_full} \
            --snp_min_af {params.snp_min_af} \
            --indel_min_af {params.indel_min_af} \
            --min_contig_size {params.min_contig_size} \
            --cmd_fn {output.cmd}
        """

rule pileup_variants:
    threads:
        1
    resources:
        mem_mb=10000
    input:
        xam="results/alignment/{sample}.cram",
        xam_index="results/alignment/{sample}.cram.crai",
        model="results/resources/r1041_e82_400bps_hac_v410",
        reference=REFERENCE,
        index="results/resources/genome.dna.homo_sapiens.GRCh38.105.fasta.mmi",
        command="results/clair3/{sample}/tmp/CMD",
    output:
        vcf="results/clair3/{sample}/tmp/pileup_output/{sample}.{contig}.{chunk_id}_of_{total_chunks}.vcf",
    params:
        snp_min_af=0.08,
        indel_min_af=0.15,
        min_mq=5,
        min_cov=2,
        gvcf_tmp_path="results/clair3/{sample}/tmp/gvcf_tmp_output",
        named_pipe="pipe/pileup_variants/{sample}.{contig}.{chunk_id}_of_{total_chunks}.pipe",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        # // note: the VCF output here is required to use the contig
        # //       name since that's parsed in the SortVcf step
        # // note: snp_min_af and indel_min_af have an impact on performance
        # // TODO Fix REF_PATH
        # mkdir -p pipe/pileup_variants
        # rm -f {params.named_pipe}
        # mkfifo {params.named_pipe}
        # cat {params.named_pipe} > {output.vcf} & 
        '''
        python $(which clair3.py) CallVariantsFromCffi \
            --chkpnt_fn {input.model}/pileup \
            --bam_fn {input.xam} \
            --ref_fn {input.reference} \
            --ctgName {wildcards.contig} \
            --chunk_id {wildcards.chunk_id} \
            --chunk_num {wildcards.total_chunks} \
            --platform ont \
            --fast_mode False \
            --snp_min_af {params.snp_min_af} \
            --indel_min_af {params.indel_min_af} \
            --minMQ {params.min_mq} \
            --minCoverage {params.min_cov} \
            --call_snp_only False \
            --gvcf True \
            --temp_file_dir {params.gvcf_tmp_path} \
            --cmd_fn {input.command} \
            --pileup \
            --sampleName {wildcards.sample} \
            --call_fn {output.vcf}
        touch {output.vcf}
        '''

def aggregate_pileup_variants_input(wc):
    with checkpoints.make_chunks.get(sample=wc.sample).output.chunks.open() as f:
        df = pd.read_csv(f"results/clair3/{wc.sample}/tmp/CHUNK_LIST", names=("contig", "chunk_id", "total_chunks"), sep=" ")
        return expand(f"results/clair3/{wc.sample}/tmp/pileup_output/{wc.sample}.{{contig}}.{{chunk_id}}_of_{{total_chunks}}.vcf", zip, contig=df["contig"], chunk_id=df["chunk_id"], total_chunks=df["total_chunks"])

rule aggregate_pileup_variants:
    # Aggregates and sorts all variants (across all chunks of all contigs)
    # from pileup network. Determines quality filter for selecting variants
    # to use for phasing.
    threads: 2
    input:
        ref=REFERENCE,
        index="results/resources/genome.dna.homo_sapiens.GRCh38.105.fasta.mmi",
        # these need to be named as original, as program uses info from
        # contigs file to filter
        vcfs=aggregate_pileup_variants_input,
        contigs="results/clair3/{sample}/tmp/CONTIGS",
        command="results/clair3/{sample}/tmp/CMD",
    output:
        pileup_vcf="results/clair3/{sample}/pileup.vcf.gz",
        pileup_vcf_index="results/clair3/{sample}/pileup.vcf.gz.tbi",
        phase_qual="results/clair3/{sample}/phase_qual",
    params:
        input_dir="results/clair3/{sample}/tmp/pileup_output",
        vcf="results/clair3/{sample}/pileup.vcf",
        phase_qual_prefix="results/clair3/{sample}",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        '''
        python $(which clair3.py) SortVcf \
            --input_dir {params.input_dir}/ \
            --output_fn {params.vcf} \
            --sampleName {wildcards.sample} \
            --ref_fn {input.ref} \
            --contigs_fn {input.contigs} \
            --cmd_fn {input.command}
        bcftools index -n {output.pileup_vcf}
        bgzip -@ {threads} -fdc {output.pileup_vcf} | python $(which clair3.py) SelectQual --phase --output_fn {params.phase_qual_prefix}
        '''

rule select_het_snps:
    threads: 2
    input:
        pileup_vcf="results/clair3/{sample}/pileup.vcf.gz",
        pileup_vcf_index="results/clair3/{sample}/pileup.vcf.gz.tbi",
        phase_qual="results/clair3/{sample}/phase_qual",
    output:
        vcf="results/clair3/{sample}/split/{contig}.vcf",
        vcf_gz="results/clair3/{sample}/split/{contig}.vcf.gz",
        vcf_gz_index="results/clair3/{sample}/split/{contig}.vcf.gz.tbi",
    params:
        split_folder="results/clair3/{sample}/split"
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        '''
        python $(which clair3.py) SelectHetSnp \
            --vcf_fn {input.pileup_vcf} \
            --split_folder {params.split_folder} \
            --ctgName {wildcards.contig} \
            --qual_fn {input.phase_qual}
        bgzip -c {output.vcf} > {output.vcf_gz}
        tabix {output.vcf_gz}
        '''

rule phase_contig_haplotag:
    # Tags reads in an input BAM from heterozygous SNPs
    # Also haplotag for those modes that need it (--str)
    threads: 4
    input:
        vcf_gz="results/clair3/{sample}/split/{contig}.vcf.gz",
        vcf_gz_index="results/clair3/{sample}/split/{contig}.vcf.gz.tbi",
        xam="results/alignment/{sample}.cram",
        xam_index="results/alignment/{sample}.cram.crai",
        ref=REFERENCE,
    output:
        bam="results/phased/{sample}/{contig}_hp.cram",
        bai="results/phased/{sample}/{contig}_hp.cram.crai",
        vcf_gz="results/phased/{sample}/{contig}.vcf.gz",
    params:
        prefix="results/phased/{sample}/{contig}",
        vcf="results/phased/{sample}/{contig}.vcf",
    conda:
        "../envs/phasing.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        /projects/humgen/pipelines/dna-seq-nanopore/workflow/tools/longphase_linux-x64 phase --ont -o {params.prefix} \
            -s {input.vcf_gz} -b {input.xam} -r {input.ref} -t {threads}
        bgzip {params.vcf}

        tabix -f -p vcf {output.vcf_gz}

        whatshap haplotag \
            --reference {input.ref} \
            --ignore-read-groups \
            --regions {wildcards.contig} \
            {output.vcf_gz} \
            {input.xam} \
        | samtools view -b -1 -@ 3 -o {output.bam} -O cram,embed_ref=2 --write-index
        """


def get_contigs(sample):
    with checkpoints.make_chunks.get(sample=sample).output.contigs.open() as f:
        return f.read().splitlines()

rule merge_haplotagged_contigs:
    threads:
        8
    input:
        xams=lambda wc: expand("results/phased/{{sample}}/{contig}_hp.cram", contig=get_contigs(wc.sample)),
        ref=REFERENCE,
    output:
        cram="results/phased/{sample}.cram",
        crai="results/phased/{sample}.cram.crai",
    conda:
        "../envs/phasing.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        samtools merge -@ {threads} -O cram,embed_ref=2 -o {output.cram} {input.xams}
        samtools index -@ {threads} {output.cram}
        """

rule get_qual_filter:
    threads: 2
    input:
        vcf_gz="results/clair3/{sample}/pileup.vcf.gz",
        vcf_gz_index="results/clair3/{sample}/pileup.vcf.gz.tbi",
    output:
        qual="results/clair3/qual/{sample}/qual"
    params:
        var_pct_full=0.7,
        ref_pct_full=0.1,
        prefix="results/clair3/qual/{sample}"
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        '''
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        bgzip -fdc {input.vcf_gz} | \
        python $(which clair3.py) SelectQual \
                --output_fn {params.prefix} \
                --var_pct_full {params.var_pct_full} \
                --ref_pct_full {params.ref_pct_full} \
                --platform ont 
        '''

checkpoint create_candidates:
    # Create BED files for candidate variants for "full alignment" network
    # from the previous full "pileup" variants across all chunks of all chroms
    # 
    # Performed per chromosome; output a list of bed files one for each chunk.
    threads: 2
    input:
        ref=REFERENCE,
        vcf_gz="results/clair3/{sample}/pileup.vcf.gz",
        vcf_gz_index="results/clair3/{sample}/pileup.vcf.gz.tbi",
        qual="results/clair3/qual/{sample}/qual",
        # // this is used implicitely by the program
        # // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectCandidates.py#L146
    output:
        full_aln_file="results/clair3/candidate_bed/{sample}/{contig}/FULL_ALN_FILE_{contig}",
    params:
        var_pct_full=0.7,
        ref_pct_full=0.1,
        candidate_bed_dir="results/clair3/candidate_bed/{sample}/{contig}",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        # // This creates BED files as candidate_bed/<ctg>.0_14 with candidates
        # // along with a file the FULL_ALN_FILE_<ctg> listing all of the BED
        # // files.  All we really want are the BEDs, the file of filenames is
        # // used for the purposes of parallel in the original workflow.
        # // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/scripts/clair3.sh#L218
        # // TODO: would be nice to control the number of BEDs produced to enable
        # // better parallelism.
        '''
        python $(which clair3.py) SelectCandidates \
            --pileup_vcf_fn {input.vcf_gz} \
            --split_folder {params.candidate_bed_dir} \
            --ref_fn {input.ref} \
            --var_pct_full {params.var_pct_full} \
            --ref_pct_full {params.ref_pct_full} \
            --platform ont \
            --ctgName {wildcards.contig} \
            --qual_fn {input.qual}
        '''

rule evaluate_candidates:
    # // Run "full alignment" network for variants in a candidate bed file.
    # // phased_bam just references the input BAM as it no longer contains phase information.
    threads: 1
    input:
        phased_xam="results/phased/{sample}/{contig}_hp.cram",
        phased_xai="results/phased/{sample}/{contig}_hp.cram.crai",
        phased_vcf_gz="results/phased/{sample}/{contig}.vcf.gz",
        full_aln_file="results/clair3/candidate_bed/{sample}/{contig}/FULL_ALN_FILE_{contig}",
        ref=REFERENCE,
        cmd="results/clair3/{sample}/tmp/CMD",
    output:
        vcf="results/clair3/full_alignment/{sample}/{contig}.{n}_{total}.vcf"
    params:
        model="results/resources/r1041_e82_400bps_hac_v410",
        min_mq=5,
        min_cov=2,
        snp_min_af=0.08,
        indel_min_af=0.15,
        candidate_bed="results/clair3/candidate_bed/{sample}/{contig}/{contig}.{n}_{total}",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        python $(which clair3.py) CallVariantsFromCffi \
            --chkpnt_fn {params.model}/full_alignment \
            --bam_fn {input.phased_xam} \
            --call_fn {output.vcf} \
            --sampleName {wildcards.sample} \
            --ref_fn {input.ref} \
            --full_aln_regions {params.candidate_bed} \
            --ctgName {wildcards.contig} \
            --add_indel_length \
            --gvcf False \
            --minMQ {params.min_mq} \
            --minCoverage {params.min_cov} \
            --snp_min_af {params.snp_min_af} \
            --indel_min_af {params.indel_min_af} \
            --platform ont \
            --cmd_fn {input.cmd} \
            --phased_vcf_fn {input.phased_vcf_gz}
        """

import os

def get_full_alignment_for_contig(sample, contig):
    with checkpoints.create_candidates.get(sample=sample, contig=contig).output.full_aln_file.open() as f:
        for part in f.read().splitlines():
            yield f"results/clair3/full_alignment/{sample}/{os.path.basename(part)}.vcf"


rule aggregate_full_align_variants_contigs:
    input:
        lambda wc: get_full_alignment_for_contig(wc.sample, wc.contig)
    output:
        "results/clair3/full_alignment/{sample}/{contig}.done"
    shell:
        """
        touch {output}
        """

rule aggregate_full_align_variants:
    # // Sort and merge all "full alignment" variants
    threads: 2
    input:
        ref=REFERENCE,
        vcfs=lambda wc: expand("results/clair3/full_alignment/{{sample}}/{contig}.done", contig=get_contigs(wc.sample)),
        # full_aln_files=expand("results/clair3/candidate_bed/{sample}/{contig}/FULL_ALN_FILE_{contig}"),
        contigs="results/clair3/{sample}/tmp/CONTIGS",
        cmd="results/clair3/{sample}/tmp/CMD",
    output:
        full_aln_vcf="results/clair3/{sample}/full_alignment.vcf.gz",
        full_alignment_vcf_index="results/clair3/{sample}/full_alignment.vcf.gz.tbi",
        # non_var = "results/clair3/{sample}/non_var.gvcf",
    params:
        full_alignment_vcf="results/clair3/{sample}/full_alignment.vcf",
        vcf_dir="results/clair3/full_alignment/{sample}",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        '''
        python $(which clair3.py) SortVcf \
            --input_dir {params.vcf_dir} \
            --output_fn {params.full_alignment_vcf} \
            --sampleName {wildcards.sample} \
            --ref_fn {input.ref} \
            --cmd_fn {input.cmd} \
            --contigs_fn {input.contigs}
        if [ "$( bcftools index -n {output.full_aln_vcf} )" -eq 0 ]; then
            echo "[INFO] Exit in full-alignment variant calling"
            exit 0
        fi
        '''

        # if [ "$( bcftools index -n full_alignment.vcf.gz )" -eq 0 ]; then
        #     echo "[INFO] Exit in full-alignment variant calling"
        #     exit 0
        # fi

rule aggregate_full_align_variants_gvcf:
    threads: 2
    input:
        ref=REFERENCE,
        gvcf_tmp_path="results/clair3/{sample}/tmp/gvcf_tmp_output",
        contigs="results/clair3/{sample}/tmp/CONTIGS",
        command="results/clair3/{sample}/tmp/CMD",
        # a hack, the vcf existence also ensures the gvcf existence
        vcfs=lambda wc: expand("results/clair3/full_alignment/{{sample}}/{contig}.done", contig=get_contigs(wc.sample)),
    output:
        non_var="results/clair3/{sample}/non_var.gvcf",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        python $(which clair3.py) SortVcf \
            --input_dir {input.gvcf_tmp_path} \
            --vcf_fn_suffix .tmp.gvcf \
            --output_fn {output.non_var} \
            --sampleName {wildcards.sample} \
            --ref_fn {input.ref} \
            --cmd_fn {input.command} \
            --contigs_fn {input.contigs}
        """

rule merge_pileup_and_full_vars:
    # // Merge VCFs
    threads: 2
    input:
        ref=REFERENCE,
        pileup_vcf="results/clair3/{sample}/pileup.vcf.gz",
        pileup_vcf_index="results/clair3/{sample}/pileup.vcf.gz.tbi",
        full_aln_vcf="results/clair3/{sample}/full_alignment.vcf.gz",
        full_aln_vcf_index="results/clair3/{sample}/full_alignment.vcf.gz.tbi",
        non_var="results/clair3/{sample}/non_var.gvcf",
        full_aln_file=lambda wc: expand("results/clair3/candidate_bed/{{sample}}/{contig}/FULL_ALN_FILE_{contig}", contig=get_contigs(wc.sample)),
    output:
        merged_vcf="results/clair3/{sample}/merged/{contig}.vcf.gz",
        merged_vcf_index="results/clair3/{sample}/merged/{contig}.vcf.gz.tbi",
        merged_gvcf="results/clair3/{sample}/merged/{contig}.gvcf",
    params:
        merged_vcf="results/clair3/{sample}/merged/{contig}.vcf",
        candidate_bed_prefix="results/clair3/candidate_bed/{sample}/{contig}",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        '''
        echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"
        python $(which clair3.py) MergeVcf \
            --pileup_vcf_fn {input.pileup_vcf} \
            --bed_fn_prefix {params.candidate_bed_prefix} \
            --full_alignment_vcf_fn {input.full_aln_vcf} \
            --output_fn {params.merged_vcf}\
            --platform ont \
            --print_ref_calls False \
            --gvcf True \
            --haploid_precise False \
            --haploid_sensitive False \
            --gvcf_fn {output.merged_gvcf} \
            --non_var_gvcf_fn non_var.gvcf \
            --ref_fn {input.ref} \
            --ctgName {wildcards.contig}
        bgzip -c {params.merged_vcf} > {output.merged_vcf}
        tabix {output.merged_vcf}
        '''


rule post_clair_phase_contig:
    # // Phase VCF for a contig
    # // CW-2383: now uses base image to allow phasing of both snps and indels
    threads: 4
    input:
        ref=REFERENCE,
        merged_vcf="results/clair3/{sample}/merged/{contig}.vcf.gz",
        merged_vcf_index="results/clair3/{sample}/merged/{contig}.vcf.gz.tbi",
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
    output:
        phased_vcf="results/phased_contig/{sample}/{contig}.phased.vcf.gz",
        phased_vcf_index="results/phased_contig/{sample}/{contig}.phased.vcf.gz.tbi",
    params:
        prefix="results/phased_contig/{sample}/{contig}.phased",
        vcf="results/phased_contig/{sample}/{contig}.phased.vcf",
        #use_longphase = true
    conda:
        "../envs/phasing.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        echo "Using longphase for phasing"
        /projects/humgen/pipelines/dna-seq-nanopore/workflow/tools/longphase_linux-x64 phase --ont -o {params.prefix} --indels \
            -s {input.merged_vcf} -b {input.xam} -r {input.ref} -t {threads}
        bgzip {params.vcf}
        tabix -f -p vcf {output.phased_vcf}
        """

rule aggregate_all_variants:
    threads: 4
    input:
        ref=REFERENCE,
        merged_vcfs=lambda wc: expand("results/clair3/{{sample}}/merged/{contig}.vcf.gz", contig=get_contigs(wc.sample)),
        contigs="results/clair3/{sample}/tmp/CONTIGS",
        command="results/clair3/{sample}/tmp/CMD",
    output:
        final_vcf="results/snps/{sample}.vcf.gz",
        final_vcf_index="results/snps/{sample}.vcf.gz.tbi",
    params:
        input_dir="results/clair3/{sample}/merged",
        prefix="results/snps/{sample}",
    conda:
        "../envs/clair3.yaml"
    group:
        lambda wc: f"{wc.sample}"
    shell:
        """
        pypy $(which clair3.py) SortVcf \
            --input_dir {params.input_dir} \
            --output_fn {params.prefix}.vcf \
            --sampleName {wildcards.sample} \
            --ref_fn {input.ref} \
            --cmd_fn {input.command} \
            --contigs_fn {input.contigs}
        """