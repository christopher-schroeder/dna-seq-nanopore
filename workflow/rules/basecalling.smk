import queue

rule link_input:
    input:
        dir=lambda wc: units[(units.sample_name == wc.sample) & (units.unit_name == wc.unit)]["fast5"]
    output:
        dir=directory("results/input/{sample}/{unit}")
    shell:
        """
        mkdir -p results/input
        ln -s -r {input} {output}
        """

def fast5_directories_input(wc):
    return units[(units.sample_name == wc.sample) & (units.unit_name == wc.unit)].iloc[0]["fast5"]

# rule basecalling:
#     threads:
#         8
#     resources:
#         gpu=1,
#     input:
#         fast5_dir=fast5_directories_input,
#     output:
#         ubam="results/basecalls/{sample}.{unit}.ubam"
#     params:
#         # cuda_device=lambda wc: f"cuda:{get_gpu()}",
#         model="results/resources/dna_r10.4.1_e8.2_400bps_hac@v4.1.0",
#         remora_args="",
#         basecaller_args="--modified-bases 5mCG_5hmCG",
#     conda:
#         "../envs/dorado.yaml"
#     script:
#         "../scripts/dorado.py"

# rule link_input:
#     input:
#         dir=lambda wc: units[(units.sample_name == wc.sample) & (units.unit_name == wc.U)]["fast5"]
#     output:
#         "results/input/{sample}/{unit}.{chunk}_of_{chunksize}.ubam"
#     shell:
#         """
#         mkdir -p results/input
#         ln -s -r {input} {output}
#         """

rule basecalling:
    threads:
        8
    resources:
        slurm_partition="GPUampere",
        slurm_extra="--gpus 1",
        gpus=1,
        runtime="7d",
    input:
        fast5_dir=fast5_directories_input,
    output:
        "results/basecalls/{sample}.{unit}.ubam",
    params:
        #cuda_device="cuda:all",
        #model="results/resources/dna_r10.4.1_e8.2_400bps_hac@v4.1.0",
        model="results/resources/dna_r9.4.1_e8_hac@v3.3",
        remora_args="",
        basecaller_args="--modified-bases 5mCG_5hmCG -c 1000",
        cuda_devices=lambda wc: os.environ.get("CUDA_VISIBLE_DEVICES")
    # conda:
    #     "../envs/dorado.yaml"
    log:
        "logs/basecalling/{sample}.{unit}.log"
    shell:
        """
        (/projects/humgen/pipelines/dna-seq-nanopore/workflow/tools/dorado-0.5.3-linux-x64/bin/dorado basecaller \
            {params.model} \
            {input.fast5_dir} \
            {params.remora_args} \
            {params.basecaller_args} \
            --device cuda:0,2,3,5 > {output}) > {log} 2>&1
        """

rule merge_unit_basecalls:
    threads:
        10
    resources:
        mem_mb=1000
    input:
        lambda wc: expand("results/basecalls/{{sample}}.{unit}.ubam", unit=get_units_for_sample(wc.sample))
    output:
        "results/basecalls_sample/{sample}.ubam"
    conda:
        "../envs/samtools.yaml"
    group:
        lambda wc: wc.sample
    shell:
        "samtools merge {output} {input} -@ {threads} -O bam"

# | samtools view -t {threads} --no-PG -b -o {output}
# def merge_input(wc):
#     sample = wc.sample
#     ret = []
#     print(sample)
#     for unit in units[units.sample_name == sample]["unit_name"].unique():
#         for chunk in range(1, (chunksize:=len(chunk_filenames[(sample, unit)]))+1):
#             ret.append(f"results/basecalls/{sample}/{unit}.{chunk}_of_{chunksize}.ubam")
#     return ret

# rule merge:
#     input:
#         merge_input
#     output:
#         "{sample}.txt"
#     shell:
#         "touch {output}"
