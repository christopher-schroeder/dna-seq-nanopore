#!/usr/bin/env python
import csv
import os
import sys
import textwrap
import glob

units = (
    pd.read_csv(
        "config/units.tsv",
        sep="\t",
        comment="#",
        )
)

samples = (
    pd.read_csv(
        "config/samples.tsv",
        sep="\t",
        dtype={"sample_name": str, "group": str},
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)

sample_names = units["sample_name"].unique()
groups = samples["group"].unique()


def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample_name"]

def get_units_for_sample(sample):
    return units[units["sample_name"] == sample]["unit_name"]

# def chunks(xs, n):
#     n = max(1, n)
#     return list(xs[i:i+n] for i in range(0, len(xs), n))

# def get_chunk_filenames(sample, unit):
#     path = units[(units.sample_name == sample) & (units.unit_name == unit)]["fast5"].iloc[0]
#     fast5_files = glob.glob(f"{path}/*.fast5", recursive=True)
#     return chunks(fast5_files, 10)

# chunk_filenames = dict()

# for _, row in units.iterrows():
#     sample_name, unit_name = row["sample_name"], row["unit_name"]
#     chunk_filenames[(sample_name, unit_name)] = get_chunk_filenames(sample_name, unit_name)


# """Map a basecalling model to a Clair3 model.

# An unknown basecalling model or basecalling model without an appropriate
# Clair3 model will explode with a large stderr message and exit non-zero.
# A happy basecalling model will print a Clair3 model to stdout and exit 0.
# """

# Delegating this to a Python script seems overkill but allows us to
# expand to more complex logic trivially in future.
# Plus I don't want to write this in Groovy right now.
def exit_obvious_error(header, error_msg, advice, basecaller, width=80):
    """Write an obvious looking error to stderr and quit."""
    line = ("-" * width) + '\n'
    msg = (
        f"The input basecaller configuration '{basecaller_cfg}' does not "
        "have a suitable Clair3 model "
    )
    sys.stderr.write(line)
    sys.stderr.write(f"[CRITICAL ERROR] {header}\n\n")
    sys.stderr.write(textwrap.fill(msg + error_msg, width))
    sys.stderr.write('\n\n')
    sys.stderr.write(textwrap.fill(advice, width))
    sys.stderr.write('\n')
    sys.stderr.write(line)
    sys.exit(os.EX_DATAERR)


def select_model(basecaller_cfg, lookup_table=os.path.join(workflow.basedir, "data/clair3_models.tsv")):
    with open(lookup_table) as tsv:
        for row in csv.DictReader(tsv, delimiter='\t'):
            if row["basecall_model_name"] == basecaller_cfg:
                model = row["clair3_model_name"]
                reason = row["clair3_nomodel_reason"]
                if model == "-" or model == "":
                    # Basecalling model valid but no Clair3 model
                    exit_obvious_error(
                        header="No appropriate Clair3 model.",
                        error_msg=f"because {reason}.",
                        advice=(
                            "It is not possible to run the SNP subworkflow "
                            "with this data.\n"
                        ),
                        basecaller_cfg=basecaller_cfg
                    )
                    break  # exit before here but this keeps my intention obvious
                else:
                    # Good model found
                    sys.stdout.write(model)
                    break
        else:
            # No model found (loop not broken)
            exit_obvious_error(
                header="Unknown basecaller configuration.",
                error_msg=(
                    "because the basecaller configuration has not been recognised."
                ),
                advice=(
                    "Check your --basecaller_cfg has been provided correctly. "
                ),
            )


rule get_all_chromosomes_bed:
    threads:
        1
    input:
        reference=REFERENCE,
    output:
        bed="results/resources/all_chromosomes.bed"
    conda:
        "../envs/pyfaidx.yaml"
    shell:
        "faidx --transform bed {input.reference} > {output}"


rule cut_all_chromosomes_bed:
    input:
        bed="results/resources/all_chromosomes.bed"
    output:
        bed="results/resources/all_chromosomes.cut.bed"
    params:
        depth_window_size=25000
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # extract first 3 columns of input BED to prevent col 4 leaking as label into outputs [CW-1702]
        # and convert them into windows of the given size [CW-2015]
        # The workflow now sort the bed input, merge overlapping intervals and then build windows
        # preventing crash in downstream tools [CW-2247]
        sort -k 1,1 -k2,2n {input.bed} | \
            bedtools merge -i - | \
            bedtools makewindows -b - -w {params.depth_window_size} > {output.bed}
        """


#NOTE grep MOSDEPTH_TUPLE if changing output tuple
rule mosdepth:
    threads:
        2
    input:
        xam="results/phased/{sample}.cram",
        xai="results/phased/{sample}.cram.crai",
        cut_bed="results/resources/all_chromosomes.cut.bed",
        reference=REFERENCE,
    output:
        regions="results/mosdepth/{sample}.regions.bed.gz",
        dist="results/mosdepth/{sample}.mosdepth.global.dist.txt",
        thresholds="results/mosdepth/{sample}.thresholds.bed.gz",
        summary="results/mosdepth/{sample}.mosdepth.summary.txt",
        # per_base="results/mosdepth/{sample}.per-base.bed.gz",
    params:
        prefix="results/mosdepth/{sample}",
        perbase_args="", #"--no-per-base"
    conda:
        "../envs/mosdepth.yaml"
    shell:
        """
        export REF_PATH={input.reference}
        export MOSDEPTH_PRECISION=3
        # Run mosdepth
        mosdepth \
        -x \
        -t {threads} \
        -b {input.cut_bed} \
        --thresholds 1,10,20,30 \
        {params.perbase_args} \
        {params.prefix} \
        {input.xam}
        """