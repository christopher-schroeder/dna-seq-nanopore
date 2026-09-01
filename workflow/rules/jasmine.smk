import pandas as pd
import os

# rule annotate_control:
#     input:
#         variants="results/sv/{group}.jasmine.vcf",
#     output:
#         variants="results/sv/{group}.control_annotated.vcf",
#     conda:
#         "../envs/control.yaml"
#     params:
#         controls=controls
#     script:
#         "../scripts/annotate_frequencies.py"

# def get_ethnicity(s):
#     for e in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
#         if e in s:
#             return e
#     assert False

# controls = pd.read_csv("config/controls.tsv", sep="\t")
# controls["base_name"] = [os.path.basename(url).removesuffix(".vcf.gz") for url in controls["url"]]
# controls["name"] = controls["base_name"]
# controls["ethnicity"] = [get_ethnicity(s) for s in controls["base_name"]]
# controls = controls.set_index('name')


controls_path = "/projects/humgen/science/depienne/ont1000g/results/hg38/raw/vcf_modified_unpacked"
controls_path_glob = os.path.join(controls_path, "{control}.vcf")
controls = glob_wildcards(controls_path_glob).control


# def get_controls(ethnicity):
#     return controls[controls["ethnicity"]==ethnicity]


rule fix_jasmine:
    input:
        variants="results/sv/{group}.control_annotated.vcf",
    output:
        variants="results/sv/{group}.jasmine_fix.vcf",
        primary_samples=temp("results/sv/{group}.primary_samples.txt"),
        real_names=temp("results/sv/{group}.real_names.txt"),
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # jasmine prefixes every sample column with the index of its input file;
        # index 0 is this group's own calls, everything else is a control
        bcftools view {input.variants} -h | tail -n 1 | cut -f 10- | tr "\\t" "\\n" | grep '^0_' > {output.primary_samples}
        if [ ! -s {output.primary_samples} ]; then
            echo "ERROR: no 0_* sample column in {input.variants}" >&2
            exit 1
        fi
        cut -c3- {output.primary_samples} > {output.real_names}
        bcftools view -S {output.primary_samples} {input.variants} |
            sed 's/Description=""/Description="blank">/g' |
            bcftools reheader --samples {output.real_names} |
            bcftools sort -o {output.variants}
        """


# rule fix_jasmine:
#     input:
#         variants="results/sv/{group}.control_annotated.vcf",
#     output:
#         variants="results/sv/{group}.jasmine_fix.vcf",
#     conda:
#         "../envs/bcftools.yaml"
#     shell:
#         """
#             bcftools sort {input} -o {output.variants}
#         """



rule annotate_control:
    input:
        variants="results/sv/{group}.jasmine.vcf",
    output:
        variants="results/sv/{group}.control_annotated.vcf",
    conda:
        "../envs/control.yaml"
    params:
        controls=controls
    script:
        "../scripts/annotate_frequencies.py"


rule jasmine:
    threads:
        64
    input:
        variants="results/sv/{group}.filtered.vcf",
        #controls=expand("results/controls/{control}.vcf", control=controls["base_name"]),
        reference=REFERENCE,
        fai=f"{REFERENCE}.fai"
    output:
        file_list=temp("results/jasmine/{group}.file_list.txt"),
        variants="results/sv/{group}.jasmine.vcf",
    log:
        "logs/jasmine/{group}.log"
    conda:
        "../envs/jasmine.yaml"
    params:
        content="\n".join(
            ["results/sv/{group}.filtered.vcf"] +
            expand(controls_path_glob, control=controls)
        )
    shell:
        """
        echo '{params.content}' > {output.file_list}
        # jasmine writes the merged VCF itself; keep it out of the final path
        # until it is complete so an aborted run leaves no half-written output
        tmp={output.variants}.tmp.vcf
        jasmine -Xmx50g file_list={output.file_list} out_file=$tmp threads={threads} --require_first_sample genome_file={input.reference} --output_genotypes --ignore_strand --dup_to_ins --normalize-chrs --centroid_merging --allow_intrasample > {log} 2>&1
        mv $tmp {output.variants}
        """ 