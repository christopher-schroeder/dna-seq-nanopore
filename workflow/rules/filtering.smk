rule snps_filter_by_maf:
    input:
        "results/{x}.bcf"
    output:
        "results/{x}.maf.{maf}.bcf"
    # wildcard_constraints:
    #     dataset="\d+"
    conda:
        "../envs/vembrane.yaml"
    benchmark:
        "benchmarks/filter_by_maf/{x}.{maf}.txt"
    resources:
        mem_mb=256
    group:
        "group"
    shell:
        """vembrane filter "INFO.get('gnomad_AF', 0) <= {wildcards.maf} and QUAL > 10" {input} -o {output}"""


# rule sv_filter_by_maf:
#     input:
#         "results/sv/{group}.gnomad_annotated.bcf"
#     output:
#         "results/sv/{group}.gnomad_annotated.maf.{maf}.bcf"
#     conda:
#         "../envs/vembrane.yaml"
#     benchmark:
#         "benchmarks/filter_by_maf_sv/{group}.{maf}.txt"
#     resources:
#         mem_mb=256
#     group:
#         "group"
#     shell:
#         """vembrane filter "INFO.get('gnomad_AF', 0) <= {wildcards.maf} and QUAL > 10" {input} -o {output}"""