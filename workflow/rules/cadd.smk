# rule make_chrom_sizes:
#     input:
#         fai=f"{REFERENCE}.fai"
#     output:
#         sizes=f"results/resources/{REFERENCE}.chrom_sizes.txt",
#     shell:
#         """
#         cut -f1,2 {input.fai} > {output.sizes}
#         """

# rule decompress:
#     input:
#         "results/snps/{group}.annotated.gnomad.bcf",
#     output:
#         temp("results/snps/{group}.annotated.gnomad.bcf"),
#     shell:
#         """
#         zcat {input} > {output}
#         """