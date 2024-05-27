rule make_mmi:
    threads:
        64
    input:
        reference=REFERENCE,
    output:
        index="results/resources/genome.dna.homo_sapiens.GRCh38.105.fasta.mmi"
    log:
        "log/minimap_index.log"
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -t {threads} -x map-ont -d {output.index} {input.reference}"


rule download_clair3_model:
    output:
        gz="results/resources/r1041_e82_400bps_hac_v410.tar.gz",
        dir=directory("results/resources/r1041_e82_400bps_hac_v410"),
    params:
        url="https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v410.tar.gz",
        prefix="results/resources/",
    shell:
        """
        wget {params.url} -O {output.gz}
        mkdir -p {output.dir}
        tar -xvf {output.gz} -C {params.prefix}
        """
