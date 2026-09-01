# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A **Snakemake pipeline** for whole-genome long-read (Oxford Nanopore) analysis: basecalling → alignment → phasing → variant calling (SNPs, SVs, STRs) → annotation → filtered TSV tables, plus QC and (currently disabled) methylation.

The entry point is `workflow/Snakefile`. All rules run from the **project root**, not from `workflow/`. `container.sif` at the root is declared via `containerized:`.

## Running the Pipeline

```bash
snakemake -n --use-conda                      # dry run — always do this first
snakemake --use-conda --slurm --jobs 50       # full run on SLURM
snakemake --use-conda --use-singularity --singularity-args "--nv"

# named targets (see "Targets" below)
snakemake only_mapping --use-conda --slurm --jobs 50
```

`rule all` currently builds only `results/tables/{group}.snps.maf.0.01.tsv` — the MAF-filtered SNP table per group. Every other product (SVs, STRs, QC) is reachable only via a named target or by editing `rule all`. Most alternative inputs are already there, commented out.

## Prerequisites not in the repo

None of these are created by a rule; the run fails without them.

- `config/config.yaml`, `config/units.tsv`, `config/samples.tsv` — **not in git and absent from this checkout**. Read by `Snakefile` / `common.smk` at parse time, so *every* Snakemake invocation fails without them.
  - `config.yaml` keys actually read by the code: `reference` (required, used to build `REFERENCE`), `basecalling_model` (`hac`/`sup`/`fast`, default `hac`), `jasmine` (bool, default `True`), `peddy` (bool, default `True` — whether the MultiQC report includes peddy), `vep` (optional GFF path for SV VEP), `methbat_regions` (required only if `methyl.smk` is re-enabled).
  - `units.tsv`: `sample_name`, `unit_name`, `fast5`. `samples.tsv`: `sample_name`, `group`, plus the *optional* pedigree columns `paternal_id`, `maternal_id`, `sex`, `phenotype` used to build peddy's PED (see QC below).
- `results/resources/<reference>.fasta` — the reference FASTA must be placed there by hand. `reference.smk` only builds `.fai` and `.mmi` from it.
- `resources/gnomad.genomes.v4.1.sites.vcf.gz` and `resources/gnomad.v4.1.sv.sites.no_chr.vcf.gz` (+ `.tbi`) — root-level `resources/`, hardcoded in `annotate.smk`. Note the SV database is the **no-chr-prefix** variant.
- `/projects/humgen/science/depienne/ont1000g/results/hg38/raw/vcf_modified_unpacked/*.vcf` — 1000G control SVs. `jasmine.smk` runs `glob_wildcards` over this path **at parse time**; if the directory is unreachable the control list silently becomes empty rather than erroring.
- `/local/work/cschroeder/snakemake-scratch/fs/` — local scratch on the execution node, hardcoded in the Clair3 rule.

## Sample vs. group

`units.tsv` drives per-unit and per-sample work (basecalling, mapping, phasing, STR, methylation, QC). `samples.tsv` maps samples to groups; SNP merging, SV joint calling, annotation and tables are **per group**. `wildcard_constraints` in the Snakefile pin `group` and `sample` to the values found in these files.

## Data flow — the CRAM chain

The three alignment directories are distinct stages and downstream rules are picky about which one they read:

```
results/alignment/{sample}.cram     minimap2 map-ont, qs>=10 filter (rejects → alignment_failed/)
  → results/phased/{sample}.cram    whatshap haplotag using that sample's own Clair3 SNPs
      → results/filtered/{sample}.cram   samtools -F 2308 (primary alignments only) → Sniffles
```

`snp_calling`, `samtools_stats` and `qualimap` read `alignment/`; `mosdepth`, `str_calling`, `methyl` and `filter_bam` read `phased/`; only Sniffles reads `filtered/`.

Full flow:

```
FAST5/POD5 (units.tsv)
  → dorado basecaller (GPU, per unit) → results/basecalls/{sample}.{unit}.ubam
  → samtools merge → results/basecalls_sample/{sample}.ubam
  → minimap2 → results/alignment/{sample}.cram
  → Clair3 (per sample, phased output) → results/snps_sample/{sample}.vcf.gz
      → whatshap haplotag → results/phased/{sample}.cram
      → bcftools merge per group → normalize → VEP → gnomAD → inhouse → vembrane filter → table
  → Sniffles2 .snf per sample → joint call per group → filter by depth-derived read support
      → Jasmine merge with 1000G controls → control freqs → sample fixup
      → transform → VEP → transform back → gnomAD → simple repeats → table
  → straglr: per-sample discovery → merged loci BED → per-sample genotyping
  → QC: NanoPlot, samtools stats/flagstat/idxstats, qualimap, mosdepth,
        peddy (per group) → multiqc
```

## Rule modules (`workflow/rules/`)

| File | Notes |
|------|-------|
| `common.smk` | Loads `units.tsv`/`samples.tsv`; `index_bcf`, `mosdepth`, chromosome BED helpers |
| `basecalling.smk` | Dorado GPU basecalling, unit merging |
| `mapping.smk` | minimap2 + Q-score filter → CRAM |
| `reference.smk` | faidx, minimap2 `.mmi`, Clair3 model download |
| `snp_calling.smk` | Clair3 via local scratch, group merge, whatshap haplotag |
| `sv_calling.smk` | `filter_bam`, Sniffles2 SNF + joint call, read-support filter |
| `str_calling.smk` | straglr discovery + genotyping |
| `jasmine.smk` | Jasmine merge with controls, `annotate_control`, `fix_jasmine` |
| `annotate.smk` | VEP cache/plugins, normalize, VEP + gnomAD for SNPs and SVs, repeats, inhouse DB |
| `filtering.smk` | `vembrane filter` on `gnomad_AF <= {maf}` and `QUAL > 10` |
| `visualization.smk` | vembrane TSV tables (SNP, SNP-MAF, SV, STR) |
| `qc.smk` | NanoPlot, samtools stats/flagstat/idxstats, qualimap, peddy, multiqc |
| `methyl.smk` | **`include:` is commented out in the Snakefile.** MethBat-based (pileup → profile → cohort → outliers) plus modkit phased pileup |
| `cadd.smk` | Included but entirely commented out — a no-op |
| `which_gpu.smk` | Not included anywhere |

## Named targets in the Snakefile

`only_basecalling_units`, `only_basecalling`, `only_mapping`, `only_mapping_avail` (globs whatever `basecalls_sample/` already has), `only_clair`, `only_snp_calling`, `only_snp_vep`, `only_snp_gnomad`, `only_inhouse`, `only_snf`, `only_sv_calling`, `only_sv_annotated`, `only_sv_table`, `only_jasmine`, `only_jasmine_fix`, `only_download_controls`, `only_qc` (the full `results/qc/multiqc.html`), `only_nanoplot`, `only_peddy`, `test` (hardcoded sample names).

## Implementation notes / gotchas

- **Dorado is invoked by absolute path**, `workflow/tools/dorado-2.0.0-linux-x64/bin/dorado`, and the device is hardcoded `--device cuda:1` (not `cuda:0`, and not derived from the SLURM allocation). Basecalling requests `slurm_partition="GPUampere"`, `--gpus 1`, `runtime=7d`.
- **Basecalling model mismatch**: the rule expects `results/resources/dna_r10.4.1_e8.2_400bps_{model}@v5.2.0`, but the models vendored in `workflow/tools/` are `@v4.2.0`. The v5.2.0 model must be fetched into `results/resources/` separately.
- **Clair3 stages through local scratch**: inputs are `cp`'d to `/local/work/cschroeder/snakemake-scratch/fs/`, Clair3 runs there, and the VCF is copied back. The scratch output dir is `rm -rf`'d on success.
- **Clair3 has hand-written sanity checks** because it exits 0 on partial failures (e.g. GNU parallel dying on a missing `$TMPDIR` inside the container). The rule fails if fewer than 1,000,000 variants are called, or if any contig in `tmp/CONTIGS` produced zero variants (chrY exempted). Lowering `min_variants` is the knob for non-human or targeted data.
- **The Clair3 model is pinned** to `r1041_e82_400bps_hac_v410` in `reference.smk` and `snp_calling.smk`; it does *not* follow `basecalling_model`. `select_model()` in `common.smk` and `workflow/data/clair3_models.tsv` implement lookup-table selection but are **not wired up**.
- **SV symbolic-allele round trip**: VEP cannot handle literal INS/DEL sequences, so `sv_transform` → `sv_annotate_vep` → `sv_transform_back` converts ALTs to `<INS>`/`<DEL>` and restores them. Keep those three adjacent when editing.
- **Jasmine on/off changes the SV graph**: `sv_transform_input()` in `annotate.smk` picks `{group}.jasmine_fix.vcf` when `config["jasmine"]` is true, else `{group}.filtered.vcf`.
- **`fix_jasmine` strips the control samples** by keeping only header columns starting with `0` (Jasmine's first-input prefix) and stripping that 2-char prefix to recover real sample names.
- **`ruleorder: annotate_snps_gnomad > index_bcf`** resolves the ambiguity on `*.bcf.csi`, and **`ruleorder: peddy_vcf > tabix`** the one on `*.vcf.gz.tbi`. Any new rule writing its own `.csi`/`.tbi` alongside its main output needs the same treatment — the generic `index_bcf`/`tabix` rules match every path.
- Absolute `/projects/humgen/pipelines/dna-seq-nanopore/...` paths appear in `str_calling.smk` (`strling_to_vcf.py`), `annotate.smk` (`simple_repeats.tsv`/`.hdr.txt`) and `basecalling.smk` — the repo is not relocatable as-is.

## QC subworkflow

`snakemake only_qc` builds `results/qc/multiqc.html`. The MultiQC modules that actually
fire are **nanostat, samtools (stats + flagstat + idxstats), qualimap, mosdepth and peddy**;
`workflow/data/multiqc_config.yaml` sets the report title and the `extra_fn_clean_exts`
that collapse `<sample>.NanoStats` back to `<sample>` so every tool shares one row.

- **Nanopore read stats come from NanoPlot, not NanoStat.** NanoStat itself only takes
  `--bam`, while the pipeline produces CRAM; NanoPlot accepts `--cram` and writes the same
  `NanoStats.txt` that MultiQC parses. `--tsv_stats` selects the `Metrics<TAB>dataset`
  flavour MultiQC prefers, and `--no_static` skips the matplotlib renderings, which are
  slow and need `kaleido` at WGS read counts. Output names are
  `results/qc/nanoplot/{sample}/{sample}.NanoStats.txt` — NanoPlot's `--prefix` is
  concatenated, not joined, so the trailing `.` in the prefix is load-bearing.
- **`samtools flagstat`/`idxstats` take no `--reference`** (only `samtools stats` does) and
  need none, since neither decodes sequence. Don't "fix" them by adding one — it is a hard
  CLI error.
- **`only_qc` transitively requires Clair3.** mosdepth reads the haplotagged
  `results/phased/` CRAM and peddy needs called SNPs, so QC is not a cheap post-mapping step.
  Set `peddy: False` in the config to drop peddy from the report.
- **peddy's PED is generated from `samples.tsv`** by `rule peddy_ped`: family = `group`, and
  the optional `paternal_id` / `maternal_id` / `sex` / `phenotype` columns are used when
  present (sex accepts `1`/`m`/`male` and `2`/`f`/`female`, anything else is 0 = unknown).
  With only the two required columns you still get inferred sex and relatedness, just
  nothing to check them against.
- **peddy gets its own re-merged VCF** (`rule peddy_vcf`) rather than reusing
  `results/snps/{group}.bcf`: plain `bcftools merge` leaves `./.` wherever a sample was not
  called, which peddy reads as missing rather than hom-ref and which wrecks relatedness.
  `peddy_vcf` merges the per-sample calls with `-0` and keeps biallelic SNPs only.
- **`--sites hg38` matches this reference.** peddy's bundled `GRCH38.sites` is not
  chr-prefixed, like the reference here. Those sites are autosome-only; the sex check reads
  X straight from the VCF, so `sex_check.csv` is *not* written when the VCF has no X
  variants — the rule then fails on a missing output, which is intended (a silently skipped
  sex check is worse).
- **`workflow/envs/peddy.post-deploy.sh` patches peddy after install.** peddy 0.4.8 (2021,
  final release) calls `np.fromstring(<bytes>)`, which NumPy 2.0 removed, so `het_check`,
  `sex_check`, the PCA and the HTML report all die while `ped_check.csv` is written first
  and looks fine. The script rewrites that one call to `np.frombuffer`. Pinning `numpy <2`
  instead does **not** work: every available cyvcf2 build is compiled against the NumPy 2
  ABI and fails with `numpy.dtype size changed`.

## Stale files — do not treat as live code

`workflow/rules/snp_calling.smk.bak`, `workflow/scripts/*.bak`/`*.bak2`, `workflow/rules/cadd.smk` (all commented), `which_gpu.smk` + `which_gpu_fabian.py` (not included), `scripts/dorado.py` (rule commented out), and the Borealis stack (`build_borealis_model.R`, `run_borealis_cohort.R`, `score_borealis_sample.R`, `modkit_to_bismark.py`, `workflow/containers/borealis.*`) — no rule references the Borealis/bismark path any more; methylation moved to MethBat. `workflow/tools/` vendors seven Dorado versions; only `dorado-2.0.0-linux-x64` is used.

Large chunks of `snp_calling.smk` are a commented-out hand-rolled reimplementation of Clair3's internal stages (chunking, pileup, full-alignment, longphase). The live path is the single `clair3_call_variants` rule calling `run_clair3.sh`.
