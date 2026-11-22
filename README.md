# somatic-tn

A Snakemake workflow for tumour-normal somatic variant analysis of targeted DNA panels.
It goes from paired FASTQ files to an annotated list of somatic SNVs and indels,
copy-number segments, QC metrics and one HTML report per tumour-normal pair.

I wrote it to work through the full secondary and tertiary analysis of an oncology
panel with standard, well-documented tools, and to check the calls against a known
truth set before trusting them.

**Status:** research use only. Not validated for diagnostic use.

## What it does

| Step | Tool |
|---|---|
| Read QC and adapter trimming | fastp |
| Mapping (read groups per sample) | BWA-MEM, samtools sort |
| Duplicate marking | GATK MarkDuplicates |
| Base quality recalibration | GATK BaseRecalibrator / ApplyBQSR (dbSNP, Mills indels) |
| Coverage on targets | mosdepth (mean depth, % of target bases at 100x / 250x / 500x) |
| Somatic SNVs and indels | GATK Mutect2 tumour-normal mode, germline resource (gnomAD), optional panel of normals |
| Artefact and contamination filtering | LearnReadOrientationModel, GetPileupSummaries, CalculateContamination, FilterMutectCalls |
| Normalisation | bcftools norm (left-align, split multi-allelic), PASS only |
| Annotation | Ensembl VEP offline from GTF + FASTA (gene, consequence, HGVSc / HGVSp) |
| Copy number | CNVkit (tumour vs matched normal, CBS segmentation, scatter plot, gene metrics) |
| QC summary | MultiQC (fastp, MarkDuplicates, mosdepth, samtools, bcftools, VEP) |
| Report | per-pair HTML: QC with PASS / CHECK flags, contamination, variants, CNV, benchmark, tool versions, input checksums |

Gene fusions are not covered. They need RNA (Arriba / STAR-Fusion) or a structural
variant caller on DNA, which is the next module I plan to add.

## Validation on a simulated tumour-normal pair

`test/simulate_tumour_normal.py` builds a tumour-normal pair from the 40 kb human
chr22 test reference with a known answer:

- 12 heterozygous germline SNVs in both samples (must not be called somatic)
- 65 somatic variants in the tumour only (50 SNVs, 15 indels) in five groups at
  VAF 30%, 10%, 5%, 2% and 1%
- one copy-number gain (3 copies) and one loss (1 copy) in the tumour
- tumour about 790x, normal about 400x, 2 x 150 bp reads

Results of the clean test run (`results/benchmark/`):

| Group | Truth | Detected | Recall |
|---|---|---|---|
| VAF 30% | 13 | 13 | 100% |
| VAF 10% | 13 | 13 | 100% |
| VAF 5% | 13 | 13 | 100% |
| VAF 2% | 13 | 11 | 85% |
| VAF 1% | 13 | 3 | 23% |

- germline SNVs called as somatic: 0 of 12
- false positives: 1, at the breakpoint of the simulated 9.4 kb deletion
  (the deletion junction is reported by Mutect2 as a long insertion). This is a
  real structural event read as a small variant, which is a good reminder to look
  at calls near CNV breakpoints in IGV.
- CNV: gain called at copy number 3 (log2 0.41), loss called at copy number 1
  (log2 -0.87). Both match the truth.

Under these conditions the detection limit sits between 2% and 5% VAF at about
800x. Simulated reads are cleaner than real ones (no GC bias, no FFPE damage, no
capture bias), so a real validation needs reference material with known variants
(for example Seraseq or Horizon standards) and replicate runs.

The second test pair uses the public nf-core test reads (about 9x on the target
region). It only checks that the workflow runs end to end on real reads; its QC
correctly comes out as CHECK because the depth is far too low for somatic calling.

## Install (macOS or Linux)

1. Install Miniforge: https://github.com/conda-forge/miniforge
2. Create the base environment:

```bash
mamba env create -f environment.yaml
mamba activate somatic-tn
```

On an Apple Silicon Mac some bioconda tools have no ARM build yet. Create the
tool environments as Intel builds (they run through Rosetta):

```bash
softwareupdate --install-rosetta --agree-to-license   # once
export CONDA_SUBDIR=osx-64
```

## Run the test

```bash
bash test/download_test_data.sh
python test/simulate_tumour_normal.py test/data/genome.fasta test/targets_test.bed test/sim
snakemake --use-conda --cores 4 --configfile test/config_test.yaml
```

The first run builds the conda environments and takes longer. After that the
test runs in about 5 minutes on a laptop. Open
`results/report/SIM_T_vs_SIM_N.report.html` and `results/multiqc/multiqc_report.html`.

## Run your own data

1. Copy `config/config.yaml` and `config/samples.tsv` and fill in your reference
   bundle (GRCh38 FASTA, GTF, panel BED, dbSNP, Mills, gnomAD AF-only, common SNPs).
   The GATK resource bundle has all the VCFs.
2. One normal per patient; every tumour of that patient is paired with it.
3. Set the QC thresholds in the config from your own validation data.

```bash
snakemake --use-conda --cores 8 --configfile config/my_run.yaml
```

## Outputs

```
results/
  bam/                  recalibrated BAM + index per sample
  qc/                   fastp, MarkDuplicates, mosdepth, flagstat, bcftools stats
  somatic/<pair>/       Mutect2 raw and filtered VCF, PASS VCF, variants.tsv
  annotation/<pair>/    VEP VCF and VEP summary
  cnv/<pair>/           CNVkit .cnr / .cns / .call.cns, scatter plot, gene metrics, segment table
  benchmark/            recall / precision and CNV match (only with a truth set)
  report/<pair>.report.html
  multiqc/multiqc_report.html
  pipeline_info/        software versions, input MD5 checksums, config used
```

## Layout

```
workflow/Snakefile
workflow/rules/       reference, align, qc, versions, somatic, annotate, cnv, report
workflow/scripts/     variant_table.py, report.py, benchmark.py, cnv_benchmark.py, checksums.py
workflow/envs/        one conda environment per tool group, versions pinned
config/               config and sample sheet templates
test/                 test data download, simulation, test config
```

## Author

Chahrazed Belhout
