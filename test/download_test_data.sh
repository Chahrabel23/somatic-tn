#!/usr/bin/env bash
# Download the public nf-core test data (human chr22 slice, 40 kb) used by the test run.
set -euo pipefail
BASE=https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens
OUT=$(dirname "$0")/data
mkdir -p "$OUT"
files=(
  genome/genome.fasta genome/genome.gtf
  genome/vcf/gnomAD.r2.1.1.vcf.gz genome/vcf/gnomAD.r2.1.1.vcf.gz.tbi
  genome/vcf/dbsnp_146.hg38.vcf.gz genome/vcf/dbsnp_146.hg38.vcf.gz.tbi
  genome/vcf/mills_and_1000G.indels.vcf.gz genome/vcf/mills_and_1000G.indels.vcf.gz.tbi
  illumina/fastq/test_1.fastq.gz illumina/fastq/test_2.fastq.gz
  illumina/fastq/test2_1.fastq.gz illumina/fastq/test2_2.fastq.gz
)
for f in "${files[@]}"; do
  [ -s "$OUT/$(basename $f)" ] || curl -fsSL -o "$OUT/$(basename $f)" "$BASE/$f"
  echo "ok  $(basename $f)"
done
