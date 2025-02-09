#!/bin/bash
source activate ldsc

ldsc_dir=""  # path to ldsc software
sldsc_dir=""  # path to 1000G SNPs (available here: https://doi.org/10.5281/zenodo.7768714)
out_dir=""  # output directory
gwas_dir=""  # path to GWAS summary statistics
trait=""  # GWAS trait

python ${ldsc_dir}/munge_sumstats.py \
--sumstats ${gwas_dir}/${trait}_gwas_before_munge.txt.gz \
--a1 a1 \
--a2 a0 \
--signed-sumstats zscore,0 \
--chunksize 500000 \
--merge-alleles ${sldsc_dir}/w_hm3.snplist \
--out ${out_dir}/${trait}_gwas_after_munge
