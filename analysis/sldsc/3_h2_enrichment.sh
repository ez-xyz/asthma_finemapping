#!/bin/bash
source activate ldsc

ldsc_dir=""  # path to ldsc software
sumstats_dir=""  # path to munged GWAS sumstats
sldsc_dir=""  # path to 1000G SNPs (available here: https://doi.org/10.5281/zenodo.7768714)
ldscore_dir=""  # path to ld scores
out_dir=""  # output directory
trait=""  # GWAS trait

mkdir -p ${out_dir}

python ${ldsc_dir}/ldsc.py \
--h2 ${sumstats_dir}/${trait}_gwas_after_munge.sumstats.gz \
--ref-ld-chr ${sldsc_dir}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ldscore_dir}/annot1.,${ldscore_dir}/annot2. \
--frqfile-chr ${sldsc_dir}/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${sldsc_dir}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out ${out_dir}/${trait}_h2_enrichment
