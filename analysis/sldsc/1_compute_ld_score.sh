#!/bin/bash
source activate ldsc

ldsc_dir=""  # path to ldsc software
sldsc_dir=""  # path to 1000G SNPs (available here: https://doi.org/10.5281/zenodo.7768714)
bed_dir=""  # path to bed files
ldscore_out_dir=""  # output directory

echo "Computing chromosome ${CHROM} LD score for ${BNAME}..."

## Step 1: Creating an annot file
echo "Step 1: creating annot file..."

python ${ldsc_dir}/make_annot.py \
--bed-file ${bed_dir}/${BNAME}.bed \
--bimfile ${sldsc_dir}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHROM}.bim \
--annot-file ${ldscore_out_dir}/${BNAME}.${CHROM}.annot.gz

## Step 2: Computing LD scores with an annot file
echo "Step 2: computing LD scores with annot file..."

python ${ldsc_dir}/ldsc.py \
--l2 \
--bfile ${sldsc_dir}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHROM} \
--print-snps ${sldsc_dir}/listHM3.txt \
--ld-wind-cm 1 \
--annot ${ldscore_out_dir}/${BNAME}.${CHROM}.annot.gz \
--thin-annot \
--out ${ldscore_out_dir}/${BNAME}.${CHROM}
