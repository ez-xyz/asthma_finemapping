#!/usr/bin/env nextflow

process finemapping {
  publishDir "${params.outdir}", mode: 'copy'

  input:
    val gwas
    val prior
    val locus
    val ldr
    val rvar
    val n
    val L
    val purity

  output:
    path "*_susie_res.RDS"
    path "*_susie_sumstats.RDS"

  script:
  """
  module unload java
  module load gsl
  module load gcc
  module load R/4.2.0
  Rscript 1_finemapping.R \
    --gwas $gwas \
    --prior $prior \
    --locus $locus \
    --ldr $ldr \
    --rvar $rvar \
    --n $n \
    --L $L \
    --purity $purity \
    --outname ${params.outname} \
    --outdir .
  """
}

process get_cs_info {
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path susie_res
    path susie_sumstats
    val geneinfo

  output:
    path "*_cs_list.RDS"
    path "*_cs_summary.txt"

  script:
  """
  module unload java
  module load R/4.2.0
  Rscript 2_get_cs_info.R \
    --res $susie_res \
    --sumstats $susie_sumstats \
    --gene-info $geneinfo \
    --outname ${params.outname} \
    --outdir .
  """
}

process pip_attribution {
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path cs_list
    val annot
    val weight
    val pt
    val palette

  output:
    path "*_pip_attribution.pdf"

  script:
  """
  module unload java
  module load R/4.2.0
  Rscript 3_pip_attribution.R \
    --cs $cs_list \
    --annot $annot \
    --weight $weight \
    --pt $pt \
    --palette $palette \
    --outname ${params.outname} \
    --outdir .
  """
}

workflow {
  (susie_res_file, susie_sumstats_file) = finemapping(
    params.gwas,
    params.prior,
    params.locus,
    params.ldr,
    params.rvar,
    params.n,
    params.L,
    params.purity
  )

  (cs_list_file, cs_summary_file) = get_cs_info(
    susie_res_file,
    susie_sumstats_file,
    params.geneinfo
  )

  pip_attribution(
    cs_list_file,
    params.annot,
    params.weight,
    params.pt,
    params.palette
  )
}
