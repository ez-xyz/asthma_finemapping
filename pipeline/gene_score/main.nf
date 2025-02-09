#!/usr/bin/env nextflow

process intersect_cre_snp {
  publishDir "${params.outdir}", mode: 'copy'

  input:
    val bedtools_dir
    val cre
    val snp
    
  output:
    path "*_snp_cre_pairs.RDS"

  script:
  """
  Rscript 1_intersect_cre_snp.R \
    --bedtools-dir $bedtools_dir \
    --cre $cre \
    --snp $snp \
    --outname ${params.outname} \
    --outdir .
  """
}

process calculate_gene_score {
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path snp_cre_pairs
    val pip_cutoff
    val score_cutoff
    val contact
    val qtl
    val exonic
    val nearest
    val plot_params
    val gene_info
    
  output:
    path "*_gene_score_table.txt"
    path "*_high_confidence_genes.txt"
    path "*_gene_score_heatmap.pdf"

  script:
  """
  module unload java
  module load R/4.2.0
  Rscript 2_calculate_gene_score.R \
    --pip-cutoff $pip_cutoff \
    --score-cutoff $score_cutoff \
    --contact $contact \
    --qtl $qtl \
    --exonic $exonic \
    --nearest $nearest \
    --cre $snp_cre_pairs \
    --plot-params $plot_params \
    --gene-info $gene_info \
    --outname ${params.outname} \
    --outdir .
  """
}

workflow {
  snp_cre_pairs = intersect_cre_snp(
    params.bedtools_dir,
    params.cre,
    params.snp
  )

  calculate_gene_score(
    snp_cre_pairs,
    params.pip_cutoff,
    params.score_cutoff,
    params.contact,
    params.qtl,
    params.exonic,
    params.nearest,
    params.plot_params,
    params.gene_info
  )
}
