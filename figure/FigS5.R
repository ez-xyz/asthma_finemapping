rm(list = ls())
library(data.table)
library(ggplot2)

susie.dir <- ""  # path to fine-mapping result
score.dir <- ""  # path to gene score result
plot.dir <- ""  # output directory
fill.values <- c("#A16AE8", "#FD49A0", "#B4FEE7")  # fill colors

for (trait in c("aoa", "coa")) {
  score.df <- fread(file.path(score.dir, paste0(trait, "_gene_score_table.txt")), header = T)
  susie.sumstats <- readRDS(file.path(susie.dir, paste0(trait, "_functional_susie_gwas.RDS")))
  susie.sumstats$locus_cs <- paste0(susie.sumstats$locus, "_", susie.sumstats$CS)
  score.df$locus_cs <- susie.sumstats$locus_cs[match(score.df$snp, susie.sumstats$snp)]
  
  gg.df <- score.df |>
    dplyr::group_by(gene, locus_cs) |>
    dplyr::mutate(cs_score = sum(link_score * pip)) |>
    dplyr::ungroup() |>
    dplyr::mutate(score_prop = cs_score / gene_score) |>
    dplyr::distinct(gene, locus_cs, .keep_all = T) |>
    dplyr::select(gene, locus_cs, gene_score, score_prop) |>
    dplyr::arrange(desc(gene_score), desc(score_prop)) |>
    dplyr::group_by(gene) |>
    dplyr::mutate(cs = paste0("credible set ", row_number())) |>
    dplyr::mutate(n_cs = sum(score_prop > 0.1)) |>
    dplyr::ungroup() |>
    dplyr::filter(gene_score > 0.95 & n_cs > 1)
  
  gg.df$gene <- factor(gg.df$gene, levels = unique(gg.df$gene))
  
  gg <- ggplot(gg.df, aes(x = fct_rev(gene), y = score_prop)) + 
    geom_col(position = "fill", aes(fill = cs), alpha = 0.9, color = "black") + 
    theme_bw() + 
    labs(x = "", y = "Proportion of gene score") + 
    theme(axis.text.x = element_text(size = 25), 
          axis.text.y = element_text(size = 25, face = "italic"), 
          axis.title = element_text(size = 20), 
          legend.text = element_text(size = 25), 
          legend.title = element_blank(), 
          legend.spacing.y = unit(0.3, 'cm')) + 
    guides(fill = guide_legend(byrow = TRUE)) + 
    coord_flip() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) + 
    scale_fill_manual(values = fill.values)
  
  pdf(file.path(plot.dir, paste0(trait, "_gene_score_by_cs_score_prop.pdf")), 
      height = 8, width = 12)
  print(gg)
  dev.off()
}
