#### Compare PIP under different prior ####
rm(list = ls())
library(ggplot2)
library(patchwork)
library(data.table)

susie.dir <- ""  # path to fine-mapping summary statistics
out.dir <- ""  # output directory
annot.dir <- ""  # torus annotation input file directory

gg_list <- list()

for (trait in c("aoa", "coa")) {
  uniform <- readRDS(paste0(susie.dir, '/', trait, '_uniform_susie_sumstats.RDS'))
  functional <- readRDS(paste0(susie.dir, '/', trait, '_functional_susie_sumstats.RDS'))
  annot <- fread(paste0(annot.dir, '/', trait, '_joint_enrichment_annotations.txt.gz'), header = T)
  annot <- annot[match(uniform$snp, annot$snp),][,-1]
  annot$label <- rowSums(annot)
  annot$label <- ifelse(annot$label == 0, 0, 1)
  df <- data.frame(snp = uniform$snp, uniform = uniform$susie_pip, functional = functional$susie_pip, annot = ifelse(annot$label == 1, "With annotation", "Without annotation"))
  df$annot <- factor(df$annot, levels = c("Without annotation", "With annotation"))
  gg_list[[trait]] <- ggplot(df |> dplyr::arrange(desc(annot)) |> dplyr::mutate(annot = factor(annot)), aes(x = uniform, y = functional, color = annot)) + 
    scale_color_manual(values = c("#A9A9A9", ifelse(trait == "aoa", "#00B6BC", "#F3BC50"))) + 
    geom_point(size = 15, alpha = 0.7) + 
    theme_minimal() + 
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
    labs(x = "PIP under uniform prior", y = 'PIP under functional prior', title = toupper(trait)) + 
    theme(legend.position = "bottom",
          plot.title = element_text(size = 78, hjust = 0.5), 
          axis.title = element_text(size = 73), 
          axis.text = element_text(size = 78), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 40)) + 
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

gg <- gg_list$aoa + gg_list$coa + 
  plot_layout(axis_titles = "collect")

png(paste0(out.dir, '/uniform_vs_functional_pip.png'), 
    height = 1200, width = 2000)
gg
dev.off()

