rm(list = ls())
library(data.table)
library(ggplot2)
library(patchwork)

in.dir <- ""  # path to gene score table
out.dir <- ""  # output directory

trait <- "coa"  # trait
score.cutoff <- 0.95  # gene score cutoff for high-confidence risk genes
link.order <- c("Nearest", "ABC", "PCHi-C", "eQTL", "Exon")  # order to display CRE-gene link

score.df <- fread(file.path(in.dir, paste0(trait, "_gene_score_table.txt")), 
                  header = T)
gene.ranked <- score.df %>% 
  dplyr::arrange(desc(gene_score)) %>% 
  dplyr::filter(gene_score > score.cutoff) %>% 
  dplyr::distinct(gene)
gene.list <- list(part1 = gene.ranked$gene[1:12], 
                  part2 = gene.ranked$gene[13:24], 
                  part3 = gene.ranked$gene[25:35])

df.list <- list()
for (j in c("part1", "part2", "part3")) {
  plot.df <- score.df |>
    dplyr::filter(gene %in% gene.list[[j]]) |>
    dplyr::distinct(gene, n_snp, link, weight, gene_score) |>
    dplyr::mutate(gene = paste0(gene, " (", n_snp, ")"))
  
  # some genes may have no support from one or multiple evidence categories - add them back
  for (i in unique(plot.df$gene)) {
    zero.link <- setdiff(link.order, plot.df$link[which(plot.df$gene == i)])
    n.zero.link <- length(zero.link)
    if (n.zero.link == 0) {
      next
    }
    supp.df <- data.frame(gene = rep(i, n.zero.link), 
                          n_snp = rep(unique(plot.df$n_snp[which(plot.df$gene == i)]), n.zero.link), 
                          link = zero.link, 
                          weight = rep(0, n.zero.link), 
                          gene_score = rep(unique(plot.df$gene_score[which(plot.df$gene == i)]), n.zero.link))
    plot.df <- bind_rows(plot.df, supp.df)
  }
  plot.df <- plot.df |> dplyr::arrange(desc(gene_score))
  if (nrow(plot.df) != length(link.order) * length(unique(plot.df$gene))) {
    stop("There should be ", length(link.order) * length(unique(plot.df$gene)), " rows in the final data frame.")
  }
  
  df.list[[j]] <- plot.df
}

gg.list <- list()
for (i in c("part1", "part2", "part3")) {
  gg.list[[i]] <- ggplot(df.list[[i]], aes(x = factor(link, levels = link.order), 
                      y = forcats::fct_rev(factor(gene, levels = unique(gene))), fill = weight)) + 
    geom_point(shape = 21, stroke = 0, size = 8) + 
    scale_x_discrete(position = "top") + 
    scale_fill_gradient(low = "#f1f1f1", high = "#F3BC50", 
                        breaks = seq(0, 2, 1), 
                        limits = c(0, 2)) + 
    coord_fixed(ratio = 0.3) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          axis.title = element_blank(), 
          axis.text.x.top = element_text(size = 15, angle = 45, vjust = 0.25), 
          axis.text.y = element_text(size = 15, face = "italic"), 
          legend.text = element_text(size = 15), 
          legend.title = element_blank()) + 
    theme(aspect.ratio = 2/1) + 
    guides(fill = guide_colorbar(ticks.colour = NA, breaks = c(0, 1, 2), 
                                 title.position = "top", label.position = "right", order = 2))
}

gg <- gg.list$part1 + 
  gg.list$part2 + 
  gg.list$part3 + 
  plot_layout(guides = "collect") & theme(plot.background = element_rect(fill = 'transparent', color = NA), 
                                          legend.background = element_rect(fill = 'transparent', color = NA))

pdf(file.path(out.dir, paste0(trait, "_gene_score_heatmap.pdf")), 
    height = 5, width = 11)
print(gg)
dev.off()

