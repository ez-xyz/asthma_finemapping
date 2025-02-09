rm(list = ls())
library(ggplot2)
library(patchwork)
library(data.table)

in.dir <- ""  # path to gene score table
out.dir <- ""  # output directory

# full results
aoa <- fread(file.path(in.dir, "aoa_enrichment_results_wg_result.txt"))
coa <- fread(file.path(in.dir, "coa_enrichment_results_wg_result.txt"))

# weighted set cover
aoa_wsc <- fread(file.path(in.dir, "aoa_enriched_geneset_wsc_topsets_wg_result.txt"))
coa_wsc <- fread(file.path(in.dir, "coa_enriched_geneset_wsc_topsets_wg_result.txt"))

aoa_df <- aoa |> 
  dplyr::filter(geneSet %in% aoa_wsc$`# Coverage: 1`) |> 
  dplyr::arrange(enrichmentRatio)
aoa_df$description <- tools::toTitleCase(aoa_df$description)
aoa_df$description <- factor(aoa_df$description, levels = aoa_df$description)

coa_df <- coa |> 
  dplyr::filter(geneSet %in% coa_wsc$`# Coverage: 1`) |> 
  dplyr::arrange(enrichmentRatio)
coa_df$description <- tools::toTitleCase(coa_df$description)
coa_df$description <- factor(coa_df$description, levels = coa_df$description)

aoa_gg <- ggplot(aoa_df, aes(x = description, y = enrichmentRatio)) + 
  geom_bar(stat = "identity", fill = "#00B6BC", alpha = 0.9) + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "Enrichment ratio", x = "GO term", title = "AOA") + 
  theme(axis.text.x = element_text(size = 27), 
        axis.text.y = element_text(size = 27), 
        axis.title = element_text(size = 27), 
        plot.title = element_text(size = 27, hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5))
aoa_gg

coa_gg <- ggplot(coa_df, aes(x = description, y = enrichmentRatio)) + 
  geom_bar(stat = "identity", fill = "#F3BC50", alpha = 0.9) + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "Enrichment ratio", x = "GO term", title = "COA") + 
  theme(axis.text.x = element_text(size = 27), 
        axis.text.y = element_text(size = 27), 
        axis.title = element_text(size = 27), 
        plot.title = element_text(size = 27, hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5))
coa_gg

gg <- aoa_gg + coa_gg + plot_layout(ncol = 1, heights = c(1, 4), axis_titles = "collect")

pdf(paste0(out.dir, '/aoa_coa_wsc_GO_terms.pdf'), 
    height = 7, width = 15)
gg
dev.off()

