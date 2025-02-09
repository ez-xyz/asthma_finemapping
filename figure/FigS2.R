rm(list = ls())
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

plot.dir <- ""  # output directory
sldsc.dir <- ""  # path to S-LDSC result

# Ideal order by group
nogp.order <- c("B_lung", "B_blood", "T_lung", "CD4_blood", "CD8_blood", "NK_lung", "NK_blood", 
                "Macrophage_lung", "mDC_blood", "pDC_blood", "Mono_blood", 
                "BEC_lung", "AT1_lung", "AT2_lung", "PNEC_lung", "Basal_lung", "Ciliated_lung", "Club_lung", 
                "ASMC_lung", "MatrixFib1_lung", "MatrixFib2_lung", "Myofibroblast_lung", "Pericyte_lung", 
                "Arterial_lung", "Cap1_lung", "Cap2_lung", "Lymphatic_lung")

# S-LDSC cell type order
nogp.category <- c("B_blood", "CD4_blood", "CD8_blood", "Mono_blood", "NK_blood", "mDC_blood", 
                   "pDC_blood", "AT1_lung", "AT2_lung", "Arterial_lung", "B_lung", "Basal_lung", 
                   "Cap1_lung", "Cap2_lung", "Ciliated_lung", "Club_lung", "Lymphatic_lung", 
                   "Macrophage_lung", "MatrixFib1_lung", "MatrixFib2_lung", "Myofibroblast_lung", 
                   "NK_lung", "PNEC_lung", "Pericyte_lung", "T_lung", "ASMC_lung", "BEC_lung")

aoa.nogp <- fread(paste0(sldsc.dir, '/aoa_h2_enrichment_nogp.results'))
coa.nogp <- fread(paste0(sldsc.dir, '/coa_h2_enrichment_nogp.results'))

aoa.nogp <- aoa.nogp[aoa.nogp$Category %in% paste0("L2_", 1:length(nogp.category)),]
coa.nogp <- coa.nogp[coa.nogp$Category %in% paste0("L2_", 1:length(nogp.category)),]

aoa.nogp$Category <- coa.nogp$Category <- nogp.category

assign_group <- function(x) {
  if (x %in% c("B_lung", "B_blood", "T_lung", "CD4_blood", "CD8_blood", "NK_lung", "NK_blood")) {
    return("Lymphoid")
  }
  if (x %in% c("Macrophage_lung", "mDC_blood", "pDC_blood", "Mono_blood")) {
    return("Myeloid")
  }
  if (x %in% c("BEC_lung", "AT1_lung", "AT2_lung", "PNEC_lung", "Basal_lung", "Ciliated_lung", "Club_lung")) {
    return("Epithelial")
  }
  if (x %in% c("ASMC_lung", "MatrixFib1_lung", "MatrixFib2_lung", "Myofibroblast_lung", "Pericyte_lung")) {
    return("Mesenchymal")
  }
  if (x %in% c("Arterial_lung", "Cap1_lung", "Cap2_lung", "Lymphatic_lung")) {
    return("Endothelial")
  }
}

aoa.nogp$group <- sapply(aoa.nogp$Category, assign_group)
coa.nogp$group <- sapply(coa.nogp$Category, assign_group)
aoa.nogp$group <- factor(aoa.nogp$group, levels = c("Lymphoid", "Myeloid", "Epithelial", "Mesenchymal", "Endothelial"))
coa.nogp$group <- factor(coa.nogp$group, levels = c("Lymphoid", "Myeloid", "Epithelial", "Mesenchymal", "Endothelial"))

fill.values <- c("#7FE5F0", "#6395ee", "#e0afff", "#9D33E2", "#541675")

gg_aoa <- ggplot(data = aoa.nogp, aes(x = factor(Category, levels = nogp.order), y = Enrichment, fill = group)) + 
  geom_bar(stat = "identity", alpha = 0.9) + 
  geom_errorbar(aes(x = factor(Category, levels = nogp.order), 
                    ymin = Enrichment - 2 * Enrichment_std_error, 
                    ymax = Enrichment + 2 * Enrichment_std_error, 
                    width = 0.2), 
                color = "black") + 
  scale_y_continuous(limits = c(-20, 20)) + 
  scale_x_discrete(limits = rev) + 
  scale_fill_manual(values = fill.values) + 
  coord_flip() + 
  labs(title = "AOA") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5), 
        legend.title = element_blank(), 
        axis.title = element_text(size = 13), 
        axis.text = element_text(size = 13), 
        legend.spacing.y = unit(0.2, 'cm'), 
        legend.text = element_text(size = 13)) + 
  labs(x = "Cell type")

gg_coa <- ggplot(data = coa.nogp, aes(x = factor(Category, levels = nogp.order), y = Enrichment, fill = group)) + 
  geom_bar(stat = "identity", alpha = 0.9) + 
  geom_errorbar(aes(x = factor(Category, levels = nogp.order), 
                    ymin = Enrichment - 2 * Enrichment_std_error, 
                    ymax = Enrichment + 2 * Enrichment_std_error, 
                    width = 0.2), 
                color = "black") + 
  scale_y_continuous(limits = c(-40, 40)) + 
  scale_x_discrete(limits = rev) + 
  scale_fill_manual(values = fill.values) + 
  coord_flip() + 
  labs(title = "COA") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5), 
        legend.title = element_blank(), 
        axis.title = element_text(size = 13), 
        axis.text = element_text(size = 13), 
        legend.spacing.y = unit(0.2, 'cm'), 
        legend.text = element_text(size = 13)) + 
  labs(x = "Cell type")

gg <- gg_aoa + gg_coa + plot_layout(axis_titles = "collect", guides = "collect")

pdf(paste0(plot.dir, '/sldsc_no_gp_enrichment_barplot.pdf'), height = 6, width = 12)
gg
dev.off()

# Find significantly enriched cell types
aoa_nogp_sig_0.05 <- aoa.nogp %>% 
  dplyr::filter(Enrichment_p < 0.05 & Enrichment > 0)

coa_nogp_sig_0.05 <- coa.nogp %>% 
  dplyr::filter(Enrichment_p < 0.05 & Enrichment > 0)

aoa_nogp_sig_0.1 <- aoa.nogp %>% 
  dplyr::filter(Enrichment_p < 0.1 & Enrichment > 0)

coa_nogp_sig_0.1 <- coa.nogp %>% 
  dplyr::filter(Enrichment_p < 0.1 & Enrichment > 0)


