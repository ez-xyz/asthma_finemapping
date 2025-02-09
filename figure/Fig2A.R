rm(list = ls())
library(data.table)
library(ggplot2)
library(forcats)

sldsc.dir <- ""  # path to S-LDSC h2 enrichment results
plot.dir <- ""  # output directory

l.order <- c("Lymphoid", "Myeloid", "Epithelial", "Mesenchymal", "Endothelial")  # order of cell lineages in S-LDSC result file

aoa.res <- fread(paste0(sldsc.dir, '/aoa_h2_enrichment.results'))
coa.res <- fread(paste0(sldsc.dir, '/coa_h2_enrichment.results'))

aoa.res <- aoa.res[aoa.res$Category %in% paste0("L2_", 1:length(l.order)),]
coa.res <- coa.res[coa.res$Category %in% paste0("L2_", 1:length(l.order)),]

aoa.coa.df <- data.frame(lineage = c(l.order, l.order), 
                         trait = c(rep("AOA", length(l.order)), rep("COA", length(l.order))), 
                         beta = c(aoa.res$Enrichment, coa.res$Enrichment), 
                         se = c(aoa.res$Enrichment_std_error, coa.res$Enrichment_std_error))

aoa.coa.df$lineage <- factor(aoa.coa.df$lineage, levels = l.order)
aoa.coa.df$trait <- factor(aoa.coa.df$trait, levels = c("AOA", "COA"))

gg <- ggplot(data = aoa.coa.df, aes(x = fct_rev(lineage), y = beta, fill = fct_rev(trait))) + 
  geom_col(position = "dodge", alpha = 0.8) + 
  geom_errorbar(aes(ymin = beta - 2 * se, ymax = beta + 2 * se, width = 0.2), 
                position = position_dodge(width = 0.9)) + 
  scale_y_continuous(limits = c(-10, 20)) + 
  scale_fill_manual(values = c("#00B6BC", "#F3BC50"), limits = c("AOA", "COA")) + 
  coord_flip() + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 30), 
        axis.title = element_text(size = 25), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 30), 
        legend.spacing.y = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE)) + 
  labs(x = "", y = "Heritability enrichment")
gg

pdf(paste0(plot.dir, '/sldsc_enrichment_barplot.pdf'), 
    height = 6, width = 8)
gg
dev.off()

