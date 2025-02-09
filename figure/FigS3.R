rm(list = ls())
library(data.table)
library(tidyverse)
library(susieR)

res.dir <- ""  # path to fine-mapping result
out.dir <- ""  # output directory

gp <- 3

n_cs_list <- list()

for (trait in c("aoa", "coa")) {
  res <- readRDS(paste0(res.dir, '/', trait, '_functional_susie_res.RDS'))
  n_cs <- c()
  for (i in 1:length(res)) {
    n_cs[i] <- length(res[[i]]$sets$cs)
  }
  n_cs_df <- data.frame(table(n_cs))
  colnames(n_cs_df) <- c("n_cs", "count")
  n_cs_df$trait <- trait
  n_cs_df <- n_cs_df[which(n_cs_df$n_cs != 0),]
  n_cs_list[[trait]] <- n_cs_df
}

df <- bind_rows(n_cs_list)
df$n_cs <- as.factor(df$n_cs)
df$trait <- toupper(df$trait)

gg <- ggplot(df, aes(x = n_cs, y = count, fill = trait)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) + 
  theme_bw() + 
  scale_fill_manual(values = c("#00B6BC", "#F3BC50")) + 
  scale_y_continuous(limits = c(0, 37), breaks = seq(0, 35, 5), expand = c(0, 0)) + 
  labs(x = "Number of credible sets at a LD block", y = "Count") + 
  theme(axis.text = element_text(size = 25), 
        axis.title = element_text(size = 20), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 25), 
        legend.spacing.y = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE))

pdf(paste0(out.dir, '/susie_gp', gp, '_cs_count_barplot.pdf'), 
    height = 8, width = 8)
print(gg)
dev.off()
