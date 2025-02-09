rm(list = ls())
library(ggplot2)
library(data.table)
library(dplyr)

susie.dir <- ""  # path to credible set summary
plot.dir <- ""  # output directory

df_list <- list()

for (trait in c("aoa", "coa")) {
  summary <- fread(paste0(susie.dir, '/', trait, '_functional_cs_summary.txt'))
  
  df <- data.frame(size = c("1", "2-5", "6-10", "11-20", "20+"), 
                   count = c(length(which(summary$size == 1)), 
                             length(which(summary$size >= 2 & summary$size <= 5)), 
                             length(which(summary$size >= 6 & summary$size <= 10)), 
                             length(which(summary$size >= 11 & summary$size <= 20)), 
                             length(which(summary$size > 20))), 
                   trait = trait)
  
  df$size <- factor(df$size, levels = c("1", "2-5", "6-10", "11-20", "20+"))
  
  df_list[[trait]] <- df
}

gg_df <- bind_rows(df_list)
gg_df$size <- factor(gg_df$size, levels = unique(gg_df$size))
gg_df$trait <- toupper(gg_df$trait)
gg_df <- gg_df |> 
  dplyr::group_by(trait) |> 
  dplyr::mutate(total_count = sum(count)) |> 
  dplyr::mutate(prop = count / total_count) |> 
  ungroup()

gg <- ggplot(gg_df, aes(x = size, y = count, fill = trait)) + 
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) + 
  labs(x = "Number of SNPs in credible set", y = "Number of credible sets") + 
  theme_minimal() + 
  scale_fill_manual(values = c("#00B6BC", "#F3BC50")) + 
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) + 
  theme(axis.text = element_text(size = 30), 
        axis.text.x = element_text(angle = 30), 
        axis.title = element_text(size = 27), 
        legend.text = element_text(size = 30), 
        legend.title = element_blank(), 
        legend.spacing.y = unit(0.6, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE))

gg

pdf(paste0(plot.dir, '/susie_cs_size_barplot.pdf'), 
    height = 6, width = 8)
print(gg)
dev.off()

