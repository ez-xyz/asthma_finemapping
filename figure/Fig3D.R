rm(list = ls())
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

cre.dir <- ""  # path to snp_cre_pairs file 
out.dir <- ""  # output directory

aoa.cre <- readRDS(file.path(cre.dir, "aoa_snp_cre_pairs.RDS"))
coa.cre <- readRDS(file.path(cre.dir, "coa_snp_cre_pairs.RDS"))
cre <- bind_rows(aoa.cre, coa.cre) |>
  dplyr::distinct(cre_id, cre_context)
context <- str_split(cre$cre_context, ",")
lineage <- lapply(context, function(x) x[!(x == "blood" | x == "lung")])  # convert context to lineage

df <- data.frame(table(lengths(lineage)))
df$prop <- df$Freq / sum(df$Freq)
colnames(df) <- c("n_lineage", "count", "proportion")

fill.values <- c("#7FE5F0", "#6395ee", "#e0afff", "#9D33E2", "#541675")

gg <- ggplot(df, aes(x = "", y = proportion, fill = n_lineage)) + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar(theta = "y") + 
  theme_void() + 
  scale_fill_manual(values = fill.values, 
                    name = "Number of lineages") + 
  theme(legend.title = element_text(size = 22), 
        legend.text = element_text(size = 22), 
        legend.spacing.y = unit(0.2, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE)) + 
  geom_text(aes(x = 1.1, label = paste0(round(proportion * 100), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 7, color = "white", fontface = "bold")
gg

pdf(file.path(out.dir, "cre_n_lineage_piechart.pdf"), 
    height = 6, width = 6)
gg
dev.off()

