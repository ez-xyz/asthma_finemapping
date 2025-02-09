rm(list = ls())
library(data.table)
library(ggplot2)

susie.dir <- ""  # path to credible set summary
plot.dir <- ""  # output directory

aoa_summary <- fread(paste0(susie.dir, '/aoa_functional_cs_summary.txt'))
coa_summary <- fread(paste0(susie.dir, '/coa_functional_cs_summary.txt'))

n.aoa.specific <- sum(aoa_summary$shared_or_specific == "specific")
n.coa.specific <- sum(coa_summary$shared_or_specific == "specific")
# sanity check: the number of shared credible sets should be the same in two tables
if (sum(aoa_summary$shared_or_specific != "specific") != sum(coa_summary$shared_or_specific != "specific")) {
  stop("The number of shared credible sets are different between AOA and COA.")
}
n.shared <- sum(aoa_summary$shared_or_specific != "specific")

df <- data.frame(type = c("AOA-specific", "COA-specific", "shared"), 
                 count = c(n.aoa.specific, n.coa.specific, n.shared))

df <- df |> dplyr::mutate(proportion = count / sum(count))

fill.values <- c("#00B6BC", "#F3BC50", "#C1A78C")

gg <- ggplot(df, aes(x = "", y = proportion, fill = type)) + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar(theta = "y") + 
  theme_void() + 
  scale_fill_manual(values = fill.values, 
                    name = "") + 
  geom_text(aes(x = 1.1, label = paste0(round(proportion * 100), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 12, color = "white", fontface = "bold") + 
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 30), 
        legend.spacing.y = unit(1.5, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE))
gg

pdf(file.path(plot.dir, "susie_cs_shared_specific_piechart.pdf"), 
    height = 6, width = 8)
gg
dev.off()

