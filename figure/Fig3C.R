rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)

cre.dir <- ""  # path to snp_cre_pairs file 
out.dir <- ""  # output directory

label.list <- list()  # label each CRE by PIP range
for (trait in c("aoa", "coa")) {
  temp <- readRDS(file.path(cre.dir, paste0(trait, '_snp_cre_pairs.RDS')))
  temp <- temp |>
    dplyr::group_by(cre_id) |>
    dplyr::mutate(cre_epip = sum(snp_pip)) |>  # calculate ePIP
    dplyr::ungroup() |>
    dplyr::distinct(cre_id, cre_epip)
  label.list[[trait]] <- data.frame(label = cut(temp$cre_epip, 
                                                breaks = c(0, 0.1, 0.5, 0.8, 10), 
                                                labels = c("<= 0.1", "0.1-0.5", "0.5-0.8", "> 0.8"), 
                                                include.lowest = F, right = T))
}

df <- bind_rows(label.list$aoa |> dplyr::group_by(label) |> dplyr::count() |> dplyr::mutate(trait = "AOA"), 
                label.list$coa |> dplyr::group_by(label) |> dplyr::count() |> dplyr::mutate(trait = "COA"))

# Define label mapping to expressions
label_map <- c("<= 0.1" = expression(""<=0.1),
               "0.1-0.5" = expression("0.1-0.5"),
               "0.5-0.8" = expression("0.5-0.8"),
               "> 0.8" = expression("">0.8))

gg <- ggplot(df) + 
  geom_col(aes(x = label, y = n, fill = trait), 
           alpha = 0.8, position = "dodge", width = 0.8) + 
  labs(x = "ePIP", y = "Count") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size = 27), 
        axis.text.x = element_text(angle = 30), 
        legend.title = element_blank(), 
        legend.spacing.y = unit(0.5, 'cm'), 
        legend.text = element_text(size = 27)) + 
  guides(fill = guide_legend(byrow = TRUE)) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_fill_manual(values = c("#00B6BC", "#F3BC50")) + 
  scale_x_discrete(labels = label_map)
gg

pdf(file.path(out.dir, "cre_epip_barplot.pdf"), 
    height = 6, width = 8)
print(gg)
dev.off()

