rm(list = ls())
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)

data.dir <- ""  # path to CRE summary
out.dir <- ""  # output directory

aoa.cre <- fread(file.path(data.dir, "aoa_cre_summary.txt"))
coa.cre <- fread(file.path(data.dir, "coa_cre_summary.txt"))

aoa.genes <- unique(unlist(str_split(unlist(aoa.cre %>% dplyr::select(nearest_gene, ABC_max_genes, PCHiC_genes, eQTL_genes)), ",")))
aoa.genes <- aoa.genes[-which(aoa.genes == "")]
coa.genes <- unique(unlist(str_split(unlist(coa.cre %>% dplyr::select(nearest_gene, ABC_max_genes, PCHiC_genes, eQTL_genes)), ",")))
coa.genes <- coa.genes[-which(coa.genes == "")]

aoa.count <- aoa.cre %>%
  dplyr::select(nearest_gene, ABC_max_genes, PCHiC_genes, eQTL_genes) %>%
  rowwise() %>%
  mutate(
    n_linked_genes = length(unique(unlist(strsplit(paste(c_across(everything()), collapse = ","), ","))) %>% 
                              .[. != ""])) %>%
  ungroup() %>%
  dplyr::select(n_linked_genes) %>%
  mutate(trait = "AOA") %>% 
  mutate(cre = aoa.cre$id)

coa.count <- coa.cre %>%
  dplyr::select(nearest_gene, ABC_max_genes, PCHiC_genes, eQTL_genes) %>%
  rowwise() %>%
  mutate(
    n_linked_genes = length(unique(unlist(strsplit(paste(c_across(everything()), collapse = ","), ","))) %>% 
                              .[. != ""])) %>%
  ungroup() %>%
  dplyr::select(n_linked_genes) %>%
  mutate(trait = "COA") %>% 
  mutate(cre = coa.cre$id)

df <- dplyr::bind_rows(aoa.count, coa.count) |>
  dplyr::mutate(label = case_when(
    n_linked_genes == 0 ~ "0", 
    n_linked_genes == 1 ~ "1", 
    n_linked_genes >= 2 & n_linked_genes < 5 ~ "2-5", 
    n_linked_genes >= 5 & n_linked_genes < 10 ~ "5-10", 
    n_linked_genes >= 10 ~ ">=10"
  ))

table(df$label[which(df$trait == "AOA")])
table(df$label[which(df$trait == "COA")])
df$label <- factor(df$label, levels = c("1", "2-5", "5-10", ">=10"))

# Define label mapping to expressions
label_map <- c("1" = expression("1"),
               "2-5" = expression("2-5"),
               "5-10" = expression("5-10"),
               ">=10" = expression("">=10))

gg <- ggplot(df, aes(x = label, fill = trait)) + 
  geom_bar(position = "dodge", color = "black", alpha = 0.8) + 
  theme_bw() + 
  labs(x = "Number of genes linked to a CRE", y = "Count") + 
  scale_fill_manual(values = c("#00B6BC", "#F3BC50")) + 
  scale_y_continuous(limits = c(0, 92), breaks = seq(0, 90, 15), expand = c(0, 0)) + 
  theme(axis.text = element_text(size = 25), 
        axis.title = element_text(size = 20), 
        legend.text = element_text(size = 25), 
        legend.title = element_blank(), 
        legend.spacing.y = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE)) + 
  scale_x_discrete(labels = label_map)

sum(df$trait == "AOA" & df$n_linked_genes > 1)  # 53 AOA CREs linked to more than 1 gene
sum(df$trait == "COA" & df$n_linked_genes > 1)  # 118 COA CREs linked to more than 1 gene

pdf(file.path(out.dir, "cre_n_linked_genes_barplot.pdf"), 
    height = 7, width = 7)
gg
dev.off()
