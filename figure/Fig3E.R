rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(stringr)

mpra.dir <- ""  # path to MPRA data
cre.dir <- ""  # path to snp_cre_pairs file 
plot.dir <- ""  # output directory

#### helper function ####
check_overlap <- function(x, y) {
  s <- unlist(stringr::str_split(x, pattern = ','))
  return(ifelse(length(intersect(s, y)) == 0, "No", "Yes"))
}

#### process MPRA data ####
bc <- fread(file.path(mpra.dir, "barcodes.txt"), header = T)
nrow(bc)  # 96327 barcodes
length(unique(bc$SNP))  # 2034 SNPs tested
bc$id <- paste0(bc$SNP, "_", bc$allele)
length(unique(bc$id))  # 4587 alleles tested, 1558 with 2 alleles, 433 with 3 alleles, 43 with 4 alleles
ts <- unique(bc$SNP)  # all tested SNPs

en3 <- fread(file.path(mpra.dir, "enhancer-3reps-full.txt"), header = T)  # constructs significant in 3/4 replicates
en4 <- fread(file.path(mpra.dir, "enhancer-4reps-full.txt"), header = T)  # constructs significant in 4/4 replicates
es <- unique(c(sub("\\.[ACGT]$", "", en3$allele), 
               sub("\\.[ACGT]$", "", en4$allele)))  # SNPs in MPRA+ enhancers

#### main function ####
df.list <- list()
for (trait in c("aoa", "coa")) {
  cre <- readRDS(file.path(cre.dir, paste0(trait, "_snp_cre_pairs.RDS")))
  cre <- cre |>
    dplyr::group_by(cre_id) |>
    dplyr::mutate(cre_epip = sum(snp_pip)) |>
    dplyr::mutate(overlapped_snp = paste(snp_id, collapse = ",")) |>
    dplyr::ungroup() |>
    dplyr::distinct(cre_id, overlapped_snp, cre_epip)
  cre$mpra_tested <- sapply(cre$overlapped_snp, check_overlap, y = ts)
  cre$mpra_positive <- sapply(cre$overlapped_snp, check_overlap, y = es)
  df.list[[trait]] <- data.frame(epip = c(cre$cre_epip[which(cre$mpra_tested == "Yes" & cre$mpra_positive == "No")], 
                                          cre$cre_epip[which(cre$mpra_tested == "Yes" & cre$mpra_positive == "Yes")]), 
                                 label = c(rep("MPRA-", sum(cre$mpra_tested == "Yes" & cre$mpra_positive == "No")), 
                                           rep("MPRA+", sum(cre$mpra_tested == "Yes" & cre$mpra_positive == "Yes"))), 
                                 trait = toupper(trait))
}

wilcox.test(epip ~ label, data = df.list$aoa, correct = T)  # AOA p = 0.21
wilcox.test(epip ~ label, data = df.list$coa, correct = T)  # COA p = 5.74e-5

aoa.wilcox.p <- "p == 0.21"
coa.wilcox.p <- "p == 5.74 %*% 10^-5"

df <- bind_rows(df.list$aoa, df.list$coa)
gg <- ggplot(df) + 
  geom_boxplot(aes(x = label, y = epip, fill = trait), color = "black", alpha = 0.8) + 
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) + 
  geom_signif(aes(x = label, y = epip), 
              annotation = c(aoa.wilcox.p, coa.wilcox.p), parse = T, 
              y_position = c(1.5,1.75), xmin = c(0.8,1.2), xmax = c(1.8,2.2),
              tip_length = c(0,0), textsize = 11, vjust = -0.1, size = 1) + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 27), 
        axis.text = element_text(size = 32), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 32), 
        legend.spacing.y = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE)) + 
  labs(x = "Candidate CREs", y = "ePIP") + 
  scale_fill_manual(values = c("#00B6BC", "#F3BC50"))
gg

pdf(file.path(plot.dir, "mpra_boxplot.pdf"), height = 8 ,width = 8)
gg
dev.off()

