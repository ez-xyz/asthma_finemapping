rm(list = ls())
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(forcats)
library(patchwork)

trait <- "aoa"
plot.dir <- ""  # output directory
susie.dir <- ""  # path to susie cs list and cs summary (from 2_get_cs_info.R)
annot.dir <- ""  # path to annotations (i.e., peaks in different cell lineages)
sldsc.dir <- ""  # path to S-LDSC h2 enrichment result (for calculating attribution weight)
sldsc.order <- c("lymphoid", "myeloid", "epithelial", "mesenchymal", "endothelial")  # order of cell lineages in S-LDSC result file

cs.list <- readRDS(file.path(susie.dir, paste0(trait, "_cs_list.RDS")))
cs.summary <- fread(file.path(susie.dir, paste0(trait, "_functional_cs_summary.txt")))
annot.list <- readRDS(file.path(annot.dir, "annot_list.RDS"))
sldsc.res <- fread(file.path(sldsc.dir, paste0(trait, "_h2_enrichment.results")))

## Helper functions ##
divide_pip <- function(x) {
  if (sum(x != 0) > 1) {  # for ambiguous SNP
    x <- x / sum((x != 0))  # evenly distribute the PIP among cell types
  }
  return(x)
}

weight_pip <- function(x, wt) {
  if (sum(x != 0) > 1) {  # for ambiguous SNP
    nz <- which(x != 0)  # non-zero entries (i.e., overlapping cell types)
    if (!all(wt[nz] == 0)) {  # if all non-zero entries have 0 weight, evenly divide the PIP
      x[nz] <- sum(x) * (wt[nz] / sum(wt[nz]))
    } else {
      x[nz] <- x[nz] / length(nz)
    }
  }
  return(x)
}

## Main function ##
# calculate weight from S-LDSC h2 explained
prop.h2 <- sldsc.res |>
  dplyr::filter(Category %in% paste0("L2_", c(1:length(sldsc.order)))) |>
  dplyr::select(`Prop._h2`) |>
  unlist() |>
  unname()
prop.h2 <- ifelse(prop.h2 > 0, prop.h2, 0)
weight <- prop.h2 / sum(prop.h2)
names(weight) <- sldsc.order

annot.name <- names(annot.list)
weight <- weight[annot.name]  # align weights with annotations
pip.prop <- list()  # a list of data frames; each data frame stores the PIP attribution results for a credible set

for (i in 1:length(cs.list)) {
  cs <- cs.list[[i]]
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = cs$chr, 
                                                           start = cs$pos, 
                                                           end = cs$pos))
  m <- matrix(0, nrow = length(cs.list[[i]]$snp), ncol = length(annot.name))
  colnames(m) <- annot.name
  for (j in annot.name) {
    qh <- findOverlaps(gr, annot.list[[j]])@from  # get the cell type(s) each SNP overlaps with
    m[qh,j] <- cs$pip[qh]
  }
  m <- t(apply(m, 1, divide_pip))
  m <- t(apply(m, 1, weight_pip, wt = weight))
  
  cs.id <- ifelse(is.null(names(cs.list)[i]), paste0("cs", i), names(cs.list)[i])
  
  temp <- data.frame(cs = rep(cs.id, length(annot.name) + 1), 
                     annot = c(annot.name, "none"), 
                     prop = c(unname(colSums(m)), sum(cs$pip) - sum(m)), 
                     cs_ng = paste0(cs.summary$id[i], " (", cs.summary$nearest_gene[i], ")"))
  temp$prop <- ifelse(temp$prop >= 0, temp$prop, 0)  # sometimes "none" can be negative due to numerical issue, and we change it to zero
  temp$prop <- temp$prop / sum(temp$prop)  # ensure the proportions sum to 1
  
  pip.prop[[i]] <- temp
}

df <- dplyr::bind_rows(pip.prop)
df$annot <- factor(df$annot, levels = c("none", annot.name))
df$cs_ng <- factor(df$cs_ng, levels = unique(df$cs_ng))

fill.values <- c("#a9a9a9", "#7FE5F0", "#6395ee", "#e0afff", "#9D33E2", "#541675")

if (trait == "aoa") {
  text.size <- 36
  gg <- ggplot(df, aes(x = forcats::fct_rev(cs_ng), y = prop, fill = annot)) + 
    geom_col(position = "fill", color = "black", alpha = 0.9) + 
    coord_flip() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) + 
    scale_fill_manual(values = fill.values) + 
    labs(x = "Credible set", y = "Proportion of total PIP") + 
    theme_bw() + 
    theme(legend.title = element_blank(), 
          # axis.text.x = element_text(angle = 45, vjust = 0.5), aspect.ratio = 1/3, 
          axis.text.x = element_text(size = text.size, angle = 30), axis.text.y = element_text(size = text.size, face = "italic"), 
          axis.title = element_text(size = text.size), 
          plot.title = element_text(size = 25, hjust = 0.5), 
          legend.text = element_text(size = text.size), 
          legend.spacing.y = unit(0.5, 'cm')) + 
    guides(fill = guide_legend(byrow = TRUE))
}

if (trait == "coa") {
  text.size <- 37
  coa_part1 <- unique(df$cs)[1:22]
  coa_part2 <- unique(df$cs)[23:44]
  coa_part3 <- unique(df$cs)[45:67]
  
  gg_coa_1 <- ggplot(df |> dplyr::filter(cs %in% coa_part1), aes(x = forcats::fct_rev(cs_ng), y = prop, fill = annot)) + 
    geom_col(position = "fill", color = "black", alpha = 0.9) + 
    coord_flip() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0,0)) + 
    scale_fill_manual(values = fill.values) + 
    labs(x = "Credible set", y = "") + 
    theme_bw() + 
    theme(legend.title = element_blank(), 
          # axis.text.x = element_text(angle = 45, vjust = 0.5), aspect.ratio = 1/3, 
          axis.text.x = element_text(size = text.size, angle = 30), axis.text.y = element_text(size = text.size, face = "italic"), 
          axis.title = element_text(size = text.size), 
          plot.title = element_text(size = 25, hjust = 0.5), 
          legend.position = "none") + 
    guides(fill = guide_legend(byrow = TRUE))
  
  gg_coa_2 <- ggplot(df |> dplyr::filter(cs %in% coa_part2), aes(x = forcats::fct_rev(cs_ng), y = prop, fill = annot)) + 
    geom_col(position = "fill", color = "black", alpha = 0.9) + 
    coord_flip() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0,0)) + 
    scale_fill_manual(values = fill.values) + 
    labs(x = "", y = "Proportion of total PIP") + 
    theme_bw() + 
    theme(legend.title = element_blank(), 
          # axis.text.x = element_text(angle = 45, vjust = 0.5), aspect.ratio = 1/3, 
          axis.text.x = element_text(size = text.size, angle = 30), axis.text.y = element_text(size = text.size, face = "italic"), 
          axis.title = element_text(size = text.size), 
          plot.title = element_text(size = 25, hjust = 0.5), 
          legend.position = "none") + 
    guides(fill = guide_legend(byrow = TRUE))
  
  gg_coa_3 <- ggplot(df |> dplyr::filter(cs %in% coa_part3), aes(x = forcats::fct_rev(cs_ng), y = prop, fill = annot)) + 
    geom_col(position = "fill", color = "black", alpha = 0.9) + 
    coord_flip() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0,0)) + 
    scale_fill_manual(values = fill.values) + 
    labs(x = "", y = "") + 
    theme_bw() + 
    theme(legend.title = element_blank(), 
          # axis.text.x = element_text(angle = 45, vjust = 0.5), aspect.ratio = 1/3, 
          axis.text.x = element_text(size = text.size, angle = 30), axis.text.y = element_text(size = text.size, face = "italic"), 
          axis.title = element_text(size = text.size), 
          plot.title = element_text(size = 25, hjust = 0.5), 
          legend.text = element_text(size = text.size), 
          legend.spacing.y = unit(0.5, 'cm')) + 
    guides(fill = guide_legend(byrow = TRUE))
  
  gg <- (gg_coa_1 + theme(plot.margin = unit(c(0,25,0,0), "pt"))) + 
    (gg_coa_2 + theme(plot.margin = unit(c(0,25,0,0), "pt"))) + 
    (gg_coa_3 + theme(plot.margin = unit(c(0,0,0,0), "pt")))
}

if (trait == "aoa") {
  pdf(file.path(plot.dir, paste0(trait, "_cs_cell_type_distribution.pdf")), 
      height = 12, width = 12)
  print(gg)
  dev.off()
}

if (trait == "coa") {
  pdf(file.path(plot.dir, paste0(trait, "_cs_cell_type_distribution.pdf")), 
      height = 14, width = 30)
  print(gg)
  dev.off()
}

