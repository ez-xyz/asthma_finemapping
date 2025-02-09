rm(list = ls())
## Load packages ##
load_package <- function(pkg) {
  if (!suppressPackageStartupMessages(require(pkg, character.only = T, quietly = T))) {
    stop(sprintf("The '%s' package is not installed.", pkg), call. = F)
  }
}
load_package("optparse")
load_package("dplyr")
load_package("stringr")
load_package("ggplot2")
load_package("forcats")
load_package("data.table")
load_package("GenomicRanges")

## Process command-line arguments ##
parser <- OptionParser()
parser <- add_option(parser, c("--pip-cutoff"), type = "double", default = 0.1)  # PIP cutoff for SNPs used to calculate gene score
parser <- add_option(parser, c("--score-cutoff"), type = "double", default = 0.95)  # score cutoff for genes in the heatmap
parser <- add_option(parser, c("--contact"), type = "character", default = "chromatin_contact_linked_genes_files.RDS")  # linked genes by chromatin contacts
parser <- add_option(parser, c("--qtl"), type = "character", default = "qtl_linked_genes_files.RDS")  # linked genes by QTLs
parser <- add_option(parser, c("--exonic"), type = "character", default = "exonic_snp.RDS")  # linked genes by exonic variants
parser <- add_option(parser, c("--nearest"), type = "character", default = "nearest_gene.RDS")  # linked genes by distance (nearest gene)
parser <- add_option(parser, c("--cre"), type = "character", default = "snp_cre_pairs.RDS")  # all snp-cre pairs (from step 1)
parser <- add_option(parser, c("--plot-params"), type = "character", default = "plot_params.RDS")  # ggplot heatmap parameters (for scale_fill_gradient: high, breaks, limits)
parser <- add_option(parser, c("--gene-info"), type = "character", default = "gene_info.RDS")  # gene information
parser <- add_option(parser, c("--outname"), type = "character", default = "aoa")  # output file name prefix
parser <- add_option(parser, c("--outdir"), type = "character", default = "")  # output directory
out <- parse_args(parser)

cat("PIP cutoff for SNPs used to calculate gene score:", out$cutoff, "\n")
cat("Linked genes by chromatin contacts:", out$contact, "\n")
cat("Linked genes by QTLs:", out$qtl, "\n")
cat("Linked genes by exonic variants:", out$exonic, "\n")
cat("Linked genes by distance (nearest gene):", out$nearest, "\n")
cat("CRE-SNP pairs:", out$cre, "\n")
cat("Heatmap parameters:", out$`plot-params`, "\n")
cat("Gene information:", out$`gene-info`, "\n")
cat("Output file name prefix:", out$outname, "\n")
cat("Output directory:", out$outdir, "\n")

## Input files ##
pip.cutoff <- out$`pip-cutoff`
score.cutoff <- out$`score-cutoff`
contact.files <- readRDS(out$contact)
qtl.files <- readRDS(out$qtl)
exonic.link <- readRDS(out$exonic)
nearest.link <- readRDS(out$nearest)
cre.snp <- readRDS(out$cre)
heatmap.params <- readRDS(out$`plot-params`)
gene.info <- readRDS(out$`gene-info`)
out.name <- out$outname
out.dir <- out$outdir
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = T)
}

rm(parser,out)

# read in different types of chromatin contacts and QTLs
if (!is.null(contact.files)) {
  contact.link <- list()
  for (i in 1:nrow(contact.files)) {
    contact.link[[contact.files$data_type[i]]] <- readRDS(contact.files$data_path[i])
  }
}
if (!is.null(qtl.files)) {
  qtl.link <- list()
  for (i in 1:nrow(qtl.files)) {
    qtl.link[[qtl.files$data_type[i]]] <- readRDS(qtl.files$data_path[i])
  }
}

## Main function ##
df.list <- list()

# link non-exonic snp to gene
for (i in cre.snp$snp_id) {
  if (i %in% exonic.link$snp) {  # if the snp is an exonic variant, skip
    next
  }
  snp.pip <- cre.snp$snp_pip[which(cre.snp$snp_id == i)]
  if (snp.pip <= pip.cutoff) {  # if the snp is below PIP cutoff, skip
    next
  }
  snp.gr <- GRanges(seqnames = cre.snp$snp_chr[which(cre.snp$snp_id == i)], 
                    ranges = IRanges(start = cre.snp$snp_end[which(cre.snp$snp_id == i)], 
                                     end = cre.snp$snp_end[which(cre.snp$snp_id == i)]))
  seqlevelsStyle(snp.gr) <- "UCSC"
  gene.overlap <- gene.info$gene_name[subjectHits(findOverlaps(snp.gr, gene.info))]  # overlapped gene(s)
  cre.overlap <- cre.snp |> dplyr::filter(snp_id == i)  # overlapped CRE(s)
  gene.list <- list()  # list of linked genes
  # nearest gene to overlapped CRE(s)
  if (!is.null(nearest.link)) {
    gene.list$Nearest <- nearest.link$nearest_gene[which(nearest.link$cre %in% cre.overlap$cre_id)]
    if (length(gene.overlap) > 0) {
      gene.list$Nearest <- unique(c(gene.list$Nearest, gene.overlap))  # for intragenic variant, also consider the overlapped gene
    }
  }
  # linked genes by different types of chromatin contacts
  if (!is.null(contact.files)) {
    for (j in names(contact.link)) {
      contact.genes <- character(0)
      # identify linked genes separately for each CRE
      for (k in 1:nrow(cre.overlap)) {
        contact.genes <- c(contact.genes, 
                           contact.link[[j]] |>
                             dplyr::filter(cre == cre.overlap$cre_id[k]) |>
                             dplyr::filter(context %in% unlist(str_split(cre.overlap$cre_context[k], ','))) |>
                             dplyr::distinct(linked_genes) |>
                             unlist() |>
                             unname())
      }
      gene.list[[j]] <- unique(contact.genes)
    }
  }
  # linked genes by different types of QTLs
  if (!is.null(qtl.files)) {
    for (j in names(qtl.link)) {
      gene.list[[j]] <- character(0)
      gene.list[[j]] <- qtl.link[[j]] |>
        dplyr::filter(snp == i) |>
        dplyr::filter(context %in% unlist(str_split(cre.overlap$cre_context, ','))) |>
        dplyr::distinct(linked_genes) |>
        unlist() |>
        unname()
    }
  }
  df <- data.frame(snp = rep(i, sum(lengths(gene.list))), 
                   link = rep(names(gene.list), times = lengths(gene.list)), 
                   gene = unname(unlist(gene.list)), 
                   link_score = 1, 
                   pip = snp.pip)
  # update link score
  df.list[[i]] <- df |>
    dplyr::group_by(link) |>
    dplyr::mutate(link_score = link_score / length(link_score)) |>
    dplyr::ungroup()
}

# link exonic snp to gene
if (!is.null(exonic.link)) {
  exonic.link$link <- "Exon"
  exonic.link$link_score <- 1
  colnames(exonic.link) <- c("snp", "pip", "gene", "link", "link_score")
  # update link score
  exonic.link <- exonic.link |>
    dplyr::filter(pip > pip.cutoff) |>
    dplyr::group_by(snp) |>
    dplyr::mutate(link_score = link_score / length(link_score)) |>
    dplyr::ungroup() |>
    dplyr::relocate(snp, link, gene, link_score, pip)
}

# combine all snp-gene links and calculate gene score
score.df <- bind_rows(bind_rows(df.list), exonic.link) |>
  dplyr::group_by(gene) |>
  dplyr::mutate(gene_score = sum(link_score * pip)) |>
  dplyr::mutate(n_snp = length(unique(snp))) |>
  dplyr::group_by(gene, link) |>
  dplyr::mutate(weight = sum(pip * link_score)) |>
  dplyr::ungroup() |>
  dplyr::arrange(desc(gene_score))
gene.df <- score.df |>
  dplyr::filter(gene_score > score.cutoff) |>
  dplyr::distinct(gene)
plot.df <- score.df |>
  dplyr::filter(gene_score > score.cutoff) |>
  dplyr::distinct(gene, n_snp, link, weight, gene_score) |>
  dplyr::mutate(gene = paste0(gene, " (", n_snp, ")"))

# order of links in heatmap
link.order <- c()
if (!is.null(nearest.link)) {
  link.order <- c(link.order, "Nearest")
}
if (!is.null(contact.link)) {
  link.order <- c(link.order, names(contact.link))
}
if (!is.null(qtl.link)) {
  link.order <- c(link.order, names(qtl.link))
}
if (!is.null(exonic.link)) {
  link.order <- c(link.order, "Exon")
}

# some genes may have no support from one or multiple evidence categories - add them back
for (i in unique(plot.df$gene)) {
  zero.link <- setdiff(link.order, plot.df$link[which(plot.df$gene == i)])
  n.zero.link <- length(zero.link)
  if (n.zero.link == 0) {
    next
  }
  supp.df <- data.frame(gene = rep(i, n.zero.link), 
                        n_snp = rep(unique(plot.df$n_snp[which(plot.df$gene == i)]), n.zero.link), 
                        link = zero.link, 
                        weight = rep(0, n.zero.link), 
                        gene_score = rep(unique(plot.df$gene_score[which(plot.df$gene == i)]), n.zero.link))
  plot.df <- bind_rows(plot.df, supp.df)
}
plot.df <- plot.df |> dplyr::arrange(desc(gene_score))
if (nrow(plot.df) != length(link.order) * length(unique(plot.df$gene))) {
  stop("There should be ", length(link.order) * length(unique(plot.df$gene)), " rows in the final data frame.")
}

gg <- ggplot(plot.df, aes(x = factor(link, levels = link.order), 
                     y = forcats::fct_rev(factor(gene, levels = unique(gene))), fill = weight)) + 
  geom_point(shape = 21, stroke = 0, size = 8) + 
  scale_x_discrete(position = "top") + 
  scale_fill_gradient(low = "#f1f1f1", high = heatmap.params$high, 
                      breaks = heatmap.params$breaks, 
                      limits = heatmap.params$limits) + 
  coord_fixed(ratio = 0.3) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        axis.title = element_blank(), 
        axis.text.x.top = element_text(size = 13, angle = 45, vjust = 0.25), 
        axis.text.y = element_text(size = 13, face = "italic"), 
        legend.text = element_text(size = 13), 
        legend.title = element_blank()) + 
  theme(aspect.ratio = 2/1) + 
  guides(fill = guide_colorbar(ticks.colour = NA, breaks = c(0, 1, 2), 
                               title.position = "top", label.position = "right", order = 2))

fwrite(score.df, file.path(out.dir, paste0(out.name, "_gene_score_table.txt")), 
       col.names = T, row.names = F, quote = F, sep = '\t')
fwrite(gene.df, file.path(out.dir, paste0(out.name, "_high_confidence_genes.txt")), 
       col.names = F, row.names = F, quote = F, sep = '\t')
pdf(file.path(out.dir, paste0(out.name, "_gene_score_heatmap.pdf")), 
    height = 8, width = 8)
gg
dev.off()

