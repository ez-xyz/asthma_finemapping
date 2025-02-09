rm(list = ls())
## Load packages ##
load_package <- function(pkg) {
  if (!require(pkg, character.only = T, quietly = T)) {
    stop(sprintf("The '%s' package is not installed.", pkg), call. = F)
  }
}
load_package("GenomicRanges")
load_package("optparse")
load_package("tidyverse")
load_package("ggsci")

## Process command-line arguments ##
parser <- OptionParser()
parser <- add_option(parser, c("--cs"), type = "character")  # credible sets
parser <- add_option(parser, c("--annot"), type = "character")  # annotations
parser <- add_option(parser, c("--weight"), type = "character")  # weights
parser <- add_option(parser, c("--pt"), type = "integer")  # font size
parser <- add_option(parser, c("--palette"), type = "character")  # palette (for ggsci)
parser <- add_option(parser, c("--outname"), type = "character")  # output file name prefix
parser <- add_option(parser, c("--outdir"), type = "character")  # output directory
out <- parse_args(parser)

cat("Credible sets:", out$cs, "\n")
cat("Annotations:", out$annot, "\n")
cat("Weights:", out$weight, "\n")
cat("Font size:", out$pt, "\n")
cat("Palette:", out$palette, "\n")
cat("Output file name prefix:", out$outname, "\n")
cat("Output directory:", out$outdir, "\n")

## Input files ##
cs.list <- readRDS(out$cs)
annot.list <- readRDS(out$annot)
weight <- readRDS(out$weight)
text.size <- out$pt
palette.name <- out$palette
out.name <- out$outname
out.dir <- out$outdir

rm(parser,out)

## Check input ##
if (length(annot.list) != length(weight)) {
  stop("Every cell type should have one and only one weight.")
}
if (is.null(names(annot.list))) {
  stop("Peaks should be stored in a named list.")
}
if (is.null(names(weight))) {
  stop("Weights should be stored in a named vector.")
}

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
    x[nz] <- sum(x) * (wt[nz] / sum(wt[nz]))
  }
  return(x)
}

# Functions to assist plotting; more details here: https://nanx.me/blog/post/ggplot2-color-interpolation/
pal_ramp <- function(values) {
  force(values)
  function(n) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}

pal_adaptive <- function(name, palette, alpha = 1) {
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db[[name]][[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  pal_ramp(unname(alpha_cols))
}

scale_fill_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("fill", name, pal_adaptive(name, palette, alpha), ...)
}

## Main function ##
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
                     prop = c(unname(colSums(m)), sum(cs$pip) - sum(m)))
  temp$prop <- ifelse(temp$prop >= 0, temp$prop, 0)  # sometimes "none" can be negative due to numerical issue, and we change it to zero
  temp$prop <- temp$prop / sum(temp$prop)  # ensure the proportions sum to 1
  
  pip.prop[[i]] <- temp
}

df <- dplyr::bind_rows(pip.prop)
df$annot <- factor(df$annot, levels = c(annot.name, "none"))
df$cs <- factor(df$cs, levels = unique(df$cs))

gg <- ggplot(df, aes(x = fct_rev(cs), y = prop, fill = annot)) + 
  geom_col(position = "fill", color = "black", alpha = 0.9) + 
  coord_flip() + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) + 
  # scale_fill_manual(values = fill.values) + 
  scale_fill_adaptive(name = palette.name, palette = "default") + 
  labs(x = "Credible set", y = "Proportion of total PIP") + 
  theme_bw() + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = text.size), 
        axis.text.y = element_text(size = text.size), 
        axis.title = element_text(size = text.size), 
        legend.text = element_text(size = text.size), 
        legend.spacing.y = unit(0.5, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE))

pdf(file.path(out.dir, paste0(out.name, "_pip_attribution.pdf")), 
    height = 8, width = 8)
gg
dev.off()

