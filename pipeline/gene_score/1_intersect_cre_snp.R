rm(list = ls())
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

## Process command-line arguments ##
parser <- OptionParser()
parser <- add_option(parser, c("--bedtools-dir"), type = "character")  # path to bedtools executable
parser <- add_option(parser, c("--cre"), type = "character")  # path to CREs .bed file
parser <- add_option(parser, c("--snp"), type = "character")  # path to credible set SNPs .bed file
parser <- add_option(parser, c("--outname"), type = "character")  # output file name prefix
parser <- add_option(parser, c("--outdir"), type = "character")  # output directory
out <- parse_args(parser)
cat("bedtools executable location:", out$`bedtools-dir`, "\n")
cat("CRE bed file:", out$cre, "\n")
cat("SNP bed file:", out$snp, "\n")
cat("Output file name prefix:", out$outname, "\n")
cat("Output directory:", out$outdir, "\n")
bedtools.dir <- out$`bedtools-dir`
cre.bed <- out$cre
snp.bed <- out$snp
out.name <- out$outname
out.dir <- out$outdir
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = T)
}

rm(parser,out)

system(paste0(file.path(bedtools.dir, "bedtools"), 
              " intersect -a ", cre.bed, 
              " -b ", snp.bed, " -loj > ", 
              file.path(out.dir, "cre_snp_intersect.txt")))  # run bedtools

df <- fread(file.path(out.dir, "cre_snp_intersect.txt"), header = F)
if (ncol(df) != 10) {
  stop("There should be exactly 10 columns in cre_snp_intersect.txt")
}
colnames(df) <- c("cre_chr", "cre_start", "cre_end", "cre_celltype", "cre_context", 
                  "snp_chr", "snp_start", "snp_end", "snp_id", "snp_pip")
snp_cre_df <- df |>
  dplyr::filter(snp_chr != ".") |>  # all SNP-CRE pairs
  dplyr::mutate(cre_id = paste0(cre_chr, "_", cre_start, "_", cre_end)) |>
  dplyr::relocate(cre_id, .before = cre_context)
cre_df <- snp_cre_df |>
  dplyr::distinct(cre_chr, cre_start, cre_end)  # get a list of unique CREs

saveRDS(snp_cre_df, file.path(out.dir, "snp_cre_pairs.RDS"))

