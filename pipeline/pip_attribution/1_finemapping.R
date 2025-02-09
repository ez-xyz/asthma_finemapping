rm(list = ls())
## Load packages ##
load_package <- function(pkg) {
  if (!require(pkg, character.only = T, quietly = T)) {
    stop(sprintf("The '%s' package is not installed.", pkg), call. = F)
  }
}
load_package("tidyverse")
load_package("data.table")
load_package("susieR")
load_package("Rfast")
load_package("optparse")

## Process command-line arguments ##
parser <- OptionParser()
parser <- add_option(parser, c("--gwas"), type = "character", default = "gwas.RDS")  # GWAS sum stats
parser <- add_option(parser, c("--prior"), type = "character", default = "prior.RDS")  # prior
parser <- add_option(parser, c("--locus"), type = "character", default = "finemap_locus.RDS")  # fine-map locus
parser <- add_option(parser, c("--ldr"), type = "character", default = "LD_matrix.RDS")  # LD matrix
parser <- add_option(parser, c("--rvar"), type = "character", default = "LD_snp_info.RDS")  # LD snp information
parser <- add_option(parser, c("--n"), type = "integer", default = 100000)  # GWAS sample size
parser <- add_option(parser, c("--L"), type = "integer", default = 5)  # SuSiE L
parser <- add_option(parser, c("--purity"), type = "double", default = 0.1)  # purity
parser <- add_option(parser, c("--outname"), type = "character", default = "aoa")  # output file name prefix
parser <- add_option(parser, c("--outdir"), type = "character", default = "")  # output directory
out <- parse_args(parser)

cat("GWAS file:", out$gwas, "\n")
cat("Prior file:", out$prior, "\n")
cat("Fine-map locus file:", out$locus, "\n")
cat("LD matrix file:", out$ldr, "\n")
cat("LD snp information file:", out$rvar, "\n")
cat("GWAS sample size:", out$n, "\n")
cat("L:", out$L, "\n")
cat("Purity:", out$purity, "\n")
cat("Output file name prefix:", out$outname, "\n")
cat("Output directory:", out$outdir, "\n")

## Input files ##
gwas <- readRDS(out$gwas)
prior <- readRDS(out$prior)
finemap_locus <- readRDS(out$locus)
ldr <- readRDS(out$ldr)
rvar <- readRDS(out$rvar)
n <- out$n
L <- out$L
purity <- out$purity
out.name <- out$outname
out.dir <- out$outdir

rm(parser,out)

if (!is.null(prior)) {
  gwas$prior <- unname(prior[match(gwas$snp, names(prior))])
}

## Main function ##
susie_res <- list()  # store fine-mapping results
susie_sumstats <- list()  # store fine-mapping summary statistics

for (i in finemap_locus) {
  cat(sprintf("Finemapping locus %s...\t", i))
  
  # subset summary statistics based on LD blocks
  sub.sumstats <- gwas[gwas$locus == i,]
  
  # load ld matrix from .RDS and variant info from .Rvar
  # the order of variants match between them
  ld_matrix <- readRDS(ldr[i])  # LD matrix
  vars <- fread(rvar[i], header = T)  # SNP info (for snp_match)
  colnames(vars) <- c("chr", "id", "pos", "a1", "a0")
  
  if (nrow(sub.sumstats) > 1) {
    # use bigsnpr to index each SNP within the corresponding LD block
    # NOTE: this pipeline doesn't account for strand flip - be cautious when using out of sample LD reference panel
    sumstats.matched <- bigsnpr::snp_match(sub.sumstats, info_snp = vars, strand_flip = F)  # match by chr, a0, a1, and pos
    
    # use SNP index from the `_NUM_ID_` column to extract LD matrix that match with GWAS SNPs
    R <- ld_matrix[sumstats.matched$`_NUM_ID_`, sumstats.matched$`_NUM_ID_`]  # _NUM_ID_ refers to index in vars file
    
    # compute z scores 
    # original z scores from sumstats should not be used as their signs are not corrected when bigsnpr reverses the alleles to match with ones in LD reference
    zhat <- sumstats.matched$beta / sumstats.matched$se
    
    # run susie
    if (!is.null(prior)) {
      res <- susieR::susie_rss(z = zhat, R = R, n = n, L = L, prior_weights = sumstats.matched$prior, verbose = F, min_abs_corr = purity)
    } else {
      res <- susieR::susie_rss(z = zhat, R = R, n = n, L = L, verbose = F, min_abs_corr = purity)
    }
    
    susie_res[[as.character(i)]] <- res
    susie_sumstats[[as.character(i)]] <- sumstats.matched[,c("chr", "pos", "a0", "a1", 
                                                             "beta", "se", "snp", "pval", 
                                                             "zscore", "locus")]
  }
  cat(sprintf("%.0f%% completed.\n", length(susie_res) / length(finemap_locus) * 100))
}

susie_sumstats <- dplyr::bind_rows(susie_sumstats)
susie_sumstats$susie_pip <- 0
susie_sumstats$CS <- 0
loci <- names(susie_res)
for (l in loci) {
  n.snps <- length(susie_res[[l]]$pip)  # number of SNPs fine-mapped for each locus
  susie_sumstats[susie_sumstats$locus == l, "susie_pip"] <- susie_res[[l]]$pip  # PIP of each SNP
  cs.index <- rep(0, n.snps)  # credible set index of each SNP
  if (!is.null(susie_res[[l]]$sets$cs)) {
    for (j in 1:length(susie_res[[l]]$sets$cs_index)) {
      cs.index[susie_res[[l]]$sets$cs[[j]]] <- susie_res[[l]]$sets$cs_index[j]
    }
  }
  susie_sumstats[susie_sumstats$locus == l, "CS"] <- cs.index
}

saveRDS(susie_res, file.path(out.dir, paste0(out.name, '_susie_res.RDS')))
saveRDS(susie_sumstats, file.path(out.dir, paste0(out.name, '_susie_sumstats.RDS')))

