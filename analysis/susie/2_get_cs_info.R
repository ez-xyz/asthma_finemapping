rm(list = ls())
## Load packages ##
load_package <- function(pkg) {
  if (!require(pkg, character.only = T, quietly = T)) {
    stop(sprintf("The '%s' package is not installed.", pkg), call. = F)
  }
}
load_package("tidyverse")
load_package("data.table")
load_package("optparse")
load_package("GenomicRanges")

## Process command-line arguments ##
parser <- OptionParser()
parser <- add_option(parser, c("--res"), type = "character", default = "")  # path to fine-mapping result
parser <- add_option(parser, c("--sumstats"), type = "character", default = "")  # path to fine-mapping summary statistics
parser <- add_option(parser, c("--gene-info"), type = "character", default = "")  # path to gene information file
parser <- add_option(parser, c("--outname"), type = "character", default = "")  # output file name prefix
parser <- add_option(parser, c("--outdir"), type = "character", default = "")  # output directory
out <- parse_args(parser)

cat("Fine-mapping result:", out$res, "\n")
cat("Fine-mapping summary statistics:", out$sumstats, "\n")
cat("Gene information:", out$`gene-info`, "\n")
cat("Output file name prefix:", out$outname, "\n")
cat("Output directory:", out$outdir, "\n")

## Input files ##
susie.res <- readRDS(out$res)
susie.sumstats <- readRDS(out$sumstats)
gene.info <- readRDS(out$`gene-info`)
out.name <- out$outname
out.dir <- out$outdir

rm(parser,out)

## Main function ##
#### Part 1: Store credible set info in a list for PIP attribution ####
cs.list <- list()  # a list to store credible set info

for (i in 1:length(susie.res)) {
  locus.id <- names(susie.res)[i]
  dat <- susie.res[[i]]
  gwas <- susie.sumstats %>% dplyr::filter(locus == locus.id)
  if (length(dat$sets$cs) > 0) {
    for (j in 1:length(dat$sets$cs)) {
      cs.id <- paste0(locus.id, "-", dat$sets$cs_index[j])
      cs.list[[cs.id]] <- list()
      cs.list[[cs.id]][["snp"]] <- gwas$snp[dat$sets$cs[[j]]]
      cs.list[[cs.id]][["chr"]] <- gwas$chr[dat$sets$cs[[j]]]
      cs.list[[cs.id]][["pos"]] <- gwas$pos[dat$sets$cs[[j]]]
      cs.list[[cs.id]][["pip"]] <- gwas$susie_pip[dat$sets$cs[[j]]]
    }
  }
}

saveRDS(cs.list, file.path(out.dir, paste0(out.name, "_cs_list.RDS")))

#### Part 2: Summarize credible set info in a table ####
if (!is.null(gene.info)) {
  cat("Gene location file provided, start summarizing credible set information...\n")
  
  # use TSS as gene location; TSS is start for "+" strand and end for "-" strand
  gene.locations <- as.data.frame(gene.info)[, c("seqnames", "start", "end", "gene_name", "strand")]
  gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.info, width = 1))
  gene.locations.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations, 
                                                               seqnames.field = "seqnames", start.field = "start", end.field = "end", 
                                                               keep.extra.columns = T)
  gene.tss.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.locations[,c("seqnames", "tss", "gene_name")], 
                                                         seqnames.field = "seqnames", start.field = "tss", end.field = "tss", 
                                                         keep.extra.columns = T)
  
  df.header <- c("locus", "chr", "start", "end", "span", "size", "content", "pip", "nearest_gene")  # credible set information to be stored in the table
  n.cs <- c()  # number of credible sets at each locus
  for (i in 1:length(susie.res)) {
    n.cs[i] = length(susie.res[[i]]$sets$cs)  # count the number of CSs at each LD block
  }
  df <- data.frame(matrix(nrow = sum(n.cs), ncol = length(df.header)))
  colnames(df) <- df.header
  
  c <- 1  # counter
  for (i in 1:length(susie.res)) {
    if (n.cs[i] != 0) {
      for (j in 1:n.cs[i]) {  # for each credible set
        df$locus[c] = as.numeric(names(susie.res)[i])  # get locus id
        sub = susie.sumstats[which(susie.sumstats$locus == names(susie.res)[i]),]  # get fine-mapped SNPs at the locus
        df$chr[c] = unique(sub$chr)  # get chromosome
        var_id = susie.res[[i]]$sets$cs[[j]]  # get SNPs in credible set
        df$start[c] = susie.sumstats %>% dplyr::filter(snp %in% sub$snp[var_id]) %>% dplyr::select(pos) %>% dplyr::pull() %>% min()  # min credible set SNP position as start
        df$end[c] = susie.sumstats %>% dplyr::filter(snp %in% sub$snp[var_id]) %>% dplyr::select(pos) %>% dplyr::pull() %>% max()  # max credible set SNP position as end
        df$span[c] = df$end[c] - df$start[c] + 1  # length of credible set
        df$size[c] = length(var_id)  # number of SNPs in credible set
        pip_temp = susie.res[[i]]$pip[var_id]  # unsorted PIPs
        var_id_ord = var_id[order(pip_temp, decreasing = T)]  # sort SNPs in CS by PIP
        df$pip[c] = paste(round(susie.res[[i]]$pip[var_id_ord], digits = 6), collapse = ",")  # sorted PIP
        var_id_sorted = sub$snp[var_id_ord]  # sorted SNPs
        df$content[c] = paste(var_id_sorted, collapse = ",")  # SNP rsid, ordered by PIP
        top_df <- susie.sumstats %>% dplyr::filter(snp == var_id_sorted[1])
        top_gr <- GenomicRanges::makeGRangesFromDataFrame(top_df, seqnames.field = "chr", start.field = "pos", end.field = "pos")
        seqlevelsStyle(top_gr) <- seqlevelsStyle(gene.locations.gr) <- seqlevelsStyle(gene.tss.gr) <- "UCSC"
        sh <- subjectHits(findOverlaps(top_gr, gene.locations.gr))
        if (length(sh) > 0) {  # if SNP in gene body, make overlapped gene the nearest gene
          df$nearest_gene[c] <- gene.locations.gr$gene_name[sh][which.min(GenomicRanges::distance(top_gr, gene.tss.gr)[sh])]  # if overlapped multiple genes, use the nearest one
        } else {  # otherwise, find the nearest gene by distance to TSS
          df$nearest_gene[c] <- gene.tss.gr$gene_name[which.min(GenomicRanges::distance(top_gr, gene.tss.gr))]
        }
        c = c + 1
      }
    }
  }
  df$id <- paste0("cs", c(1:nrow(df)))
  fwrite(df, file.path(out.dir, paste0(out.name, "_cs_summary.txt")), 
         col.names = T, row.names = F, quote = F, sep = '\t')
}

