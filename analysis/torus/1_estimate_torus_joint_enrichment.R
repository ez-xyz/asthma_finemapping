rm(list = ls())
library(data.table)
library(optparse)

# process command-line arguments
parser <- OptionParser()
parser <- add_option(parser, c("--trait"), type = "character", default = "")  # trait
out <- parse_args(parser)
trait <- out$trait
cat(paste0("Trait: ", trait, " \n"))
rm(out)

out.dir <- ""  # output directory
torus.path <- ""  # path to TORUS executable
zscore.file <- ""  # path to GWAS z-scores
annot.file <- ""  # path to TORUS annotations
torus.args <- c("-d", zscore.file, "-annot", annot.file, 
                "--load_zval", "-est", "-dump_prior", "prior")

torus.res <- list()
res <- processx::run(command = torus.path, args = torus.args, 
                     echo_cmd = T, echo = T)
enrich <- fread(res$stdout, skip = 1, header = F)
colnames(enrich) <- c("term", "estimate", "low", "high")
torus.res$enrich <- enrich
files <- list.files(path = "prior/", pattern = "*.prior", full.names = T)
files.str <- paste0(files, collapse = " ")
system(paste("cat", files.str, "> prior/allchunks.txt"))
snp.prior <- fread("prior/allchunks.txt", header = F, sep = " ")
colnames(snp.prior) <- c("snp", "torus_prior")
system("rm -rf prior/")
torus.res$snp_prior <- snp.prior

saveRDS(torus.res, file.path(out.dir, paste0(trait, '_joint_enrichment_torus_res.RDS')))

