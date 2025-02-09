rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomicRanges))
suppressMessages(library(org.Hs.eg.db))  # match gene ID to gene symbol
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))  # plot genes
suppressMessages(library(GenomicInteractions)) # visualize Hi-C interactions
suppressMessages(library(Gviz))

helper.dir <- ""  # path to helper function
susie.dir <- ""  # path to susie fine-mapping summary statistics
annot.dir <- ""  # path to gene annotation
hic.dir <- ""  # path to PCHi-C data
ocr.dir <- ""  # path to open chromatin data
plot.dir <- ""  # output directory

source(file.path(helper.dir, "trackplot_helper.R"))  # helper function

#### fine-mapping summary statistics ####
trait <- ""
finemapstats <- readRDS(file.path(susie.dir, paste0(trait, "_functional_susie_gwas.RDS")))
finemapstats <- finemapstats |>
  dplyr::select(snp, chr, pos, susie_pip, pval, zscore, CS, locus) |>
  dplyr::rename(cs = CS, pip = susie_pip)
finemapstats.gr <- makeGRangesFromDataFrame(finemapstats, 
                                            seqnames.field = "chr", 
                                            start.field = "pos", 
                                            end.field = "pos", 
                                            keep.extra.columns = T)
finemapstats.gr$chr <- finemapstats$chr
finemapstats.gr$pos <- finemapstats$pos

#### gene track ####
genomic.annots <- readRDS(file.path(annot.dir, 'gencode.v46lift37.basic.annotation.RDS'))
gene.annots <- genomic.annots$genes

#### PCHi-C track ####
blood_hic <- file.path(hic.dir, 'PCHiC_blood.gr')
asmc_hic <- file.path(hic.dir, 'PCHiC_asmc.gr')
bec_hic <- file.path(hic.dir, 'PCHiC_bec.gr')
pcHiC_name <- c(blood_hic, asmc_hic, bec_hic)
pcHiC_list <- lapply(pcHiC_name, readRDS)
names(pcHiC_list) <- c("blood", "asmc", "bec")

for (i in 1:length(pcHiC_list)) {
  seqlevelsStyle(pcHiC_list[[i]]) <- "UCSC"
}

#### ATAC-seq track ####
lym <- file.path(ocr.dir, 'lymphoid.bed')
mye <- file.path(ocr.dir, 'myeloid.bed')
epi <- file.path(ocr.dir, 'epithelial.bed')
mes <- file.path(ocr.dir, 'mesenchymal.bed')
endo <- file.path(ocr.dir, 'endothelial.bed')
ocr_name <- c(lym, mye, epi, mes, endo)
ocr_list <- lapply(ocr_name, rtracklayer::import, format = "BED")
names(ocr_list) <- c("lymphoid", "myeloid", "epithelial", "mesenchymal", "endothelial")

for (i in 1:length(ocr_list)) {
  seqlevelsStyle(ocr_list[[i]]) <- "UCSC"
}

#### make track plot ####
highlight.snps <- ""  # highlighted SNPs in the plot
target.gene <- ""  # likely target gene of the highlighted SNPs

# define plotting region
anchor.snp <- ""
lb <- 100000
ub <- 100000
region.target <- GRanges(seqnames = finemapstats$chr[which(finemapstats$snp == anchor.snp)], IRanges(start = max(finemapstats$pos[which(finemapstats$snp == anchor.snp)] - lb, 0), end = finemapstats$pos[which(finemapstats$snp == anchor.snp)] + ub))  # target region
seqlevelsStyle(region.target) <- "UCSC"

pdf(file.path(plot.dir, paste0(target.gene, '_', anchor.snp, '_', trait, '_trackplot.pdf')), width = 12, height = 12)
finemapping_annot_trackplot(finemapstats.gr, 
                            region.target, 
                            gene.annots, 
                            genome = "hg19", 
                            data_colors = c("red", "green"), 
                            peaks = list("Lymphoid" = ocr_list$lymphoid, 
                                         "Myeloid" = ocr_list$myeloid, 
                                         "Epithelial" = ocr_list$epithelial, 
                                         "Mesenchymal" = ocr_list$mesenchymal, 
                                         "Endothelial" = ocr_list$endothelial), 
                            HiC_loops = list("Blood immune cell PCHi-C" = pcHiC_list$blood, 
                                             "ASMC PCHi-C" = pcHiC_list$asmc, 
                                             "BEC PCHi-C" = pcHiC_list$bec), 
                            highlight_snps = highlight.snps)
dev.off()

