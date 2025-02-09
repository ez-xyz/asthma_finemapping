rm(list = ls())
library(data.table)
library(stringr)

aoa_cs <- ""  # aoa credible set summary from step 2
coa_cs <- ""  # coa credible set summary from step 2
aoa_susie_gwas <- ""  # aoa fine-mapping summary statistics
coa_susie_gwas <- ""  # coa fine-mapping summary statistics
aoa_susie_res <- ""  # aoa fine-mapping result
coa_susie_res <- ""  # coa fine-mapping result

sum(table(aoa_cs$locus) > 1)  # 5/16 AOA LD blocks had more than 1 CS
sum(table(coa_cs$locus) > 1)  # 16/48 COA LD blocks had more than 1 CS

# store whether a cs is shared or specific
aoa_label <- c()
coa_label <- c()

for (i in 1:nrow(aoa_cs)) {
  coa_sub <- coa_cs[which(coa_cs$locus == aoa_cs$locus[i]),]  # check COA cs at the same locus
  aoa_snps <- unlist(str_split(aoa_cs$content[i], ","))  # get snps in the AOA cs
  aoa_pips <- as.numeric(unlist(str_split(aoa_cs$pip[i], ",")))  # get pips
  if (nrow(coa_sub) == 0) {  # if no COA cs at the same locus
    aoa_label[i] <- "specific"  # then this AOA cs must be specific
    next
  }
  aoa_label[i] <- "specific"
  for (j in 1:nrow(coa_sub)) {  # otherwise we iterate through all COA cs at the same locus
    coa_snps <- unlist(str_split(coa_sub$content[j], ","))  # get snps in the COA cs
    coa_pips <- as.numeric(unlist(str_split(coa_sub$pip[j], ",")))
    shared_snps <- intersect(aoa_snps, coa_snps)
    shared_snps_aoa_pip <- as.numeric(aoa_susie_gwas$susie_pip[match(shared_snps, aoa_susie_gwas$snp)])
    shared_snps_coa_pip <- as.numeric(coa_susie_gwas$susie_pip[match(shared_snps, coa_susie_gwas$snp)])
    if (length(shared_snps) > floor(length(aoa_snps) / 2) | length(shared_snps) > floor(length(coa_snps) / 2) | sum(shared_snps_aoa_pip) > sum(aoa_pips) / 2 | sum(shared_snps_coa_pip) > sum(coa_pips) / 2) {
      aoa_label[i] <- paste0("shared-", coa_sub$id[j])
    }
  }
}

for (i in 1:nrow(coa_cs)) {
  aoa_sub <- aoa_cs[which(aoa_cs$locus == coa_cs$locus[i]),]  # check AOA cs at the same locus
  coa_snps <- unlist(str_split(coa_cs$content[i], ","))  # get snps in the COA cs
  coa_pips <- as.numeric(unlist(str_split(coa_cs$pip[i], ",")))  # get pips
  if (nrow(aoa_sub) == 0) {  # if no AOA cs at the same locus
    coa_label[i] <- "specific"  # then this COA cs must be specific
    next
  }
  coa_label[i] <- "specific"
  for (j in 1:nrow(aoa_sub)) {  # otherwise we iterate through all AOA cs at the same locus
    aoa_snps <- unlist(str_split(aoa_sub$content[j], ","))  # get snps in the AOA cs
    aoa_pips <- as.numeric(unlist(str_split(aoa_sub$pip[j], ",")))
    shared_snps <- intersect(aoa_snps, coa_snps)
    shared_snps_aoa_pip <- as.numeric(aoa_susie_gwas$susie_pip[match(shared_snps, aoa_susie_gwas$snp)])
    shared_snps_coa_pip <- as.numeric(coa_susie_gwas$susie_pip[match(shared_snps, coa_susie_gwas$snp)])
    if (length(shared_snps) > floor(length(aoa_snps) / 2) | length(shared_snps) > floor(length(coa_snps) / 2) | sum(shared_snps_aoa_pip) > sum(aoa_pips) / 2 | sum(shared_snps_coa_pip) > sum(coa_pips) / 2) {
      coa_label[i] <- paste0("shared-", aoa_sub$id[j])
    }
  }
}

length(aoa_label)  # 21 AOA CS
length(which(aoa_label == "specific"))  # 9 AOA-specific
length(unique(aoa_cs$locus[which(aoa_label == "specific")]))  # spanning 7 loci
length(aoa_label) - length(which(aoa_label == "specific"))  # 12 shared
length(unique(aoa_cs$locus[grep("shared", aoa_label)]))  # spanning 10 loci

length(coa_label)  # 67 COA CS
length(which(coa_label == "specific"))  # 55 COA-specific
length(unique(coa_cs$locus[which(coa_label == "specific")]))  # spanning 45 loci
length(coa_label) - length(which(coa_label == "specific"))  # 12 shared
length(unique(coa_cs$locus[grep("shared", coa_label)]))  # spanning 10 loci

aoa_cs$shared_or_specific <- aoa_label
coa_cs$shared_or_specific <- coa_label

fwrite(aoa_cs, "", 
       col.names = T, row.names = F, quote = F, sep = '\t')
fwrite(coa_cs, "", 
       col.names = T, row.names = F, quote = F, sep = '\t')
