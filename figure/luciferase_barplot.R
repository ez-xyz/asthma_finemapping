rm(list = ls())
library(openxlsx)
library(ggplot2)
library(ggsignif)
library(matrixStats)

dir <- ""  # path to luciferase result
out.dir <- ""  # output directory
center_snp <- ""  # center SNP in the luciferase construct
target_gene <- ""  # likely target gene of the center SNP
pl <- c(1:3)  # experimental replicates (plates) to analyze
n_haplo <- 2  # number of haplotypes

candidate <- paste0(target_gene, "_", center_snp)  # luciferase candidate
plates <- paste0("plate.", pl)

# read in luciferase result
res <- read.xlsx(paste0(dir, "/", candidate, ".xlsx"))
res <- res |> 
  dplyr::select(1, matches(plates))

ratio <- data.frame(matrix(nrow = nrow(res) / 2, ncol = ncol(res)))  # store LUC/REN ratios
colnames(ratio) <- colnames(res)

ratio[,1] <- unique(str_remove(res[,1], " LUC| REN"))
ratio[n_haplo+1,1] <- "control"
for (i in 1:nrow(ratio)) {
  ratio[i,-1] <- res[2*i-1,-1] / res[2*i,-1]
}

for (i in plates) {
  ratio[[i]] <- rowMeans(ratio |> dplyr::select(ends_with(i)))
}
ratio[["MEAN"]] <- rowMeans(ratio |> dplyr::select(plates))
ratio[["SD"]] <- rowSds(as.matrix(ratio |> dplyr::select(plates)))
ratio[["SEM"]] <- ratio[["SD"]] / sqrt(length(plates))
ratio[["HIGH"]] <- ratio[["MEAN"]] + 2 * ratio[["SEM"]]
ratio[["LOW"]] <- ratio[["MEAN"]] - 2 * ratio[["SEM"]]

df <- ratio |> 
  dplyr::select(Plasmid, plates, MEAN, SD, SEM, HIGH, LOW) |> 
  dplyr::filter(Plasmid != "SV40")

df_norm <- cbind(Plasmid = df$Plasmid, as.data.frame(t(t(df[,-1])/unlist(df[which(df$Plasmid == "control"),-1]))))
df_norm[,2:ncol(df_norm)] <- log2(df_norm[,2:ncol(df_norm)])

df_norm[["MEAN"]] <- rowMeans(df_norm |> dplyr::select(plates))
df_norm[["SD"]] <- rowSds(as.matrix(df_norm |> dplyr::select(plates)))
df_norm[["SEM"]] <- df_norm[["SD"]] / sqrt(length(plates))
df_norm[["HIGH"]] <- df_norm[["MEAN"]] + 2 * df_norm[["SEM"]]
df_norm[["LOW"]] <- df_norm[["MEAN"]] - 2 * df_norm[["SEM"]]

df_norm <- df_norm |> dplyr::filter(Plasmid != "control")

## t-test ##
h1_h2_t <- paste0("p = ", round(t.test(as.numeric(df_norm[1,2:(length(pl)+1)]), as.numeric(df_norm[2,2:(length(pl)+1)]), paired = T)$p.value, digits = 2))

## plot ##
colors <- c("#895129", "#016236")

gg <- ggplot(df_norm) + 
  geom_col(aes(x = Plasmid, y = MEAN, fill = Plasmid), width = 0.5, color = "black", alpha = 0.8) + 
  geom_signif(aes(x = Plasmid, y = MEAN), 
              annotation = h1_h2_t, 
              y_position = 2.7, xmin = 1, xmax = 2,
              tip_length = c(0,0), textsize = 17, vjust = -0.5, size = 1) + 
  geom_errorbar(aes(x = Plasmid, ymin = LOW, ymax = HIGH), width = 0.2, linewidth = 1) + 
  geom_point(data = df_norm |> dplyr::select(Plasmid, plates) |> pivot_longer(cols = plates), 
             aes(x = Plasmid, y = value), size = 5) + 
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'black', linewidth = 1.8) + 
  scale_fill_manual(name = "", values = colors) + 
  labs(x = "Construct", y = "log2 fold change of luciferase activity") + 
  theme_minimal() + 
  theme(axis.title.x = element_text(size = 47, margin = unit(c(5, 0, 0, 0), "mm")), 
        axis.title.y = element_text(size = 47, margin = unit(c(0, 5, 0, 0), "mm")), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 53), 
        axis.line = element_line(color = "black"), 
        legend.text = element_text(size = 47), 
        legend.title = element_blank(), 
        legend.spacing.y = unit(0.7, 'cm')) + 
  guides(fill = guide_legend(byrow = TRUE)) + 
  scale_y_continuous(limits = c(-1, 4.5), expand = c(0, 0))

pdf(file.path(out.dir, paste0(candidate, "_luciferase_barplot.pdf")), 
    height = 12, width = 13)
gg
dev.off()

