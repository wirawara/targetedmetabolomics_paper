#!/usr/bin/env Rscript --vanilla



### change the path accordingly 
path <- "~/Documents/martina/Kerndl_et_al_2023/"

setwd(path)


## heatmap_AVMKe39a
library(openxlsx)
library(pheatmap)
AVMKe39a <- read.xlsx("./raw/20230116_mke_avmke39a_heatmap_calculations.xlsx")

group <-  AVMKe39a$sample

annot.avmke39a <- data.frame(Group = group)
rownames(annot.avmke39a) <- AVMKe39a$sample


rownames(AVMKe39a) <- AVMKe39a$sample
AVMKe39a <- AVMKe39a[,-c(1:2)]

pdf.options(encoding='ISOLatin2.enc')
pdf("./plots/SupplFig4D_heatmap_AVMKe39a_16Jan2023.pdf", 7, 4)
pheatmap(t(AVMKe39a),
         color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
         border_color = "gray40",
         angle_col = 45,
         fontsize_row = 9,
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 15,
         annotation_col = annot.avmke39a,
         main = "AVMMKe39a"
)
dev.off()




#### ---------------------------
#### in vitro metabolomics


setwd("/Users/andreakomljenovic/Documents/martina/invitrometabolomics/AVMKe08g/")


## supernatant
Ke08g.supernatant <- read.xls("AVMKe08g_supernatant.xlsx")
Ke08g.supernatant.new <-  Ke08g.supernatant[,-which(colnames(Ke08g.supernatant) %in% c("AMP", "NADH", "GlutathioneRED", "NAD", "Ribose5Phosphate", "UDPglucose", "ATP", "X3Phosphoglyceric.acid", "UTP", "Fructose16biphosphate", "Diphosphoglyceric.acid"))]
# because of n.q.


# only 48h
Ke08g.supernatant.new.48h.RPMI <- Ke08g.supernatant.new[9:18,]

Ke08g.supernatant.48h.RPMI <- Ke08g.supernatant.new.48h.RPMI %>%
  dplyr::select(-specimen, -sample, -group) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::select(where(is.numeric))

Ke08g.supernatant.48h.RPMI<- t(Ke08g.supernatant.48h.RPMI)
colnames(Ke08g.supernatant.48h.RPMI) <- Ke08g_supernatant.new.48h.RPMI$group
colnames(Ke08g.supernatant.48h.RPMI)[9:10] <- c(rep("RPMI", 2))

logm.ke08gsupernatant48.RPMI <- log(Ke08g.supernatant.48h.RPMI, 2)
# you can also use vsn here instead of scale
logm.ke08gsupernatant48.RPMI  <- scale(logm.ke08gsupernatant48.RPMI, center = T, scale = T)
logm.ke08gsupernatant48.RPMI[is.nan(logm.ke08gsupernatant48.RPMI)] <- 0
logm.ke08gsupernatant48.RPMI.removed <- logm.ke08gsupernatant48.RPMI[, -which(colnames(logm.ke08gsupernatant48.RPMI) %in% c("KO.2", "RPMI", "RPMI.1"))]
colnames(logm.ke08gsupernatant48.RPMI.removed) <- c(rep("KO",3), rep("WT",4))
saveRDS(logm.ke08gsupernatant48.RPMI.removed, file = "logm.ke08gsupernatant48.RPMI.removed.rds")





setwd("/Users/andreakomljenovic/Documents/martina/invitrometabolomics/AVMKe46c/")

### CSF data - all samples no exclusion
# Ke08g <- read.xls("AVMKe08g_Metabolomics data_Andrea.xlsx")

AVMKe46c.supernatant <- read.xls("AVMKe46c_Metabolomics results_Andrea.xlsx")
colnames(AVMKe46c.supernatant) <- AVMKe46c.supernatant[1,]

AVMKe46c.supernatant <- AVMKe46c.supernatant[-1,]
AVMKe46c.supernatant <- AVMKe46c.supernatant[,1:30]

AVMKe46c.supernatant.new <- AVMKe46c.supernatant[,- which(colnames(AVMKe46c.supernatant) %in% c("Cysteine", "Glutathione reduced", "Serine", "Taurine", "Tryptophan", "Phenylalanine", "Isoleucine", "Leucine"))]



AVMKe46c.supernatant.groups <- AVMKe46c.supernatant.new %>%
  dplyr::select(-Group, -Sample) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::select(where(is.numeric))

# AVMKe46c.supernatant.groups <- as.data.frame(sapply(AVMKe46c.supernatant.groups, as.numeric))
AVMKe46c.supernatant.groups <- t(as.matrix(AVMKe46c.supernatant.groups))
colnames(AVMKe46c.supernatant.groups) <- AVMKe46c.supernatant.new$Group

#normalization for the PCA
logm.supernatant.groups <- log(AVMKe46c.supernatant.groups, 2)
# you can also use vsn here instead of scale
logm.supernatant.groups  <- scale(logm.supernatant.groups, center = T, scale = T)



anova_summary <- function(df, annotation){

  # df - dataframe with the data from the knock out or wt
  # annotation - that is the original dataframe from metabodiff with the metabolites as rownames

  coeff <- sapply(1:nrow(df), function(x) aov(unlist(df[x,])~as.factor(colnames(df)))$coefficient)

  res_time <- sapply(1:nrow(df),
            function(x) summary(aov(unlist(df[x,])~as.factor(colnames(df)))))

  res_df <- data.frame(pval=as.vector(sapply(sapply(res_time,"[",i=5),"[",i=1)),
                          adj_pval=p.adjust(as.vector(sapply(sapply(res_time,"[",i=5),"[",i=1)),method = "fdr"),
                          ## check this still
                          dm=coeff[2,])

  rownames(res_df) <- annotation
  return(res_df)
}


logdf <- as.data.frame(logm.supernatant.groups )
logdf.4groups <- logdf[, c(1:4,13:16,21:24,17:20)]
#logdf.csf <- t(logdf.csf)
colnames(logdf.4groups) <- c(rep("Unstim",4), rep("LI", 4), rep("LIG", 4), rep("LG",4))
anno.logdf.4groups <- rownames(logdf.4groups)

# anova CSF
dfanova.4groups <- anova_summary(logdf.4groups, anno.logdf.4groups)

dfanova.4groups.signif <- dfanova.4groups[dfanova.4groups$adj_pval < 0.05,] 
dim(dfanova.4groups.signif) # 19


library(pheatmap)

logm.supernatant.groups.signif <- logdf.4groups[match(rownames(dfanova.4groups.signif), rownames(logdf.4groups)),]

saveRDS(logm.supernatant.groups.signif, file = "logm.supernatant.groups.signif.rds")


#### --------------------------

avmke08g <- readRDS("./raw/logm.ke08gsupernatant48.RPMI.removed.rds")
avmke46c <- readRDS("./raw/logm.supernatant.groups.signif.rds")


## remove .1 from the rownames

rownames(avmke46c) <- gsub("\\.[0-9]*$", "", rownames(avmke46c))
rownames(avmke08g)[17] <- "Glutathione oxidated"


avmke08g.new <- avmke08g[na.omit(match(rownames(avmke46c), rownames(avmke08g))),] # 11
avmke46c.new <- avmke46c[na.omit(match(rownames(avmke08g.new), rownames(avmke46c))),] # 11


avmke08g.new <- avmke08g.new[, c(4:7,1:3)]


#### pdf for hilic + rp 

colnames(avmke46c.new) <- make.names(colnames(avmke46c.new), unique=TRUE)
Group.hilic.rp <- colnames(avmke46c.new)

annotation.signif.hilic.rp<- data.frame(Group = Group.hilic.rp)
# colnames(Combined.paths.csf)<- colnames(logdf.csf[,c(1:10)])
rownames(annotation.signif.hilic.rp) <- colnames(avmke46c.new)


pdf("./plots/Fig4I_heatmap_hilicrp_signif_ANOVA_AVMKe46c_according_to_wtko.pdf",7,3)
library(pheatmap)
pheatmap(avmke46c.new,
         color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
         border_color = "gray40",
         scale = "row",
         angle_col = 45,
         fontsize_row = 9,
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 15,
         annotation_col =  annotation.signif.hilic.rp,
         # cutree_rows = 3, 
         clustering_method = "average",
         # annotation_row = F,
         # annotation_names_col = F,
         main = "AVMMKe46c according to wt/ko 08g"
)
dev.off()



colnames(avmke08g.new) <- make.names(colnames(avmke08g.new), unique=TRUE)
Group.08g <- colnames(avmke08g.new)

annotation.08g <- data.frame(Group = Group.08g)
# colnames(Combined.paths.csf)<- colnames(logdf.csf[,c(1:10)])
rownames(annotation.08g) <- colnames(avmke08g.new)

pdf("./plots/Fig4_heatmap_AVMKE08g_to_anova_signif_AVMKE46c_09Feb2023.pdf", 7,3)
library(pheatmap)
pheatmap(avmke08g.new,
         color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
         border_color = "gray40",
         scale = "row",
         angle_col = 45,
         fontsize_row = 9,
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 15,
         annotation_col =  annotation.08g,
         # cutree_rows = 3, 
         clustering_method = "average",
         # annotation_row = F,
         # annotation_names_col = F,
         main = "AVMMKe08g to ANOVA signif AVMKe46c"
)
dev.off()















