#!/usr/bin/env Rscript --vanilla


# R < scriptName.R --no-save  

######## -------------------------------------------------------
######## -------------------------------------------------------


## load the CSF, SC and SC wt-ko for average effect


library(gdata)
library(reshape2)
library(grid)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(EnrichmentBrowser)
library(mixOmics)
library(dplyr)
library(stringr)
library(org.Mm.eg.db)
library(ggrepel)
library(radiant.data)
library(plyr)
library(pheatmap)
library(factoextra)
library(openxlsx)
library(edgeR)
options(max.print=999999)

### change the path accordingly 
path <- "~/Documents/martina/Kerndl_et_al_2023/"

setwd(path)


### --------------------------------------------
#### 1. Load the data

## ---------
#### CSF

# setwd("/Users/andreakomljenovic/Documents/martina/roko/Arginase_Metabolomics/AVMKe33_healthy_peak_remission_raw_CSF/")

### CSF data - all samples no exclusion
CSF <- read.xls("./raw/AVMKe33_healthy_peak_remission_raw_CSF.xls")
CSF <- CSF[, -which(colnames(CSF) == "Quinolinic.Acid")] # remove duplicates

CSF <- CSF %>%
  arrange(GroupName)

###  cleaning data
CSFNum <- CSF %>%
  dplyr::select(-X, -GroupNum, -GroupName, -UserSampleName, -Replicate) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::select(where(is.numeric))



## --------
#### SC 
SC <- read.xls("./raw/AVMKe33_healthy_peak_remission_raw_spinal_cord_tissue_renamed.xls")
SC <- SC[, -which(colnames(SC) %in% c("Quinolinic.Acid", "SAH", "Taurocholic.acid"))]

SC <- SC %>%
  arrange(GroupName)

SCNum <- SC %>%
  dplyr::select(-X,-GroupNum, -GroupName, -UserSampleName, -Replicate) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::select(where(is.numeric))

rownames(SCNum) <- SC$Replicate
SCNum <- SCNum[-which(rownames(SCNum)=="Peak_3"),]
# apply(as.data.frame(SCNum),2, function(x) is.na(x))
# IP data
IP.SC <- read.xlsx("./raw/AVMKe33a_IP_SpTis_PoolSize.xlsx") 
# NIP data 
NIP.SC <- read.xlsx("./raw/AVMKe33a_NIP_SpTis_PoolSize.xlsx")


###### -------------------------------------
####  2. Preprocessing


### """"""""""""""""""""""
##### CSF data
## arrange according to the health, peak, remission
### log data
logm.csf <- log(CSFNum, 2)

## transform the data
logm.csf <- scale(logm.csf, center = T, scale = T)
rownames(logm.csf) <- CSF$Replicate
t.logm.csf <- t(logm.csf)
# remove NAs
t.logm.csf <- t.logm.csf[complete.cases(t.logm.csf),]
dim(t.logm.csf) 
# 111


### """"""""""""""""""""""
###### SC data
### clean up the ip matrix
mouse.ip.sc <- IP.SC[1:12,]
### clean up the nip matrix
mouse.nip.sc <- NIP.SC

mouse.sc.eae.healthy <- cbind(mouse.ip.sc, mouse.nip.sc)
mouse.sc.eae.healthy <- mouse.sc.eae.healthy[,!duplicated(colnames(mouse.sc.eae.healthy))]

rownames(mouse.sc.eae.healthy) <- mouse.sc.eae.healthy$sample


mouse.sc.eae.healthy <- mouse.sc.eae.healthy %>%
  dplyr::select(-mouse.ID, -sample, -group) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::select(where(is.numeric))


dim(mouse.sc.eae.healthy)
# [1]  12 171
dim(SCNum)
# 12 176

mouse.sc.eae.healthy.matched <- mouse.sc.eae.healthy[, na.omit(match(colnames(SCNum), colnames(mouse.sc.eae.healthy)))]
mouse.sc.eae.healthy.notmatched <- mouse.sc.eae.healthy[, -na.omit(match(colnames(SCNum), colnames(mouse.sc.eae.healthy)))]

SCNum.matched <- SCNum[, na.omit(colnames(mouse.sc.eae.healthy.matched), colnames(SCNum))]

       
mergedSC <- rbind(mouse.sc.eae.healthy.matched, SCNum.matched)
#  24 155 # because I removed 2

### batch correction
library(edgeR)
batch <- c(rep("b1", 12), rep("b2", 12))

#### log the matrix
logm.merged <- log(mergedSC, 2)

## scale the matrix
logm.merged <- scale(logm.merged, center = T, scale = T)

# remove the batch
log.merged.corrected.sc <- removeBatchEffect(t(logm.merged), batch)



####### ----------------------------------------
#### 3. analysis



#####  PCA analysis 

# ------------
# CSF data
logm.csf <- t.logm.csf
colnames(logm.csf) <- gsub("\\_.*","",colnames(logm.csf))


# ------------
# SC data
log.merged.corr.sc <- log.merged.corrected.sc
colnames(log.merged.corr.sc) <- c(rep("Peak",6), rep("Healthy",11), rep("Peak",2), rep("Remission", 5))



# x - logged and scaled matrix of the metabolites
pca.analysis <- function(x, pdf.name, title){
        
        pca.merged <- prcomp(t(x),center = F)
        score.data.merged <- data.frame(pca.merged$x)[,1:2]
        eigs.merged <- pca.merged$sdev^2
        perc.vars.merged <- round(eigs.merged[1:2] / sum(eigs.merged),digits=2)*100
        # colnames of the conditions
        score.data.merged$grouping <- as.vector(colnames(x))


        pdf(pdf.name, 10,10)

        p.merged <- ggplot(score.data.merged, aes(x = PC1, y = PC2)) +
        geom_point(aes(shape = grouping, colour = grouping)) + 
        geom_label_repel(label = score.data.merged$grouping) +
        stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill =grouping)) +
        xlab(paste("PC1:", perc.vars.merged[1], "%")) + 
        ylab(paste("PC2:", perc.vars.merged[2], "%")) +
        ggtitle(title) + theme_classic()

        print(p.merged)

       dev.off()

}


## CSF
pca.analysis(logm.csf, "./plots/metabolomics/SupplFig2B_PCA_CSF.pdf", "PCA Cerebrospinal Fluid")
### spinal cord
pca.analysis(log.merged.corr.sc, "./plots/metabolomics/SupplFig2A_PCA_SC.pdf", "PCA Spinal Cord")




#### density / barplots

#### -----------------------------
## CSF
## density plots 
metabolomics.conditions.csf <- as.data.frame(sapply(seq(4, ncol(t.logm.csf), 5), function(j) rowMeans(t.logm.csf[, j+(-3:1)])))     
colnames(metabolomics.conditions.csf) <- c("Healthy", "Peak", "Remission")

metabolomics.data.avg.csf <- melt(metabolomics.conditions.csf)

#create overlaying density plots
g.final.csf <- ggplot(metabolomics.data.avg.csf, aes(x=value, fill= variable )) +
  geom_density(alpha=.25) + xlab("Metabolites z-scores") + 
           ylab("Density")

pdf("./plots/metabolomics/Fig1B_density_plots_CSF.pdf", 10, 7)
g.final.csf + guides(fill=guide_legend(title="Conditions")) + theme_bw() + ylim(0, 1.5)
dev.off()

metabolomics.data.avg.csf$metabolites <- rownames(metabolomics.conditions.csf)
metabolomics.data.avg.csf$color <- c(rep("#E6DFCA", 111), rep("#7985B2", 111), rep("#ADBDCA", 111))

pdf("./plots/metabolomics/Fig1B_barplots_densityplots_CSF.pdf", 10, 15)
ggplot(metabolomics.data.avg.csf, aes( metabolites,value)) + 
      facet_wrap(~variable, nrow=1) + 
      geom_bar(stat="identity") + 
      scale_fill_manual(values=c("#7985B2", "#ADBDCA", "#E6DFCA")) +
      coord_flip() +   
      theme_bw(base_size=10)+ 
      geom_density(stat = "identity", alpha = 0.3, aes(group = variable, fill = color)) + ylab("z-scores") + guides(fill=guide_legend(title="Conditions"))


dev.off()




### -------------------------------
##### spinal cord
log.merged.corrected.reordered <- log.merged.corrected[, c(7:17, 18:19, 1:6, 20:24)]
###
metabolomics.conditions.healthy <- rowMeans(log.merged.corrected.reordered[, c(1:11)])
metabolomics.conditions.peak <- rowMeans(log.merged.corrected.reordered[,c(12:19)])
metabolomics.conditions.remission <- rowMeans(log.merged.corrected.reordered[,c(20:24)])
metabolomics.conditions.sc <- cbind(metabolomics.conditions.healthy, metabolomics.conditions.peak, metabolomics.conditions.remission) 
colnames(metabolomics.conditions.sc) <- c("Healthy", "Peak", "Remission")
 
metabolomics.data.avg.sc <- melt(metabolomics.conditions.sc)

#create overlaying density plots
g.final <- ggplot(metabolomics.data.avg.sc, aes(x=value, fill= Var2 )) +
  geom_density(alpha=.25) + xlab("Metabolites") + 
           ylab("Density")

pdf("./plots/metabolomics/Fig1B_density_plots_SC.pdf", 10, 7)
g.final + guides(fill=guide_legend(title="Conditions")) + theme_bw()
dev.off()

## with barplots

# adding metabolites
metabolomics.data.avg.sc$metabolites <- rownames(metabolomics.conditions.sc)
metabolomics.data.avg.sc$color <- c(rep("#E6DFCA", 155), rep("#7985B2", 155), rep("#ADBDCA", 155))


pdf("./plots/metabolomics/Fig1B_barplots_densityplots_SC.pdf", 10, 15)
ggplot(metabolomics.data.avg.sc, aes( metabolites,value)) + 
      facet_wrap(~Var2, nrow=1) + 
      geom_bar(stat="identity") + 
      scale_fill_manual(values=c("#7985B2", "#ADBDCA", "#E6DFCA")) +
      coord_flip() + 
      theme_bw(base_size=10)+ 
      geom_density(stat = "identity", alpha = 0.3, aes(group = Var2, fill = color)) + ylab("z-scores") + guides(fill=guide_legend(title="Conditions"))

dev.off()




#### Peak vs Healthy


#### ------------------
## CSF

## T-TEST
pvalPeakvsHealthy.CSF <- apply(t.logm.csf, 1, function(x) {
  t.test(x[6:10], x[1:5])$p.value
})

## adjusted p values from the t-tests pairwise comparison

pvalHealthyPeak.adjusted.csf <- p.adjust(pvalPeakvsHealthy.CSF, method = "BH")
pvalHealthyPeak.adjusted.signif.csf <- pvalHealthyPeak.adjusted.csf[pvalHealthyPeak.adjusted.csf < 0.15]
# 47

write.table(names(pvalHealthyPeak.adjusted.signif.csf), "./results/significant_metabolites_CSF_PeakVSHealthy.txt", row.names= F, quote = F)





##### -----------------
# spinal cord
pvalPeakvsHealthy.SC <- apply(log.merged.corrected.sc, 1, function(x) {
  t.test(x[c(1:6, 18:19)], x[7:17])$p.value
})

## adjusted p values from the t-tests pairwise comparison

pvalHealthyPeak.adjusted.sc <- p.adjust(pvalPeakvsHealthy.SC, method = "BH")
pvalHealthyPeak.adjusted.sc.signif <- pvalHealthyPeak.adjusted.sc[pvalHealthyPeak.adjusted.sc < 0.15]
## 117

write.table(names(pvalHealthyPeak.adjusted.sc.signif), "./results/significant_metabolites_SC_corrected_PeakVSHealthy.txt", row.names= F, quote = F)




### MSEA

MSEAprep <- function(data) {
  
  data$Total <- as.numeric(data$Total)
  data$Hits <- as.numeric(data$Hits)
  data$FDR <- as.numeric(data$FDR)
  
  data <- data %>%
    filter(FDR < 0.1) %>%
    mutate(EnrichmentRatio = Hits / Total)
  
}


#Arginine biosynthesis
arginine <- read.xlsx('./raw/Arginine_metabolism.xlsx')
arginine$ID <- str_trim(arginine$ID, side = 'both')




### -------
### CSF

msea.HealthyPeak.csf <- read.csv("./raw/CSF_HealthyPeak_Enrich_MetaboAnalyst/pathway_results.csv", header = TRUE)
colnames(msea.HealthyPeak.csf)[1] <- "Metabolite.Set" 

msea.HealthyPeak.calcs.csf <- MSEAprep(msea.HealthyPeak.csf)
msea.HealthyPeak.calcs.csf <- msea.HealthyPeak.calcs.csf[order(msea.HealthyPeak.calcs.csf$EnrichmentRatio, decreasing = TRUE),]
msea.HealthyPeak.calcs.csf$Order <- c(1,2,3)

# write.csv(msea.HealthyPeak.calcs.csf, file = "./raw/CSF_HealthyPeak_Enrich_MetaboAnalyst/pathway_results_HealthyPeak_Signif.csv")


pdf("./plots/metabolomics/Fig2A_MSEA_CSF_HealthyPeak.pdf", 6,3)
ggplot(msea.HealthyPeak.calcs.csf, aes(x = EnrichmentRatio, y = reorder(Metabolite.Set, -Order), fill = FDR)) +
   geom_segment(aes(xend=0, yend=Metabolite.Set)) +
   geom_point(shape=21,, size = 4) +
    scale_fill_gradient(
     low =  "#D6604D",
     high = "#4393C3", limits = c(0, 0.10),
     aesthetics = "fill") +
      labs(y = "") 
dev.off()


hits.msea.healthypeak.csf <- read.xlsx('./raw/Signif_CSF_HeatlhyPeak_Metaboanalyst_Metabolites.xlsx')
logdf.csf <- as.data.frame(t.logm.csf)
logdf.csf.healthypeak <- logdf.csf[,c(1:10)]
colnames(logdf.csf.healthypeak) <- c(rep("Healthy",5), rep("Peak", 5))
anno.logdf.csf.healthypeak  <- rownames(logdf.csf.healthypeak)

PeakHealthy.arginineMSEA.CSF <- hits.msea.healthypeak.csf[hits.msea.healthypeak.csf$ID %in% arginine$ID,]
PeakHealthy.arginine.logm.csf <- logdf.csf.healthypeak[rownames(logdf.csf.healthypeak) %in% PeakHealthy.arginineMSEA.CSF$Metabolites_original,]


Combined.paths.csf <- PeakHealthy.arginine.logm.csf

annotation_row <- data.frame(
  Pathway = factor(rep(c( "Arginine biosynthesis"), 3))
)
nams <- rownames(Combined.paths.csf)
rownames(annotation_row) <- make.names(nams, unique=TRUE)
rownames(Combined.paths.csf) <- make.names(nams, unique=TRUE)
rownames(annotation_row) <- rownames(Combined.paths.csf)

Group <- colnames(logdf.csf.healthypeak)

annot <- data.frame(Group = Group)
colnames(Combined.paths.csf)<- colnames(logdf.csf[,c(1:10)])
rownames(annot) <- colnames(logdf.csf[,c(1:10)])

pdf("./plots/metabolomics/Fig2B_MSEA_heatmaps_CSF_ArginineBiosynthesis.pdf", 15,1.75)
pheatmap(Combined.paths.csf,
         color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
         border_color = "gray40",
         scale = "row",
         angle_col = 45,
         fontsize_row = 9,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 15,
         annotation_col =  annot,
         annotation_row = annotation_row,
         annotation_names_col = F,
         main = "AVMMKe33 Healthy vs Peak CSF KEGG Significant Pathways Enrichment Arginine Biosynth"
)
dev.off()





#### -----
# spinal cord
msea.HealthyPeak.sc <- read.csv("./raw/SC_HealthyPeak_Enrich_MetaboAnalyst/pathway_results.csv", header = TRUE)

colnames(msea.HealthyPeak.sc)[1] <- "Metabolite.Set" 

msea.HealthyPeak.calcs.sc<- MSEAprep(msea.HealthyPeak.sc)
msea.HealthyPeak.calcs.sc <- msea.HealthyPeak.calcs.sc[order(msea.HealthyPeak.calcs.sc$EnrichmentRatio, decreasing = TRUE),]
msea.HealthyPeak.calcs.sc$Order <- c(1,2,3,4,5,6,7,8,9,10,11)

# write.csv(msea.HealthyPeak.calcs.sc, file = "./raw/SC_HealthyPeak_Enrich_MetaboAnalyst/pathway_results_HealthyPeak_Signif_SC.csv")


### lollipop plot
pdf("./plots/metabolomics/Fig2A_MSEA_SC_HealthyPeak.pdf", 6,3)

ggplot(msea.HealthyPeak.calcs.sc, aes(x = EnrichmentRatio, y = reorder(Metabolite.Set, -Order), fill = FDR)) +
   geom_segment(aes(xend=0, yend=Metabolite.Set)) +
   geom_point(shape=21,, size = 4) +
    scale_fill_gradient(
     low =  "#D6604D",
     high = "#4393C3", limits = c(0, 0.10),
     aesthetics = "fill") +
      labs(y = "") 
  dev.off()




####### heatmaps of the arginine biosynthesis

library(stringr)
library(openxlsx)


hits.msea.healthypeak.sc <- read.xlsx('./raw/Signif_SC_HealthPeak_Metaboanalyst_Metabolites.xlsx')


#Arginine biosynthesis

PeakHealthy.arginineMSEA.SC <- hits.msea.healthypeak.sc[hits.msea.healthypeak.sc$ID %in% arginine$ID,]


## prepare the log scale metabolomics matrix
logdf.sc.corr <- as.data.frame(log.merged.corrected.sc)
logdf.sc.corr.healthypeak <- logdf.sc.corr[,c(7:17,1:6,18:19)]
colnames(logdf.sc.corr.healthypeak ) <- c(rep("Healthy",11), rep("Peak", 8))
anno.logdf.sc.corr.healthypeak  <- rownames(logdf.sc.corr.healthypeak)

PeakHealthy.arging.logm.sc <- logdf.sc.corr.healthypeak[rownames(logdf.sc.corr.healthypeak) %in% PeakHealthy.arginineMSEA.SC$Metabolites_original,]
Combined.paths.sc <- PeakHealthy.arging.logm.sc

annotation.row.sc <- data.frame(
  Pathway = factor(rep("Arginine biosynthesis", 6))
)
nams.sc <- rownames(Combined.paths.sc)
rownames(annotation.row.sc) <- make.names(nams.sc, unique=TRUE)
rownames(Combined.paths.sc) <- make.names(nams.sc, unique=TRUE)
rownames(annotation.row.sc) <- rownames(Combined.paths.sc)

Group <- colnames(logdf.sc.corr.healthypeak)

annot.sc <- data.frame(Group = Group)
colnames(Combined.paths.sc)<- colnames(logdf.sc.corr[,c(7:17,1:6,18:19)])
rownames(annot.sc) <- colnames(logdf.sc.corr[,c(7:17,1:6,18:19)])

pdf("./plots/metabolomics/Fig2B_MSEA_heatmaps_SC_ArginineBiosynthesis.pdf", 15,3)
pheatmap(Combined.paths.sc,
         color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
         border_color = "gray40",
         scale = "row",
         angle_col = 45,
         fontsize_row = 9,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 15,
         annotation_col =  annot.sc,
         annotation_row = annotation.row.sc,
         annotation_names_col = F,
         main = "AVMMKe33 Healthy vs Peak SC KEGG Significant Pathway Arginine Enrichment"
)
dev.off()





#### clusters



cluster.analysis <- function(averaged.conditions.metabolomics, tissue, pdf.file.motherline.name){

          library(Mfuzz)
          library(reshape2)
          library(tidyr)
          library(dplyr)
          library(patchwork)
          library(openxlsx)


          # setwd(path)

          metabolomics.eset <- ExpressionSet(assayData = averaged.conditions.metabolomics)
          # standardize matrix of metabolomics data
          z.mat <- standardise(metabolomics.eset)


          mestimate <- function(df){
                N <-  dim(df)[[1]]
                D <- dim(df)[[2]]
                m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
                return(m.sj)
            }


          mmm <- mestimate(z.mat)

          # 3 clusters
          clustering <- mfuzz(z.mat, centers=3, m=mmm)
 

          #get the centroids into a long dataframe:
          fcm.centroids <- clustering$centers
          fcm.centroids.df <- data.frame(fcm.centroids)
          fcm.centroids.df$cluster <- row.names(fcm.centroids.df)
          centroids.long <- tidyr::gather(fcm.centroids.df,"sample",'value',1:3)


          #start with the input data
          fcm.plotting.df <- data.frame(z.mat)
          fcm.plotting.df <- t(fcm.plotting.df)
          fcm.plotting.df <- as.data.frame(fcm.plotting.df)


          #add genes
          fcm.plotting.df$gene <- rownames(fcm.plotting.df)

          #bind cluster assignment
          fcm.plotting.df$cluster <- clustering$cluster
          #fetch the membership for each gene/top scoring cluster
          fcm.plotting.df$membership <- sapply(1:length(fcm.plotting.df$cluster),function(row){
              clust <- fcm.plotting.df$cluster[row]
              clustering$membership[row,clust] 
            })


          k.to.plot <- 1

              #subset the dataframe by the cluster and get it into long form
              #using a little tidyr action
              cluster.plot.df <- dplyr::filter(fcm.plotting.df, cluster == k.to.plot) %>%
                  dplyr::select(.,1:3,membership,gene) %>%
                  tidyr::gather(.,"sample",'value',1:3)

            #order the dataframe by score
            cluster.plot.df <- cluster.plot.df[order(cluster.plot.df$membership),]
            #set the order by setting the factors using forcats
            cluster.plot.df$gene <- forcats::fct_inorder(cluster.plot.df$gene)

            #subset the cores by cluster
            core <- dplyr::filter(centroids.long, cluster == k.to.plot)

            p1 <- ggplot(cluster.plot.df, aes(x=sample,y=value)) + 
                geom_line(aes(colour=membership, group=gene)) +
                # scale_colour_gradientn(colours=c('blue1','red2')) +
                scale_colour_gradient(low = "white", high = "gray80", limits =c(0,1.00)) +
                #this adds the core 
                geom_line(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
                xlab("Time") +
                ylab("Standardized metabolite z-score") +
                labs(title= paste0("Cluster ",k.to.plot,"Pattern by Time"),color = "Degree") +
                theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


            k.to.plot2 <- 2

              #subset the dataframe by the cluster and get it into long form
              #using a little tidyr action
              cluster.plot.df2 <- dplyr::filter(fcm.plotting.df, cluster == k.to.plot2) %>%
                dplyr::select(.,1:3,membership,gene) %>%
                tidyr::gather(.,"sample",'value',1:3)

              #order the dataframe by score
              cluster.plot.df2 <- cluster.plot.df2[order(cluster.plot.df2$membership),]
              #set the order by setting the factors using forcats
              cluster.plot.df2$gene <- forcats::fct_inorder(cluster.plot.df2$gene)

              #subset the cores by cluster
              core2 <- dplyr::filter(centroids.long, cluster == k.to.plot2)

              p2 <- ggplot(cluster.plot.df2, aes(x=sample,y=value)) + 
                  geom_line(aes(colour=membership, group=gene)) +
                  # scale_colour_gradientn(colours=c('blue1','red2')) +
                  scale_colour_gradient(low = "white", high = "gray80",limits =c(0,1.00)) +
                  #this adds the core 
                  geom_line(data=core2, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
                  xlab("Time") +
                  ylab("Standardized metabolite z-score") +
                  labs(title= paste0("Cluster ",k.to.plot2,"Pattern by Time"),color = "Degree") +
                  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


              k.to.plot3 <- 3

              #subset the dataframe by the cluster and get it into long form
              #using a little tidyr action
              cluster.plot.df3 <- dplyr::filter(fcm.plotting.df, cluster == k.to.plot3) %>%
                dplyr::select(.,1:3,membership,gene) %>%
                tidyr::gather(.,"sample",'value',1:3)

              #order the dataframe by score
              cluster.plot.df3 <- cluster.plot.df3[order(cluster.plot.df3$membership),]
              #set the order by setting the factors using forcats
              cluster.plot.df3$gene <- forcats::fct_inorder(cluster.plot.df3$gene)

              #subset the cores by cluster
              core3 <- dplyr::filter(centroids.long, cluster == k.to.plot3)

              p3 <- ggplot(cluster.plot.df3 , aes(x=sample,y=value)) + 
                  geom_line(aes(colour=membership, group=gene)) +
                  # scale_colour_gradientn(colours=c('blue1','red2')) +
                  scale_colour_gradient(low = "white", high = "gray80",limits =c(0,1.00)) +
                  #this adds the core 
                  geom_line(data=core3, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
                  xlab("Time") +
                  ylab("Standardized metabolite z-score") +
                  labs(title= paste0("Cluster ",k.to.plot3,"Pattern by Time"),color = "Degree") +
                  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




                  ### plotting the clusters
                  pdf(pdf.file.motherline.name, 20, 7)
                  print(p1)+print(p2)+print(p3)
                  dev.off()


                  #### write the metabolites within the cluster
                  for(i in 1:3){
                      write.table(names(clustering$cluster[which(clustering$cluster == i)]), paste0("./results/", "cluster", tissue, i, "_metabolites_fuzzy_averaged.txt"), append = FALSE, sep = " ", 
                              row.names = FALSE, quote = FALSE, col.names = FALSE)
                  }



                  ## print out the membership of metabolites (gene = metabolites for the annotation only) 
                  printing.cl1 <- cluster.plot.df[,c("membership", "gene", "sample")]
                  write.xlsx(x = printing.cl1, paste0("./results/cluster1", tissue, "_info_averaged.xlsx"),rowNames = TRUE)

                  printing.cl2 <- cluster.plot.df2[,c("membership", "gene", "sample")]
                  write.xlsx(x = printing.cl2, paste0("./results/cluster2", tissue, "_info_averaged.xlsx"),rowNames = TRUE)

                  printing.cl3 <- cluster.plot.df3[,c("membership", "gene", "sample")]
                  write.xlsx(x = printing.cl3, paste0("./results/cluster3", tissue, "_info_averaged.xlsx"),rowNames = TRUE)



    return(clustering)

}

### CSF 
t.logm.csf.healthy <- rowMeans(t.logm.csf[,1:5])
t.logm.csf.peak <- rowMeans(t.logm.csf[,6:10])
t.logm.csf.remission <- rowMeans(t.logm.csf[,11:15])

t.logm.csf.av <- cbind(t.logm.csf.healthy, t.logm.csf.peak , t.logm.csf.remission)

## clustering analysis
cluster.analysis(t.logm.csf.av, "CSF", "./plots/metabolomics/Fig1C_CSF_clusters.pdf")



##### SC
log.merged.corrected.reordered <- log.merged.corrected.sc[, c(7:17, 18:19, 1:6, 20:24)]

metabolomics.conditions.healthy <- rowMeans(log.merged.corrected.reordered[, c(1:11)])
metabolomics.conditions.peak <- rowMeans(log.merged.corrected.reordered[,c(12:19)])
metabolomics.conditions.remission <- rowMeans(log.merged.corrected.reordered[,c(20:24)])
metabolomics.conditions.sc <- cbind(metabolomics.conditions.healthy, metabolomics.conditions.peak, metabolomics.conditions.remission) 

# clustering analysis
cluster.analysis(metabolomics.conditions.sc, "SC","./plots/metabolomics/Fig1C_SC_clusters.pdf")




#### DAMs, HAMs and RAMs

statistical.framework.dams.hams.rams <- function(logdf.tissue, anno.logdf.tissue, fdr, tissue){

      library(dplyr)


        ### function for anova analyysis
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



          ### anova calculations
          dfanova.tissue <- anova_summary(logdf.tissue, anno.logdf.tissue)
          print(dfanova.tissue)
          dfanova.tissue.signif.v3 <- dfanova.tissue[dfanova.tissue$adj_pval < fdr,] 
          print(dim(dfanova.tissue.signif.v3) )



          #### tukey p-values
          tukeyhsdf.tissue <- sapply(1:nrow(logdf.tissue),
                function(x) TukeyHSD(aov(unlist(logdf.tissue[x,])~as.factor(colnames(logdf.tissue)))))


          names(tukeyhsdf.tissue) <- rownames(dfanova.tissue)
          ## only significant ones
          tukeyhsdf.tissue.signif <- tukeyhsdf.tissue[na.omit(match(rownames(dfanova.tissue.signif.v3), names(tukeyhsdf.tissue)))] 





          ###### DISEASE METABOLITES


          ## if significant at healthy-peak, and healthy-remission = DISEASE
          disease.metabolites <- lapply(tukeyhsdf.tissue.signif, function(x) which(x[1,4] < 0.05 && x[2,4] < 0.05))
          diseases.metab <- lapply(disease.metabolites, function(x) length(x) == 1 )
          disease.metab <- as.data.frame(do.call(rbind, diseases.metab ))

          dis.metab <- filter(disease.metab, V1 == TRUE) # 57 out of 113


          ### take log matrix and select disease metabolites
          disease.metabolites.logm <- logdf.tissue[na.omit(match(rownames(dis.metab), rownames(logdf.tissue))),] # 55
          colnames(disease.metabolites.logm) <- make.names(colnames(disease.metabolites.logm), unique=TRUE)

          group.tissue.disease <- colnames(disease.metabolites.logm )

          annot.tissue.disease <- data.frame(Group = group.tissue.disease)
          rownames(annot.tissue.disease ) <- colnames(disease.metabolites.logm )

          # setwd("/Users/andreakomljenovic/Documents/martina/new_data_spinal_cord/")


          # pdf("disease_metabolites_55_metab_AVMKE33_HealthyPeakRemission_anova_SC_29March2023_v3.pdf", 10, 13.5)
          pdf(paste0("./plots/metabolomics/Fig1",tissue, "_disease_metabolites.pdf"), 10, 30)

          print(pheatmap(disease.metabolites.logm,
              color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
              border_color = "gray40",
              scale = "row",
              angle_col = 45,
              fontsize_row = 9,
              cluster_cols = F,
              cluster_rows = T,
              cellwidth = 15,
              annotation_col =  annot.tissue.disease,
              # annotation_row = F,
              # annotation_names_col = F,
         main = "AVMMKe33 Healthy vs Peak vs Remission ANOVA DISEASE METABOLITES"
        ))
        dev.off()




          ###### HOMEOSTASIS METABOLITES

          homeostasis.metabolites <- lapply(tukeyhsdf.tissue.signif, function(x) which(x[1,4] < 0.05 && x[2,4] > 0.05 && x[3,4] < 0.05))
          homeostasis.metab <- lapply(homeostasis.metabolites, function(x) length(x) == 1 )
          homeostasis.metab <- as.data.frame(do.call(rbind, homeostasis.metab ))
          
          homeo.metab <- filter(homeostasis.metab, V1 == TRUE) 


          homeostasis.metabolites.logm <- logdf.tissue[na.omit(match(rownames(homeo.metab), rownames(logdf.tissue))),] # 52 
          colnames(homeostasis.metabolites.logm) <- make.names(colnames(homeostasis.metabolites.logm), unique=TRUE)

          group.tissue.homeostasis <- colnames(homeostasis.metabolites.logm)

          annot.tissue.homeostasis <- data.frame(Group = group.tissue.disease)
          rownames(annot.tissue.homeostasis) <- colnames(homeostasis.metabolites.logm)

        pdf(paste0("./plots/metabolomics/Fig1",tissue, "_homeostasis_metabolites.pdf"), 20, 3)

          print(pheatmap(homeostasis.metabolites.logm,
                color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
                border_color = "gray40",
                scale = "row",
                angle_col = 45,
                fontsize_row = 9,
                cluster_cols = F,
                cluster_rows = T,
                cellwidth = 15,
                annotation_col =  annot.tissue.homeostasis,
                # annotation_row = F,
                # annotation_names_col = F,
          main = "AVMMKe33 Healthy vs Peak vs Remission ANOVA HOMEOST METABOLITES"
        ))
        dev.off()



        #### RECOVERY METABOLITES
        recovery.metabolites <- lapply(tukeyhsdf.tissue.signif, function(x) which(x[1,4] > 0.05 && x[2,4] < 0.05 && x[3,4] < 0.05))
        recovery.metab <- lapply(recovery.metabolites, function(x) length(x) == 1 )
        recovery.metab <- as.data.frame(do.call(rbind, recovery.metab ))
  
        reco.metab <- filter(recovery.metab, V1 == TRUE) 

        recovery.metabolites.logm <- logdf.tissue[na.omit(match(rownames(reco.metab), rownames(logdf.tissue))),] # 52
        colnames(recovery.metabolites.logm) <- make.names(colnames(recovery.metabolites.logm), unique=TRUE)

        group.tissue.recovery <- colnames(recovery.metabolites.logm)      

        annot.tissue.recovery <- data.frame(Group = group.tissue.recovery)
        # colnames(Combined.paths.csf)<- colnames(logdf.csf[,c(1:10)])
        rownames(annot.tissue.recovery) <- colnames(recovery.metabolites.logm)

        # pdf("recovery_metabolites_1_metab_AVMKE33_HealthyPeakRemission_anova_SC_29March2023_v3.pdf",15, 1.5)
        pdf(paste0("./plots/metabolomics/Fig1",tissue, "_recovery_metabolites.pdf"),15, 3)

        print(pheatmap(recovery.metabolites.logm ,
                 color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
                 border_color = "gray40",
                 scale = "row",
                 angle_col = 45,
                 fontsize_row = 9,
                 cluster_cols = F,
                 cluster_rows = F,
                 cellwidth = 15,
                 annotation_col =  annot.tissue.recovery,
                 # annotation_row = F,
                 # annotation_names_col = F,
                 main = "AVMMKe33 Healthy vs Peak vs Remission ANOVA RECOVERY METABOLITES",
                #legend=T, legend_breaks = -3:3, legend_labels = c(-3,-2,-1,0,1,2,3)  
        ))
         dev.off()


        ### WEAK METABOLITES

        names.80 <- c(rownames(reco.metab), rownames(homeo.metab), rownames(dis.metab))
        weak.metabolites <- tukeyhsdf.tissue.signif[!names(tukeyhsdf.tissue.signif) %in% names.80]


      weak.metabolites <- lapply(weak.metabolites, function(x) as.data.frame(x))
      weak.metabolites<- lapply(weak.metabolites, function(x) cbind(" "=rownames(x), x))
      # writexl::write_xlsx(weak.metabolites, "weak_metabolites_SC_16April2023_removedmetabol.xlsx")


      weak.metabolites.logm <- logdf.tissue[na.omit(match(names(weak.metabolites), rownames(logdf.tissue))),] # 47
      colnames(weak.metabolites.logm) <- make.names(colnames(weak.metabolites.logm), unique=TRUE)

      group.tissue.weak <- colnames(weak.metabolites.logm)

      annot.tissue.weak <- data.frame(Group = group.tissue.weak)
      rownames(annot.tissue.weak) <- colnames(weak.metabolites.logm)


      pdf(paste0("./plots/metabolomics/Fig1",tissue, "_weak_metabolites.pdf"), 10, 10)

      print(pheatmap(weak.metabolites.logm ,
         color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
         border_color = "gray40",
         scale = "row",
         angle_col = 45,
         fontsize_row = 9,
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 15,
         annotation_col =  annot.tissue.weak,
         # annotation_row = F,s
         # annotation_names_col = F,
         main = "AVMMKe33 Healthy vs Peak vs Remission ANOVA WEAK METABOLITES"
      ))
    dev.off()


     return(list()) 


}

### CSF
colnames(logdf.csf) <- c(rep("Healthy",5), rep("Peak", 5), rep("Remission", 5))
anno.logdf.csf <- rownames(logdf.csf)
statistical.framework.dams.hams.rams(logdf.csf, anno.logdf.csf, 0.15, "CSF")




##### SC
logdf.sc <- as.data.frame(log.merged.corrected.reordered)
colnames(logdf.sc) <- c(rep("Healthy",11), rep("Peak", 8), rep("Remission", 5))
anno.logdf.sc <- rownames(logdf.sc)
### finding disease and weak metabolites
statistical.framework.dams.hams.rams(logdf.sc, anno.logdf.sc, 0.15, "SC")



# -------------------------------------
### Average Effect analysis

### calculating the differences in means for the Healthy - Peak
means.sc <- as.data.frame(rowMeans(log.merged.corrected.reordered[, 12:19])) # Peak
means.sc[,2] <- rowMeans(log.merged.corrected.reordered[, 1:11]) # Healthy
colnames(means.sc) <- c("Peak", "Healthy")

means.sc[,3] <- means.sc["Peak"] - means.sc["Healthy"]
colnames(means.sc)[3] <- "DiffInMeans"


#### WT - KO 

IP.SC.wtko <- read.xls("./raw/AVMKe33a_IP_SpTis_PoolSize_wtko.xlsx")

### exclusion of 465, 464, 491 sample
IP.SC.reduced <- IP.SC.wtko[-which(IP.SC.wtko$core.ID %in% c("465", "464", "491")),]


# 6 replicaets WT, 6 replicates KO
NIP.SC.wtko <- read.xls("./raw/AVMKe33a_NIP_SpTis_PoolSize_wtko.xlsx")

NIP.SC.reduced <- NIP.SC.wtko[-which(NIP.SC.wtko$Core.ID %in% c("465", "464", "491")),]


## combined IP and NIP
CombinedSC <- cbind(IP.SC.reduced,NIP.SC.reduced)

CombinedSC <- CombinedSC[, !duplicated(colnames(CombinedSC))]
# 9 175

CombinedNumSC <- CombinedSC %>%
  dplyr::select(-core.ID, -Core.ID, -sample, -group) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::select(where(is.numeric))


rownames(CombinedNumSC) <- CombinedSC$sample 

CombinedNumSC <- t(CombinedNumSC)
colnames(CombinedNumSC) <-  c(rep("WT", 4), rep("KO", 5))

logm.combined.sc <- log(CombinedNumSC, 2)

logm.combined.sc <- scale(logm.combined.sc, center = T, scale = T)


####### 
means.combined.sc <- as.data.frame(rowMeans(logm.combined.sc[,1:4]))
means.combined.sc[,2] <- rowMeans(logm.combined.sc[,5:9])
colnames(means.combined.sc) <- c("WT", "KO")


means.combined.sc[,3] <- means.combined.sc["KO"] - means.combined.sc["WT"]
colnames(means.combined.sc)[3] <- "DiffInMeans"


means.sc.reduced <- means.sc[na.omit(match(rownames(means.combined.sc), rownames(means.sc))),]
nrow(means.sc.reduced) # 155

means.combined.sc.reduced <- means.combined.sc[na.omit(match(rownames(means.sc.reduced), rownames(means.combined.sc))),]
nrow(means.combined.sc.reduced) #155


pdf("./plots/metabolomics/Fig5C_AverageEffect_spinalCord.pdf", 20,20)
plot(means.combined.sc.reduced[,3], means.sc.reduced[,3], pch = 16, ylab = "Difference in means (Peak - Healthy) (Average effect)", xlab = "Difference in means (ko-wt) AVMKe33a (Average Effect)", xlim = c(-0.3, 0.4))
abline(h = 0, v= 0, lty = 2)
text(means.combined.sc.reduced[,3], means.sc.reduced[,3], row.names(means.sc.reduced), cex=0.6, pos=4, col="red")
dev.off()



