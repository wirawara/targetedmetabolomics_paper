#!/usr/bin/env Rscript --vanilla

#### scRNA seq analysis



library(Seurat)
library(ggplot2)

### change the path accordingly 
path <- "~/Documents/martina/Kerndl_et_al_2023/"

setwd(path)



rds.file <- readRDS("./raw/Mike.B6.EAE.all.cells.figure.1.data.rds")

pdf("./plots/scRNAseq/SupplFgi4A_scRNA_seq_all.pdf", 10, 10)
DimPlot(rds.file,label.size = 4,pt.size =0.8, repel = T,label = T)
dev.off()

## specific genes
pdf("./plots/scRNAseq/Fig4A_scRNA_seq_all_Arg1.pdf", 5, 5)
# Arg1
FeaturePlot(rds.file, features = c("Arg1"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()

pdf("./plots/scRNAseq/Fig4A_scRNA_seq_all_Lyz2.pdf", 5, 5)
# Lyz2
FeaturePlot(rds.file, features = c("Lyz2"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()

pdf("./plots/scRNAseq/Fig4A_scRNA_seq_all_Ly6c2.pdf", 5, 5)
# Ly6c2
FeaturePlot(rds.file, features = c("Ly6c2"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()

pdf("./plots/scRNAseq/Fig4B_scRNA_seq_all_Ccr2.pdf", 5, 5)
# Ccr2
FeaturePlot(rds.file, features = c("Ccr2"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()

pdf("./plots/scRNAseq/Fig4A_scRNA_seq__all_Ly6i.pdf", 5, 5)
# Ly6i
FeaturePlot(rds.file, features = c("Ly6i"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()


#### ------------------------------------------
### only myeloid clusters

rds.file.subset2 <- subset(x = rds.file, idents = c("2", "0", "1","3","4","18"))

##### scRNAseq only myeloid clusters extraction
my_cols2 <- c("2"="#E41A1C","0"="#A43E55","1"="#64638E","3"="#3983AC", "4"="#419583","18"="#c51b8a")
pdf("./plots/scRNAseq/Fig4B_scRNAseq_with_onlymyeloidclusters.pdf", 10,10)
DimPlot(rds.file.subset2, cols = my_cols2,pt.size = 0.8, label = T)
dev.off()


### dot plot of differentially 
pdf("./plots/scRNAseq/Fig4C_dotplot_Arg1_Nos2_Fabp5_Spp1_Ccl6.pdf")
DotPlot(rds.file.subset2, features=c("Arg1","Nos2", "Fabp5", "Spp1", "Ccl6"), dot.scale=10, idents=c('2', '0', '1', '3', '4', '18'), cols="RdBu")
dev.off()



#### ---------------------------------
# subclustering

arg_subset <- subset(rds.file, idents = c("2"))
# An object of class Seurat 
# 37873 features across 2675 samples within 2 assays 
# Active assay: RNA (35873 features, 0 variable features)
#  1 other assay present: integrated
#  2 dimensional reductions calculated: pca, tsne
DefaultAssay(arg_subset) <- "integrated"
arg_subset <- FindNeighbors(arg_subset, dims = 1:10)
arg_subset <- FindClusters(arg_subset,resolution = 0.1)
# Generate a new column called sub_cluster in the metadata
arg_subset$sub_cluster <- as.character(Idents(arg_subset))
# Change the information of cells containing sub-cluster information
rds.file$sub_cluster <- as.character(Idents(rds.file))
rds.file$sub_cluster[Cells(arg_subset)] <- paste("c2",Idents(arg_subset))

Idents(rds.file) <- rds.file$sub_cluster

rds.file.subset <- subset(x = rds.file, idents = c("c2 0", "c2 1", "0", "1","3","4","18"))

pdf("./plots/scRNAseq/Fig4E_scRNAseq_with_onlymyeloidclusters_subclustering_with_colors.pdf", 10,10)
my_cols <- c("c2 0"="#bcbddc", "c2 1"="#E41A1C","0"="#A43E55","1"="#64638E","3"="#3983AC", "4"="#419583","18"="#c51b8a")
DimPlot(rds.file.subset, cols = my_cols,pt.size = 0.8, label = T)
dev.off()

pdf("./plots/scRNAseq/SupplFig4C_scRNAseq_with_onlymyeloidclusters_Csf2ra.pdf", 5,5)
FeaturePlot(rds.file.subset, features = c("Csf2ra"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()

pdf("./plots/scRNAseq/SupplFig4C_scRNAseq_with_onlymyeloidclusters_Csf2rb.pdf", 5,5)
FeaturePlot(rds.file.subset, features = c( "Csf2rb"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()


pdf("./plots/scRNAseq/Fig4D_scRNAseq_with_onlymyeloidclusters_Slc7a2.pdf", 5,5)
FeaturePlot(rds.file.subset, features = c( "Slc7a2"),cols = c("lightgrey", "maroon"), pt.size = 0.8, ncol = 1) 
dev.off()




#### -----------------------------------------------------------
# fgsea
markers_cluster_cl2_0_vs_myeloid <- FindMarkers(rds.file, ident.1 = "c2 0", ident.2 = c(0,1,3,4,18))

#### excel file with those genes
markers_cluster_cl2_0_vs_myeloid_signif <- markers_cluster_cl2_0_vs_myeloid[markers_cluster_cl2_0_vs_myeloid$p_val_adj < 0.05,]
# 1096

library(biomaRt)
library(fgsea)
# convert ensembl id to entrez

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = "http://www.ensembl.org")

top.Scl7a2vsmyelo <- markers_cluster_cl2_0_vs_myeloid

genes.Scl7a2vsmyelo.cluster <- getBM(filters = "external_gene_name",
                        attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                        values = rownames(top.Scl7a2vsmyelo), 
                        mart = mart)

genes.Scl7a2vsmyelo.cluster$external_gene_name <- make.names(genes.Scl7a2vsmyelo.cluster$external_gene_name, unique = TRUE)

topgene.list.entrez.Scl7a2vsmyelo <- top.Scl7a2vsmyelo[na.omit(match(genes.Scl7a2vsmyelo.cluster$external_gene_name, rownames(top.Scl7a2vsmyelo))),]
genes.Scl7a2vsmyelo.cluster.ordered <- genes.Scl7a2vsmyelo.cluster[na.omit(match(rownames(topgene.list.entrez.Scl7a2vsmyelo),genes.Scl7a2vsmyelo.cluster$external_gene_name)),]


topgene.list.entrez.Scl7a2vsmyelo$entrez <- genes.Scl7a2vsmyelo.cluster.ordered$entrezgene_id
topgene.list.entrez.ranked.Scl7a2vsmyelo <- topgene.list.entrez.Scl7a2vsmyelo[order(topgene.list.entrez.Scl7a2vsmyelo$p_val_adj),]


# we want the log2 fold change 
original.gene.list.avglogFC.Scl7a2vsmyelo <- topgene.list.entrez.ranked.Scl7a2vsmyelo$avg_log2FC

# name the vector
names(original.gene.list.avglogFC.Scl7a2vsmyelo) <- topgene.list.entrez.ranked.Scl7a2vsmyelo$entrez

# omit any NA values 
gene.list.Scl7a2vsmyelo <- na.omit(original.gene.list.avglogFC.Scl7a2vsmyelo)

# sort the list in decreasing order (required for clusterProfiler)
gene.list.Scl7a2vsmyelo <- sort(gene.list.Scl7a2vsmyelo, decreasing = TRUE)




gobp.up.highlighted <- c("GO:1905952", "GO:0033993", "GO:0032496", "GO:0071346", "GO:0032368", "GO:0010876", "GO:0015711", "GO:0001816", "GO:0042886", "GO:0072593")
gobp.down.highlighted <- c("GO:0090304", "GO:1902679", "GO:0051252", "GO:0006366", "GO:0016070", "GO:0032774", "GO:0045934", "GO:0060284", "GO:0009966", "GO:0010558")


library(GO.db)
goTerms <- toTable(GOTERM)
goTerms.unique <- goTerms[!duplicated(goTerms$go_id),]
goTermsBP <- goTerms.unique[which(goTerms.unique$Ontology == "BP"),]

## gene to GO mappings for an organism Mus musculus
library(org.Mm.eg.db)
xx2 <- as.list(org.Mm.egGO2ALLEGS)

xx.gobp <- xx2[match(goTermsBP$go_id,names(xx2))]
xx.gobp <- xx.gobp[!is.na(names(xx.gobp))]

fgsea_go_bp <- fgsea(pathways = xx.gobp, 
                         stats = gene.list.Scl7a2vsmyelo,
                          minSize=15,
                          maxSize=500,
                          nperm=100000)

sum(fgsea_go_bp[, padj < 0.10]) # 463
head(fgsea_go_bp[order(padj), ])

ffbp <- fgsea_go_bp

###
ffbp <- ffbp[match(goTermsBP$go_id, ffbp$pathway),] 
ffbp$Term <- goTermsBP$Term
ffbp <- na.omit(ffbp)

# order according to NES and adj p values

ffbp_signif <- ffbp[which(ffbp$padj < 0.05),]
ffbp_signif$Term <- factor(ffbp_signif$Term,                                    # Factor levels in decreasing order
                  levels = ffbp_signif$Term[order(ffbp_signif$NES, decreasing = FALSE)])

#### padj
head(as.data.frame(ffbp_signif[order(padj), ]), 60)


#### take only highlighted


ff_bp_highlighted_up <- ffbp_signif[ffbp_signif$pathway %in% gobp.up.highlighted,]
ff_bp_highlighted_down <- ffbp_signif[ffbp_signif$pathway %in% gobp.down.highlighted,]


ff_gobp_highlighted_slc7a2highvsothers <- rbind(ff_bp_highlighted_up, ff_bp_highlighted_down)


c4 <- c(rep("indianred1", 10), rep("gray80", 10))
ff_gobp_highlighted_slc7a2highvsothers <- cbind(ff_gobp_highlighted_slc7a2highvsothers, c4)

ff_gobp_highlighted_slc7a2highvsothers$Term <- factor(ff_gobp_highlighted_slc7a2highvsothers$Term,                                    # Factor levels in decreasing order
                  levels = ff_gobp_highlighted_slc7a2highvsothers$Term[order(ff_gobp_highlighted_slc7a2highvsothers$NES, decreasing = FALSE)])

ff_gobp_highlighted_slc7a2highvsothers$Term <- paste(ff_gobp_highlighted_slc7a2highvsothers$pathway, ff_gobp_highlighted_slc7a2highvsothers$Term, sep=", ")
ff_gobp_highlighted_slc7a2highvsothers$Term <- factor(ff_gobp_highlighted_slc7a2highvsothers$Term, levels=rev(ff_gobp_highlighted_slc7a2highvsothers$Term))

ff_gobp_highlighted_slc7a2highvsothers$Term <- factor(ff_gobp_highlighted_slc7a2highvsothers$Term,                                    # Factor levels in decreasing order
                  levels = ff_gobp_highlighted_slc7a2highvsothers$Term[order(ff_gobp_highlighted_slc7a2highvsothers$NES, decreasing = FALSE)])



library(ggplot2)
pdf("./plots/scRNAseq/Fig4F_fgsea_GOBP_highlighted_Slc7a2_high_vs_othermyeloid.pdf", 15,15)
### this works, there is no negative enrichment scores
ggplot(ff_gobp_highlighted_slc7a2highvsothers, aes(Term, NES)) + theme_classic() +
    # scale_fill_brewer("", palette="Paired") +
    geom_bar(position= "identity", stat="identity", width=0.8, fill = c4) +
    coord_flip() +
    scale_y_continuous(position = "right") +
    labs( y = "Normalized Enrichment Score GO BP") +
    theme(axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    text = element_text(size=14))
dev.off()



