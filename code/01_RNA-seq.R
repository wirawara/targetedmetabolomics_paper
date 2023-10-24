#!/usr/bin/env Rscript --vanilla

#### RNA seq analysis

library(utils)
library(edgeR)

### change the path accordingly 
path <- "~/Documents/martina/Kerndl_et_al_2023/"

setwd(path)

# load data
annotation <- read.table("./raw/rnaseq_deseq_AVMKe27f_annotation_gene.tsv", header = TRUE)
gene.counts <- read.table("./raw/rnaseq_deseq_AVMKe27f_counts_raw.tsv", header = TRUE)    

gene.counts <- gene.counts[gene.counts$gene_biotype == "protein_coding",]
dim(gene.counts)
# 21936    22
rownames(gene.counts) <- gene.counts$gene_id
# remove not necessary columns
gene.counts <- gene.counts[,c(13:22)]
dim(gene.counts)


# expression matrix with negative and positive arginase - 5vs 5
countTable <- gene.counts[, c(1,7,8,9,10,2,3,4,5,6)]

## design matrix
treatment <- c(rep("negative", 5), rep("positive", 5))
dge <- DGEList(counts=countTable, group=treatment)
design <- model.matrix(~treatment)
rownames(design) <- colnames(dge)
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)

## limma - voom
library(limma)
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v,design)
fit <- eBayes(fit)
top <- topTable(fit,coef=2,number=Inf,sort.by="P")
sum(top$adj.P.Val<0.05)
# 1892
top.orig <- top

library(ggrepel)

## rename the gene names

top.gene.names <- annotation[match(rownames(top), annotation$gene_id),]
rownames(top) <- top.gene.names$gene_name

top$Gene <- top.gene.names$gene_name


top$genelabels <- ""
top$genelabels <- ifelse(top$Gene == "Mrc1"
                              | top$Gene == "Fabp4"
                              | top$Gene == "Fabp5"
                              | top$Gene == "Scd4"
                              | top$Gene == "Lipn"
                              | top$Gene == "Plpp3"
                              | top$Gene == "Cyp2j6"
                              | top$Gene == "Slc6a8"
                              | top$Gene == "Slc36a2"
                              | top$Gene == "Cdo1"
                              | top$Gene == "Akr1b8"
                              | top$Gene == "Mgst2"
                              | top$Gene == "Ccr2"
                              | top$Gene == "Flt3"
                              | top$Gene == "Ccr7"
                              | top$Gene == "Cd83"
                              | top$Gene == "Zbtb46"
                              | top$Gene == "Ccl22"
                              | top$Gene == "Ccl17"
                              | top$Gene == "Ccl6"
                              | top$Gene == "Spp1", TRUE, FALSE)



volcano_plot_summary_update <- function(results, title){

    #results are dataframe with log2FC and adj_pval

    library(ggplot2)
    library(gridExtra)
    library(grid)

    # add a column of NAs
    results$diffexpressed <- "NOT SIG"

    # if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
    results$diffexpressed[results$logFC > 1 & results$adj.P.Val < 0.05] <- "UP"

    # if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
    results$diffexpressed[results$logFC < -1 & results$adj.P.Val < 0.05] <- "DOWN"

    ### to change the colors of the legend
    # mycolors <- c("firebrick1", "dodgerblue2", "gray80")
    # names(mycolors) <- c("DOWN", "UP", "NO")
    print(levels(as.factor(results$diffexpressed)))

    if(length(levels(as.factor(results$diffexpressed)))==3){
        g <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val),label = rownames(results), col=diffexpressed)) +
            geom_point() +  scale_color_manual(values=c("blue", "gray", "red")) +
            geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
            labs(x = bquote(~Log[2]~ "fold change"), y = bquote(~-Log[10]~italic(FDR))) +
            theme_bw() + geom_vline(xintercept=c(-1, 1), linetype = "dashed") +
            ggtitle(title) + 
            geom_text_repel(aes(logFC, -log10(adj.P.Val)),label = ifelse(results$genelabels == TRUE, as.character(results$Gene),""), 
            box.padding = unit(0.45, "lines"),hjust=1)
            
         g2 <- grid.arrange(g,
                bottom = textGrob(title,
                x = 0.40, y = 1, gp = gpar(fontsize = 9)))

        } else {
      
          g <- ggplot(results, aes(x = logFC, y = -log10(adj_pval), col=diffexpressed)) +
                  geom_point() +  scale_color_manual(values=c( "gray")) +
                  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
                  labs(x = bquote(~Log[2]~ "fold change"), y = bquote(~-Log[10]~italic(FDR))) +
                  theme_bw() + geom_vline(xintercept=c(-1, 1), linetype = "dashed") +
                  ggtitle(title)

          g2 <- grid.arrange(g,
                      bottom = textGrob(title,
                      x = 0.45, y = 1, gp = gpar(fontsize = 9)))

        }
    # g3 <- g2 + scale_colour_manual(values = mycolors)

    return(g2)

}


pdf("./plots/RNAseq/Fig3A_volcano_plot_positiveVSnegative_highlights.pdf", 25,25)
volcano_plot_summary_update(top, "Positive vs Negative")
dev.off()



#### fgsea REACTOME / BP
library(biomaRt)
library(org.Mm.eg.db)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = "http://www.ensembl.org")

genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","entrezgene_id"),
               values = rownames(top.orig), 
               mart = mart)

# remove duplicated
genes <- genes[!duplicated(genes$ensembl_gene_id),]
# remove NAs
genes <- genes[!is.na(genes$entrezgene_id),]


topgene_list_entrez_ranked <- top.orig[match(genes$ensembl_gene_id, rownames(top.orig)),]
dim(topgene_list_entrez_ranked)
# [1] 10900     6
topgene_list_entrez_ranked$entrez <- genes$entrezgene_id

### taking the stats
original_gene_list_t <- topgene_list_entrez_ranked$t

# name the vector
names(original_gene_list_t) <- topgene_list_entrez_ranked$entrez

# omit any NA values 
gene_list_t<-na.omit(original_gene_list_t)
length(gene_list_t)
# [1] 10900
# sort the list in decreasing order (required for clusterProfiler)
gene_list_t <- sort(gene_list_t, decreasing = TRUE)
gene_list_t_clean <- gene_list_t[!duplicated(names(gene_list_t))]


### REACTOME analysis
library(reactome.db) ## there are names here
library(ReactomePA)
library(fgsea)
library(org.Mm.eg.db)

my_pathways <- reactomePathways(names(gene_list_t_clean))
fgsea_reactome_t <- fgsea(pathways = my_pathways, 
                        stats = gene_list_t_clean,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)

sum(fgsea_reactome_t[, padj < 0.01]) # 91
head(fgsea_reactome_t[order(padj), ])

## convert the entrez gene names in the leading edge
ff <- fgsea_reactome_t[, leadingEdge := mapIdsList(org.Mm.eg.db, keys=leadingEdge, column="SYMBOL", keytype="ENTREZID")]
ff_reactome_pathways <- ff$leadingEdge
names(ff_reactome_pathways) <- fgsea_reactome_t$pathway


ff_reactome <- fgsea_reactome_t[which(fgsea_reactome_t$padj < 0.01),]
ff_reactome$pathway <- factor(ff_reactome$pathway,                                    # Factor levels in decreasing order
                  levels = ff_reactome$pathway[order(ff_reactome$NES, decreasing = FALSE)])


reactome_down_highlighted_v2 <- c("Formation of a pool of free 40S subunits",
"GTP hydrolysis and joining of the 60S ribosomal subunit",
"Nonsense-Mediated Decay (NMD)",
"Cap-dependent Translation Initiation",
"Eukaryotic Translation Initiation",
"SRP-dependent cotranslational protein targeting to membrane",
"Major pathway of rRNA processing in the nucleolus and cytosol",
"rRNA processing",
"rRNA processing in the nucleus and cytosol",
"Ribosomal scanning and start codon recognition")

reactome_up_highlighted_v2 <- c(
"Organelle biogenesis and maintenance",
"Metabolism of vitamins and cofactors",
"Metabolism of carbohydrates",
"Metabolism of lipids",
"Glycosaminoglycan metabolism",
"Heparan sulfate/heparin (HS-GAG) metabolism",
"Fatty acid metabolism",
"Mitochondrial Fatty Acid Beta-Oxidation",
"Glycosphingolipid metabolism",
"Sphingolipid metabolism")


ff_reactome_pathways_up_v2 <- ff_reactome_pathways[match(as.character(reactome_up_highlighted_v2), names(ff_reactome_pathways))]
ff_reactome_pathways_down_v2 <- ff_reactome_pathways[match(as.character(reactome_down_highlighted_v2), names(ff_reactome_pathways))]


ff_reactome_pathways_up_allinfo_v2 <- ff_reactome[match(as.character(reactome_up_highlighted_v2), ff_reactome$pathway),]
ff_reactome_pathways_down_allinfo_v2 <- ff_reactome[match(as.character(reactome_down_highlighted_v2), ff_reactome$pathway),]

ff_reactome_highlighted_v2 <- rbind(ff_reactome_pathways_up_allinfo_v2, ff_reactome_pathways_down_allinfo_v2)


ff_reactome_highlighted_v2$pathway <- factor(ff_reactome_highlighted_v2$pathway,                                    # Factor levels in decreasing order
                  levels = ff_reactome_highlighted_v2$pathway[order(ff_reactome_highlighted_v2$NES, decreasing = FALSE)])


c4 <- c(rep("indianred1", 10), rep("gray80", 10))
ff_reactome_highlighted_v2 <- cbind(ff_reactome_highlighted_v2, c4)

pdf("./plots/RNAseq/Fig3B_fgsea_REACTOME_highlighted.pdf", 10,15)
ggplot(ff_reactome_highlighted_v2, aes(pathway, NES)) + theme_classic() +
    # scale_fill_brewer("", palette="Paired") +
    geom_bar(position= "identity", stat="identity", width=0.8, fill = c4) +
    coord_flip() +
    scale_y_continuous(position = "right") +
    labs( y = "Normalized Enrichment Score Reactome") +
    theme(axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    text = element_text(size=14))
dev.off()



#### GO BP
### GO terms 
library(GO.db)

goTerms <- toTable(GOTERM)
goTerms.unique <- goTerms[!duplicated(goTerms$go_id),]
goTermsBP <- goTerms.unique[which(goTerms.unique$Ontology == "BP"),]
goTermsMF <- goTerms.unique[which(goTerms.unique$Ontology == "MF"),]
goTermsCC <- goTerms.unique[which(goTerms.unique$Ontology == "CC"),]


library(org.Mm.eg.db)

xx <- as.list(org.Mm.egGO2ALLEGS)

xx.gobp <- xx[match(goTermsBP$go_id,names(xx))]
xx.gobp <- xx.gobp[!is.na(names(xx.gobp))]

xx.gomf <- xx[match(goTermsMF$go_id,names(xx))]
xx.gomf <- xx.gomf[!is.na(names(xx.gomf))]


xx.gocc <- xx[match(goTermsCC$go_id,names(xx))]
xx.gocc <- xx.gocc[!is.na(names(xx.gocc))]


## perform the fgsea
## bp
library(fgsea)
fgsea_go_bp_t <- fgsea(pathways = xx.gobp, 
                        stats = gene_list_t_clean,
                         minSize=15,
                         maxSize=500)

sum(fgsea_go_bp_t[, padj < 0.01])
head(fgsea_go_bp_t[order(padj), ])

ffbp <- fgsea_go_bp_t

###
ffbp <- ffbp[match(goTermsBP$go_id, ffbp$pathway),] 
ffbp$Term <- goTermsBP$Term
ffbp <- na.omit(ffbp)

# order according to NES and adj p values

ffbp_signif <- ffbp[which(ffbp$padj < 0.01),]
ffbp_signif$Term <- factor(ffbp_signif$Term,                                    # Factor levels in decreasing order
                  levels = ffbp_signif$Term[order(ffbp_signif$NES, decreasing = FALSE)])



ffgobp_tops <- ffbp_signif[order(ffbp_signif$NES, decreasing = FALSE),]
# library(data.table)
# fwrite(ffgobp_tops , "ffgobp_tops_fwrite_16Jan2023.csv")


ff_go_bp_top10 <- tail(ffgobp_tops, 10)
ff_go_bo_low10 <- head(ffgobp_tops, 10)


ff_gobp_highlighted <- rbind(ff_go_bp_top10, ff_go_bo_low10)


ff_gobp_highlighted$Term <- factor(ff_gobp_highlighted$Term,                                    # Factor levels in decreasing order
                  levels = ff_gobp_highlighted$Term[order(ff_gobp_highlighted$NES, decreasing = FALSE)])


c4 <- c(rep("indianred1", 10), rep("gray80", 10))
ff_gobp_highlighted <- cbind(ff_gobp_highlighted, c4)

pdf("./plots/RNAseq/Fig3E_fgsea_GOBP_signif_ordered_16Jan2023_AVMKe27f.pdf", 20,23)

ggplot(ff_gobp_highlighted, aes(Term, NES)) + theme_classic() +
    # scale_fill_brewer("", palette="Paired") +
    geom_bar(position= "identity", stat="identity", width=0.8, fill = c4) +
    coord_flip() +
    scale_y_continuous(position = "right") +
    labs( y = "GO BP Normalized Enrichment Score") +
    theme(axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    text = element_text(size=14))
dev.off()


###### analysis for the leading edge

#### run this function before plotting the GSEA plot
leading_edge <- function(observed_info) {
    core_enrichment <- lapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            leading_gene <- runningES$gene[1:i]
        } else {
            i <- which.min(runningES$runningScore)
            leading_gene <- runningES$gene[-c(1:(i-1))]
        }
        return(leading_gene)
    })

    rank <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) {
            rr <- which.max(runningES$runningScore)
        } else {
            i <- which.min(runningES$runningScore)
            rr <- nrow(runningES) - i + 1
        }
        return(rr)
    })

    tags <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        runningES <- runningES[runningES$position == 1,]
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else {
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    ll <- sapply(observed_info, function(x) {
        runningES <- x$runningES
        ES <- x$ES
        if (ES >= 0) {
            i <- which.max(runningES$runningScore)
            res <- i/nrow(runningES)
        } else {
            i <- which.min(runningES$runningScore)
            res <- (nrow(runningES) - i + 1)/nrow(runningES)
        }
        return(res)
    })

    N <- nrow(observed_info[[1]]$runningES)
    setSize <- sapply(observed_info, function(x) sum(x$runningES$position))
    signal <- tags * (1-ll) * (N / (N - setSize))

    tags <- paste0(round(tags * 100), "%")
    ll <- paste0(round(ll * 100), "%")
    signal <- paste0(round(signal * 100), "%")
    leading_edge <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)

    res <- list(rank = rank,
                tags = tags,
                list = ll,
                signal = signal,
                leading_edge = leading_edge,
                core_enrichment = core_enrichment)
    return(res)
}



#### analysis for the leading edge

library(reactome.db) #
library(ReactomePA)
library(fgsea)

library(clusterProfiler)
library(enrichplot)

my_pathways <- reactomePathways(names(gene_list_t_clean))

fgsea_reactome_t <- fgsea(pathways = my_pathways, 
                        stats = gene_list_t_clean,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)


p.adj <- p.adjust( fgsea_reactome_t$pval, method="BH")
   
Description <- fgsea_reactome_t$pathway
    

res <- data.frame(
        ID = as.character(fgsea_reactome_t$pathway),
        Description = unname(Description),
        setSize = fgsea_reactome_t$size,
        enrichmentScore = fgsea_reactome_t$ES,
        NES = fgsea_reactome_t$NES,
        pvalue = fgsea_reactome_t$pval,
        p.adjust = p.adj,
        # qvalue = qvalues,
        stringsAsFactors = FALSE
    )

res <- res[!is.na(res$pvalue),]
res <- res[ res$pvalue <= 0.05, ]
res <- res[ res$p.adjust <= 0.05, ]
idx <- order(res$p.adjust, -abs(res$NES), decreasing = FALSE)
res <- res[idx, ]

row.names(res) <- res$ID
gseaScores <- getFromNamespace("gseaScores", "DOSE")
observed_info <- lapply(my_pathways[res$ID], function(gs)
        gseaScores(geneSet=gs,
                   geneList=gene_list_t_clean,
                   exponent=1)
        )

ledge <- leading_edge(observed_info)

res$rank <- ledge$rank
res$leading_edge <- ledge$leading_edge
res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')

params <- list(pvalueCutoff = 0.05,
                       nPerm = 100000,
                        pAdjustMethod = "BH",
                        exponent = 1,
                        minGSSize = 15,
                        maxGSSize = 500
                    )

results <- new("gseaResult",
        result     = res,
        geneSets   = my_pathways,
        geneList   = gene_list_t_clean,
        params     = params,
        readable   = FALSE
        )

results@organism <- "mmu"
results@setType <- "KEGG"
results@keytype <- "UNKNOWN"
 



pdf("./plots/RNAseq/Fig3D_Formation_of_a_pool_of_free_40S_subunits.pdf",5,5)
gseaplot2(results, geneSetID = 1, title = results$Description[1])
dev.off()


pdf("./plots/RNAseq/Fig3C_Sphingolipid_metabolism.pdf",5,5)
gseaplot2(results, geneSetID = 15, title = results$Description[15])
dev.off()






