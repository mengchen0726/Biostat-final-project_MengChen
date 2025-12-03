library(GEOquery)
library(oligo)
library(affy)
library(limma)
library(sva)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(patchwork)
library(pheatmap)
library(grDevices)
library(ggplotify)
library(VennDiagram)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(WGCNA)
library(glmnet)
library(randomForest)


options(stringsAsFactors = F) 

dir.create("./plots")
dir.create("./data")

#PD dataset
gse1 <- getGEO(filename = "./data/GSE16134_series_matrix.txt.gz", getGPL = F)
#pSS dataset
gse2 <- getGEO(filename = "./data/GSE40611_series_matrix.txt.gz", getGPL = F)

#clinical information
pd1 <- pData(gse1)
pd2 <- pData(gse2)

untar("./data/GSE16134_RAW.tar", exdir = "./data/GSE16134_RAW")
untar("./data/GSE40611_RAW.tar", exdir = "./data/GSE40611_RAW")

raw1 <- ReadAffy(celfile.path = "./data/GSE16134_RAW")
#background correction and normalization
exp1 <- rma(raw1)
exp1 <- exprs(exp1)


raw2 <- ReadAffy(celfile.path = "./data/GSE40611_RAW")
exp2 <- rma(raw2)
exp2 <- exprs(exp2)

save(exp1, exp2, file = "./data/step0_Affy_raw.Rdata")


###quality control
load("./data/step0_Affy_raw.Rdata")
range(exp1)
range(exp2)

#expression distribution and variance homogeneity
boxplot(exp1)
boxplot(exp2)

#normalization
#exp1 <- normalizeBetweenArrays(exp1)
#exp2 <- normalizeBetweenArrays(exp2)

#Convert probe IDs to gene symbols
gpl <- read.csv("./data/GPL570-55999.txt", sep = "\t", header = T, comment.char = "#")
id2s <- data.frame(ID = gpl$ID, SYMBOL = gpl$Gene.Symbol)
id2s[1:4,]
id2s <- id2s[!str_detect(id2s$SYMBOL, "///"), ]

#ID conversion
IDtrans <- function(id2s, exp){
  
  exp <- exp[match(id2s$ID, rownames(exp)),]
  rownames(exp) <- id2s$SYMBOL
  exp <- as.data.frame(exp)
  
  exp$RowName <- id2s$SYMBOL
  
  exp <- exp %>%
    filter(!is.na(RowName) & RowName != "")
  
  result <- exp %>%
    group_by(RowName) %>%
    summarise(across(everything(), mean))
  
  result <- as.data.frame(result)
  rownames(result) <- result$RowName
  result <- result[,-1]
  
  return(result)
}

exp1 <- IDtrans(id2s, exp1)
exp2 <- IDtrans(id2s, exp2)

dim(exp1); dim(exp2)

save(exp1, exp2, file = "./data/Step1_transid_exp.RData")

###Clinical information
colnames(exp1)
rownames(pd1)

colnames(exp1) <- str_split(colnames(exp1), "\\.", simplify = T)[,1]

colnames(exp2)
rownames(pd2)

colnames(exp2) <- str_split(colnames(exp2), "_", simplify = T)[,1]

#Check consistency
identical(rownames(pd1), colnames(exp1))
identical(rownames(pd2), colnames(exp2))

table(pd1$`gingival tissue phenotype:ch1`)
group1 <- ifelse(str_detect(pd1$`gingival tissue phenotype:ch1`, "Diseased"), "Diseased", "Healthy")

table(pd2$source_name_ch1)
group2 <- pd2$`disease status:ch1`


#save only pss & control
filt <- group2 != "Sicca"
group2 <- group2[filt]
exp2 <- exp2[,filt]
length(group2); ncol(exp2)


table(group1)
group1 <- factor(group1, levels = c("Diseased", "Healthy"))
table(group2)
group2 <- factor(group2, levels = c("pSS", "Control"))

save(exp1, group1, pd1,
     exp2, group2, pd2,
     file = "./data/step2_exp&group.RData")


###PCA
qc_PCA <- function(exp, group){
  cg = names(tail(sort(apply(exp,1,sd)),1000))
  dat = t(exp[cg,])
  dat = as.data.frame(dat)
  
  dat.pca <- PCA(dat , graph = FALSE)
  p <- fviz_pca_ind(dat.pca,
                    geom.ind = "point", # show points only (nbut not "text")
                    col.ind =  group, # color by groups 
                    addEllipses = T,
                    legend.title = "Groups"
  )
  return(p)
}




p1 <- qc_PCA(as.matrix(exp1), group1) + ggtitle("PD_PCA")
p2 <- qc_PCA(exp2, group2) + ggtitle("pSS_PCA")

p1 + p2

ggsave('./plots/QC_all_samples_PCA.png', p1 + p2, height = 5, width = 8)


###Differential expression analysis
doDEG <- function(dat, group, logFC_t = NULL, padj_t = 0.05, vs = NULL){
  design <- model.matrix(~0 + group)
  colnames(design)=levels(group)
  rownames(design)=colnames(dat)
  if(is.null(vs)){
    vs <- paste0(levels(group)[1], "-", levels(group)[2])
  }
  
  contrast.matrix<-makeContrasts(contrasts = vs,levels = design)
  
  fit <- lmFit(dat,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  deg = topTable(fit2, coef=1, n=Inf)
  deg = na.omit(deg)
  if(is.null(logFC_t)){
    logFC_t <- with(deg,mean(abs(deg$logFC)) + 2*sd(abs(deg$logFC)) )
  }
  print(paste0("The threshold of logFoldChange is ", logFC_t))
  deg$v= -log10(deg$adj.P.Val)
  deg$change <- ifelse(deg$adj.P.Val > padj_t,'stable',
                       ifelse(deg$logFC > logFC_t,'up',
                              ifelse(deg$logFC < -logFC_t,'down','stable'))
  )
  deg$name=rownames(deg)
  return(deg)
}

# |logFoldChange|= 0.5, p.adjust = 0.05
DEG1 <- doDEG(exp1, group1, logFC_t = 0.5)
table(DEG1$change)

DEG2 <- doDEG(exp2, group2, logFC_t = 0.5)
table(DEG2$change)
save(DEG1, DEG2, file = "./data/step3_deg.Rdata")

###volcano plot
DEG_volcano <- function(deg,
                        volcano_tile = F,
                        logFC_cut_off = 1,
                        padj_cut_off = 0.05,
                        max_x=10,
                        min_x=-10,
                        max_y=0,
                        min_y=50,
                        target_gene=NULL,
                        color=NULL)
{
  label <-  deg %>%
    dplyr::filter(deg$name %in% target_gene)
  if(is.null(color)){
    c("#354587","#A0A0A0","#FC574C")->color
  }
  p <- ggplot(data = deg,
              aes(x = logFC,
                  y = v)) +
    geom_point(alpha=0.6, size=1.5,
               aes(color=change)) +
    geom_point(size = 3, shape = 1 , data = label) +
    ggrepel::geom_text_repel(aes(label = name),data = label,
                             max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                             segment.size = 0.2, segment.color = "black",
                             color="black") +
    ylab("-log10(P.adjust)")+
    scale_color_manual(values=color)+
    geom_vline(xintercept= 0,lty=4,col="grey",lwd=0.8) +
    xlim(min_x, max_x)+
    ylim(min_y, max_y)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(size=8,hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=16))
  
  if(T) p<-p+ggtitle(paste0("Threshold of logFC is ", logFC_cut_off,"\nThe number of up gene is ",
                            nrow(deg[deg$change == "up", ]), "\nThe number of down gene is ", nrow(deg[deg$change == "down", ])))+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=16,face="bold"))+
    theme(plot.title = element_text(size = 16, face = "bold"))+
    theme(legend.key.size = unit(40, "pt"))
  
  if(T){
    p<-p+
      geom_vline(xintercept=c(-logFC_cut_off,logFC_cut_off),lty=4,col="black",lwd=0.4) +
      geom_hline(yintercept = -log10(padj_cut_off),lty=4,col="black",lwd=0.4)
  }
  
  return(p)
}

target_gene <- c("CXCR4", "LYN", "CSF2RB")
#target_gene <- c("CSF2RB")

p1 <- DEG_volcano(DEG1,
                  logFC_cut_off = 0.5,
                  padj_cut_off = 0.05,
                  max_x = 4,
                  min_x = -4,
                  max_y = 40,
                  min_y = 0,
                  target_gene = target_gene)
p2 <- DEG_volcano(DEG2,
                  logFC_cut_off = 0.5,
                  padj_cut_off = 0.05,
                  max_x = 5,
                  min_x = -5,
                  max_y = 8,
                  min_y = 0,
                  target_gene = target_gene)
ggsave(p1 + p2,
       filename = "./plots/Fig2_DEG_volcano.png",
       width = 10,       height = 5)

###Heatmap
DEG_heatmap <- function(dat, group, target_gene) 
{
  n = dat
  ac = data.frame(group = group)
  rownames(ac) = colnames(n)
  color <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", 
             "#FF9896FF", "#9467BDFF", "#C5B0D5FF", "#8C564BFF", "#98DF8AFF", 
             "#C49C94FF", "#E377C2FF", "#FFBB78FF", "#F7B6D2FF", "#7F7F7FFF", 
             "#C7C7C7FF")[1:length(levels(group))]
  names(color) <- names(table(group))
  
  p <- pheatmap(n,
                color = colorRampPalette(c("#00008B", "#555597", "#FFFFB0", 
                                           "#B1553A", "#8B0000"))(100),
                labels_row = target_gene, 
                show_colnames = F,
                show_rownames = T, 
                cluster_cols = F,
                annotation_colors = list(group = color), 
                annotation_col = ac)
  return(p)
}

genes1 <- rownames(subset(DEG1, change != "stable"))
genes2 <- rownames(subset(DEG2, change != "stable"))

genes1 <- intersect(genes1, rownames(exp1))
genes2 <- intersect(genes2, rownames(exp2))

plot_exp1 <- exp1[genes1, , drop = FALSE]
plot_exp2 <- exp2[genes2, , drop = FALSE]

plot_exp1 <- plot_exp1[DEG1[DEG1$change!="stable", "name"],]
rownames_custom1 <- ifelse(rownames(plot_exp1) %in% target_gene,
                           rownames(plot_exp1), "")
plot_exp2 <- plot_exp2[DEG2[DEG2$change!="stable", "name"],]
rownames_custom2 <- ifelse(rownames(plot_exp2) %in% target_gene, rownames(plot_exp2), "")

plot_group1 <- group1
plot_group2 <- group2

ord1 <- order(plot_group1)
ord2 <- order(plot_group2)

plot_exp1   <- plot_exp1[, ord1, drop = FALSE]
plot_exp2   <- plot_exp2[, ord2, drop = FALSE]
plot_group1 <- droplevels(plot_group1[ord1])
plot_group2 <- droplevels(plot_group2[ord2])


p1 <- DEG_heatmap(plot_exp1, plot_group1, rownames_custom1)
p2 <- DEG_heatmap(plot_exp2, plot_group2, rownames_custom2)

ggsave(as.ggplot(p1) + as.ggplot(p2),
       filename = "./plots/Fig2_DEG_heatmap.png",
       width = 12,       height = 6)

###Venn
PD_up <- DEG1[DEG1$change == "up","name"]
PD_down <- DEG1[DEG1$change == "down","name"]
pSS_up <- DEG2[DEG2$change == "up","name"]
pSS_down <- DEG2[DEG2$change == "down","name"]

png(filename = "./plots/Fig2_venn.png", width = 480, height = 480)
venn.plot <- venn.diagram(
  x = list(Set1 = PD_up, Set2 = pSS_down, Set3 = pSS_up, Set4 = PD_down),
  category.names = c("PD_up", "pSS_down", "pSS_up", "PD_down"),
  filename = NULL,
  output = TRUE,
  cex = 3,
  fill = c("#D62728FF", "#2CA02CFF", "#FF7F0EFF", "#1F77B4FF"),
  cat.cex = 1.5,
  cat.col = c("#D62728FF", "#2CA02CFF", "#FF7F0EFF", "#1F77B4FF")
)
grid.draw(venn.plot)
dev.off()

###Enrichment analysis
#commonDEG from PD & pSS
commonDEGs <- intersect(DEG1[DEG1$change!="stable", "name"],
                        DEG2[DEG2$change!="stable", "name"])

library(org.Hs.eg.db)
library(AnnotationDbi)
keytypes(org.Hs.eg.db)
k <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID") 
e2s <- AnnotationDbi::select(org.Hs.eg.db, keys = k,
                             columns="SYMBOL", keytype = "ENTREZID")
ENTREZID_commonDEGs <- e2s[match(commonDEGs, e2s$SYMBOL),"ENTREZID"] %>% na.omit %>% as.numeric

library(clusterProfiler)
library(enrichplot)
#GO
ont <- c("BP", "MF", "CC")
for(i in ont){
  go <- enrichGO(ENTREZID_commonDEGs,
                 OrgDb = "org.Hs.eg.db",
                 ont = i,
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05)
  p <- enrichplot::dotplot(go, showCategory = 10) +
    labs(title = i) +
    theme(plot.title = element_text(size = 20))
  ggsave(p, filename = paste0("./plots/Fig3_ORA_GO_", i, ".png"),
         width    = 8,
         height   = 6,
         dpi      = 300)
}

#KEGG
kegg <- enrichKEGG(ENTREZID_commonDEGs,
                   organism = "hsa")
p <- enrichplot::dotplot(kegg, showCategory = 10) +
  labs(title = "KEGG enrichment") +
  theme(plot.title = element_text(size = 20))
ggsave(p, filename = paste0("./plots/Fig3_ORA_KEGG.png"),
       width    = 8,
       height   = 6,
       dpi      = 300)

###WGCNA
load("./data/step2_exp&group.RData")
library(WGCNA)


Pick_SFT <- function(top_sd_gene_number = 4000,
                     exp,
                     name,
                     RsquaredCut = 0.85){
  apply(exp, 1, sd) %>% sort(.,decreasing = T) %>% head(., top_sd_gene_number) %>% names %>% exp[.,] -> WGCNA_exp
  
  t_exp <- t(WGCNA_exp)
  
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(t_exp, 
                          RsquaredCut = RsquaredCut,
                          powerVector = powers, 
                          verbose = 5)
  po <- sft$powerEstimate
  assign(paste0("po_", name), po, envir = .GlobalEnv)
  print(po)
  
  png(paste0("./plots/Fig4_pick_SFT_", name, "_RsquaredCut_", RsquaredCut,".png"), width = 900, height = 500)
  
  par(mfrow = c(1,2));
  cex1 = 0.8;
  
  plot(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  
  abline(h=RsquaredCut,col="red")
  
  plot(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       xlab="Soft Threshold (power)",
       ylab="Mean Connectivity", 
       type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], 
       sft$fitIndices[,5], 
       labels=powers, cex=cex1,col="red")
  
  dev.off()
}

Pick_SFT(top_sd_gene_number = 10000, exp1, name = "exp1", RsquaredCut = 0.85)
Pick_SFT(top_sd_gene_number = 10000, exp2, name = "exp2", RsquaredCut = 0.8)





###WGCNA-PD
apply(exp1, 1, sd) %>% sort(.,decreasing = T) %>% head(., 10000) %>% names %>% exp1[.,] -> WGCNA_exp1
t_exp1 <- t(WGCNA_exp1)

allowWGCNAThreads()
adjacency = adjacency(t_exp1, power = po_exp1)
TOM = TOMsimilarity(adjacency)
save(TOM, file = paste0("./data/step6_TOM_exp1_po", po_exp1, ".RData"))
load(paste0("./data/step6_TOM_exp1_po", po_exp1, ".RData"))
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

png("geneTree.png", width = 800, height = 600)
par(mar = c(1,1,1,1))

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()


dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = 50) 
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(t_exp1, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

abline(h=0.25, col = "red")
merge = mergeCloseModules(t_exp1, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors = merge$colors
table(mergedColors)
mergedMEs = merge$newMEs

png("./plots/Fig4_dendrogram_and_module_colors_exp1.png")
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()


moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "./data/step6_networkConstruction_stepByStep_exp1.RData")

#Modules-traits relationships
nGenes = ncol(t_exp1)
nSamples = nrow(t_exp1)

design=model.matrix(~0+ group1)
colnames(design)=levels(group1)
moduleTraitCor = cor(MEs, design, use = "p")
head(moduleTraitCor)

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 
                                     nSamples)
head(moduleTraitPvalue)

textMatrix =  paste(signif(moduleTraitCor, 2), 
                    "\n(",signif(moduleTraitPvalue, 1),
                    ")", 
                    sep = "");
dim(textMatrix) = dim(moduleTraitCor)

design1 <- as.data.frame(design)
png("./plots/Fig4_module_traits_relationship_exp1.png", width = 480, height = 720)
par(mar = c(5, 10, 4, 2));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#Module membership
trait = as.data.frame(design1$Diseased);
names(trait) = "trait"

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(t_exp1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(t_exp1, trait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(trait), sep="")
names(GSPvalue) = paste("p.GS.", names(trait), sep="")

module = "turquoise"
column = match(module, modNames);
moduleGenes = mergedColors == module;

png("./plots/Fig4_GS_MM_exp1.png")
par(mar = c(4,4,4,4));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

PD_module <- colnames(t_exp1)[mergedColors == "turquoise"]
PD_module
save(PD_module, file = "./data/step6_WGCNA_PD_module.RData")

###WGCNA-pSS
apply(exp2, 1, sd) %>% sort(.,decreasing = T) %>% head(., 10000) %>% names %>% exp2[.,] -> WGCNA_exp2
t_exp2 <- t(WGCNA_exp2)

allowWGCNAThreads()
adjacency = adjacency(t_exp2, power = po_exp2)
TOM = TOMsimilarity(adjacency)
save(TOM, file = paste0("./data/step6_TOM_exp2_po", po_exp2, ".RData"))
load(paste0("./data/step6_TOM_exp2_po", po_exp2, ".RData"))
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 50) 
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(t_exp2, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

abline(h=0.25, col = "red")
merge = mergeCloseModules(t_exp2, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors = merge$colors
table(mergedColors)
mergedMEs = merge$newMEs

png("./plots/Fig4_dendrogram_and_module_colors_exp2.png")
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "./data/step6_networkConstruction_stepByStep_exp2.RData")

#Modules-traits relationships
nGenes = ncol(t_exp2)
nSamples = nrow(t_exp2)
group2 <- factor(group2)
design=model.matrix(~0+ group2)
colnames(design)=levels(group2)
moduleTraitCor = cor(MEs, design, use = "p")
head(moduleTraitCor)

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 
                                     nSamples)

head(moduleTraitPvalue)

textMatrix =  paste(signif(moduleTraitCor, 2), 
                    "\n(",signif(moduleTraitPvalue, 1),
                    ")", 
                    sep = "");
dim(textMatrix) = dim(moduleTraitCor)
textMatrix[1:4,1:2]

design2 <- as.data.frame(design)
png("./plots/Fig4_module_traits_relationship_exp2.png", width = 480, height = 720)
par(mar = c(5, 10, 4, 2));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#Module membership
trait = as.data.frame(design2$pSS);
names(trait) = "trait"

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(t_exp2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(t_exp2, trait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(trait), sep="")
names(GSPvalue) = paste("p.GS.", names(trait), sep="")

module = "salmon"
column = match(module, modNames);
moduleGenes = mergedColors == module;

png("./plots/Fig4_GS_MM_exp2.png")
par(mar = c(4,4,4,4));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

pSS_module <- colnames(t_exp2)[mergedColors == "salmon"]
pSS_module
save(pSS_module, file = "./data/step6_WGCNA_pSS_module.RData")

###Venn
library(VennDiagram)

png(filename = "./plots/Fig5_venn_WCGNA.png", width = 480, height = 480)
venn.plot <- venn.diagram(
  x = list(set1 = PD_module, set2 = pSS_module),
  category.names = c("PD_module", "pSS_module"),
  cat.cex = 1.5,
  margin = 0.1,
  cex = 1.5,
  scaled = F,
  filename = NULL,
  output = TRUE,
  fill = c("#2CA02CFF", "#D62728FF"),
  cat.col = c("#2CA02CFF", "#D62728FF")
)
grid.draw(venn.plot)
dev.off()
common_WGCNA_genes <- intersect(PD_module, pSS_module)


load("./data/step3_deg.Rdata")
common_DEG_genes <- intersect(DEG1$name[DEG1$change != "stable"], DEG2$name[DEG2$change != "stable"])
png(filename = "./plots/Fig5_venn_W&D.png", width = 480, height = 480)
venn.plot <- venn.diagram(
  x = list(set1 = common_DEG_genes, set2 = common_WGCNA_genes),
  category.names = c("common_DEG_genes", "common_WGCNA_genes"),
  cat.cex = 1.5,
  margin = 0.15,
  cex = 1.5,
  scaled = F,
  filename = NULL,
  output = TRUE,
  fill = c("#2CA02CFF", "#D62728FF"),
  cat.col = c("#2CA02CFF", "#D62728FF")
)
grid.draw(venn.plot)
dev.off()

gene_for_ML <- intersect(common_DEG_genes, common_WGCNA_genes)
gene_for_ML

###Enrichment analysis on intersecting genes
keytypes(org.Hs.eg.db)
k <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID") 
e2s <- AnnotationDbi::select(org.Hs.eg.db, keys = k,
                             columns="SYMBOL", keytype = "ENTREZID")
ENTREZID_commonWGCNA <- e2s[match(common_WGCNA_genes, e2s$SYMBOL),"ENTREZID"] %>% na.omit %>% as.numeric

ont <- c("BP", "MF", "CC")
for(i in ont){
  go <- enrichGO(ENTREZID_commonWGCNA,
                 OrgDb = "org.Hs.eg.db",
                 ont = i,
                 minGSSize = 0,
                 pvalueCutoff = 0.05)
  p <- enrichplot::dotplot(go, showCategory = 10) +
    labs(title = i) +
    theme(plot.title = element_text(size = 20))
  ggsave(p,
         filename = paste0("./plots/Fig5_ORA_GO_", i, ".png"),
         width = 6,
         height = 6)
}

kegg <- enrichKEGG(ENTREZID_commonWGCNA,
                   organism = "hsa",
                   pAdjustMethod = "none",
                   minGSSize = 0,
                   qvalueCutoff = 0.9)
p <- enrichplot::dotplot(kegg, showCategory = 12) +
  labs(title = "KEGG enrichment") +
  theme(plot.title = element_text(size = 20))
p
ggsave(p,
       filename = paste0("./plots/Fig5_ORA_KEGG.png"),
       width = 6, height = 6)





###Machine learning - PD
#LASSO
#genes for machine learning(expression matrix)
exprSet = exp1[gene_for_ML,]

x = as.matrix(t(exprSet))
#y represents healthy or disease
y = ifelse(group1 == "Healthy", 0, 1)

library(glmnet)
set.seed(21)
png("./plots/Fig6_lasso_1.png", width = 960, height = 480)
par(mfrow = c(1, 2))
fit <- glmnet(x=x, y=y)
plot(fit,xvar = "lambda")

cv_fit <- cv.glmnet(x=x, y=y)
plot(cv_fit)
dev.off()

model_lasso_min <- glmnet(x=x, y=y, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, lambda=cv_fit$lambda.1se)


idx <- model_lasso_min$beta@i + 1       
gene_names <- model_lasso_min$beta@Dimnames[[1]]
lasso_PD_min <- gene_names[idx]
lasso_PD_min

idx <- model_lasso_1se$beta@i + 1       
gene_names <- model_lasso_1se$beta@Dimnames[[1]]
lasso_PD_se <- gene_names[idx]
lasso_PD_se






#Random Forest
library(randomForest)

x = t(exprSet)
y = factor(group1)
set.seed(21) 
rf_output = randomForest(x=x, y=y, importance = TRUE, ntree = 500,
                         proximity=TRUE )

png("./plots/Fig6_RF_1.png", width = 1080, height = 480)
par(mfrow = c(1, 2))
plot(rf_output)
varImpPlot(rf_output, type=2, n.var=10, scale=FALSE,
           main="MeanDecreaseGini",
           cex =2)
dev.off()

rf_importances = importance(rf_output, scale=FALSE)
head(rf_importances, 10)

#top10
rf_PD = rev(rownames(tail(rf_importances[order(rf_importances[,4]),],10)))
rf_PD


common_ML_PD <- intersect(lasso_PD_se,rf_PD)
common_ML_PD


#Venn
library(ggvenn)
venn_plot <- ggvenn(
  list(LASSO = lasso_PD_se, Random_forest = rf_PD),
  fill_color = c("green","lightblue"),
  stroke_size = 0.5,
  text_size = 4,
  set_name_size = 4,
  show_percentage = FALSE
)

venn_plot <- venn_plot +
  ggtitle(
    paste("Common genes:", 
          paste(common_ML_PD, collapse = ", "))
  )

ggsave("./plots/Fig6_venn_exp1.png", venn_plot, width = 6, height = 6)




###Machine learning - pSS
#LASSO
set.seed(21)
exprSet = exp2[gene_for_ML,]

x = as.matrix(t(exprSet))
y = ifelse(group2 == "Control", 0, 1)

library(glmnet)
png("./plots/Fig6_lasso_2.png", width = 960, height = 480)
par(mfrow = c(1, 2))
fit <- glmnet(x=x, y=y)
plot(fit,xvar = "lambda")

cv_fit <- cv.glmnet(x=x, y=y ,nfolds = 10)
plot(cv_fit)
dev.off()

model_lasso_min <- glmnet(x=x, y=y, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, lambda=cv_fit$lambda.1se)


idx <- model_lasso_min$beta@i + 1       
gene_names <- model_lasso_min$beta@Dimnames[[1]]
lasso_pSS_min <- gene_names[idx]
lasso_pSS_min


idx <- model_lasso_1se$beta@i + 1       
gene_names <- model_lasso_1se$beta@Dimnames[[1]]
lasso_pSS_se <- gene_names[idx]
lasso_pSS_se


#Random forest
library(randomForest)

x = t(exprSet)
y = factor(group2)
set.seed(21)
rf_output = randomForest(x=x, y=y, importance = TRUE, ntree = 500,
                         proximity=TRUE )

png("./plots/Fig6_RF_2.png", width = 1080, height = 480)
par(mfrow = c(1, 2))
plot(rf_output)
varImpPlot(rf_output, type=2, n.var=10, scale=FALSE, 
           main="MeanDecreaseGini",
           cex = 2)
dev.off()

rf_importances = importance(rf_output, scale=FALSE)
head(rf_importances, 10)


rf_pSS = rev(rownames(tail(rf_importances[order(rf_importances[,4]),],10)))
rf_pSS


common_ML_pSS <- intersect(lasso_pSS_min,rf_pSS)
common_ML_pSS

#Venn
venn_plot <- ggvenn(
  list(LASSO = lasso_pSS_min, Random_forest = rf_pSS),
  fill_color = c("green","lightblue"),
  stroke_size = 0.5,
  text_size = 4,
  set_name_size = 4,
  show_percentage = FALSE
)

venn_plot <- venn_plot + 
  ggtitle(
    paste("Common genes:", 
          paste(common_ML_pSS, collapse = ", "))
  )

ggsave("./plots/Fig6_venn_exp2.png",
       venn_plot,
       width = 6, height = 6)





#final veen
venn_plot <- ggvenn(
  list(PD = common_ML_PD, pSS = common_ML_pSS),
  fill_color = c("green", "red"),
  stroke_size = 0.5,
  text_size = 4,
  set_name_size = 4,
  show_percentage = FALSE
)

venn_plot + 
  ggtitle(
    paste("Common genes:", 
          paste(intersect(common_ML_PD, common_ML_pSS), collapse = ", "))
  )
ggsave("./plots/Fig6_venn_biomarker.png")



