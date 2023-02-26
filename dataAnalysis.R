## Proteomics data analysis
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(ComplexHeatmap)
library(RColorBrewer)
library(org.Hs.eg.db)

load("notebook/data/DataCleaned.RData")

sum(paste(sampleDetail$SlideId, sampleDetail$Subarray, sep = "_") == rownames(somaSel)) == nrow(somaSel)
sum(sampleDetail$SampleId == phenotype$ID) == nrow(phenotype)

#### PCA ####

pr <- prcomp(scale(log2(somaSel), center = TRUE, scale = TRUE))
spr <- summary(pr)

dplot <- data.frame(
  pr$x,
  group = phenotype$Grp,
  id = phenotype$ID,
  sex = phenotype$Sex,
  age = phenotype$Age,
  bmi = phenotype$BMI,
  stringsAsFactors = FALSE
)

nplot <- 5
idx <- 1
gps <- list()
axisTextSize <- 8
for(i in 1:nplot){
  for(j in 1:nplot){
    if(i==j){
      gp <- ggplot(dplot) + 
        geom_density(aes_string(x = paste0("PC", i))) + 
        labs(x = paste0("PC", i, "(", format(round(spr$importance[2,i]*100, 2), nsmall = 2), "%)"),
             y = "Density") +
        theme_light() + 
        theme(axis.text = element_text(size = axisTextSize))
      gps[[idx]] <- gp
      idx <- idx + 1
    } else if(i<j){
      gp <- ggplot(dplot) + 
        geom_point(aes_string(x=paste0("PC", j), y=paste0("PC",i), color="group"))+
        labs(x=paste0("PC",j), y=paste0("PC",i)) + 
        ggsci::scale_color_nejm() + 
        theme_light() + 
        theme(legend.title = element_blank(), axis.text = element_text(size = axisTextSize))
      gps[[idx]] <- gp
      idx <- idx + 1
    } else {
      gps[[idx]] <- NULL
      idx <- idx + 1
    }
  }
}

gp <- ggpubr::ggarrange(
  plotlist = gps,
  nrow = nplot, ncol = nplot,
  common.legend = TRUE, legend = "bottom"
)


png("figures/PCAGrp.png", width = 9, height = 8, units = "in", res = 300)
print(gp)
invisible(dev.off())




nplot <- 5
idx <- 1
gps <- list()
axisTextSize <- 8
for(i in 1:nplot){
  for(j in 1:nplot){
    if(i==j){
      gp <- ggplot(dplot) + 
        geom_density(aes_string(x = paste0("PC", i))) + 
        labs(x = paste0("PC", i, "(", format(round(spr$importance[2,i]*100, 2), nsmall = 2), "%)"),
             y = "Density") +
        theme_light() + 
        theme(axis.text = element_text(size = axisTextSize))
      gps[[idx]] <- gp
      idx <- idx + 1
    } else if(i<j){
      gp <- ggplot(dplot) + 
        geom_point(aes_string(x=paste0("PC", j), y=paste0("PC",i), color="sex"))+
        labs(x=paste0("PC",j), y=paste0("PC",i)) + 
        ggsci::scale_color_nejm() + 
        theme_light() + 
        theme(legend.title = element_blank(), axis.text = element_text(size = axisTextSize))
      gps[[idx]] <- gp
      idx <- idx + 1
    } else {
      gps[[idx]] <- NULL
      idx <- idx + 1
    }
  }
}

gp <- ggpubr::ggarrange(
  plotlist = gps,
  nrow = nplot, ncol = nplot,
  common.legend = TRUE, legend = "bottom"
)


png("figures/PCASex.png", width = 9, height = 8, units = "in", res = 300)
print(gp)
invisible(dev.off())


gp <- ggplot(dplot) + 
  geom_point(aes(x=PC1, y=PC2, color = sex, shape = group)) + 
  labs(x=paste0("PC1 (", format(round(spr$importance[2,1]*100, 2), nsmall = 2), "%)"),
       y=paste0("PC2 (", format(round(spr$importance[2,2]*100, 2), nsmall = 2), "%)")) + 
  ggsci::scale_color_nejm() + 
  theme_light() + 
  theme(legend.title = element_blank(), legend.position = "bottom")

pdf("figures/PC12_Grp_Sex.pdf", width = 5, height = 5)
print(gp)
dev.off()

# p = 0.026
t.test(PC1 ~ sex, data = dplot)
# p = 0.83
t.test(PC1 ~ group, data = dplot)



#### All ####
sex <- ifelse(phenotype$Sex=="F", 0, 1)
age <- phenotype$Age
bmi <- phenotype$BMI
rls <- ifelse(phenotype$Grp == "Control", 0, 1)

pvalue <- NULL
beta <- NULL
for(i in 1:ncol(somaSel)){
  y <- log2(somaSel[,i])
  
  tmp <- lm(y ~ sex + age + bmi + rls)
  tmp <- summary(tmp)
  
  pvalue <- c(pvalue, tmp$coefficients[5,4])
  beta <- c(beta, tmp$coefficients[5,1])
}

rlt <- data.frame(
  AptName = anaInfoSel$AptName,
  UniProt = anaInfoSel$UniProt,
  EntrezGeneID = anaInfoSel$EntrezGeneID,
  EntrezGeneSymbol = anaInfoSel$EntrezGeneSymbol,
  pvalue = pvalue,
  beta = beta,
  fdr = p.adjust(pvalue, method = "fdr"),
  stringsAsFactors = FALSE
)

write.csv(rlt, file = "rlt/pvProteinAdjAll.csv", quote = TRUE, row.names = FALSE)
saveRDS(rlt, file = "rlt/pvProteinAdjAll.rds")

dplot <- data.frame(
  UniProt = rlt$UniProt,
  Symbol = rlt$EntrezGeneSymbol,
  beta = rlt$beta,
  pvalue = -log10(rlt$pvalue),
  stringsAsFactors = FALSE
)

pvCutoff <- -log10(0.05)
betaCutoff <- log2(1.5)

dplot$group <- "G0"
dplot$group[dplot$pvalue> pvCutoff & dplot$beta < -betaCutoff] <- "G1"
dplot$group[dplot$pvalue> pvCutoff & dplot$beta > betaCutoff] <- "G2"

dplotText <- dplot[dplot$group!="G0",]

gp <- ggplot(dplot) +
  geom_point(aes(x=beta, y = pvalue, color = group), size = 1) + 
  geom_vline(xintercept = c(-betaCutoff,betaCutoff), color = "grey", linetype = 2) + 
  geom_hline(yintercept = c(pvCutoff), color = "grey", linetype = 2) + 
  geom_text_repel(data = dplotText, aes(x=beta, y=pvalue, label = Symbol, color = group), size = 3)+
  scale_color_manual(values = c(G0="grey", G1="#0072B5FF", G2="#BC3C29FF")) + 
  labs(x= "Log<sub>2</sub> Fold Change", y = "-Log<sub>10</sub> <i> P</i>") + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.title.x = element_markdown(), axis.title.y = element_markdown())

pdf("figures/volcanoAdjAll.pdf", width = 5, height = 4)
print(gp)
dev.off()


#### Male ####
idxM <- phenotype$Sex == "M"
age <- phenotype$Age[idxM]
bmi <- phenotype$BMI[idxM]
rls <- ifelse(phenotype$Grp[idxM] == "Control", 0, 1)

pvalue <- NULL
beta <- NULL
for(i in 1:ncol(somaSel)){
  y <- log2(somaSel[idxM,i])
  
  tmp <- lm(y ~ age + bmi + rls)
  tmp <- summary(tmp)
  
  pvalue <- c(pvalue, tmp$coefficients[4,4])
  beta <- c(beta, tmp$coefficients[4,1])
}

rlt <- data.frame(
  AptName = anaInfoSel$AptName,
  UniProt = anaInfoSel$UniProt,
  EntrezGeneID = anaInfoSel$EntrezGeneID,
  EntrezGeneSymbol = anaInfoSel$EntrezGeneSymbol,
  pvalue = pvalue,
  beta = beta,
  fdr = p.adjust(pvalue, method = "fdr"),
  stringsAsFactors = FALSE
)

write.csv(rlt, file = "rlt/pvProteinAdjMale.csv", quote = TRUE, row.names = FALSE)
saveRDS(rlt, file = "rlt/pvProteinAdjMale.rds")

dplot <- data.frame(
  UniProt = rlt$UniProt,
  Symbol = rlt$EntrezGeneSymbol,
  beta = rlt$beta,
  pvalue = -log10(rlt$pvalue),
  stringsAsFactors = FALSE
)

pvCutoff <- -log10(0.05)
betaCutoff <- log2(1.5)

dplot$group <- "G0"
dplot$group[dplot$pvalue> pvCutoff & dplot$beta < -betaCutoff] <- "G1"
dplot$group[dplot$pvalue> pvCutoff & dplot$beta > betaCutoff] <- "G2"

dplotText <- dplot[dplot$group!="G0",]

gp <- ggplot(dplot) +
  geom_point(aes(x=beta, y = pvalue, color = group), size = 1) + 
  geom_text_repel(data = dplotText, aes(x=beta, y=pvalue, label = Symbol, color = group), size = 2, max.overlaps = 20)+
  geom_vline(xintercept = c(-betaCutoff,betaCutoff), color = "grey", linetype = 2) + 
  geom_hline(yintercept = c(pvCutoff), color = "grey", linetype = 2) + 
  scale_color_manual(values = c(G0="grey", G1="#0072B5FF", G2="#BC3C29FF")) + 
  labs(x= "Log<sub>2</sub> Fold Change", y = "-Log<sub>10</sub> <i> P</i>") + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.title.x = element_markdown(), axis.title.y = element_markdown())

pdf("figures/volcanoAdjMale.pdf", width = 5, height = 4)
print(gp)
dev.off()


#### Female ####

idxF <- phenotype$Sex == "F"
age <- phenotype$Age[idxF]
bmi <- phenotype$BMI[idxF]
rls <- ifelse(phenotype$Grp[idxF] == "Control", 0, 1)

pvalue <- NULL
beta <- NULL
for(i in 1:ncol(somaSel)){
  y <- log2(somaSel[idxF,i])
  
  tmp <- lm(y ~ age + bmi + rls)
  tmp <- summary(tmp)
  
  pvalue <- c(pvalue, tmp$coefficients[4,4])
  beta <- c(beta, tmp$coefficients[4,1])
}

rlt <- data.frame(
  AptName = anaInfoSel$AptName,
  UniProt = anaInfoSel$UniProt,
  EntrezGeneID = anaInfoSel$EntrezGeneID,
  EntrezGeneSymbol = anaInfoSel$EntrezGeneSymbol,
  pvalue = pvalue,
  beta = beta,
  fdr = p.adjust(pvalue, method = "fdr"),
  stringsAsFactors = FALSE
)

write.csv(rlt, file = "rlt/pvProteinAdjFemale.csv", quote = TRUE, row.names = FALSE)
saveRDS(rlt, file = "rlt/pvProteinAdjFemale.rds")

dplot <- data.frame(
  UniProt = rlt$UniProt,
  Symbol = rlt$EntrezGeneSymbol,
  beta = rlt$beta,
  pvalue = -log10(rlt$pvalue),
  stringsAsFactors = FALSE
)

pvCutoff <- -log10(0.05)
betaCutoff <- log2(1.5)

dplot$group <- "G0"
dplot$group[dplot$pvalue> pvCutoff & dplot$beta < -betaCutoff] <- "G1"
dplot$group[dplot$pvalue> pvCutoff & dplot$beta > betaCutoff] <- "G2"

dplotText <- dplot[dplot$group!="G0",]

gp <- ggplot(dplot) +
  geom_point(aes(x=beta, y = pvalue, color = group), size = 1) + 
  geom_text_repel(data = dplotText, aes(x=beta, y=pvalue, label = Symbol, color = group), size = 2, max.overlaps = 20)+
  geom_vline(xintercept = c(-betaCutoff,betaCutoff), color = "grey", linetype = 2) + 
  geom_hline(yintercept = c(pvCutoff), color = "grey", linetype = 2) + 
  scale_color_manual(values = c(G0="grey", G1="#0072B5FF", G2="#BC3C29FF")) + 
  labs(x= "Log<sub>2</sub> Fold Change", y = "-Log<sub>10</sub> <i> P</i>") + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.title.x = element_markdown(), axis.title.y = element_markdown())

pdf("figures/volcanoAdjFemale.pdf", width = 5, height = 4)
print(gp)
dev.off()


#### overlap ####


rlt <- readRDS("rlt/pvProteinAdjAll.rds")
rltMale <- readRDS("rlt/pvProteinAdjMale.rds")
rltFemale <- readRDS("rlt/pvProteinAdjFemale.rds")

write.csv(rlt[rlt$pvalue<0.05 & abs(rlt$beta) > log2(1.5),], "rlt/sigAdjAll.csv",
          row.names = FALSE, quote = TRUE)

write.csv(rltMale[rltMale$pvalue<0.05 & abs(rltMale$beta) > log2(1.5),], "rlt/sigAdjMale.csv",
          row.names = FALSE, quote = TRUE)

write.csv(rltFemale[rltFemale$pvalue<0.05 & abs(rltFemale$beta) > log2(1.5),], "rlt/sigAdjFemale.csv",
          row.names = FALSE, quote = TRUE)


x <- list(
  All = rlt$AptName[rlt$pvalue<0.05 & abs(rlt$beta) > log2(1.5)],
  Male = rltMale$AptName[rltMale$pvalue<0.05 & abs(rltMale$beta) > log2(1.5)],
  Female = rltFemale$AptName[rltFemale$pvalue<0.05 & abs(rltFemale$beta) > log2(1.5)]
)

pdf("figures/vennGeneAdj.pdf", width = 5, height = 4)
ggvenn::ggvenn(
  x, c("Male", "Female", "All"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4, text_size = 3,
)
dev.off()


## gene SLMAP
tmp <- intersect(x$Male, x$Female)
rltMale[rltMale$EntrezGeneSymbol == tmp,]
rltFemale[rltFemale$EntrezGeneSymbol == tmp,]

## M: RLS high
log2(mean(somaSel$seq.7191.32[phenotype$Sex == "M" & phenotype$Grp!="Control"])/
mean(somaSel$seq.7191.32[phenotype$Sex == "M" & phenotype$Grp=="Control"]))

## F: RLS low
log2(mean(somaSel[phenotype$Sex == "F" & phenotype$Grp!="Control",colnames(somaSel)=="seq.7191.32"])/
mean(somaSel[phenotype$Sex == "F" & phenotype$Grp=="Control",colnames(somaSel)=="seq.7191.32"]))



#### Heatmap ####
idxM <- phenotype$Sex == "M"
ma <- somaSel[which(idxM), which(colnames(somaSel) %in% x$Male)]
ma <- scale(as.matrix(ma))
ma <- t(ma)
rownames(ma) <- rltMale$EntrezGeneSymbol[rltMale$pvalue<0.05 & abs(rltMale$beta) > log2(1.5)]

col_fun <- circlize::colorRamp2(c(-4, 0, 4), brewer.pal(7, "RdYlBu")[c(7, 4, 1)])
colG <- c("#BC3C2999", "#0072B599")
names(colG) <- c("RLS", "Control")

ha <- columnAnnotation(
  Group = phenotype$Grp[which(idxM)],
  col = list(
    Group = colG
  ),
  annotation_legend_param = list(
    Group = list(nrow=1, direction = "horizontal")
  ),
  show_annotation_name = FALSE
)

#set.seed(123)
hm <- Heatmap(
  ma,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title = "Scaled Protein Level",
    direction = "horizontal",
    legend_width = unit(4, "cm")
  ),
  use_raster = TRUE,
  raster_device = "png",
  
  col = col_fun,
  
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  
  #column_km = 2,
  #column_km_repeats = 20,
  
  top_annotation = ha
)


pdf("figures/heatmapMale.pdf", width = 5, height = 8)
draw(
  hm,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  legend_gap = unit(1.5, "cm"),
  merge_legend = TRUE
)
dev.off()



idxF <- phenotype$Sex == "F"
ma <- somaSel[which(idxF), which(colnames(somaSel) %in% x$Female)]
ma <- scale(as.matrix(ma))
ma <- t(ma)
rownames(ma) <- rltMale$EntrezGeneSymbol[rltFemale$pvalue<0.05 & abs(rltFemale$beta) > log2(1.5)]

col_fun <- circlize::colorRamp2(c(-4, 0, 4), brewer.pal(7, "RdYlBu")[c(7, 4, 1)])
colG <- c("#BC3C2999", "#0072B599")
names(colG) <- c("RLS", "Control")

ha <- columnAnnotation(
  Group = phenotype$Grp[which(idxF)],
  col = list(
    Group = colG
  ),
  annotation_legend_param = list(
    Group = list(nrow=1, direction = "horizontal")
  ),
  show_annotation_name = FALSE
)

#set.seed(123)
hm <- Heatmap(
  ma,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title = "Scaled Protein Level",
    direction = "horizontal",
    legend_width = unit(4, "cm")
  ),
  use_raster = TRUE,
  raster_device = "png",
  
  col = col_fun,
  
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  
  #column_km = 2,
  #column_km_repeats = 20,
  
  top_annotation = ha
)


pdf("figures/heatmapFemale.pdf", width = 5, height = 4)
draw(
  hm,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  legend_gap = unit(1.5, "cm"),
  merge_legend = TRUE
)
dev.off()



#### Gene set enrichment analysis ####

rlt <- readRDS("rlt/pvProteinAdjAll.rds")
geneList <- rlt$beta
names(geneList) <- rlt$UniProt
geneList <- sort(geneList, decreasing = TRUE)

set.seed(123)
gKEGG <- clusterProfiler::gseKEGG(
  geneList     = geneList,
  organism     = 'hsa', 
  keyType = "uniprot",
  pvalueCutoff = 0.05,
  verbose      = FALSE)


gGO <- clusterProfiler::gseGO(
  geneList     = geneList,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  eps = 0,
  pvalueCutoff = 0.05,
  verbose      = FALSE)

save(gKEGG, gGO, file = "rlt/GSEAAdjAll.RData")

pdf("figures/dotPlotKEGGAdjAll.pdf", width = 7, height = 7)
#enrichplot::dotplot(gKEGG, x="Count", showCategory = gKEGG@result$Description[idx[1:10]])
enrichplot::dotplot(gKEGG)
dev.off()

pdf("figures/upsetKEGGAdjAll.pdf", width = 10, height = 7)
enrichplot::upsetplot(gKEGG)
dev.off()


pdf("figures/dotPlotGOAdjAll.pdf", width = 7, height = 7)
#enrichplot::dotplot(gKEGG, x="Count", showCategory = gKEGG@result$Description[idx[1:10]])
enrichplot::dotplot(gGO)
dev.off()

pdf("figures/upsetGOAdjAll.pdf", width = 10, height = 7)
enrichplot::upsetplot(gGO)
dev.off()





rlt <- readRDS("rlt/pvProteinAdjMale.rds")
geneList <- rlt$beta
names(geneList) <- rlt$UniProt
geneList <- sort(geneList, decreasing = TRUE)

set.seed(123)
gKEGG <- clusterProfiler::gseKEGG(
  geneList     = geneList,
  organism     = 'hsa', 
  keyType = "uniprot",
  pvalueCutoff = 0.05,
  verbose      = FALSE)


gGO <- clusterProfiler::gseGO(
  geneList     = geneList,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  eps = 0,
  pvalueCutoff = 0.05,
  verbose      = FALSE)

save(gKEGG, gGO, file = "rlt/GSEAAdjMale.RData")

# no significant for KEGG
# pdf("figures/dotPlotKEGGAdjMale.pdf", width = 7, height = 7)
# #enrichplot::dotplot(gKEGG, x="Count", showCategory = gKEGG@result$Description[idx[1:10]])
# enrichplot::dotplot(gKEGG)
# dev.off()
# 
# pdf("figures/upsetKEGGAdjMale.pdf", width = 10, height = 7)
# enrichplot::upsetplot(gKEGG)
# dev.off()

write.csv(gGO@result, "rlt/GOAdjMale.csv", row.names = FALSE, quote = TRUE)

pdf("figures/dotPlotGOAdjMale.pdf", width = 7, height = 7)
#enrichplot::dotplot(gKEGG, x="Count", showCategory = gKEGG@result$Description[idx[1:10]])
enrichplot::dotplot(gGO)
dev.off()

pdf("figures/upsetGOAdjMale.pdf", width = 15, height = 7)
enrichplot::upsetplot(gGO)
dev.off()


#enrichplot::dotplot(gGO, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)




rlt <- readRDS("rlt/pvProteinAdjFemale.rds")
geneList <- rlt$beta
names(geneList) <- rlt$UniProt
geneList <- sort(geneList, decreasing = TRUE)

set.seed(123)
gKEGG <- clusterProfiler::gseKEGG(
  geneList     = geneList,
  organism     = 'hsa', 
  keyType = "uniprot",
  pvalueCutoff = 0.05,
  verbose      = FALSE)


gGO <- clusterProfiler::gseGO(
  geneList     = geneList,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  eps = 0,
  pvalueCutoff = 0.05,
  verbose      = FALSE)

save(gKEGG, gGO, file = "rlt/GSEAAdjFemale.RData")


write.csv(gKEGG@result, file = "rlt/KEGGAdjFemale.csv", row.names = FALSE, quote = TRUE)

pdf("figures/dotPlotKEGGAdjFemale.pdf", width = 7, height = 7)
#enrichplot::dotplot(gKEGG, x="Count", showCategory = gKEGG@result$Description[idx[1:10]])
enrichplot::dotplot(gKEGG)
dev.off()

pdf("figures/upsetKEGGAdjFemale.pdf", width = 10, height = 7)
enrichplot::upsetplot(gKEGG)
dev.off()



write.csv(gGO@result, file = "rlt/GOAdjFemale.csv", row.names = FALSE, quote = TRUE)

pdf("figures/dotPlotGOAdjFemale.pdf", width = 7, height = 7)
#enrichplot::dotplot(gKEGG, x="Count", showCategory = gKEGG@result$Description[idx[1:10]])
enrichplot::dotplot(gGO)
dev.off()

pdf("figures/upsetGOAdjFemale.pdf", width = 15, height = 7)
enrichplot::upsetplot(gGO)
dev.off()










