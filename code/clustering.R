source("./code/global.R", local = FALSE)

# Load astrocyte and neuroblast sample
astro_npcs <- readRDS("./rds/raw_astros_npcs.rds")

# Normalize and find variable genes
astro_npcs <- NormalizeData(astro_npcs)

astro_npcs <- FindVariableGenes(astro_npcs, 
                                x.low.cutoff = 0.0125, 
                                x.high.cutoff = 5, 
                                y.cutoff = 0.5, 
                                do.plot = F)

length(x =astro_npcs@var.genes)

astro_npcs <- ScaleData(astro_npcs, 
                        vars.to.regress = c("percent.mito", "nUMI"))

# Compute PCA and tSNE
astro_npcs <- RunPCA(astro_npcs, 
                     pc.genes = astro_npcs@var.genes)

PCElbowPlot(astro_npcs)
PCAPlot(astro_npcs)

VizPCA(astro_npcs,
       pcs.use = 1:8, 
       nCol = 4)

astro_npcs <- RunTSNE(astro_npcs, 
                      dims.use = 1:10,
                      do.fast = F, 
                      check_duplicates = F, 
                      perplexity = 100, 
                      seed.use = 2)

astro_npcs <- FindClusters(astro_npcs, 
                           reduction.type = "pca", 
                           dims.use = 1:10, 
                           resolution = 0.2)
TSNEPlot(astro_npcs)

#saveRDS(astro_npcs, "./rds/tsne_astros.rds")

# DEA
markers <- FindAllMarkers(astro_npcs, 
                          only.pos = T)

#write.csv(markers, "markers_astro.ident.csv")

old.ident <- c("0", "1", "2", "3", "4")
new.ident <- c("AC1", "AC2", "AC4", "NB", "AC3")
astro_npcs@ident <- plyr::mapvalues(x = astro_npcs@ident, 
                                    from = old.ident, 
                                    to = new.ident)

# Plot DEGs
notch.genes <- c("Hes5", "Hes1", "Nrarp", "Hey1", "Ubc", "Tmem100", "Ascl1", "Hes6", "Tmsb4x")
neurog.genes <- c("Fzd2", "Bcl2", "Chrdl1", "Vcl", "Id4", "Btg2", "Egfr", "Olig2", "Ednrb")

VlnPlot(astro_npcs, 
        notch.genes, 
        point.size.use = 0, 
        cols.use = my.palette)

VlnPlot(astro_npcs, 
        neurog.genes, 
        point.size.use = 0, 
        cols.use = my.palette)

#saveRDS(astro_npcs, "./rds/tsne.astros.rds")

# Clear dataset of low quality cells and re-run tSNE & DEA
astros_clean <- SubsetData(astro_npcs, 
                           ident.remove = "AC4") # Possibly low quality cells

astros_clean <- FilterCells(astros_clean, 
                            subset.names = c("Pllp", "Mog", "Mag", "Vtn"), 
                            high.thresholds = c(rep(1.5, 3), 2))

astros_clean <- RunTSNE(astros_clean, 
                        dims.use = 1:10, 
                        seed.use = 7)

astros_clean@ident <- factor(x = astros_clean@ident, 
                             levels = c("AC1", "AC2", "AC3", "NB"))
TSNEPlot(astros_clean)

markers <- FindAllMarkers(astros_clean, 
                          only.pos = T)

adj.markers <- filter(markers, 
                      p_val_adj < 0.05)

top.50 <- adj.markers %>%
  group_by(cluster) %>% 
  top_n(50, avg_logFC)
rownames(top.50) <- top.50$gene

cluster.averages <- AverageExpression(astros_clean, 
                                      return.seurat = T)

averages_hm_genes <- cluster.averages@scale.data[top.50$gene,] %>% 
  as.data.frame()
averages_hm_genes$gene <- rownames(averages_hm_genes)

# Order genes by avg_fc across groups
cl.avg <- list()
for(i in levels(astros_clean@ident)) {
  
  genes <- top.50 %>% 
    filter(cluster == i)
  
  cl.avg[[i]] <- averages_hm_genes %>% 
    filter(gene %in% genes$gene) %>%
    select(i, "gene")
  
  colnames(cl.avg[[i]]) = c("AVG", "gene")
}

merge_cl_avg = do.call(rbind, cl.avg)
top.50_avg <- merge(merge_cl_avg, top.50)
ordered.50 <- top.50_avg %>% 
  arrange(desc(AVG)) %>% 
  arrange(cluster)

DoHeatmap(cluster.averages, genes.use = ordered.50$gene) +
  scale_fill_gradient2(mid = "#A52C60FF", high = "#F7D03CFF", low = "#1B0C42FF")

#pdf("expression_top50_markers.pdf", width = 5, height = 5, paper = 'special')
#genes <- top.50$gene
#for (i in 1:length(x = genes)) {
#  FeaturePlot(astros_clean, 
#              features.plot = genes[[i]])
#}
#dev.off()

markers_AC2_AC1 <- FindMarkers(astros_clean, 
                               ident.1 = "AC2", 
                               ident.2 = "AC1", 
                               only.pos = FALSE)

#saveRDS(astro_npcs, "./rds/tsne.filtered_astros.rds")






astros_clean <- RunDiffusion(astros_clean, dims.use = c(1, 3:15), 
                          q.use = 0.001, scale.clip = 3, max.dim = 3,
                           reduction.use = "pca", reduction.name = "dm")
DimPlot(astros_clean, reduction.use = "dm", dim.1 = 1, dim.2 = 2)
DimPlot(astros_clean, reduction.use = "dm", dim.1 = 1, dim.2 = 3)

FeaturePlot(astros_clean,  reduction.use = "dm", 
            cols.use = my.inferno, no.axes = T, dim.1 = 1, dim.2 = 3,
            pt.size = 2,
            c("Aqp4","Ascl1", "Slc38a1", "Chd7", "Pak3", "Nell2", "Mki67",
              "Nnat", "Traf4",
              "Dcx", "Gad1", "Htr3a",
              "Ifit1", "Iigp1", "Stat1"))

library("scatterplot3d") # load
dm1<-GetCellEmbeddings(astros_clean, reduction.type = "dm", dims.use = 1,
                        cells.use = NULL)
dm2<-GetCellEmbeddings(astros_clean, reduction.type = "dm", dims.use = 2,
                        cells.use = NULL)
dm3<-GetCellEmbeddings(astros_clean, reduction.type = "dm", dims.use = 3,
                        cells.use = NULL)

colors <- ifelse((astros_clean@ident == "AC1"), "#A52C60FF",
                 ifelse((astros_clean@ident == "AC2"),"#2274A5", 
                        ifelse((astros_clean@ident == "NB"), "#5FAD56", "#F7B32B")))

scatterplot3d(dm2, dm3, dm1, grid=TRUE, box=FALSE, color = colors, pch = 16, angle = 50)

# Look at small cluster
select.cells <- DimPlot(astros_clean, reduction.use = "dm", dim.1 = 1, dim.2 = 3, do.identify = TRUE)
astros_clean <- SetIdent(astros_clean, cells.use = select.cells, ident.use = "NewCells")

# Now, we find markers that are specific to the new cells, and find clear DC
# markers
newcells.markers_vsAC3 <- FindMarkers(astros_clean, ident.1 = "NewCells", ident.2 = "AC3", 
                                min.diff.pct = 0.3, only.pos = TRUE)
newcells.markers_vsAC3$gene <- rownames(newcells.markers_vsAC3)
newcells.markers_vsAC3 <- filter(newcells.markers_vsAC3, p_val_adj < 0.05)
newcells.markers_vsNB <- FindMarkers(astros_clean, ident.1 = "NewCells", ident.2 = "NB", 
                                      min.diff.pct = 0.3, only.pos = TRUE)
newcells.markers_vsNB$gene <- rownames(newcells.markers_vsNB)
newcells.markers_vsNB <- filter(newcells.markers_vsNB, p_val_adj < 0.05)
newcells.markersAC3_vsNB <- FindMarkers(astros_clean, ident.1 = c("NewCells", "AC3"), ident.2 = "NB", 
                                     min.diff.pct = 0.3, only.pos = TRUE)
newcells.markersAC3_vsNB$gene <- rownames(newcells.markersAC3_vsNB)
newcells.markersAC3_vsNB <- filter(newcells.markersAC3_vsNB, p_val_adj < 0.05)

pdf("newcells.markers_vsNB.pdf", width = 5, height = 5, paper = 'special')
genes <- newcells.markers_vsNB$gene
features <- genes[[1]]
for (i in 1:length(x = genes)) {
  FeaturePlot(astros_clean, features.plot = genes[[i]], reduction.use = "dm", dim.1 = 1, dim.2 = 3)}
dev.off()

pdf("newcells.markersAC3_vsNB.pdf", width = 5, height = 5, paper = 'special')
genes <- newcells.markersAC3_vsNB$gene
features <- genes[[1]]
for (i in 1:length(x = genes)) {
  FeaturePlot(astros_clean, features.plot = genes[[i]], reduction.use = "dm", dim.1 = 1, dim.2 = 3)}
dev.off()


newcells.markers <- FindMarkers(astros_clean, ident.1 = "NewCells", 
                                        min.diff.pct = 0.3, only.pos = TRUE)
newcells.markers$gene <- rownames(newcells.markers)
newcells.markers <- filter(newcells.markers, p_val_adj < 0.05)
pdf("newcells.markers.pdf", width = 5, height = 5, paper = 'special')
genes <- newcells.markers$gene
features <- genes[[1]]
for (i in 1:length(x = genes)) {
  FeaturePlot(astros_clean, features.plot = genes[[i]], reduction.use = "dm", dim.1 = 1, dim.2 = 3)}
dev.off()

astros_clean <- SetAllIdent(astros_clean, "res.0.2")
saveRDS(astros_clean, "astros_clean_filtered_tsne_dm.rds")

pdf("all_AC1AC2.markers.pdf", width = 5, height = 5, paper = 'special')
genes <- rownames(markers_AC1vsAC2)
features <- genes[[1]]
for (i in 1:length(x = genes)) {
  FeaturePlot(astros_clean, features.plot = genes[[i]], reduction.use = "tsne", cols.use = my.inferno)}
dev.off()


averages_hm_genes <- cluster.averages@scale.data[top.50$gene,]

top.50$cl.average <- c(averages_hm_genes[1:50,1], 
                       averages_hm_genes[51:100,2],
                       averages_hm_genes[101:150,3],
                       averages_hm_genes[151:200,4])

ordered.50 <- top.50 %>% arrange(desc(cl.average)) %>% arrange(cluster)

DoHeatmap(cluster.averages, genes.use = ordered.50$gene, group.label.rot = TRUE, group.cex = 0, 
          col.mid = "#A52C60FF", col.high = "#F7D03CFF", col.low = "#1B0C42FF")

pdf("nsc_markers.pdf", width = 6, height = 6, paper = 'special')
genes <- as.character(nsc.genes[,1])
features <- genes[[1]]
for (i in 1:length(x = genes)) {
  FeaturePlot(astros_clean, 
              features.plot = genes[[i]], 
              cols.use = my.inferno, 
              reduction.use = "dm", dim.1 = 1, dim.2 = 3)}
dev.off()
