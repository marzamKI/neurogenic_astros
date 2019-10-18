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
