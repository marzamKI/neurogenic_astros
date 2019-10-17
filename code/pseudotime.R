source("./code/global.R", local = FALSE)

#astros_clean <- readRDS("astros_clean_filtered_tsne_dm.rds")
astros_clean <- RunDiffusion(astros_clean, 
                             dims.use = c(1, 3:10), 
                             q.use = 0.001, 
                             scale.clip = 10, 
                             max.dim = 3,
                             reduction.use = "pca", 
                             reduction.name = "dm")

FeaturePlot(astros_clean,  
            reduction.use = "dm", 
            cols.use = my.inferno, 
            dim.1 = 1, 
            dim.2 = 3,
            c("Aqp4","Ascl1", "Mki67",
              "Dcx", "Gad1", "Htr3a",
              "Ifit1", "Iigp1", "Stat1"))

# 3D plot
dm.emb <- list()
for(i in 1:3){
  dm.emb[[i]]<-GetCellEmbeddings(astros_clean, 
                                 reduction.type = "dm", 
                                 dims.use = i,
                                 cells.use = NULL)
  
}

colors <- ifelse((astros_clean@ident == "AC1"), 
                 "#A52C60FF",
                 ifelse((astros_clean@ident == "AC2"),
                        "#2274A5", 
                        ifelse((astros_clean@ident == "NB"), 
                               "#5FAD56", 
                               "#F7B32B")
                 )
)

scatterplot3d(dm.emb[[1]], 
              dm.emb[[2]], 
              dm.emb[[3]], 
              grid=TRUE, 
              box=FALSE, 
              color = colors, 
              pch = 16, 
              angle = 20)

# Plot ribo score
ribo.genes <- c(grep(pattern = "^Rpl", 
                     x = rownames(x = astros_clean@data), 
                     value = TRUE), 
                grep(pattern = "^Rps", 
                     x = rownames(x = astros_clean@data), 
                     value = TRUE))

percent.ribo <- Matrix::colSums(astros_clean@raw.data[ribo.genes,])/Matrix::colSums(astros_clean@raw.data)

astros_clean <- AddMetaData(object = astros_clean, 
                            metadata = percent.ribo, 
                            col.name = "percent.ribo")

dm_emb <- cbind(astros_clean@dr$dm@cell.embeddings, 
                astros_clean@meta.data)

# Plot cell cycle and ribosome score 
p <- list()
for(i in c("S.Score", "G2M.Score", "percent.ribo")){
  p[[i]] <- ggplot(dm_emb, 
              aes_string(x="DM1", 
                  y="DM3", 
                  color = i)) +
    geom_point(aes_string(fill=i), 
               colour = alpha("black", 0.2),
               pch=21, 
               size=2.5, 
               lwd=0) +
    scale_fill_viridis(option = "inferno", 
                       begin = 0.1, 
                       end = 1) +
    theme +
    ggtitle(i)
  
}

plot_grid(plotlist = p)

## SNN graph demonstrates small NB cluster is not connected with AC3.
astros <- astros_clean
astros@dr$tsne@cell.embeddings <- cbind(dm_emb$DM1, dm_emb$DM3) #Substitute tsne coordinates with dm
astros_clean <- BuildSNN(astros, 
                         genes.use = astros@var.genes, 
                         reduction.type = "pca", 
                         dims.use = 1:10, 
                         plot.SNN = T,
                         force.recalc = T)

# Run Monocle standard pipeline
astros_clean@meta.data$cluster <- astros_clean@ident
astros_clean@meta.data$dm1 <- dm_emb$DM1

astro_clean_no_reactive <- SubsetData(astros_clean, 
                                      ident.remove = "AC3", 
                                      do.clean = T)
 
astro_clean_no_reactive <- FindVariableGenes(astro_clean_no_reactive, 
                                             mean.function = ExpMean, 
                                             dispersion.function = LogVMR, 
                                             x.low.cutoff = 0.0125, 
                                             x.high.cutoff = 5, 
                                             y.cutoff = 0.05, 
                                             do.plot = F)
length(astro_clean_no_reactive@var.genes)

astro_clean_no_reactive <- ScaleData(astro_clean_no_reactive, 
                                     vars.to.regress = c("percent.mito", "nUMI"))
astro_clean_no_reactive <- RunPCA(astro_clean_no_reactive, 
                                  pc.genes = astro_clean_no_reactive@var.genes)

monocle <- importCDS(astro_clean_no_reactive, 
                     import_all = TRUE)
monocle <- estimateSizeFactors(monocle)

monocle <- suppressWarnings(estimateDispersions(monocle))

monocle <- detectGenes(monocle, 
                       min_expr = 0.01)

monocle <- setOrderingFilter(monocle, 
                             astro_clean_no_reactive@var.genes)
monocle <- reduceDimension(monocle, 
                           max_components=2, 
                           method = 'DDRTree') #residualModelFormulaStr

monocle <- orderCells(monocle)
plot_cell_trajectory(monocle, color_by = "Pseudotime")

# Use DM1 as pseudotime score
monocle@phenoData@data$Pseudotime <- astro_clean_no_reactive@meta.data$dm1

diff_test_res <- differentialGeneTest(monocle[astro_clean_no_reactive@var.genes,], 
                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
                                      reducedModelFormulaStr = "~1", 
                                      relative_expr = TRUE)

sig_gene_names <- row.names (subset(diff_test_res, 
                                    qval < 0.01))

top_gene_names <- row.names(diff_test_res[head(order(diff_test_res$qval),1000),])

monocle_hm <- plot_pseudotime_heatmap(monocle[top_gene_names,], 
                                      return_heatmap=T, 
                                      num_clusters = 5,
                                      hmcols = viridis(1000, 
                                                       option = "inferno", 
                                                       alpha = 1, 
                                                       begin = 0))

#saveRDS(monocle, "./rds/monocle2.rds")

top.genes <- as.data.frame(cutree(monocle_hm$tree_row, k=5))

colnames(top.genes) <- "Cluster"

top.genes$Gene <- rownames(top.genes)

#write.csv(top.genes, "pseudotime.top1000genes_k5.csv")

