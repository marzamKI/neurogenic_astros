source("./code/global.R", local = FALSE)

#svz <- readRDS("./rds/neurogenic_svz.rds")
#astros_clean <- readRDS(".rds/tsne.filtered_astros.rds")

astros_clean <- StashIdent(astros_clean, "ctx_cells")

## CC correction
# Filter cells
other_cells <- c("Pdgfra", "Mog", "Klk6", "Foxj1", "Flt1", "Vtn", "Acta2", "P2ry12", "Tac1", "Penk")

svz <- FilterCells(svz, 
                   other_cells, 
                   high.thresholds = rep(0.4, length(other_cells)))

svz <- SubsetData(svz, ident.remove = "6")

svz <- FindVariableGenes(svz, mean.function = ExpMean,
                         dispersion.function = LogVMR, 
                         x.low.cutoff = 0.05,
                         x.high.cutoff = 3.5, 
                         y.cutoff = 0.5, 
                         do.plot = F)

astros_clean <- FindVariableGenes(astros_clean, mean.function = ExpMean,
                                  dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.05,
                                  x.high.cutoff = 3.5, 
                                  y.cutoff = 0.5, 
                                  do.plot = F)

# CCA integration
hvg.svz <- svz@var.genes
hvg.ctx <- astros_clean@var.genes
hvg.union <- union(x = hvg.svz, y = hvg.ctx)

svz@meta.data[, "protocol"] <- "svz"
astros_clean@meta.data[, "protocol"] <- "ctx"

obj <- RunCCA(object = astros_clean, object2 = svz, genes.use = hvg.union)

p1 <- DimPlot(object = obj, 
              reduction.use = "cca", 
              group.by = "protocol", 
              pt.size = 0.5,
              do.return = TRUE)

p2 <- VlnPlot(object = obj, 
              features.plot = "CC1", 
              group.by = "protocol", 
              do.return = TRUE)

plot_grid(p1, p2)

obj <- CalcVarExpRatio(object = obj, 
                       reduction.type = "pca", 
                       grouping.var = "protocol",
                       dims.use = 1:13)

obj.all.save <- obj

obj <- SubsetData(object = obj, 
                  subset.name = "var.ratio.pca", 
                  accept.low = 0.5)

obj.discard <- SubsetData(object = obj.all.save, 
                          subset.name = "var.ratio.pca",
                          accept.high = 0.5)

median(x = obj@meta.data[, "nGene"])
median(x = obj.discard@meta.data[, "nGene"])

obj <- AlignSubspace(object = obj, 
                     reduction.type = "cca", 
                     grouping.var = "protocol",
                     dims.align = 1:13)

p1 <- VlnPlot(obj, 
              features.plot = "ACC1", 
              group.by = "protocol",
              do.return = TRUE)

p2 <- VlnPlot(obj, features.plot = "ACC2", 
              group.by = "protocol",
              do.return = TRUE)

plot_grid(p1, p2)

obj <- RunTSNE(object = obj, 
               reduction.use = "cca.aligned", 
               dims.use = 1:13,
               do.fast = TRUE, 
               check_duplicates = F)

obj <- FindClusters(object = obj, 
                    reduction.type = "cca.aligned", 
                    dims.use = 1:13,
                    save.SNN = TRUE, 
                    force.recalc = T, 
                    resolution = 0.5)

obj@meta.data$svz_cells[is.na(obj@meta.data$svz_cells)] <- "ctx"
my_levels <- c("AC", "NSC", "TAP", "NB", "ctx")
obj@meta.data$svz_cells <- factor(x = obj@meta.data$svz_cells, levels = my_levels)

obj@meta.data$ctx_cells[is.na(obj@meta.data$ctx_cells)] <- "svz"
my_levels <- c("AC1", "AC2", "AC3", "NB", "svz")
obj@meta.data$ctx_cells <- factor(x = obj@meta.data$ctx_cells, levels = my_levels)

p <- list()
for( i in c("protocol", "res.0.5", "svz_cells", "ctx_cells")){
  p[[i]] <- TSNEPlot(object = obj, 
                     group.by = i,
                     colors.use = my.palette,
                     pt.size = 0.5)
}

plot_grid(plotlist = p)


obj <- RunUMAP(object = obj, 
               reduction.use = "cca.aligned", 
               dims.use = 1:13,
               do.fast = TRUE, 
               seed.use = 2)

#saveRDS(obj, "cca_svz_ctx.rds")

current.cluster.ids <- c(0, 1, 2, 3, 4, 5)
new.cluster.ids <- c("AC1", "AC2", "NB", "NSC", "TAP", "AC3")
obj@ident <- plyr::mapvalues(x = obj@ident, from = current.cluster.ids, to = new.cluster.ids)
obj <- StashIdent(obj, "clusters")

# DEA
markers <- FindAllMarkers(obj, 
                          only.pos = T)

adj.markers <- filter(markers, 
                      p_val_adj < 0.05)

top.50 <- adj.markers %>%
  group_by(cluster) %>% 
  top_n(50, avg_logFC)

cluster.averages <- AverageExpression(obj, 
                                      return.seurat = T)

averages_hm_genes <- cluster.averages@scale.data[top.50$gene,] %>% 
  as.data.frame()
averages_hm_genes$gene <- rownames(averages_hm_genes)

# Order genes by avg_fc across groups
cl.avg <- list()
for(i in levels(obj@ident)) {
  
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

### Compare neurogenic cells across tissues
obj@meta.data$celltype_protocol <- paste(obj@ident, 
                                         obj@meta.data$protocol, 
                                         sep = "_")

obj <- SetAllIdent(obj, "celltype_protocol")

DimPlot(obj, "umap", group.by = "celltype_protocol")

degs_obj <- list()
obj_subset <- list()
avg_expr <- list()
plots <- list()
clusters <- c("AC1", "AC2", "AC3", "NSC", "TAP", "NB")

mito.genes <- grep(pattern = "^mt-", x = rownames(obj@raw.data), value = TRUE)
ribo.genes <- c(grep(pattern = "^Rps", x = rownames(obj@raw.data), value = TRUE),
                grep(pattern = "^Rpl", x = rownames(obj@raw.data), value = TRUE))

for(i in clusters){
  # Degs
  ident1 = paste(i, "ctx", sep = "_")
  ident2 = paste(i, "svz", sep = "_")
  markers_obj <- FindMarkers(obj, logfc.threshold = 0, min.pct = 0, 
                             ident.1 = ident1, ident.2 = ident2) 
  degs_obj[[ident1]] <- markers_obj
}

palette <- as.data.frame(cbind(c("AC1", "AC2", "AC3", "NSC", "TAP", "NB"),
                               c(purple, blue, green, light_blue, red, yellow)))
colnames(palette) = c("cluster", "color")
rownames(palette) = palette$cluster

obj <- SetAllIdent(obj, "clusters")

for(i in clusters){
  # Subset data
  obj_subset[[i]] <- SubsetData(obj, ident.use = i)
  obj_subset[[i]] <- SetAllIdent(obj_subset[[i]], "protocol")
  
  # Average expression
  avg_expr[[i]] <- log1p(x = AverageExpression(object = obj_subset[[i]]))
  ident1 = paste(i, "ctx", sep = "_")
  avg_expr[[i]] <- avg_expr[[i]][!(rownames(avg_expr[[i]]) %in% 
                                     c(mito.genes, "tdTomatoWRPE", ribo.genes)),]
  avg_expr[[i]] <- avg_expr[[i]][rownames(degs_obj[[ident1]]),]
  
  # Plots
  plots[[i]] <- ggplot(avg_expr[[i]], aes(ctx, svz)) +
    geom_point(color = palette[i,2], alpha = 0.4) + 
    ggtitle(i) + 
    xlab("Ctx") + 
    ylab("SVZ") +
    geom_abline(intercept = 1.09, 
                linetype="dashed", 
                colour = "black") +
    geom_abline(intercept = -1.09, 
                linetype="dashed", 
                colour = "black") +
    annotate("text", y = 5, x = 1, 
             label = paste("R^2 =", 
                           round(summary(
               lm(avg_expr[[i]]$ctx ~ avg_expr[[i]]$svz, 
                  data=avg_expr[[i]]))$adj.r.squared, 
               digits = 2),sep = " ")) +
    theme
  
}

pdf("degs_svz_ctx.pdf", width = 6, height = 5)
plot_grid(plotlist = plots)
dev.off()

# Representation of samples across clusters
obj@meta.data$orig.ident <- ifelse(obj@meta.data$orig.ident %in% c("an002", "an003F", "an003L"),
                                   "svz", 
                                   obj@meta.data$orig.ident)

counts = table(obj@meta.data$orig.ident, obj@ident)
prop <- prop.table(counts, 2)
barplot(prop,
        col= palette$color,
        cex.names=1, horiz = T)

# Sankey plot
all_df <- cbind(as.character(obj@meta.data$ctx_cells),
                as.character(obj@meta.data$svz_cells),
                as.character(obj@ident),
                obj@meta.data$protocol)
all_df <- as.data.frame(all_df)
colnames(all_df) <- c("ctx_cells", "svz_cells", "ident", "protocol")

all_df$clusters <- ifelse((all_df$ctx_cells != "svz"), 
                          paste(all_df$ctx_cells, all_df$protocol, sep = "_"),
                          paste(all_df$svz_cells, all_df$protocol, sep = "_"))

all_df %>% select(c(ident, clusters)) -> all_df
all_df %>% table() -> tit2d

# Specify id.vars: the variables to keep but not split apart on
tit <- melt(tit2d, id.vars=c("clusters", "V3"))

tit$ident <- ordered(tit$ident, 
                     levels = c("NB", "TAP", "NSC", "AC3", "AC2", "AC1"))

tit$clusters <- ordered(tit$clusters, 
                        levels = c("NB_svz", "TAP_svz", "NSC_svz", 
                                                 "AC_svz","NB_ctx", "AC3_ctx", 
                                                 "AC2_ctx", "AC1_ctx"))

alluvial(tit[,2:1], 
         freq=tit$value,
         col = ifelse(tit$clusters %in% c("NB_svz", "TAP_svz", "NSC_svz", "AC_svz"), 
                      "#63C6C0", 
                      "#D1482E"),
         border = ifelse(tit$clusters %in% c("NB_svz", "TAP_svz", "NSC_svz", "AC_svz"), 
                         "#63C6C0", 
                         "#D1482E"),
         hide = tit$value < 10,
         cex = 0.7)


