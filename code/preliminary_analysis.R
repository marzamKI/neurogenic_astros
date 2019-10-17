source("./code/global.R", local = FALSE)

# Load samples
samples <- list.dirs(path = "./data/raw_internal", full.names = F)[-1]

# Make Seurat objects
gbm <- list()
obj <- list()

for(i in samples) {
  gbm[[i]] <- Read10X(file.path(i))
  
  obj[[i]] <- CreateSeuratObject(gbm[[i]], 
                                 min.cells = 3, 
                                 min.genes = 200, 
                                 project = i)
  
  obj[[i]] <- RenameCells(obj[[i]], 
                          add.cell.id = obj[[i]]@meta.data$orig.ident)
}

all.cells <- obj[["healthy_wt"]]
for(i in samples[-1]) {
  all.cells <- MergeSeurat(object1 = all.cells, 
                           object2 = obj[[i]], 
                           add.cell.id2 = i)
}

#saveRDS(all.cells, "./rds/raw.all.cells.rds")

# Filter genes and cells
table(all.cells@meta.data$orig.ident)

# Mitochondrial content
mito.genes <- grep(pattern = "^mt-", 
                   x = rownames(x = all.cells@data), 
                   value = TRUE)

percent.mito <- Matrix::colSums(all.cells@raw.data[mito.genes,])/Matrix::colSums(all.cells@raw.data)
all.cells <- AddMetaData(object = all.cells, 
                         metadata = percent.mito, 
                         col.name = "percent.mito")

VlnPlot(all.cells, 
        c("nGene", "nUMI", "percent.mito"), 
        nCol = 3, 
        group.by = "orig.ident", 
        point.size.use = 0.02)

print(paste("Filter", 
            length(which(all.cells@meta.data$percent.mito > .2)), 
            "cells with mitochondrial content > 20%"))

all.cells <- FilterCells(object = all.cells, 
                         subset.names = "percent.mito", 
                         low.thresholds = -Inf, 
                         high.thresholds = 0.2)

# Cell cycle score
cc.genes <- readLines(con = "./data/meta/regev_lab_cell_cycle_genes.txt")

cc.genes <- CaseMatch(search = cc.genes, 
                     match = rownames(x = all.cells@raw.data))

all.cells <- CellCycleScoring(all.cells, 
                              s.genes = cc.genes[1:43], 
                              g2m.genes = cc.genes[44:97], 
                              set.ident = F)

GenePlot(all.cells,"nGene","S.Score",cex.use = 1)
GenePlot(all.cells,"nGene","G2M.Score",cex.use = 1)

# Filter cells with low cc.score (<0.2) nGenes >3sd from the mean.
outlier_ngenes <- 3*sd(all.cells@meta.data$nGene) + mean(all.cells@meta.data$nGene)
print(paste("Filter", 
            length(which(all.cells@meta.data$filter == "FILTER")), 
            "cells with nGenes >",
            round(outlier_ngenes, 1), 
            "& cell cycle score < 0.2"))   

VlnPlot(all.cells, 
        c("nGene", "S.Score", "G2M.Score"), 
        group.by = "filter")
all.cells <- SetAllIdent(all.cells, "filter")
all.cells <- SubsetData(all.cells, 
                        ident.remove = "FILTER")
# Scale
all.cells <- SetAllIdent(all.cells, "orig.ident")
all.cells <- ScaleData(all.cells, 
                       vars.to.regress = c("nUMI", "percent.mito"))

all.cells <- FindVariableGenes(all.cells, 
                               mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, 
                               x.high.cutoff = 5, 
                               y.cutoff = 0.5, 
                               do.plot = F)

length(x = all.cells@var.genes)

# Run PCA
all.cells <- RunPCA(all.cells, 
                    pc.genes = all.cells@var.genes)

PCElbowPlot(all.cells)

VizPCA(all.cells, 
       pcs.use = 1:10, 
       nCol = 5)

all.cells <- RunTSNE(all.cells, 
                     dims.use = 1:10,
                     do.fast = TRUE, 
                     check_duplicates =F)

all.cells <- FindClusters(all.cells, 
                          reduction.type = "pca", 
                          dims.use = 1:10, 
                          resolution = 0.1)

TSNEPlot(all.cells, group.by = "res.0.1", colors.use = my.palette)

cell.markers <- c("Aqp4", "Dcx", "Dcn", "Pllp", "Csf1r")
FeaturePlot(all.cells, 
            features.plot = cell.markers, 
            reduction.use = "tsne", 
            cols.use = my.inferno)

old.ident <- c("0", "1", "2", "3", "4", "5")
new.ident <- c("AC1", "VLMC", "MG", "NB", "AC2", "OL")
all.cells@ident <- plyr::mapvalues(x = all.cells@ident, 
                                   from = old.ident, to = new.ident)

VlnPlot(all.cells, 
        features.plot = c("nGene", "Xist", "Gjb6", "tdTomatoWRPE"), 
        point.size.use = 0.2,
        cols.use = my.palette)

#saveRDS(all.cells, "./rds/all.cells.rds")

# Filter cells lacking Xist and Tomato
all.cells_filtered <- FilterCells(object = all.cells, 
                                  subset.names = c("Xist", "tdTomatoWRPE"), 
                                  low.thresholds = c(1, 1))

#saveRDS(all.cells, "./rds/filtered.all.cells.rds")

# Remove contaminating cells, select for astrocytes and neurogenic progeny
astro_npcs <- SubsetData(all.cells, 
                         ident.use = c("AC1", "AC2", "NB"))

# Filter cells expressing Car4, Eng (contaminating endothelial cells) and Cspg4 (contaminating OPCs)
astro_npcs <- FilterCells(astro_npcs, 
                          subset.names = c("Car4", "Eng", "Cspg4", "Plp1"), 
                          high.thresholds = c(0.2, 0.2, 0.2, 4))

all.cells@meta.data$selected <- ifelse((all.cells_filtered@cell.names %in% astro_npcs@cell.names),
                                       "SELECTED", 
                                       "FILTERED") 
TSNEPlot(all.cells, group.by = "filter", colors.use = c("grey85", "#A52C60FF"), no.axes = T)

#saveRDS(astro_npcs, "./rds/raw_astros_npcs.rds")
