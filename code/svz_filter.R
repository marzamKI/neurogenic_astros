source("./code/global.R", local = FALSE)

# Load gbm and create objects
sample_files <- list.files(path = "./data/raw_external/", full.names = F)

samples <- strsplit(sample_files, "_") %>% 
  sapply(extract2, 2)

# Make Seurat objects
gbm <- list()
obj <- list()

for(i in sample_files) {
  for(j in samples) {
    gbm[[j]] <- read.delim(file.path("./data/raw_external", i), header = T)
    
    rownames(gbm[[j]]) <- gbm[[j]]$GENE
    gbm[[j]]  <- gbm[[j]][,-1]
    gbm[[j]]  <- as.matrix(gbm[[j]] )
    
    obj[[j]] <- CreateSeuratObject(gbm[[j]], 
                                   min.cells = 3, 
                                   min.genes = 200, 
                                   project = j)
    
    obj[[j]] <- RenameCells(obj[[j]], 
                            add.cell.id = obj[[j]]@meta.data$orig.ident)
  }
  
}

svz <- obj[["an002"]]
for(i in samples[-1]) {
  svz <- MergeSeurat(object1 = svz, 
                     object2 = obj[[i]], 
                     add.cell.id2 = i)
}

mito.genes <- grep(pattern = "^mt-", 
                   x = rownames(x = svz@data), 
                   value = TRUE)
percent.mito <- Matrix::colSums(svz@raw.data[mito.genes,])/Matrix::colSums(svz@raw.data)

obj_prefilter <- AddMetaData(svz, 
                             metadata = percent.mito, 
                             col.name = "percent.mito")

VlnPlot(obj_prefilter, 
        c("nGene", "nUMI", "percent.mito"),
        nCol = 3, 
        group.by = "orig.ident",
        point.size.use = 0.02, 
        size.x.use = 0)

svz <- FilterCells(obj_prefilter, 
                   subset.names = "percent.mito", 
                   low.thresholds = -Inf, 
                   high.thresholds = 0.2)

svz <- ScaleData(object = svz, 
                 vars.to.regress = c("nUMI", "percent.mito", "ident"))

svz <- FindVariableGenes(svz, 
                         mean.function = ExpMean,
                         dispersion.function = LogVMR, 
                         x.low.cutoff = 0.05,
                         x.high.cutoff = 3.5, 
                         y.cutoff = 0.5, 
                         do.plot = F)

length(x = svz@var.genes)

# Run PCA and tSNE for clustering cell types
svz <- RunPCA(svz, pc.genes = svz@var.genes)

# look at the distribution of samples on PCA plot
PCAPlot(svz, group.by = "ident")

# determine number of PCs to use for downstream analysis
PCElbowPlot(svz)

svz <- RunTSNE(svz, 
               dims.use = 1:13, 
               do.fast = TRUE, 
               check_duplicates =F)

markers <- c("Aqp4", "Thbs4", "Mki67", "Dcx", 
             "Pdgfra", "Mog", "Flt1", "Vtn", 
             "P2ry12", "Mrc1", "Penk")

FeaturePlot(svz, 
            markers, 
            cols = my.inferno, 
            nCol = 3)

svz <- FindClusters(svz, 
                    reduction.type = "pca", 
                    dims.use = 1:13, 
                    resolution = 0.3)

TSNEPlot(svz, group.by = "ident")

#saveRDS(svz, "./rds/full_svz.rds")
svz <- SubsetData(svz, ident.use = c("1", "7", "6", "2"))

other_cells <- c("Pdgfra", "Mog", "Klk6", 
                 "Foxj1", "Flt1", "Vtn", 
                 "Acta2", "P2ry12", "Tac1", "Penk")

svz <- FilterCells(svz, other_cells, 
                   high.thresholds = rep(0.4, length(other_cells)))

svz <- RunTSNE(svz, dims.use = 1:13, 
               do.fast = TRUE, 
               check_duplicates =F)

TSNEPlot(svz, group.by = "ident")

markers <- c("Slc1a3", "Thbs4", "Ascl1", "Mki67", "Dcx")
FeaturePlot(svz, 
            markers, 
            cols.use = my.inferno, 
            nCol = 3)

current.cluster.ids <- c("1", "2", "6", "7")
new.cluster.ids <- c("AC", "NB", "TAP", "NSC")
svz@ident <- plyr::mapvalues(svz@ident, from = current.cluster.ids, to = new.cluster.ids)
svz <- StashIdent(svz, "svz_cells")

#saveRDS(obj, "./rds/neurogenic_svz.rds")
