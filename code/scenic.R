# Follow standard workflow by Aerts lab (https://github.com/aertslab/SCENIC)

## Run on Server
suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SingleCellExperiment)
  library(Seurat)
})

options(stringsAsFactors=FALSE) # Do not convert strings to factors

astros_npcs <- readRDS("astros_clean_filtered_tsne_dm.rds")
exprMat <- as.matrix(astros_npcs@data)

cellInfo <- as.data.frame(astros_npcs@meta.data)
saveRDS(cellInfo, file="int/cellInfo.Rds")

org="mgi" # or hgnc, or dmel
dbDir="cisTarget_databases" # RcisTarget databases location
myDatasetTitle="NeurogAstros" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, 
                                  dbDir=dbDir, 
                                  dbs=dbs, 
                                  datasetTitle=myDatasetTitle, 
                                  nCores=20) 

scenicOptions@inputDatasetInfo$datasetTitle <- "NeurogAstros"
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

genesKept <- geneFiltering(exprMat, 
                           scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)

runGenie3(exprMat_filtered, scenicOptions)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 2
scenicOptions@settings$seed <- 123



## Run on local
source("./code/global.R", local = FALSE)

scenicOptions <- readRDS("int/scenicOptions.Rds")

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

runSCENIC_4_aucell_binarize(scenicOptions)

regulons <- loadInt(scenicOptions, "regulons")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

## Add SCENIC output as metadata
astro_scenic <- AddMetaData(astros_npcs, as.data.frame(cbind(cellInfo,t(getAUC(regulonAUC)))))

## Identify AC4
dm <- cbind(astros_clean@dr$dm@cell.embeddings, 
            as.character(astros_clean@ident))

dm <- as.data.frame(dm)

colnames(dm) <- c("DM1", "DM2", "DM3", "clusters")

dm$DM1 <- as.numeric(as.character(dm$DM1))
dm$DM3 <- as.numeric(as.character(dm$DM3))
dm$clusters <- as.character(dm$clusters)

dm %>% filter(clusters == "AC1") %>% select(DM1) %>% max() -> max_AC1

astros_clean@meta.data$new_clusters <- ifelse((dm$DM1 >= max_AC1 & 
                                                 dm$clusters == "AC2"), 
                                              "AC4", 
                                              dm$clusters)

# Adjust threshold for binary regulon activity
logMat <- astros_clean@scale.data

aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds

scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")   

# Heatmap plotting per cluster percentage of cells with active regulons
minPerc <- .4
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")

cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo) %in% 
                                            colnames(binaryRegulonActivity)),, drop=FALSE]

regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), 
                                                     cellInfo_binarizedCells$ident_neurog), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))

binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5, 
                   cluster_cols=T,
                   color = viridis(n = 100, begin = 0.3, end = 1, option = "inferno"), 
                   breaks=seq(0, 1, length.out = 100),
                   treeheight_row=20, 
                   treeheight_col=30, 
                   border_color=NA)

FeaturePlot(astro_scenic, 
            reduction.use ="dm", 
            dim.1 = 1, 
            dim.2 = 3, 
            features.plot = c("Olig2_extended (193g)","Ascl1 (72g)", 
                              "Sox11 (185g)","Sp8_extended (204g)"), 
            cols.use = my.inferno, 
            no.axes = T, 
            pt.size = 2)
