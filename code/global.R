# Load pkgs and color palettes
packages <- c("Seurat", "tidyverse", "viridis", 
              "gridExtra", "grid", "ggpubr", "ggplot2",
              "scatterplot3d", "alluvial", "reshape2")

if (!requireNamespace(packages, quietly = TRUE))
    install.packages(packages)

if (!requireNamespace("BiocManager", "monocle", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")

library(Seurat)
library(dplyr)
library(viridis)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(monocle)
library(scatterplot3d)
library(alluvial)
library(reshape2)


my.inferno = c(alpha("grey85", 
                     alpha = 0.4),
               viridis(10, 
                       option = "inferno", 
                       alpha = 1, 
                       begin = 0))

my.palette = c("#A52C60FF", "#2274A5", "#D1462F", "#5FAD56","#F7B32B", "#61C6C0")
purple = "#A52C60FF"; blue = "#2274A5"; red = "#D1462F"; green = "#5FAD56"; yellow = "#F7B32B"; light_blue = "#61C6C0"

# ggplot theme
theme = theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(), 
              legend.position = "none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())

# Load cc genes
cc.genes <- readLines(con = "./data/meta/regev_lab_cell_cycle_genes.txt")
cc.genes <- CaseMatch(search = cc.genes, 
                      match = rownames(x = all.cells@raw.data))



