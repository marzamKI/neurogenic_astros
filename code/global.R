# Load pkgs and color palettes
library(Seurat)
library(dplyr)
library(viridis)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(monocle)
library("scatterplot3d")

my.inferno = c(alpha("grey85", 
                     alpha = 0.4),
               viridis(10, 
                       option = "inferno", 
                       alpha = 1, 
                       begin = 0))

my.palette = c("#A52C60FF", "#2274A5", "#D1462F", "#5FAD56","#F7B32B", "#39A0ED")
