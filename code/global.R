# Load pkgs and color palettes
library(Seurat)
library(dplyr)
library(viridis)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(monocle)
library(scatterplot3d)

my.inferno = c(alpha("grey85", 
                     alpha = 0.4),
               viridis(10, 
                       option = "inferno", 
                       alpha = 1, 
                       begin = 0))

my.palette = c("#A52C60FF", "#2274A5", "#D1462F", "#5FAD56","#F7B32B", "#39A0ED")

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
