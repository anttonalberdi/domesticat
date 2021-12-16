library(ggtree)
library(ggtreeExtra)
library(phyloseq)
library(dplyr)
library(ape)

#https://yulab-smu.top/treedata-book/chapter10.html

setwd("/Users/anttonalberdi/github/domesticat/")

MAGtree <- read.tree("data/MAG.tree")
MAGinfo <- read.csv("data/MAG_info.csv")
rownames(MAGinfo) <- MAGinfo[,1]

phylumcolors <- c("#5b828e","#5e6668","#bbcfd7","#ba9a88","#ac7e62","#aca69f","#adab76","#666b3a","#0f211a","#012e67")

base <- ggtree(tree, layout="circular", width=0.5) %<+% info +
  geom_tippoint(aes(color=Phylum)) +
  scale_color_manual(values=phylumcolors)

qualitycolors <- c("#79ad9f","#193439","#d3b78f","#f4f4f4","#f4f4f4","#f4f4f4","#b37d5f","#503143")

pdf("figures/MAGtree.pdf",width=5,height=5)
gheatmap(base, MAGinfo[,c("Completeness","Contamination")], offset = 0, width=.1,
        colnames_position="top",
        colnames_angle=90, colnames_offset_y = 1,
        hjust=0, font.size=2) +
        scale_fill_gradientn(colours = qualitycolors) +
        theme(legend.position="none")
dev.off()
