## 0. Prepare R environment
## ******************************************

#Set working directory (specific to each user)
setwd("/Users/anttonalberdi/github/domesticat/")

#Load required libraries
library(ggplot2)

## 1. Load and modify data
## ************************************

MAGtree <- read.tree("data/MAG.tree")
MAGinfo <- read.csv("data/MAG_info.csv", row.names=1)
MAGannot <- read.delim("data/DRAM_product.tsv",sep='\t', row.names=1)

# Remove MAGs absent in MAG_info
MAGannot <- MAGannot[rownames(MAGinfo),]

# Convert True/False character to 1/0 binary, numeric data
for(i in 1:ncol(MAGannot)){
  if(is.character(MAGannot[,i])){
    MAGannot[,i][MAGannot[,i]=="True"]=1
    MAGannot[,i][MAGannot[,i]=="False"]=0
    MAGannot[,i]=as.numeric(MAGannot[,i])
  }
}

# Discard completely absent functions from further exploration
colnames(MAGannot[,colSums(MAGannot)==0])
MAGannot <- MAGannot[,colSums(MAGannot)>0]

## 2. Exploration of MAG completeness and redundancy
## *************************************************

phylumcolors <- c("#5b828e","#5e6668","#bbcfd7","#ba9a88","#ac7e62","#aca69f","#adab76","#666b3a","#0f211a","#012e67")
phylumcolorsalpha <- c("#5b828e70","#5e666870","#bbcfd770","#ba9a8870","#ac7e6270","#aca69f70","#adab7670","#666b3a70","#0f211a70","#012e6770")

#Completeness vs Contamination plot
pdf("figures/MAG_comp-cont.pdf",width=8,height=4)
ggplot(MAGinfo,aes(y=Contamination,x=Completeness,fill=Phylum, color=Phylum)) +
  geom_point() +
  scale_color_manual(values=phylumcolors) +
  theme_classic()
dev.off()

#Completeness & Contamination per Phylum
pdf("figures/MAG_comp-cont_phylum.pdf",width=8,height=4)
ggplot(MAGinfo) +
  geom_boxplot(aes(y=Completeness,x=Phylum,fill=Phylum, color=Phylum)) +
  geom_boxplot(aes(y=Contamination,x=Phylum,fill=Phylum, color=Phylum)) +
  ylim(0, 100) +
  scale_fill_manual(values=phylumcolorsalpha) +
  scale_color_manual(values=phylumcolors) +
  theme_classic() +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")
dev.off()

## 3. Exploration of ETH complexes
## *************************************************

library(ggtree)
library(ggtreeExtra)

phylumcolors <- c("#5b828e","#5e6668","#bbcfd7","#ba9a88","#ac7e62","#aca69f","#adab76","#666b3a","#0f211a","#012e67")

base <- ggtree(MAGtree, layout="rectangular") %<+% MAGinfo +
  geom_tippoint(aes(color=Phylum)) +
  scale_color_manual(values=phylumcolors)

gheatmap(base, MAGannot[,c(13:25)], offset = 0, width=.1,
          colnames_position="top",
          colnames_angle=90, colnames_offset_y = 1,
          hjust=0, font.size=2) +
          scale_fill_gradientn(colours = c("#ffffff","#000000")) +
          theme(legend.position="none")
