library(ggtree)
library(ggtreeExtra)
library(ape)

#https://yulab-smu.top/treedata-book/chapter10.html

setwd("/Users/anttonalberdi/github/domesticat/")

tree <- read.tree("data/MAG.tree")
info <- read.csv("data/MAG_info.csv",row.names=1)
