## 0. Prepare R environment
## ******************************************

#Load required libraries
library(vegan)
library(data.table)
library(hilldiv)
library(ggplot2)
library(RColorBrewer)
library(ape)
library(phytools)
library(UpSetR)
library(edgeR)
library(phyloseq)
library(ANCOMBC)

#Set working directory (specific to each user)
setwd("/Users/anttonalberdi/github/domesticat/")

## 1. Import data, transform and make arrays conformable
## ******************************************

## 1.1- Import and prepare metadata
cat_metadata=read.csv("data/DomestiCAT_metadata.csv")
head(cat_metadata)
cat_metadata=cat_metadata[-c(96:98),]
cat_metadata=cat_metadata[cat_metadata$ExtractionID!="Blank",]
cat_metadata$CombinedID <- gsub("-","_",cat_metadata$CombinedID , fixed = TRUE)

# Get general statistics
table(cat_metadata$Origin)
table(cat_metadata$Location)
table(cat_metadata$Origin,cat_metadata$Location)

## 1.2- Import MAG data
MAGcounts <- read.csv("data/MAG_counts.csv", row.names=1, header=TRUE)
MAGinfo <- read.csv("data/MAG_info.csv", row.names=1, header=TRUE)

#Check whether rows and columns are identical
identical(sort(rownames(MAGcounts)),sort(rownames(MAGinfo)))

#Sort all matrices
MAGcounts <- MAGcounts[sort(rownames(MAGcounts)),]
MAGinfo <- MAGinfo[sort(rownames(MAGinfo)),]

#Visualise
MAGcounts[1:6,1:6]
MAGinfo[1:6,]

# Filter MAGs below coverage threshold (threshold: 0.3 coverage)
MAGcov <- sweep((MAGcounts * 300), 1, MAGinfo$Length, FUN = '/')
MAGcov[MAGcov <= 0.3] <- 0
MAGcov[MAGcov > 0.3] <- 1
MAGcountsFilt <- MAGcounts * MAGcov

# Relative abundance normalisation

# EBI style
# MAGra =  [Reads mapped to MAG] / ([MAG length]/1000000 * [Total reads in sample] / 1000000)
# MAGra =  MAGcountsFilt / (Length/1000000 * colSums(MAGcountsFilt) / 1000000) # EBI style
#MAGra = sweep(sweep(MAGcountsFilt, 1, (MAGinfo$Length / 1000000), FUN = '/'),2,(colSums(MAGcountsFilt) / 1000000),FUN = '/')

# Antton (normalisation accounting for relative MAG length)
# MAGra =  tss([Reads mapped to MAG] / ([MAG length]/[All MAGs length])
tss <- function(x){sweep(x, 2, colSums(x), FUN="/")}
MAGrel = tss(sweep(MAGcountsFilt, 1, (MAGinfo$Length / sum(MAGinfo$Length)), FUN = '/'))
MAGrel[1:6,1:6]

#Remove incorrect samples
MAGrel = MAGrel[,!is.na(colSums(MAGrel))]

# Transpose the MAG table to have MAGs as columns and samples as rows
MAGrel_t=data.frame(t(MAGrel))

# Import phylo tree
#phylo_tree=read.newick("data/gtdbtk.bac120.classify.tree")
#MAGtree=keep.tip(phylo_tree,colnames(MAGcounts_relL_rel_t))
#is.ultrametric(MAGtree)
#MAGtree=force.ultrametric(MAGtree,method = "nnls")
#is.ultrametric(MAGtree)
#write.tree(MAGtree,"data/MAG.tree")
#plot(MAGtree)
MAGtree=read.tree("data/MAG.tree")

# Detect samples in the count table not present in metadata
colnames(MAGrel)[which(!colnames(MAGrel)%in%cat_metadata$CombinedID)]
# Detect samples in metadata not present in the count table
cat_metadata$CombinedID[which(!cat_metadata$CombinedID%in%colnames(MAGrel))]

#Intersect MAG counts and metadata and validate
cat_metadata_red=cat_metadata[which(cat_metadata$CombinedID%in%colnames(MAGrel)),]
MAGrel_red=MAGrel[,which(colnames(MAGrel)%in%cat_metadata_red$CombinedID)]
dim(cat_metadata_red)[1]==dim(MAGrel_red)[2]
colnames(MAGrel_red)%in%cat_metadata_red$CombinedID
cat_metadata_red=cat_metadata_red[match(colnames(MAGrel_red),cat_metadata_red$CombinedID),]
colnames(MAGrel_red)==cat_metadata_red$CombinedID

## 2. MAG tree with quality annotations
## ************************************

library(ggtree)
library(ggtreeExtra)
library(phyloseq)
library(dplyr)
library(ape)

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

## 3. Barplots clustered by sites
## ************************************

MAGrel_melt <- melt(as.matrix(MAGrel))
MAGrel_melt <- merge(MAGrel_melt,cat_metadata,by.x="Var2",by.y="CombinedID")
MAGrel_melt <- merge(MAGrel_melt,MAGinfo,by.x="Var1",by.y="row.names")
colnames(MAGrel_melt)[1:3] <- c("MAG","Sample","Value")

phylumcolors <- c("#5b828e","#5e6668","#bbcfd7","#ba9a88","#ac7e62","#aca69f","#adab76","#666b3a","#0f211a","#012e67")

pdf("figures/MAGbarplot.pdf",width=8,height=4)
ggplot(MAGrel_melt,aes(y=Value,x=Sample,fill=Phylum)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=phylumcolors) +
  facet_grid(~ Location,
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  theme_minimal() +
  theme(axis.text.x=element_blank(), panel.grid.major.x = element_blank())
dev.off()



## 3. Analysis of MAG structure in cats
## ************************************

# Jaccard type overlap metric, Neutral, q=0
MAG_dis_Neu_q0=pair_dis(MAGrel_red, qvalue=0)
saveRDS(MAG_dis_Neu_q0,file="results/MAG_dis_Neu_q0.rds")
MAG_dis_Neu_q0=readRDS(file = "results/MAG_dis_Neu_q0.rds")
MAG_dis_Neu_q0$L1_UqN[upper.tri(MAG_dis_Neu_q0$L1_UqN)]=t(MAG_dis_Neu_q0$L1_UqN)[upper.tri(MAG_dis_Neu_q0$L1_UqN)]
diag(MAG_dis_Neu_q0$L1_UqN)=1

set.seed(1)
MAG_dist_Neu_q0_Location_MHV<-betadisper(as.dist(MAG_dis_Neu_q0$L1_UqN),cat_metadata_red$Location)
permutest(MAG_dist_Neu_q0_Location_MHV)
MAG_dist_Neu_q0_Origin_MHV<-betadisper(as.dist(MAG_dis_Neu_q0$L1_UqN),cat_metadata_red$Origin)
permutest(MAG_dist_Neu_q0_Origin_MHV)

adonis(as.dist(MAG_dis_Neu_q0$L1_UqN)~Location*Origin,data = cat_metadata_red)

# Jaccard tpe overlap metric, Neutral, q=1
MAG_dis_Neu_q1=pair_dis(MAGrel_red,qvalue = 1)
saveRDS(MAG_dis_Neu_q1,file="results/MAG_dis_Neu_q1.rds")
MAG_dis_Neu_q1=readRDS(file="results/MAG_dis_Neu_q1.rds")
MAG_dis_Neu_q1$L1_UqN[upper.tri(MAG_dis_Neu_q1$L1_UqN)]=t(MAG_dis_Neu_q1$L1_UqN)[upper.tri(MAG_dis_Neu_q1$L1_UqN)]
diag(MAG_dis_Neu_q1$L1_UqN)=1

set.seed(1)
MAG_dist_Neu_q1_Location_MHV<-betadisper(as.dist(MAG_dis_Neu_q1$L1_UqN),cat_metadata_red$Location)
anova(MAG_dist_Neu_q1_Location_MHV)
MAG_dist_Neu_q1_Origin_MHV<-betadisper(as.dist(MAG_dis_Neu_q1$L1_UqN),cat_metadata_red$Origin)
anova(MAG_dist_Neu_q1_Origin_MHV)

adonis(as.dist(MAG_dis_Neu_q1$L1_UqN)~Location*Origin,data = cat_metadata_red)


# Jaccard tpe overlap metric, Phylo, q=0
MAG_dis_Phy_q0=pair_dis(MAGrel_red,qvalue = 0,tree = MAGtree)
saveRDS(MAG_dis_Phy_q0,file="results/MAG_dis_Phy_q0.rds")
MAG_dis_Phy_q0=readRDS(file="results/MAG_dis_Phy_q0.rds")
MAG_dis_Phy_q0$L1_UqN[upper.tri(MAG_dis_Phy_q0$L1_UqN)]=t(MAG_dis_Phy_q0$L1_UqN)[upper.tri(MAG_dis_Phy_q0$L1_UqN)]
diag(MAG_dis_Phy_q0$L1_UqN)=1

set.seed(1)
MAG_dist_Phy_q0_Location_MHV<-betadisper(as.dist(MAG_dis_Phy_q0$L1_UqN),cat_metadata_red$Location)
anova(MAG_dist_Phy_q0_Location_MHV)
MAG_dist_Phy_q0_Origin_MHV<-betadisper(as.dist(MAG_dis_Phy_q0$L1_UqN),cat_metadata_red$Origin)
anova(MAG_dist_Phy_q0_Origin_MHV)

adonis(as.dist(MAG_dis_Phy_q0$L1_UqN)~Location*Origin,data = cat_metadata_red)

# Jaccard tpe overlap metric, Phylo, q=1
MAG_dis_Phy_q1=pair_dis(MAGrel_red,qvalue = 1,tree=MAGtree)
saveRDS(MAG_dis_Phy_q1,file="results/MAG_dis_Phy_q1.rds")
MAG_dis_Phy_q1=readRDS(file="results/MAG_dis_Phy_q1.rds")
MAG_dis_Phy_q1$L1_UqN[upper.tri(MAG_dis_Phy_q1$L1_UqN)]=t(MAG_dis_Phy_q1$L1_UqN)[upper.tri(MAG_dis_Phy_q1$L1_UqN)]
diag(MAG_dis_Phy_q1$L1_UqN)=1

set.seed(1)
MAG_dist_Phy_q1_Location_MHV<-betadisper(as.dist(MAG_dis_Phy_q1$L1_UqN),cat_metadata_red$Location)
anova(MAG_dist_Phy_q1_Location_MHV)
MAG_dist_Phy_q1_Origin_MHV<-betadisper(as.dist(MAG_dis_Phy_q1$L1_UqN),cat_metadata_red$Origin)
anova(MAG_dist_Phy_q1_Origin_MHV)

adonis(as.dist(MAG_dis_Phy_q1$L1_UqN)~Location*Origin,data = cat_metadata_red)

# Neutral and q = 1
set.seed(1)
MAG_nmds=metaMDS(MAG_dis_Neu_q1$L1_UqN,k=2)
MAG_nmds$stress # stress = 0.15
MAG_nmds_scores=data.frame(MAG_nmds$points,group=factor(cat_metadata_red$Location))
MAG_nmds_mean=aggregate(MAG_nmds_scores[,1:2],list(MAG_nmds_scores$group),mean)
## Function to make ellipses in plot
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(MAG_nmds_scores$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(MAG_nmds_scores[MAG_nmds_scores$group==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}
myColors=c('#c4d7d1','#408892','#2d3749','#c04062','#6b3a59','#e08683')
names(myColors)=levels(cat_metadata_red$Location)
colScale=scale_colour_manual(name = "Location",values = myColors)
shapeScale=scale_shape_manual(name="Origin",values = c(1,16))
windows(h=8,w=10)
ggplot(MAG_nmds_scores,
       aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=as.factor(cat_metadata_red$Origin),
                 color=as.factor(cat_metadata_red$Location)),size=3) +
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,color=group), size=0.5, linetype=2) +
  colScale +
  shapeScale +
  theme_bw()

## 3. MAGs in different location and shared among locations
## ********************************************************

MAGcounts_pa=1*(MAGcounts_relL_rel_red>0)
MAGcounts_pa[1:6,1:6]
table_upset_analysis_cont=t(aggregate(t(MAGcounts_pa),by=list(cat_metadata_red$Location),FUN=sum)[,-1])
colnames(table_upset_analysis_cont)=levels(as.factor(cat_metadata_red$Location))
table_upset_analysis=(table_upset_analysis_cont>0)*1
table_upset_analysis=data.frame(table_upset_analysis)
table_upset_analysis=apply(table_upset_analysis,2,as.integer)
rownames(table_upset_analysis) <- rownames(MAGcounts_pa)

upset(as.data.frame(table_upset_analysis), sets = c("Aruba","Brazil","CaboVerde","Denmark","Malaysia","Spain"), mb.ratio = c(0.55, 0.45), order.by = "freq")


## 4. Enrichment analysis between Domestic/Feral
## *****************************************

## MAG info
MAGinfo <- read.csv("data/MAG_info.csv", row.names=1, header=TRUE)

## EdgeR
#Reference: Domestic (Tame)
#Alternative: Feral
MAGcounts_edger=MAGcounts
MAGcounts_edger=MAGcounts_edger[,which(colnames(MAGcounts_edger)%in%cat_metadata_red$CombinedID)]
dim(cat_metadata_red)[1]==dim(MAGcounts_edger)[2]
colnames(MAGcounts_edger)==cat_metadata_red$CombinedID

#Run EdgeR
design=model.matrix(~Location+Origin,data = cat_metadata_red)
dge=DGEList(counts=MAGcounts_edger)
dge=estimateDisp(dge,design,robust = TRUE)
fit=glmQLFit(dge, design)
qlf=glmQLFTest(fit)

#Plot
plotMD(qlf)
abline(h=c(-1, 1), col="grey")

results=topTags(qlf, n=20)
summary(decideTests(qlf))
#13 MAGs were were found to be of higher abundance in feral than domestic cats
#1 MAG was found to be of lower abundance in feral than domestic cats
difMAGs <- results[results$table$FDR<0.05,]
difMAGs <- cbind(as.data.frame(results[results$table$FDR<0.05,]),table_upset_analysis_cont[rownames(difMAGs),], taxonomy[rownames(difMAGs),], quality[rownames(difMAGs),])
write.csv(difMAGs,"results/enrichment.csv")
