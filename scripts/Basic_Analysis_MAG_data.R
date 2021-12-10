library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ape)
library(phytools)
library(UpSetR)
library(edgeR)

## 1. Import data and make arrays conformable
## ******************************************

# Import metadata
cat_metadata=read.csv("data/DomestiCAT metadata.csv")
head(cat_metadata)
cat_metadata=cat_metadata[-c(96:98),]
cat_metadata$CombinedID=gsub("_",".",cat_metadata$CombinedID)
cat_metadata$CombinedID=gsub("-",".",cat_metadata$CombinedID)
cat_metadata$CombinedID=gsub(" ","",cat_metadata$CombinedID)
cat_metadata=cat_metadata[cat_metadata$Extraction.ID!="Blank",]
cat_metadata$Tame.Feral[grepl("Tame",cat_metadata$Tame.Feral)]="Tame"
table(cat_metadata$Tame.Feral)
table(cat_metadata$Location)
table(cat_metadata$Tame.Feral,cat_metadata$Location)

# Import MAG data
MAGcounts <- read.delim("data/DomestiCat_all_MAG_counts.txt", sep = '\t')
genomelength <- read.delim("data/DomestiCat_all_MAG_length.txt", sep = '\t')
MAGcounts[1:6,1:6]
genomelength[1:6,1:6]
# Proportion of GenomeID in two tables in same order
mean(MAGcounts$Genome==genomelength$Genome)
Genome_ID=MAGcounts$Genome
# Sequence counts per genome relative to total sequences in genome
MAGcounts_relLength=MAGcounts[,-1]/genomelength[,-1]
dim(MAGcounts_relLength)
MAGcounts_relLength[1:6,1:6]
# Sequence counts per genome relative to total sequences in genome and in sample
MAGcounts_relL_rel=t(t(MAGcounts_relLength)/colSums(MAGcounts_relLength))

colnames(MAGcounts_relL_rel)=gsub(".Read.Count","",colnames(MAGcounts_relL_rel))
colnames(MAGcounts_relL_rel)=gsub("_",".",colnames(MAGcounts_relL_rel))
rownames(MAGcounts_relL_rel)=Genome_ID

# Transpose the MAG table to have MAGs as columns and samples as rows
MAGcounts_relL_rel_t=data.frame(t(MAGcounts_relL_rel))


# Import phylo tree
phylo_tree=read.newick("data/gtdbtk.bac120.classify.tree")
phylo.tree=drop.tip(phylo_tree,phylo_tree$tip.label[-match(colnames(MAGcounts_relL_hel), phylo_tree$tip.label)])
is.ultrametric(phylo.tree)
phylo.tree=force.ultrametric(phylo.tree,method = "nnls")
is.ultrametric(phylo.tree)
plot(phylo.tree)

# Sort and trim data to make them conformable
colnames(MAGcounts_relL_rel)[which(!colnames(MAGcounts_relL_rel)%in%cat_metadata$CombinedID)]
cat_metadata$CombinedID[which(!cat_metadata$CombinedID%in%colnames(MAGcounts_relL_rel))]

cat_metadata_red=cat_metadata[which(cat_metadata$CombinedID%in%colnames(MAGcounts_relL_rel)),]
MAGcounts_relL_rel_red=MAGcounts_relL_rel[,which(colnames(MAGcounts_relL_rel)%in%cat_metadata_red$CombinedID)]
dim(cat_metadata_red)[1]==dim(MAGcounts_relL_rel_red)[2]
colnames(MAGcounts_relL_rel_red)%in%cat_metadata_red$CombinedID

cat_metadata_red=cat_metadata_red[match(colnames(MAGcounts_relL_rel_red),cat_metadata_red$CombinedID),]
colnames(MAGcounts_relL_rel_red)==cat_metadata_red$CombinedID

## 2. Analysis of MAG structure in cats
## ************************************

# Jaccard type overlap metric, Neutral, q=0
#MAG_dis_Neu_q0=pair_dis(MAGcounts_relL_rel_red,qvalue = 0)
#saveRDS(MAG_dis_Neu_q0,file="data/MAG_dis_Neu_q0.rds")
MAG_dis_Neu_q0=readRDS(file = "data/MAG_dis_Neu_q0.rds")
MAG_dis_Neu_q0$L1_UqN[upper.tri(MAG_dis_Neu_q0$L1_UqN)]=t(MAG_dis_Neu_q0$L1_UqN)[upper.tri(MAG_dis_Neu_q0$L1_UqN)]
diag(MAG_dis_Neu_q0$L1_UqN)=1

set.seed(1)
MAG_dist_Neu_q0_Location_MHV<-betadisper(as.dist(MAG_dis_Neu_q0$L1_UqN),cat_metadata_red$Location)
permutest(MAG_dist_Neu_q0_Location_MHV)
MAG_dist_Neu_q0_TameFeral_MHV<-betadisper(as.dist(MAG_dis_Neu_q0$L1_UqN),cat_metadata_red$Tame.Feral)
permutest(MAG_dist_Neu_q0_TameFeral_MHV)

adonis(as.dist(MAG_dis_Neu_q0$L1_UqN)~Location*Tame.Feral,data = cat_metadata_red)

# Jaccard tpe overlap metric, Neutral, q=1
#MAG_dis_Neu_q1=pair_dis(MAGcounts_relL_rel_red,qvalue = 1)
#saveRDS(MAG_dis_Neu_q1,file="data/MAG_dis_Neu_q1.rds")
MAG_dis_Neu_q1=readRDS(file="data/MAG_dis_Neu_q1.rds")
MAG_dis_Neu_q1$L1_UqN[upper.tri(MAG_dis_Neu_q1$L1_UqN)]=t(MAG_dis_Neu_q1$L1_UqN)[upper.tri(MAG_dis_Neu_q1$L1_UqN)]
diag(MAG_dis_Neu_q1$L1_UqN)=1

set.seed(1)
MAG_dist_Neu_q1_Location_MHV<-betadisper(as.dist(MAG_dis_Neu_q1$L1_UqN),cat_metadata_red$Location)
anova(MAG_dist_Neu_q1_Location_MHV)
MAG_dist_Neu_q1_TameFeral_MHV<-betadisper(as.dist(MAG_dis_Neu_q1$L1_UqN),cat_metadata_red$Tame.Feral)
anova(MAG_dist_Neu_q1_TameFeral_MHV)

adonis(as.dist(MAG_dis_Neu_q1$L1_UqN)~Location*Tame.Feral,data = cat_metadata_red)

# Jaccard tpe overlap metric, Phylo, q=0
#MAG_dis_Phy_q0=pair_dis(MAGcounts_relL_rel_red,qvalue = 0,tree = phylo.tree)
#saveRDS(MAG_dis_Phy_q0,file="data/MAG_dis_Phy_q0.rds")
MAG_dis_Phy_q0=readRDS(file="data/MAG_dis_Phy_q0.rds")
MAG_dis_Phy_q0$L1_UqN[upper.tri(MAG_dis_Phy_q0$L1_UqN)]=t(MAG_dis_Phy_q0$L1_UqN)[upper.tri(MAG_dis_Phy_q0$L1_UqN)]
diag(MAG_dis_Phy_q0$L1_UqN)=1

set.seed(1)
MAG_dist_Phy_q0_Location_MHV<-betadisper(as.dist(MAG_dis_Phy_q0$L1_UqN),cat_metadata_red$Location)
anova(MAG_dist_Phy_q0_Location_MHV)
MAG_dist_Phy_q0_TameFeral_MHV<-betadisper(as.dist(MAG_dis_Phy_q0$L1_UqN),cat_metadata_red$Tame.Feral)
anova(MAG_dist_Phy_q0_TameFeral_MHV)

adonis(as.dist(MAG_dis_Phy_q0$L1_UqN)~Location*Tame.Feral,data = cat_metadata_red)

# Jaccard tpe overlap metric, Phylo, q=1
#MAG_dis_Phy_q1=pair_dis(MAGcounts_relL_rel_red,qvalue = 1,tree=phylo.tree)
#saveRDS(MAG_dis_Phy_q1,file="data/MAG_dis_Phy_q1.rds")
MAG_dis_Phy_q1=readRDS(file="data/MAG_dis_Phy_q1.rds")
MAG_dis_Phy_q1$L1_UqN[upper.tri(MAG_dis_Phy_q1$L1_UqN)]=t(MAG_dis_Phy_q1$L1_UqN)[upper.tri(MAG_dis_Phy_q1$L1_UqN)]
diag(MAG_dis_Phy_q1$L1_UqN)=1

set.seed(1)
MAG_dist_Phy_q1_Location_MHV<-betadisper(as.dist(MAG_dis_Phy_q1$L1_UqN),cat_metadata_red$Location)
anova(MAG_dist_Phy_q1_Location_MHV)
MAG_dist_Phy_q1_TameFeral_MHV<-betadisper(as.dist(MAG_dis_Phy_q1$L1_UqN),cat_metadata_red$Tame.Feral)
anova(MAG_dist_Phy_q1_TameFeral_MHV)

adonis(as.dist(MAG_dis_Phy_q1$L1_UqN)~Location*Tame.Feral,data = cat_metadata_red)

## Using q-value = 1 instead of 0 increases the differences observed between locations
## as well as between Tame/Feral. Differences between Tame/Feral in microbiome composition
## appear to be more marked in some locations than in others, as indicated by significant
## interactions for Neutral/q=1 PERMANOVA (and near sig for Phylo/q=1). The results when inclduing
## vs. excluding phylogeny are very similar, but location effects are slightly stronger
## when including phylogeny and Tame/Feral is slightly stronger when excluding it.

## NMDS ordinations of Neutral/q=1 and Phylo/q=1.
## The ordinations of q=0 are heavily influenced by some extreme observations and are
## not meaningful.

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
myColors=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')
names(myColors)=levels(cat_metadata_red$Location)
colScale=scale_colour_manual(name = "Location",values = myColors)
shapeScale=scale_shape_manual(name="Tame/Feral",values = c(0,15))
windows(h=8,w=10)
ggplot(MAG_nmds_scores,
       aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=as.factor(cat_metadata_red$Tame.Feral),
                 color=as.factor(cat_metadata_red$Location)),size=4) + 
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,color=group), size=1, linetype=2) +
  colScale + 
  shapeScale + 
  theme_bw()

# Phylo and q = 1

set.seed(1)
MAG_nmds=metaMDS(MAG_dis_Phy_q1$L1_UqN,k=2)
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
myColors=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')
names(myColors)=levels(cat_metadata_red$Location)
colScale=scale_colour_manual(name = "Location",values = myColors)
shapeScale=scale_shape_manual(name="Tame/Feral",values = c(0,15))
windows(h=8,w=10)
ggplot(MAG_nmds_scores,
       aes(x=MDS1,y=MDS2))+
  geom_point(aes(shape=as.factor(cat_metadata_red$Tame.Feral),
                 color=as.factor(cat_metadata_red$Location)),size=4) + 
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,color=group), size=1, linetype=2) +
  colScale + 
  shapeScale + 
  theme_bw()

set.seed(1)
pairwise.adonis(MAG_dis_Neu_q1$L1_UqN,cat_metadata_red$Location,p.adjust.m = "holm",perm = 99999)
pairwise.adonis(MAG_dis_Phy_q1$L1_UqN,cat_metadata_red$Location,p.adjust.m = "holm",perm=99999)

## Based on Phylo/q=1, there are "two groups" of locations.
## Based on Neutral/q=1 the "two groups" get connected through Spain-Malaysia.

## 3. MAGs in different location and shared among locations
## ********************************************************

MAGcounts_pa=1*(MAGcounts_relL_rel_red>0)
MAGcounts_pa[1:6,1:6]
table_upset_analysis=t(aggregate(t(MAGcounts_pa),by=list(cat_metadata_red$Location),FUN=sum)[,-1])
colnames(table_upset_analysis)=levels(as.factor(cat_metadata_red$Location))
table_upset_analysis=(table_upset_analysis>0)*1
table_upset_analysis=data.frame(table_upset_analysis)
table_upset_analysis=apply(table_upset_analysis,2,as.integer)

upset(table_upset_analysis, sets = c("Aruba","Brazil","Cabo.Verde","Denmark",
                       "Malaysia","Spain"), mb.ratio = c(0.55, 0.45), 
      order.by = "freq")

## All MAGs are globally widespread with 227 present in the 6 countries and
## the other 2 in 5 countries each.


## 4. Enrichment analysis between Tame/Feral
## *****************************************

## EdgeR

#Prepare MAG counts for edgeR
rownames(MAGcounts)=Genome_ID
MAGcounts_edger=MAGcounts[,-1]
colnames(MAGcounts_edger)=gsub(".Read.Count","",colnames(MAGcounts_edger))
colnames(MAGcounts_edger)=gsub("_",".",colnames(MAGcounts_edger))

MAGcounts_edger=MAGcounts_edger[,which(colnames(MAGcounts_edger)%in%cat_metadata_red$CombinedID)]
dim(cat_metadata_red)[1]==dim(MAGcounts_edger)[2]
colnames(MAGcounts_edger)==cat_metadata_red$CombinedID

#Run EdgeR
design=model.matrix(~Location+Tame.Feral,data = cat_metadata_red)
dge=DGEList(counts=MAGcounts_edger)
dge=estimateDisp(dge, design,robust = TRUE)
fit=glmQLFit(dge, design)
qlf=glmQLFTest(fit)

results=topTags(qlf, n=nrow(lrt))
summary(decideTests(qlf))
#13 MAGs were were found to be of lower abundance in tame cats
#1 MAG was found to be of higher abundance in feral cats
results[results$table$FDR<0.05,]

## ANCOM

rownames(cat_metadata_red)=cat_metadata_red$CombinedID
MAGCounts_phyloseq_MAG_metadata=phyloseq(otu_table(MAGcounts_edger,taxa_are_rows = TRUE),sample_data(cat_metadata_red))
out = ancombc(phyloseq = MAGCounts_phyloseq_MAG_metadata, formula = "Location + Tame.Feral", 
              p_adj_method = "fdr",group = "Tame.Feral", tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res = out$res
res_global = out$res_global

#Coefficients
tab_coef = res$beta
col_name = c("Brazil", "Cabo_Verde", "Denmark", "Malaysia", "Spain", 
             "Tame")
colnames(tab_coef) = col_name
tab_coef$Tame

#SE-s
tab_se = res$se
colnames(tab_se) = col_name
tab_se$Tame

#Test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w $Tame

#p-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p$Tame

#Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q$Tame

#Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
sum(tab_diff$Tame)

## There are no differentially abundant taxa after FDR correction based on ANCOM

## Among the taxa found to be DA between Tame.Feral in edgeR, significant p-value in ANCOM
## (unadjusted p-val).
intersect(rownames(results[results$table$FDR<0.05,]),rownames(tab_p)[tab_p$Tame<0.05])

## Among the taxa found to be DA between Tame.Feral in edgeR, non significant p-value in ANCOM
## (unadjusted p-val).
setdiff(rownames(results[results$table$FDR<0.05,]),rownames(tab_p)[tab_p$Tame<0.05])
