library(data.table)
library(vegan)
library(dplyr)
library(tidyr)
library(stringr)
library(ggcorrplot)
library(factoextra)
library(ape)
library(logisticPCA)
library(vegan)
library(psych)

taxonomy <- read.delim("data/gtdbtk.bac120.summary.tsv", sep = '\t')
taxonomyclean <- taxonomy %>% 
  mutate(classification = str_replace_all(classification, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
taxonomyclean1 <- taxonomyclean[,-1]
# Set MAG name as row names of the dataframe
rownames(taxonomyclean1) <- taxonomyclean[,1]

# Get the necessary columns, convert to matrix
taxmatrix <- as.matrix(taxonomyclean1[1:7])

# Pool Firmicutes A,B C into Firmicutes
taxmatrix_mod=data.frame(taxmatrix)
# taxmatrix_mod$Phylum[grepl("Firmi",taxmatrix_mod$Phylum)]="Firmicutes"

# Load phylogenetic tree
phylo_tree=read.newick("data/gtdbtk.bac120.classify.tree")
phylo.tree=drop.tip(phylo_tree,phylo_tree$tip.label[-match(colnames(MAGcounts_relL_hel), phylo_tree$tip.label)])
is.ultrametric(phylo.tree)
phylo.tree=force.ultrametric(phylo.tree,method = "nnls")
is.ultrametric(phylo.tree)
plot(phylo.tree)
# Load DRAM data
DRAM_data=read.delim("data/DRAM_product.tsv",sep='\t')
dim(DRAM_data)

# Discard the "fasta" genome
DRAM_data=DRAM_data[DRAM_data$genome!="fasta",]
rownames(DRAM_data)=DRAM_data$genome
DRAM=DRAM_data[,-1]
# Convert True/False character to 1/0 binary, numeric data
for(i in 1:ncol(DRAM)){
  if(is.character(DRAM[,i])){
    DRAM[,i][DRAM[,i]=="True"]=1
    DRAM[,i][DRAM[,i]=="False"]=0
    DRAM[,i]=as.numeric(DRAM[,i])
  }
}
str(DRAM)

# List of completely absent functions
colnames(DRAM[,colSums(DRAM)==0])
# [1] "X3.Hydroxypropionate.bi.cycle"                                                  
# [2] "Complex.I..NAD.P.H.quinone.oxidoreductase..chloroplasts.and.cyanobacteria"      
# [3] "Complex.I..NADH.dehydrogenase..ubiquinone..1.alpha.subcomplex"                  
# [4] "Complex.II..Succinate.dehydrogenase..ubiquinone."                               
# [5] "Complex.IV.Low.affinity..Cytochrome.aa3.600.menaquinol.oxidase"                 
# [6] "Complex.V..F.type.ATPase..eukaryotes"                                           
# [7] "Complex.V..V.type.ATPase..eukaryotes"                                           
# [8] "Methanogenesis.and.methanotrophy..Key.functional.gene"                          
# [9] "Methanogenesis.and.methanotrophy..dimethylamine....monomethylamine"             
# [10] "Methanogenesis.and.methanotrophy..methane....methanol..with.oxygen..mmo."       
# [11] "Methanogenesis.and.methanotrophy..methane....methanol..with.oxygen..pmo."       
# [12] "Methanogenesis.and.methanotrophy..methanol....methane"                          
# [13] "Methanogenesis.and.methanotrophy..monomethylamine....ammonia"                   
# [14] "Methanogenesis.and.methanotrophy..putative.but.not.defining.CO2....methane"     
# [15] "Methanogenesis.and.methanotrophy..trimethylamine....dimethylamine"              
# [16] "Nitrogen.metabolism..Bacterial..aerobic.specific..ammonia.oxidation"            
# [17] "Nitrogen.metabolism..Bacterial..anaerobic.specific..ammonia.oxidation"          
# [18] "Nitrogen.metabolism..Bacterial.Archaeal.ammonia.oxidation"                      
# [19] "Nitrogen.metabolism..Nitrogen.fixation.altennative"                             
# [20] "Nitrogen.metabolism..ammonia....nitrite"                                        
# [21] "Nitrogen.metabolism..nitric.oxide....nitrous.oxide"                             
# [22] "Nitrogen.metabolism..nitrous.oxide....nitrogen"                                 
# [23] "Other.Reductases..TMAO.reductase"                                               
# [24] "Other.Reductases..arsenate.reduction..pt.2"                                     
# [25] "Other.Reductases..selenate.Chlorate.reduction"                                  
# [26] "Photosynthesis..Photosystem.I"                                                  
# [27] "Photosynthesis..Photosystem.II"                                                 
# [28] "SCFA.and.alcohol.conversions..Propionate..pt.1"                                 
# [29] "SCFA.and.alcohol.conversions..acetate..pt.3"                                    
# [30] "Sulfur.metabolism..Thiosulfate.oxidation.by.SOX.complex..thiosulfate....sulfate"

# Discard completely absent functions from further exploration
dram=DRAM[,colSums(DRAM)>0]

# MAGs in taxonomy table and dram table in same order.
mean(rownames(dram)==rownames(taxmatrix_mod))

# Load bin quality data
bin_q=read.csv("data/All_drep_bin_quality.csv")
bin_q$genome=gsub(".fa","",bin_q$genome)
bin_q=bin_q[bin_q$genome%in%rownames(dram),]
bin_q=bin_q[match(rownames(dram),bin_q$genome),]

# Barchart of bin qualities per Phylum
data.frame(MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),
           Phylum=taxmatrix_mod$Phylum,
           completeness=bin_q$completeness)%>%
  ggplot(aes(y=completeness,x=MAG,fill=Phylum,color=Phylum)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# group DRAM functions
module=names(dram)[1:12]
etc=names(dram)[13:25]
cazy=names(dram)[26:44]
metab=names(dram)[c(45:54,66:67)]
scfa=names(dram)[55:65]

### Exploration of module variables ###

# Histograms of variables
data.frame(dram[module])%>%
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  coord_cartesian(xlim=c(0,1))+
  facet_wrap(~name, scale = "free") +
  theme_bw()

# Barchart MAGs sorted by Phyla 
data.frame(dram[module],MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),Phylum=taxmatrix_mod$Phylum)%>%
  pivot_longer(cols = 1:ncol(dram[module])) %>% 
  ggplot(aes(y=value,x=reorder(MAG,Phylum),fill=Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap(~name, scale = "free") +
  coord_cartesian(ylim = c(0,1))+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Module coverage vs MAG completeness (by Phylum) 
data.frame(dram[module],
           MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),
           Phylum=taxmatrix_mod$Phylum,
           completeness=bin_q$completeness)%>%
  pivot_longer(cols = 1:ncol(dram[module])) %>% 
  ggplot(aes(y=value,x=completeness,fill=Phylum,color=Phylum)) +
  geom_smooth(method = "lm",se=FALSE)+
  facet_wrap(~name, scale = "free") +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Correlations between variables
p_mat=cor_pmat(dram[module],  method = "spearman", use = "complete.obs")
cor_mat=cor(dram[module],method = "spearman")
ggcorrplot(cor_mat,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           p.mat = p_mat,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())

# PCA
module_pca=prcomp(dram[module],scale = TRUE)
# variance explained by dimensions
fviz_screeplot(module_pca)
fviz_pca_biplot(module_pca, 
                repel = TRUE,
                label = "var",
                # geom.ind = shapes_ind, 
                pointshape = 21,
                pointsize = 3,
                col.ind = "black",
                fill.ind = taxmatrix_mod$Phylum,  # Individuals color
                palette=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"),
                # col.var = "blue", # Variables color
                mean.point = FALSE,
                ggtheme = theme_bw()
)


### Exploration of etc variables ###

# Histograms of variables
data.frame(dram[etc])%>%
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  coord_cartesian(xlim=c(0,1))+
  facet_wrap(~name, scale = "free") +
  theme_bw()

# Barchart MAGs sorted by Phyla 
data.frame(dram[etc],MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),Phylum=taxmatrix_mod$Phylum)%>%
  pivot_longer(cols = 1:ncol(dram[etc])) %>% 
  ggplot(aes(y=value,x=reorder(MAG,Phylum),fill=Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap(~name, scale = "free") +
  coord_cartesian(ylim = c(0,1))+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Module coverage vs MAG completeness (by Phylum) 
data.frame(dram[etc],
           MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),
           Phylum=taxmatrix_mod$Phylum,
           completeness=bin_q$completeness)%>%
  pivot_longer(cols = 1:ncol(dram[etc])) %>% 
  ggplot(aes(y=value,x=completeness,fill=Phylum,color=Phylum)) +
  geom_smooth(method = "lm",se=FALSE)+
  facet_wrap(~name, scale = "free") +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Correlations between variables
p_mat=cor_pmat(dram[etc],  method = "spearman", use = "complete.obs")
cor_mat=cor(dram[etc],method = "spearman")
ggcorrplot(cor_mat,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           p.mat = p_mat,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())

# PCA
module_pca=prcomp(dram[etc],scale = TRUE)
# variance explained by dimensions
fviz_screeplot(module_pca)
fviz_pca_biplot(module_pca, 
                repel = TRUE,
                label = "var",
                # geom.ind = shapes_ind, 
                pointshape = 21,
                pointsize = 3,
                col.ind = "black",
                fill.ind = taxmatrix_mod$Phylum,  # Individuals color
                palette=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"),
                # col.var = "blue", # Variables color
                mean.point = FALSE,
                ggtheme = theme_bw()
)


### Exploration of cazy variables ###

# Frequency plots
dram[cazy] %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name, value) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = value, y = count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~name, scale = "free_x") +
  theme_bw()

# Barchart of number of cazy-s in MAGs, sorted by Phyla 
data.frame(dram[cazy],MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),Phylum=taxmatrix_mod$Phylum)%>%
  pivot_longer(cols = 1:ncol(dram[cazy])) %>% 
  ggplot(aes(y=value,x=reorder(MAG,Phylum),fill=Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap(~Phylum, scale = "free") +
  coord_cartesian(ylim = c(0,15))+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Number of CAZYs vs MAG completeness (by Phylum) 
data.frame(dram[cazy],
           MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),
           Phylum=taxmatrix_mod$Phylum,
           completeness=bin_q$completeness)%>%
  pivot_longer(cols = 1:ncol(dram[cazy])) %>% 
  ggplot(aes(y=value,x=completeness,fill=Phylum,color=Phylum)) +
  geom_smooth(method = "lm",se=FALSE)+
  facet_wrap(~name, scale = "free") +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Correlations between CAZYs (binary data)
# Pearson correlation and tetrachoric correlation, an alternative that might be
# more meaningful for binary data.
cor_mat=cor(dram[cazy],method = "pearson")
cor_tetra=tetrachoric(dram[cazy])$rho
ggcorrplot(cor_mat,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())
ggcorrplot(cor_tetra,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())

# Logistic PCA for binary data
logpca_cv = cv.lpca(dram[cazy], ks = 2, ms = 1:10)
plot(logpca_cv)
logpca_model = logisticPCA(dram[cazy], k = 2, m = which.min(logpca_cv))
plot(logpca_model,"scores")+geom_point(aes(colour = taxmatrix_mod$Phylum),size=4)+
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()


### Exploration of metab variables ###

# Frequency plots
dram[metab] %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name, value) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = value, y = count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~name, scale = "free_x") +
  theme_bw()

# Barchart of number of metab pathways in MAGs, sorted by Phyla
data.frame(dram[metab],MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),Phylum=taxmatrix_mod$Phylum)%>%
  pivot_longer(cols = 1:ncol(dram[metab])) %>% 
  ggplot(aes(y=value,x=reorder(MAG,Phylum),fill=Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap(~Phylum, scale = "free") +
  coord_cartesian(ylim = c(0,8))+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Number of metab genes vs MAG completeness (by Phylum) 
data.frame(dram[metab],
           MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),
           Phylum=taxmatrix_mod$Phylum,
           completeness=bin_q$completeness)%>%
  pivot_longer(cols = 1:ncol(dram[metab])) %>% 
  ggplot(aes(y=value,x=completeness,fill=Phylum,color=Phylum)) +
  geom_smooth(method = "lm",se=FALSE)+
  facet_wrap(~name, scale = "free") +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Correlations between metab (binary data)
# Pearson correlation and tetrachoric correlation, an alternative that might be
# more meaningful for binary data.
cor_mat=cor(dram[metab],method = "pearson")
cor_tetra=tetrachoric(dram[metab])$rho
ggcorrplot(cor_mat,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())
ggcorrplot(cor_tetra,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())

# Logistic PCA for binary data
logpca_cv = cv.lpca(dram[metab], ks = 2, ms = 1:10)
plot(logpca_cv)
logpca_model = logisticPCA(dram[metab], k = 2, m = which.min(logpca_cv))
plot(logpca_model,"scores")+geom_point(aes(colour = taxmatrix_mod$Phylum),size=4)+
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()


### Exploration of scfa variables ###

# Frequency plots
dram[scfa] %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name, value) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = value, y = count)) +
  geom_bar(stat = "identity") +
  facet_wrap(~name, scale = "free_x") +
  theme_bw()

# Barchart of number of scfa-s in MAGs, sorted by Phyla 
data.frame(dram[scfa],MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),Phylum=taxmatrix_mod$Phylum)%>%
  pivot_longer(cols = 1:ncol(dram[scfa])) %>% 
  ggplot(aes(y=value,x=reorder(MAG,Phylum),fill=Phylum)) +
  geom_bar(stat = "identity")+
  facet_wrap(~Phylum, scale = "free") +
  coord_cartesian(ylim = c(0,9))+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()

# Number of SCFAs vs MAG completeness (by Phylum) 
data.frame(dram[scfa],
           MAG=factor(rownames(dram),levels = rownames(dram)[order(taxmatrix_mod$Phylum)]),
           Phylum=taxmatrix_mod$Phylum,
           completeness=bin_q$completeness)%>%
  pivot_longer(cols = 1:ncol(dram[scfa])) %>% 
  ggplot(aes(y=value,x=completeness,fill=Phylum,color=Phylum)) +
  geom_smooth(method = "lm",se=FALSE)+
  facet_wrap(~name, scale = "free") +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()


# Correlations between SCFAs (binary data)
# Pearson correlation and tetrachoric correlation, an alternative that might be
# more meaningful for binary data.
cor_mat=cor(dram[scfa],method = "pearson")
cor_tetra=tetrachoric(dram[scfa])$rho
ggcorrplot(cor_mat,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())
ggcorrplot(cor_tetra,
           outline.color = "black",
           show.diag  = F,
           hc.order = TRUE,
           type = "upper",
           lab = T,
           digits = 1,
           insig = "blank",
           ggtheme = theme_bw())

# Logistic PCA for binary data
logpca_cv = cv.lpca(dram[scfa], ks = 2, ms = 1:10)
plot(logpca_cv)
logpca_model = logisticPCA(dram[scfa], k = 2, m = which.min(logpca_cv))
plot(logpca_model,"scores")+geom_point(aes(colour = taxmatrix_mod$Phylum),size=4)+
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  theme_bw()
