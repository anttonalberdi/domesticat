library(vegan)

# Import metadata
cat_metadata=read.csv("data/DomestiCAT metadata.csv")
head(cat_metadata)
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
# Transpose the MAG table to have MAGs as columns and samples as rows
MAGcounts_relL=data.frame(t(MAGcounts_relLength))
# Sequence counts per genome relative to total sequences in genome and in sample
MAGcounts_relL_rel=MAGcounts_relL/rowSums(MAGcounts_relL)
# Square-root relative abundance table to obtain Hellinger transformed table.
MAGcounts_relL_hel=sqrt(MAGcounts_relL_rel)

rownames(MAGcounts_relL_hel)%in%cat_metadata$CombinedID
