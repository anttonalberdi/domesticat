library(dplyr)
library(phyloseq)
library(tidyr)

## Import data
MAGcounts <- read.delim("data/DomestiCat_all_MAG_counts.txt", sep = '\t')
genomelength <- read.delim("data/DomestiCat_all_MAG_length.txt", sep = '\t')
taxonomy <- read.delim("data/gtdbtk.bac120.summary.tsv", sep = '\t')

## Coerce count table into a phylseq OTU table
# Convert to dataframe
MAGdf <- as.data.frame(MAGcounts)
# Remove first column (MAG names)
MAGdfrows <- MAGdf[,-1]
# Set MAG names as row names of the dataframe
rownames(MAGdfrows) <- MAGdf[,1] 

MAGotu <- otu_table(MAGdfrows, taxa_are_rows = T)


## Coerce gtdb-tk output to phyloseq taxonomy table
# Remove the "p__" etc. designations from the classification string
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

taxtable <- tax_table(taxmatrix) 


## Merge into a phyloseq object
physeq <- phyloseq(MAGotu, taxtable)


