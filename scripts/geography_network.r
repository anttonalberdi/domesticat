#https://datascience.blog.wzb.eu/2018/05/31/three-ways-of-visualizing-a-graph-on-a-map/

setwd("/Users/anttonalberdi/github/domesticat/")

library(assertthat)
library(dplyr)
library(purrr)
library(igraph)
library(ggplot2)
library(ggraph)
library(ggmap)

#Aruba: 12.5,-70.1
#Brazil: -10.7,-51.6
#CaboVerde: 15.8,-25.4
#Spain: 40.7,-3.5
#Denmark: 55.5,11.5
#Malaysia: 4.2, 102.3

#Load data
MAGrelative <- read.csv("data/MAG_relative.csv", row.names=1, header=TRUE)
cat_metadata=read.csv("data/DomestiCAT_metadata.csv")
head(cat_metadata)
cat_metadata=cat_metadata[-c(96:98),]
cat_metadata=cat_metadata[cat_metadata$ExtractionID!="Blank",]
cat_metadata$CombinedID <- gsub("-","_",cat_metadata$CombinedID , fixed = TRUE)

#Get location dissimilarity matrix
MAG_divpart_q0=pair_dis(MAGrelative,qvalue = 0,hierarchy=cat_metadata[,c(4,1)])
location_dis <- as.data.frame(1-MAG_divpart_q0$L2_UqN)
location_dis_table <- data.frame(t(combn(names(location_dis),2)), dist=location_dis[lower.tri(location_dis)])
location_dis_table$category <- round(location_dis_table$dist*10)
colnames(location_dis_table) <- c("from","to","weight","category")
edges <- as_tibble(location_dis_table)

#Create nodes
nodes <- as.data.frame(matrix(
  c(12.5,-70.1,-10.7,-51.6,15.8,-25.4,40.7,-3.5,55.5,11.5,4.2, 102.3),
  nrow = 6,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(c("Aruba","Brazil","CaboVerde","Spain","Denmark","Malaysia"), c("lat","lon"))
))
locations$name <- c("Aruba","Brazil","CaboVerde","Spain","Denmark","Malaysia")
locations$id <- c("Aruba","Brazil","CaboVerde","Spain","Denmark","Malaysia")
rownames(locations) <- c(1:nrow(locations))
nodes <- locations[,c(4,1,2,3)]

edges <- edges %>% mutate(category = as.factor(category))

# create the igraph graph object
g <- graph_from_data_frame(edges, directed = F, vertices = nodes)

# --------------------------------------------------------------------- #
# Common data structures and ggplot objects for all the following plots #
# --------------------------------------------------------------------- #

# create a data frame for plotting the edges
# join with nodes to get start and end positions for each
# edge (x, y and xend, yend)

edges_for_plot <- edges %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('from' = 'id')) %>%
  rename(x = lon, y = lat) %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('to' = 'id')) %>%
  rename(xend = lon, yend = lat)

assert_that(nrow(edges_for_plot) == nrow(edges))

# use the node degree for scaling the node sizes
nodes$weight = degree(g)

# common plot theme
maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#ffffff")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))

# common polygon geom for plotting the country shapes
country_shapes <- geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),
                               fill = "#f4f4f4", color = "#CECECE", size = 0.15)

# common coordinate system for all the following plots
mapcoords <- coord_fixed(xlim = c(-150, 180), ylim = c(-55, 80))

# Results in warning: "Scale for 'size' is already present. Adding another scale for
# 'size', which will replace the existing scale."

# now a plot with static node size:
ggplot(nodes) + country_shapes +
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend,     # draw edges as arcs
                 color = category, size = weight),
             data = edges_for_plot, curvature = 0.33, alpha = 0.5) +
  scale_color_manual(values=c("#bcb8ce", "#917898","#4c394f","#2e1a1e")) +
  scale_size_continuous(guide = FALSE, range = c(1, 4)) + # scale for edge widths
  geom_point(aes(x = lon, y = lat),                          # draw nodes
             shape = 21, size = 3,
             fill = 'white', color = 'black', stroke = 0.5) +
  geom_text(aes(x = lon, y = lat, label = name),             # draw text labels
            hjust = 0, nudge_x = 1, nudge_y = 4,
            size = 3, color = "black", fontface = "bold") +
  mapcoords + maptheme
