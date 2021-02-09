library(igraph)

################################################################################
#read networks ----
################################################################################

#read a matrix
x <- read.csv(file = "data/ZacharyKarateNetwork.txt") 
x <- as.matrix(x) #convert to matrix 
#get network from matrix 
g1 <- graph_from_adjacency_matrix(x, mode = "undirected")
#quick viz
plot(g1)


#read an edge list ----
g2 <- read_graph(file = "data/ZacharyKarateNetwork_edgelist.txt", 
                format = "edgelist", 
                directed=F)
g2
plot(g2)

#or read the data frame ----
df <- read.table(file = "data/ZacharyKarateNetwork_edgelist.txt")
g3 <- graph_from_data_frame(d = df, directed = F)
plot(g3)


#read a graphml ----

g4 <- read_graph(file = "data/ZacharyKarateNetwork.graphml", "graphml")
plot(g4)

################################################################################
#write networks ----
################################################################################

write.graph(graph = g4, file = "results/zachary_result.graphml", "graphml")

################################################################################
#write networks ----
################################################################################

write.graph(graph = g4, file = "results/zachary_result.graphml", "graphml")

matrix_out <- get.adjacency(graph = g4, sparse = F)
write.csv(x = matrix_out, 
          file = "results/zachary_matrix.txt", 
          quote = F, row.names = F)


################################################################################
#analyze paths ----
################################################################################

g <- read_graph(file = "data/ZacharyKarateNetwork.graphml", "graphml")

path_matrix <- shortest.paths(g)
path_matrix[1:5, 1:5]

avg_shortest_path <- average.path.length(g)

my_diameter <- diameter(g)

################################################################################
#connected components ----
################################################################################

components_g <- components(g)
components_g$membership
components_g$csize
components_g$no

#get a directed graph to compare 

directed_g <-read_graph(file = "https://raw.githubusercontent.com/guillermodeandajauregui/WorkshopAdvancedBioinformatics2021/main/data/kegg_nw.graphml", 
                        format = "graphml")

components_dg <- components(directed_g, mode = "weak")
components_dg$no

components_dg.strong <- components(directed_g, mode = "strong")
components_dg.strong$no

# degree ----

V(g) 
degree(g)

g <- set.vertex.attribute(graph = g, name = "degree", value = degree(g))
V(g)$degree_2 <- degree(g)

head(get.data.frame(g, what = "vertices"))

V(directed_g)$all_degree <-  degree(directed_g, mode = "all")
V(directed_g)$in_degree <-  degree(directed_g, mode = "in")
V(directed_g)$out_degree <-  degree(directed_g, mode = "out")

head(get.data.frame(directed_g, what = "vertices"))


#strength ---- 

#we don't have good weighted network example so let's add some random weights 

g_copy <- g

set.seed(2021)
some_weights <- rnorm(n = ecount(g_copy))
E(g_copy)$weight <- some_weights

V(g_copy)$strength <- strength(g_copy)
head(get.data.frame(g_copy, what = "vertices"))

# betweenness centrality ----

V(g)$betweenness_centrality <- betweenness(g)
V(g)$betweenness_centrality.estimate <- betweenness.estimate(graph = g, cutoff = 5)
V(g)$betweenness_centrality.estimate.badchoice <- betweenness.estimate(graph = g, cutoff = 3)
head(get.data.frame(g, what = "vertices"))

#edge betweenness ----
E(g)$edge_betweenness <- edge.betweenness(g)
head(get.data.frame(g, what = "edges"))

#clustering coefficient ---- 

V(g)$cc   <- transitivity(graph = g, type = "local", isolates = "zero")
cc_global <- transitivity(graph = g, type = "global")
head(get.data.frame(g, what = "vertices"))
#module detection ---- 

comm.louvain <- cluster_louvain(g)
comm.louvain

#add membership data ---- 
V(g)$comm.louvain <- membership(comm.louvain)
head(get.data.frame(g, what = "vertices"))

comm.infomap <- cluster_infomap(g)
comm.walktrp <- cluster_walktrap(g)
comm.fstgred <- cluster_fast_greedy(g)

#try different algorithms 
V(g)$comm.infomap <- membership(comm.infomap) #random walker based
V(g)$comm.walktrp <- membership(comm.walktrp) #random walker based
V(g)$comm.fstgred <- membership(comm.fstgred) #modularity maximization


head(get.data.frame(g, what = "vertices"))


set.seed(1)
plot(comm.louvain, g, layout = layout_nicely)
set.seed(1)
plot(comm.infomap, g, layout = layout_nicely)
set.seed(1)
plot(comm.walktrp, g, layout = layout_nicely)
set.seed(1)
plot(comm.fstgred, g, layout = layout_nicely)

#is edge inter or intra community? ---- 
E(g)$crossing_louvain <- crossing(communities = comm.louvain, graph = g)

#some nicer looking plotting 
library(ggraph)

set.seed(1)
ggraph(graph = g, layout = "kk") + 
  geom_edge_link(aes(color = crossing_louvain)) + 
  geom_node_point(aes(color = as.factor(comm.louvain))) + 
  theme_graph() + 
  scale_edge_color_manual("crossing", values = c("grey", "black")) + 
  scale_color_discrete("Community (Louvain)")

set.seed(1)
ggraph(graph = g, layout = "circle") + 
  geom_edge_link(aes(color = crossing_louvain)) + 
  geom_node_point(aes(color = as.factor(comm.louvain))) + 
  theme_graph() + 
  scale_edge_color_manual("crossing", values = c("grey", "black")) + 
  #scale_color_discrete("Community (Louvain)") + 
  scale_color_viridis(name="Community (Louvain)", discrete = T, option = "C")


set.seed(1)
ggraph(graph = g, layout = "fr") + 
  geom_edge_link(aes(color = crossing_louvain)) + 
  geom_node_point(aes(color = as.factor(comm.louvain)), size = 10 ) + 
  theme_graph() + 
  scale_edge_color_manual("crossing", values = c("grey", "black")) + 
  #scale_color_discrete("Community (Louvain)") + 
  scale_color_viridis(name="Community (Louvain)", discrete = T, option = "C")
