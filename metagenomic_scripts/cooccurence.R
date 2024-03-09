require(SpiecEasi)
require(igraph)
require(Matrix)

meta <- fread("data/metadata/parsed_patient_metadata.filt.csv")
patients <- meta %>% filter(hap_vap_cap != "Water control")
controls <- meta %>% filter(hap_vap_cap == "Water control")

df <- fread("results/metagenomic_out/read_counts.G.zeroed.csv") %>%
  filter(run_id %in% meta$run_id)

X_patient <- df %>% 
  filter(run_id %in% patients$run_id) %>%
  column_to_rownames("run_id")

X_controls <- df %>% 
  filter(run_id %in% controls$run_id) %>%
  column_to_rownames("run_id")

get_graph <- function(X) {
  set.seed(10010)
  
  # Run SparCC
  sparcc.amgut <- sparcc(X)
  
  ## Define arbitrary threshold for SparCC correlation matrix for the graph
  sparcc.graph <- sparcc.amgut$Cor
  sparcc.graph[sparcc.graph < 0.2] <- 0
  diag(sparcc.graph) <- 0
  sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
  
  ## Create igraph object
  vertex.names <- setNames(seq(ncol(X)), colnames(X))
  ig.sparcc <- adj2igraph(sparcc.graph, vertex.attr = vertex.names)
  V(ig.sparcc)$name <- colnames(X)
  
  bad.vs <- V(ig.sparcc)[degree(ig.sparcc) == 0] 
  sparcc.filt <- delete.vertices(ig.sparcc, bad.vs)
  
  return(sparcc.filt)
}

g_patients <- get_graph(X_patient)
g_controls <- get_graph(X_controls)

plot(g_patients)

patient_edges <- E(g_patients)
to_remove <- graph.intersection(g_patients, g_controls) # Get boolean index from septic that are in healthy
idx <- which(as_ids(patient_edges) %in% as_ids(E(to_remove)))
corrected_g <- delete.edges(g_patients, idx)

# Remove edgeless vertices
bad.vs <- V(corrected_g)[degree(corrected_g) == 0] 
corrected_g <- delete.vertices(corrected_g, bad.vs)


# Annotate vertex colors
V(corrected_g)$color <- "#999999"
confirmed_v <- V(corrected_g)$name %in% confirmed
CR_v <- V(corrected_g)$name %in% CR
CR_confirmed_v <- V(corrected_g)$name %in% intersect(confirmed, CR)
V(corrected_g)[confirmed_v]$color <- "red"
V(corrected_g)[CR_v]$color <- "blue"
V(corrected_g)[CR_confirmed_v]$color <- "purple"


# Plot
png(file = "results/sparCC_networks_karius.png", 
    width = 10, 
    height = 8,
    units = 'in',
    res = 300)

clust <- cluster_edge_betweenness(corrected_g)
clust$membership[clust$membership != 1] <- NA
clust$membership[clust$names %in% c("Helicobacter", "Haemophilus")] <- 1
plot(corrected_g, 
     margin = c(0, -0.5, 0.5, 0),
     layout = layout.fruchterman.reingold(corrected_g), 
     vertex.size = 5,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     vertex.label.dist = 1,
     edge.width = E(corrected_g)$weight * 10)

legend("topright", legend=c("Confirmed pathogen", 
                            "In CR feature space", 
                            "Confirmed pathogen and in CR feature space"),
       col=c("red", "blue", "purple"),
       border = "black",
       pch = 19, cex = 1.2)

text(-0.4, 0.2, "Oral commensals", font = 2)

dev.off()