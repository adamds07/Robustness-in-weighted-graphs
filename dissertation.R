#Install packages necessary for analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gprege")
install.packages("robin")
library(robin)
install.packages("perturbR")
library(perturbR)
library(ggplot2)

#Read in data from a CSV file
openflights2 <- read.table("openflights2.txt", header=TRUE)
#Convert the data from a data frame to a graph
openflights2_graph <- graph.data.frame(openflights2, FALSE)
#Remove multiple edges
of2 <- simplify(openflights2_graph, remove.multiple = TRUE, edge.attr.comb = "sum")
#Check that the graph has weights
is_weighted(openflights2_graph)
#Visualise the finished network
plot.igraph(of2, vertex.size = 10, vertex.label = NA)

#############Describing the data#########################

#Finding the number of vertices in the network
vertices <- get.vertex.attribute(of2, name = "name")
length(vertices)

#Finding the number of edges in the network
edges <- get.edgelist(of2)
dim(edges)

#Calculating the global transitivity of the network
transitivity(of2, "global")

#Calculating the betweenness centrality of each node and printing top 5
ofbetween <- betweenness(of2, v=V(of2), weights = of2$weight)
ofbetsort <- sort(ofbetween, decreasing = TRUE)
head(ofbetsort, 5)

#Calculating the eigenvector centrality of each node and printing top 5
ofeigen <- eigen_centrality(of2, weights = of2$weight)
ofeigsort <- sort(ofeigen$vector, decreasing = TRUE)
head(ofeigsort, 5)

#Calculating the weighted degree centrality of each node and printing the top 5
ofstre <- graph.strength(of2)
ofstresort <- sort(ofstre, decreasing = TRUE)
head(ofstresort, 5)

#Calculating the eccentricity of the network, the diameter and the radius
ofecc <- eccentricity(of2, V(of2))
max(ofecc)
min(ofecc)

#Using Fast Greedy and Louvain to find a clustering of the network
clustfg <- cluster_fast_greedy(of2)
clustlv <- cluster_louvain(of2)

#Using VI to compare the clusterings
compare(clustfg, clustlv, method = "VI")

#Finding the modularity of both clusterings
modularity(of2,clustfg$membership)
modularity(of2,clustlv$membership)

#Plotting the network colour coded by the Louvain clustering
my_collv <- as.numeric(as.factor(clustlv$membership))
plot.igraph(of2, vertex.color = my_collv, vertex.label = NA)

#Plotting the network colour coded by the Fast Greedy clustering
my_colfg <- as.numeric(as.factor(clustfg$membership))
plot.igraph(of2, vertex.color = my_colfg, vertex.label = NA)

##################Testing the robustness#######################

#Creating the adjacency matrix of of2
of2adj <- as_adjacency_matrix(of2, type = "both", attr = "weight", edges = FALSE, names = FALSE, sparse = FALSE)

#Rewiring the adjacency matrix to create the weighted configuration model
of2rewire <- rewireR(of2adj,4317391,dist="NegBinom")

#Creating a graph of the weighted configuration model
nullmatrix <- as.matrix(of2rewire)
nullrandommodel <- graph_from_adjacency_matrix(nullmatrix, mode = "undirected", weighted = TRUE)

#Storing the edge weights of the original network and null model
of2w<- E(of2)$weight
nullw<- E(nullrandommodel)$weight

#Using Louvain algorithm to cluster of2 and null model
of2c <- cluster_louvain(of2, weights = of2w)
nullc <- cluster_louvain(nullrandommodel, weights = nullw)

#Setting up vectors to store the VI at each perturbation level
VIC0 <- c(0,0,0,0,0,0,0,0,0,0,0)
VIC5 <- c()
VIC10 <- c()
VIC15 <- c()
VIC20 <- c()
VIC25 <- c()
VIC30 <- c()
VIC35 <- c()
VIC40 <- c()
VIC45 <- c()
VIC50 <- c()
VIC55 <- c()
VIC60 <- c()

VICrandom0 <- c(0,0,0,0,0,0,0,0,0,0,0)
VICrandom5 <- c()
VICrandom10 <- c()
VICrandom15 <- c()
VICrandom20 <- c()
VICrandom25 <- c()
VICrandom30 <- c()
VICrandom35 <- c()
VICrandom40 <- c()
VICrandom45 <- c()
VICrandom50 <- c()
VICrandom55 <- c()
VICrandom60 <- c()

#################5% Perturbation level#####################

#First trial at 5% perturbation level
#Perturbing 5% of edges in the original network
of2perturb5a <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5am <- as.matrix(of2perturb5a)
of2p5a <- graph_from_adjacency_matrix(of2p5am, mode = "undirected", weighted = TRUE)
#Perturbing 5% of edges in the null random model
nullperturb5a <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5am <- as.matrix(nullperturb5a)
nullp5a <- graph_from_adjacency_matrix(nullp5am, mode = "undirected", weighted = TRUE)

#Clustering the perturbed networks with the Louvain algorithm
of2c5a <- cluster_louvain(of2p5a, weights = E(of2p5a)$weight)
nullc5a <- cluster_louvain(nullp5a, weights = E(nullp5a)$weight)

#Calculating the VI between the original networks and the perturbed networks
VIC5[2] <- compare(of2c,of2c5a,method = "vi")
VICrandom5[2] <- compare(nullc, nullc5a, method = "vi")

#Removing the variables used for this part of the analysis from the global environment to save space
rm(of2perturb5a,of2p5a,of2p5am,nullperturb5a,nullp5am,nullp5a,of2c5a,nullc5a)

#We now repeat this same process ten times for each perturbation level

#SECOND TRIAL
of2perturb5b <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5bm <- as.matrix(of2perturb5b)
of2p5b <- graph_from_adjacency_matrix(of2p5bm, mode = "undirected", weighted = TRUE)
nullperturb5b <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5bm <- as.matrix(nullperturb5b)
nullp5b <- graph_from_adjacency_matrix(nullp5bm, mode = "undirected", weighted = TRUE)

of2c5b <- cluster_louvain(of2p5b, weights = E(of2p5b)$weight)
nullc5b <- cluster_louvain(nullp5b, weights = E(nullp5b)$weight)

VIC5[3] <- compare(of2c,of2c5b,method = "vi")
VICrandom5[3] <- compare(nullc, nullc5b, method = "vi")

rm(of2perturb5b,of2p5b,of2p5bm,nullperturb5b,nullp5bm,nullp5b,of2c5b,nullc5b)

#THIRD TRIAL
of2perturb5c <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5cm <- as.matrix(of2perturb5c)
of2p5c <- graph_from_adjacency_matrix(of2p5cm, mode = "undirected", weighted = TRUE)
nullperturb5c <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5cm <- as.matrix(nullperturb5c)
nullp5c <- graph_from_adjacency_matrix(nullp5cm, mode = "undirected", weighted = TRUE)

of2c5c <- cluster_louvain(of2p5c, weights = E(of2p5c)$weight)
nullc5c <- cluster_louvain(nullp5c, weights = E(nullp5c)$weight)

VIC5[4] <- compare(of2c,of2c5c,method = "vi")
VICrandom5[4] <- compare(nullc, nullc5c, method = "vi")

rm(of2perturb5c,of2p5c,of2p5cm,nullperturb5c,nullp5cm,nullp5c,of2c5c,nullc5c)

#FOURTH TRIAL
of2perturb5d <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5dm <- as.matrix(of2perturb5d)
of2p5d <- graph_from_adjacency_matrix(of2p5dm, mode = "undirected", weighted = TRUE)
nullperturb5d <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5dm <- as.matrix(nullperturb5d)
nullp5d <- graph_from_adjacency_matrix(nullp5dm, mode = "undirected", weighted = TRUE)

of2c5d <- cluster_louvain(of2p5d, weights = E(of2p5d)$weight)
nullc5d <- cluster_louvain(nullp5d, weights = E(nullp5d)$weight)

VIC5[5] <- compare(of2c,of2c5d,method = "vi")
VICrandom5[5] <- compare(nullc, nullc5d, method = "vi")

rm(of2perturb5d,of2p5d,of2p5dm,nullperturb5d,nullp5dm,nullp5d,of2c5d,nullc5d)

#FIFTH TRIAL
of2perturb5e <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5em <- as.matrix(of2perturb5e)
of2p5e <- graph_from_adjacency_matrix(of2p5em, mode = "undirected", weighted = TRUE)
nullperturb5e <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5em <- as.matrix(nullperturb5e)
nullp5e <- graph_from_adjacency_matrix(nullp5em, mode = "undirected", weighted = TRUE)

of2c5e <- cluster_louvain(of2p5e, weights = E(of2p5e)$weight)
nullc5e <- cluster_louvain(nullp5e, weights = E(nullp5e)$weight)

VIC5[6] <- compare(of2c,of2c5e,method = "vi")
VICrandom5[6] <- compare(nullc, nullc5e, method = "vi")

rm(of2perturb5e,of2p5e,of2p5em,nullperturb5e,nullp5em,nullp5e,of2c5e,nullc5e)

#SIXTH TRIAL
of2perturb5f <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5fm <- as.matrix(of2perturb5f)
of2p5f <- graph_from_adjacency_matrix(of2p5fm, mode = "undirected", weighted = TRUE)
nullperturb5f <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5fm <- as.matrix(nullperturb5f)
nullp5f <- graph_from_adjacency_matrix(nullp5fm, mode = "undirected", weighted = TRUE)

of2c5f <- cluster_louvain(of2p5f, weights = E(of2p5f)$weight)
nullc5f <- cluster_louvain(nullp5f, weights = E(nullp5f)$weight)

VIC5[7] <- compare(of2c,of2c5f,method = "vi")
VICrandom5[7] <- compare(nullc, nullc5f, method = "vi")

rm(of2perturb5f,of2p5f,of2p5fm,nullperturb5f,nullp5fm,nullp5f,of2c5f,nullc5f)

#SEVENTH TRIAL
of2perturb5g <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5gm <- as.matrix(of2perturb5g)
of2p5g <- graph_from_adjacency_matrix(of2p5gm, mode = "undirected", weighted = TRUE)
nullperturb5g <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5gm <- as.matrix(nullperturb5g)
nullp5g <- graph_from_adjacency_matrix(nullp5gm, mode = "undirected", weighted = TRUE)

of2c5g <- cluster_louvain(of2p5g, weights = E(of2p5g)$weight)
nullc5g <- cluster_louvain(nullp5g, weights = E(nullp5g)$weight)

VIC5[8] <- compare(of2c,of2c5g,method = "vi")
VICrandom5[8] <- compare(nullc, nullc5g, method = "vi")

rm(of2perturb5g,of2p5g,of2p5gm,nullperturb5g,nullp5gm,nullp5g,of2c5g,nullc5g)

#EIGHTH TRIAL
of2perturb5h <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5hm <- as.matrix(of2perturb5h)
of2p5h <- graph_from_adjacency_matrix(of2p5hm, mode = "undirected", weighted = TRUE)
nullperturb5h <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5hm <- as.matrix(nullperturb5h)
nullp5h <- graph_from_adjacency_matrix(nullp5hm, mode = "undirected", weighted = TRUE)

of2c5h <- cluster_louvain(of2p5h, weights = E(of2p5h)$weight)
nullc5h <- cluster_louvain(nullp5h, weights = E(nullp5h)$weight)

VIC5[9] <- compare(of2c,of2c5h,method = "vi")
VICrandom5[9] <- compare(nullc, nullc5h, method = "vi")

rm(of2perturb5h,of2p5h,of2p5hm,nullperturb5h,nullp5hm,nullp5h,of2c5h,nullc5h)

#NINTH TRIAL
of2perturb5i <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5im <- as.matrix(of2perturb5i)
of2p5i <- graph_from_adjacency_matrix(of2p5im, mode = "undirected", weighted = TRUE)
nullperturb5i <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5im <- as.matrix(nullperturb5i)
nullp5i <- graph_from_adjacency_matrix(nullp5im, mode = "undirected", weighted = TRUE)

of2c5i <- cluster_louvain(of2p5i, weights = E(of2p5i)$weight)
nullc5i <- cluster_louvain(nullp5i, weights = E(nullp5i)$weight)

VIC5[10] <- compare(of2c,of2c5i,method = "vi")
VICrandom5[10] <- compare(nullc, nullc5i, method = "vi")

rm(of2perturb5i,of2p5i,of2p5im,nullperturb5i,nullp5im,nullp5i,of2c5i,nullc5i)

#TENTH TRIAL
of2perturb5j <- rewireR(of2adj, 215869, dist = "NegBinom")
of2p5jm <- as.matrix(of2perturb5j)
of2p5j <- graph_from_adjacency_matrix(of2p5jm, mode = "undirected", weighted = TRUE)
nullperturb5j <- rewireR(nullmatrix, 215869, dist = "NegBinom")
nullp5jm <- as.matrix(nullperturb5j)
nullp5j <- graph_from_adjacency_matrix(nullp5jm, mode = "undirected", weighted = TRUE)

of2c5j <- cluster_louvain(of2p5j, weights = E(of2p5j)$weight)
nullc5j <- cluster_louvain(nullp5j, weights = E(nullp5j)$weight)

VIC5[11] <- compare(of2c,of2c5j,method = "vi")
VICrandom5[11] <- compare(nullc, nullc5j, method = "vi")

rm(of2perturb5j,of2p5j,of2p5jm,nullperturb5j,nullp5jm,nullp5j,of2c5j,nullc5j)

###############10% Perturbation level#####################

#FIRST TRIAL
of2perturb10a <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10am <- as.matrix(of2perturb10a)
of2p10a <- graph_from_adjacency_matrix(of2p10am, mode = "undirected", weighted = TRUE)
nullperturb10a <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10am <- as.matrix(nullperturb10a)
nullp10a <- graph_from_adjacency_matrix(nullp10am, mode = "undirected", weighted = TRUE)

of2c10a <- cluster_louvain(of2p10a, weights = E(of2p10a)$weight)
nullc10a <- cluster_louvain(nullp10a, weights = E(nullp10a)$weight)

VIC10[2] <- compare(of2c,of2c10a,method = "vi")
VICrandom10[2] <- compare(nullc, nullc10a, method = "vi")

rm(of2perturb10a,of2p10a,of2p10am,nullperturb10a,nullp10am,nullp10a,of2c10a,nullc10a)

#SECOND TRIAL
of2perturb10b <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10bm <- as.matrix(of2perturb10b)
of2p10b <- graph_from_adjacency_matrix(of2p10bm, mode = "undirected", weighted = TRUE)
nullperturb10b <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10bm <- as.matrix(nullperturb10b)
nullp10b <- graph_from_adjacency_matrix(nullp10bm, mode = "undirected", weighted = TRUE)

of2c10b <- cluster_louvain(of2p10b, weights = E(of2p10b)$weight)
nullc10b <- cluster_louvain(nullp10b, weights = E(nullp10b)$weight)

VIC10[3] <- compare(of2c,of2c10b,method = "vi")
VICrandom10[3] <- compare(nullc, nullc10b, method = "vi")

rm(of2perturb10b,of2p10b,of2p10bm,nullperturb10b,nullp10bm,nullp10b,of2c10b,nullc10b)

#THIRD TRIAL
of2perturb10c <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10cm <- as.matrix(of2perturb10c)
of2p10c <- graph_from_adjacency_matrix(of2p10cm, mode = "undirected", weighted = TRUE)
nullperturb10c <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10cm <- as.matrix(nullperturb10c)
nullp10c <- graph_from_adjacency_matrix(nullp10cm, mode = "undirected", weighted = TRUE)

of2c10c <- cluster_louvain(of2p10c, weights = E(of2p10c)$weight)
nullc10c <- cluster_louvain(nullp10c, weights = E(nullp10c)$weight)

VIC10[4] <- compare(of2c,of2c10c,method = "vi")
VICrandom10[4] <- compare(nullc, nullc10c, method = "vi")

rm(of2perturb10c,of2p10c,of2p10cm,nullperturb10c,nullp10cm,nullp10c,of2c10c,nullc10c)

#FOURTH TRIAL
of2perturb10d <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10dm <- as.matrix(of2perturb10d)
of2p10d <- graph_from_adjacency_matrix(of2p10dm, mode = "undirected", weighted = TRUE)
nullperturb10d <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10dm <- as.matrix(nullperturb10d)
nullp10d <- graph_from_adjacency_matrix(nullp10dm, mode = "undirected", weighted = TRUE)

of2c10d <- cluster_louvain(of2p10d, weights = E(of2p10d)$weight)
nullc10d <- cluster_louvain(nullp10d, weights = E(nullp10d)$weight)

VIC10[5] <- compare(of2c,of2c10d,method = "vi")
VICrandom10[5] <- compare(nullc, nullc10d, method = "vi")

rm(of2perturb10d,of2p10d,of2p10dm,nullperturb10d,nullp10dm,nullp10d,of2c10d,nullc10d)

#FIFTH TRIAL
of2perturb10e <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10em <- as.matrix(of2perturb10e)
of2p10e <- graph_from_adjacency_matrix(of2p10em, mode = "undirected", weighted = TRUE)
nullperturb10e <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10em <- as.matrix(nullperturb10e)
nullp10e <- graph_from_adjacency_matrix(nullp10em, mode = "undirected", weighted = TRUE)

of2c10e <- cluster_louvain(of2p10e, weights = E(of2p10e)$weight)
nullc10e <- cluster_louvain(nullp10e, weights = E(nullp10e)$weight)

VIC10[6] <- compare(of2c,of2c10e,method = "vi")
VICrandom10[6] <- compare(nullc, nullc10e, method = "vi")

rm(of2perturb10e,of2p10e,of2p10em,nullperturb10e,nullp10em,nullp10e,of2c10e,nullc10e)

#SIXTH TRIAL
of2perturb10f <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10fm <- as.matrix(of2perturb10f)
of2p10f <- graph_from_adjacency_matrix(of2p10fm, mode = "undirected", weighted = TRUE)
nullperturb10f <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10fm <- as.matrix(nullperturb10f)
nullp10f <- graph_from_adjacency_matrix(nullp10fm, mode = "undirected", weighted = TRUE)

of2c10f <- cluster_louvain(of2p10f, weights = E(of2p10f)$weight)
nullc10f <- cluster_louvain(nullp10f, weights = E(nullp10f)$weight)

VIC10[7] <- compare(of2c,of2c10f,method = "vi")
VICrandom10[7] <- compare(nullc, nullc10f, method = "vi")

rm(of2perturb10f,of2p10f,of2p10fm,nullperturb10f,nullp10fm,nullp10f,of2c10f,nullc10f)

#SEVENTH TRIAL
of2perturb10g <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10gm <- as.matrix(of2perturb10g)
of2p10g <- graph_from_adjacency_matrix(of2p10gm, mode = "undirected", weighted = TRUE)
nullperturb10g <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10gm <- as.matrix(nullperturb10g)
nullp10g <- graph_from_adjacency_matrix(nullp10gm, mode = "undirected", weighted = TRUE)

of2c10g <- cluster_louvain(of2p10g, weights = E(of2p10g)$weight)
nullc10g <- cluster_louvain(nullp10g, weights = E(nullp10g)$weight)

VIC10[8] <- compare(of2c,of2c10g,method = "vi")
VICrandom10[8] <- compare(nullc, nullc10g, method = "vi")

rm(of2perturb10g,of2p10g,of2p10gm,nullperturb10g,nullp10gm,nullp10g,of2c10g,nullc10g)

#EIGHTH TRIAL
of2perturb10h <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10hm <- as.matrix(of2perturb10h)
of2p10h <- graph_from_adjacency_matrix(of2p10hm, mode = "undirected", weighted = TRUE)
nullperturb10h <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10hm <- as.matrix(nullperturb10h)
nullp10h <- graph_from_adjacency_matrix(nullp10hm, mode = "undirected", weighted = TRUE)

of2c10h <- cluster_louvain(of2p10h, weights = E(of2p10h)$weight)
nullc10h <- cluster_louvain(nullp10h, weights = E(nullp10h)$weight)

VIC10[9] <- compare(of2c,of2c10h,method = "vi")
VICrandom10[9] <- compare(nullc, nullc10h, method = "vi")

rm(of2perturb10h,of2p10h,of2p10hm,nullperturb10h,nullp10hm,nullp10h,of2c10h,nullc10h)

#NINTH TRIAL
of2perturb10i <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10im <- as.matrix(of2perturb10i)
of2p10i <- graph_from_adjacency_matrix(of2p10im, mode = "undirected", weighted = TRUE)
nullperturb10i <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10im <- as.matrix(nullperturb10i)
nullp10i <- graph_from_adjacency_matrix(nullp10im, mode = "undirected", weighted = TRUE)

of2c10i <- cluster_louvain(of2p10i, weights = E(of2p10i)$weight)
nullc10i <- cluster_louvain(nullp10i, weights = E(nullp10i)$weight)

VIC10[10] <- compare(of2c,of2c10i,method = "vi")
VICrandom10[10] <- compare(nullc, nullc10i, method = "vi")

rm(of2perturb10i,of2p10i,of2p10im,nullperturb10i,nullp10im,nullp10i,of2c10i,nullc10i)

#TENTH TRIAL
of2perturb10j <- rewireR(of2adj, 431739, dist = "NegBinom")
of2p10jm <- as.matrix(of2perturb10j)
of2p10j <- graph_from_adjacency_matrix(of2p10jm, mode = "undirected", weighted = TRUE)
nullperturb10j <- rewireR(nullmatrix, 431739, dist = "NegBinom")
nullp10jm <- as.matrix(nullperturb10j)
nullp10j <- graph_from_adjacency_matrix(nullp10jm, mode = "undirected", weighted = TRUE)

of2c10j <- cluster_louvain(of2p10j, weights = E(of2p10j)$weight)
nullc10j <- cluster_louvain(nullp10j, weights = E(nullp10j)$weight)

VIC10[11] <- compare(of2c,of2c10j,method = "vi")
VICrandom10[11] <- compare(nullc, nullc10j, method = "vi")

rm(of2perturb10j,of2p10j,of2p10jm,nullperturb10j,nullp10jm,nullp10j,of2c10j,nullc10j)

################15% Perturbation Level#####################

#FIRST TRIAL
of2perturb15a <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15am <- as.matrix(of2perturb15a)
of2p15a <- graph_from_adjacency_matrix(of2p15am, mode = "undirected", weighted = TRUE)
nullperturb15a <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15am <- as.matrix(nullperturb15a)
nullp15a <- graph_from_adjacency_matrix(nullp15am, mode = "undirected", weighted = TRUE)

of2c15a <- cluster_louvain(of2p15a, weights = E(of2p15a)$weight)
nullc15a <- cluster_louvain(nullp15a, weights = E(nullp15a)$weight)

VIC15[2] <- compare(of2c,of2c15a,method = "vi")
VICrandom15[2] <- compare(nullc, nullc15a, method = "vi")

rm(of2perturb15a,of2p15a,of2p15am,nullperturb15a,nullp15am,nullp15a,of2c15a,nullc15a)

#SECOND TRIAL
of2perturb15b <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15bm <- as.matrix(of2perturb15b)
of2p15b <- graph_from_adjacency_matrix(of2p15bm, mode = "undirected", weighted = TRUE)
nullperturb15b <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15bm <- as.matrix(nullperturb15b)
nullp15b <- graph_from_adjacency_matrix(nullp15bm, mode = "undirected", weighted = TRUE)

of2c15b <- cluster_louvain(of2p15b, weights = E(of2p15b)$weight)
nullc15b <- cluster_louvain(nullp15b, weights = E(nullp15b)$weight)

VIC15[3] <- compare(of2c,of2c15b,method = "vi")
VICrandom15[3] <- compare(nullc, nullc15b, method = "vi")

rm(of2perturb15b,of2p15b,of2p15bm,nullperturb15b,nullp15bm,nullp15b,of2c15b,nullc15b)

#THIRD TRIAL
of2perturb15c <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15cm <- as.matrix(of2perturb15c)
of2p15c <- graph_from_adjacency_matrix(of2p15cm, mode = "undirected", weighted = TRUE)
nullperturb15c <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15cm <- as.matrix(nullperturb15c)
nullp15c <- graph_from_adjacency_matrix(nullp15cm, mode = "undirected", weighted = TRUE)

of2c15c <- cluster_louvain(of2p15c, weights = E(of2p15c)$weight)
nullc15c <- cluster_louvain(nullp15c, weights = E(nullp15c)$weight)

VIC15[4] <- compare(of2c,of2c15c,method = "vi")
VICrandom15[4] <- compare(nullc, nullc15c, method = "vi")

rm(of2perturb15c,of2p15c,of2p15cm,nullperturb15c,nullp15cm,nullp15c,of2c15c,nullc15c)

#FOURTH TRIAL
of2perturb15d <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15dm <- as.matrix(of2perturb15d)
of2p15d <- graph_from_adjacency_matrix(of2p15dm, mode = "undirected", weighted = TRUE)
nullperturb15d <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15dm <- as.matrix(nullperturb15d)
nullp15d <- graph_from_adjacency_matrix(nullp15dm, mode = "undirected", weighted = TRUE)

of2c15d <- cluster_louvain(of2p15d, weights = E(of2p15d)$weight)
nullc15d <- cluster_louvain(nullp15d, weights = E(nullp15d)$weight)

VIC15[5] <- compare(of2c,of2c15d,method = "vi")
VICrandom15[5] <- compare(nullc, nullc15d, method = "vi")

rm(of2perturb15d,of2p15d,of2p15dm,nullperturb15d,nullp15dm,nullp15d,of2c15d,nullc15d)

#FIFTH TRIAL
of2perturb15e <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15em <- as.matrix(of2perturb15e)
of2p15e <- graph_from_adjacency_matrix(of2p15em, mode = "undirected", weighted = TRUE)
nullperturb15e <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15em <- as.matrix(nullperturb15e)
nullp15e <- graph_from_adjacency_matrix(nullp15em, mode = "undirected", weighted = TRUE)

of2c15e <- cluster_louvain(of2p15e, weights = E(of2p15e)$weight)
nullc15e <- cluster_louvain(nullp15e, weights = E(nullp15e)$weight)

VIC15[6] <- compare(of2c,of2c15e,method = "vi")
VICrandom15[6] <- compare(nullc, nullc15e, method = "vi")

rm(of2perturb15e,of2p15e,of2p15em,nullperturb15e,nullp15em,nullp15e,of2c15e,nullc15e)

#SIXTH TRIAL
of2perturb15f <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15fm <- as.matrix(of2perturb15f)
of2p15f <- graph_from_adjacency_matrix(of2p15fm, mode = "undirected", weighted = TRUE)
nullperturb15f <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15fm <- as.matrix(nullperturb15f)
nullp15f <- graph_from_adjacency_matrix(nullp15fm, mode = "undirected", weighted = TRUE)

of2c15f <- cluster_louvain(of2p15f, weights = E(of2p15f)$weight)
nullc15f <- cluster_louvain(nullp15f, weights = E(nullp15f)$weight)

VIC15[7] <- compare(of2c,of2c15f,method = "vi")
VICrandom15[7] <- compare(nullc, nullc15f, method = "vi")

rm(of2perturb15f,of2p15f,of2p15fm,nullperturb15f,nullp15fm,nullp15f,of2c15f,nullc15f)

#SEVENTH TRIAL
of2perturb15g <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15gm <- as.matrix(of2perturb15g)
of2p15g <- graph_from_adjacency_matrix(of2p15gm, mode = "undirected", weighted = TRUE)
nullperturb15g <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15gm <- as.matrix(nullperturb15g)
nullp15g <- graph_from_adjacency_matrix(nullp15gm, mode = "undirected", weighted = TRUE)

of2c15g <- cluster_louvain(of2p15g, weights = E(of2p15g)$weight)
nullc15g <- cluster_louvain(nullp15g, weights = E(nullp15g)$weight)

VIC15[8] <- compare(of2c,of2c15g,method = "vi")
VICrandom15[8] <- compare(nullc, nullc15g, method = "vi")

rm(of2perturb15g,of2p15g,of2p15gm,nullperturb15g,nullp15gm,nullp15g,of2c15g,nullc15g)

#EIGHTH TRIAL
of2perturb15h <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15hm <- as.matrix(of2perturb15h)
of2p15h <- graph_from_adjacency_matrix(of2p15hm, mode = "undirected", weighted = TRUE)
nullperturb15h <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15hm <- as.matrix(nullperturb15h)
nullp15h <- graph_from_adjacency_matrix(nullp15hm, mode = "undirected", weighted = TRUE)

of2c15h <- cluster_louvain(of2p15h, weights = E(of2p15h)$weight)
nullc15h <- cluster_louvain(nullp15h, weights = E(nullp15h)$weight)

VIC15[9] <- compare(of2c,of2c15h,method = "vi")
VICrandom15[9] <- compare(nullc, nullc15h, method = "vi")

rm(of2perturb15h,of2p15h,of2p15hm,nullperturb15h,nullp15hm,nullp15h,of2c15h,nullc15h)

#NINTH TRIAL
of2perturb15i <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15im <- as.matrix(of2perturb15i)
of2p15i <- graph_from_adjacency_matrix(of2p15im, mode = "undirected", weighted = TRUE)
nullperturb15i <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15im <- as.matrix(nullperturb15i)
nullp15i <- graph_from_adjacency_matrix(nullp15im, mode = "undirected", weighted = TRUE)

of2c15i <- cluster_louvain(of2p15i, weights = E(of2p15i)$weight)
nullc15i <- cluster_louvain(nullp15i, weights = E(nullp15i)$weight)

VIC15[10] <- compare(of2c,of2c15i,method = "vi")
VICrandom15[10] <- compare(nullc, nullc15i, method = "vi")

rm(of2perturb15i,of2p15i,of2p15im,nullperturb15i,nullp15im,nullp15i,of2c15i,nullc15i)

#TENTH TRIAL
of2perturb15j <- rewireR(of2adj, 647608, dist = "NegBinom")
of2p15jm <- as.matrix(of2perturb15j)
of2p15j <- graph_from_adjacency_matrix(of2p15jm, mode = "undirected", weighted = TRUE)
nullperturb15j <- rewireR(nullmatrix, 647608, dist = "NegBinom")
nullp15jm <- as.matrix(nullperturb15j)
nullp15j <- graph_from_adjacency_matrix(nullp15jm, mode = "undirected", weighted = TRUE)

of2c15j <- cluster_louvain(of2p15j, weights = E(of2p15j)$weight)
nullc15j <- cluster_louvain(nullp15j, weights = E(nullp15j)$weight)

VIC15[11] <- compare(of2c,of2c15j,method = "vi")
VICrandom15[11] <- compare(nullc, nullc15j, method = "vi")

rm(of2perturb15j,of2p15j,of2p15jm,nullperturb15j,nullp15jm,nullp15j,of2c15j,nullc15j)

################20% Perturbation Level#######################

#FIRST TRIAL
of2perturb20a <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20am <- as.matrix(of2perturb20a)
of2p20a <- graph_from_adjacency_matrix(of2p20am, mode = "undirected", weighted = TRUE)
nullperturb20a <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20am <- as.matrix(nullperturb20a)
nullp20a <- graph_from_adjacency_matrix(nullp20am, mode = "undirected", weighted = TRUE)

of2c20a <- cluster_louvain(of2p20a, weights = E(of2p20a)$weight)
nullc20a <- cluster_louvain(nullp20a, weights = E(nullp20a)$weight)

VIC20[2] <- compare(of2c,of2c20a,method = "vi")
VICrandom20[2] <- compare(nullc, nullc20a, method = "vi")

rm(of2perturb20a,of2p20a,of2p20am,nullperturb20a,nullp20am,nullp20a,of2c20a,nullc20a)

#SECOND TRIAL
of2perturb20b <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20bm <- as.matrix(of2perturb20b)
of2p20b <- graph_from_adjacency_matrix(of2p20bm, mode = "undirected", weighted = TRUE)
nullperturb20b <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20bm <- as.matrix(nullperturb20b)
nullp20b <- graph_from_adjacency_matrix(nullp20bm, mode = "undirected", weighted = TRUE)

of2c20b <- cluster_louvain(of2p20b, weights = E(of2p20b)$weight)
nullc20b <- cluster_louvain(nullp20b, weights = E(nullp20b)$weight)

VIC20[3] <- compare(of2c,of2c20b,method = "vi")
VICrandom20[3] <- compare(nullc, nullc20b, method = "vi")

rm(of2perturb20b,of2p20b,of2p20bm,nullperturb20b,nullp20bm,nullp20b,of2c20b,nullc20b)

#THIRD TRIAL
of2perturb20c <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20cm <- as.matrix(of2perturb20c)
of2p20c <- graph_from_adjacency_matrix(of2p20cm, mode = "undirected", weighted = TRUE)
nullperturb20c <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20cm <- as.matrix(nullperturb20c)
nullp20c <- graph_from_adjacency_matrix(nullp20cm, mode = "undirected", weighted = TRUE)

of2c20c <- cluster_louvain(of2p20c, weights = E(of2p20c)$weight)
nullc20c <- cluster_louvain(nullp20c, weights = E(nullp20c)$weight)

VIC20[4] <- compare(of2c,of2c20c,method = "vi")
VICrandom20[4] <- compare(nullc, nullc20c, method = "vi")

rm(of2perturb20c,of2p20c,of2p20cm,nullperturb20c,nullp20cm,nullp20c,of2c20c,nullc20c)

#FOURTH TRIAL
of2perturb20d <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20dm <- as.matrix(of2perturb20d)
of2p20d <- graph_from_adjacency_matrix(of2p20dm, mode = "undirected", weighted = TRUE)
nullperturb20d <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20dm <- as.matrix(nullperturb20d)
nullp20d <- graph_from_adjacency_matrix(nullp20dm, mode = "undirected", weighted = TRUE)

of2c20d <- cluster_louvain(of2p20d, weights = E(of2p20d)$weight)
nullc20d <- cluster_louvain(nullp20d, weights = E(nullp20d)$weight)

VIC20[5] <- compare(of2c,of2c20d,method = "vi")
VICrandom20[5] <- compare(nullc, nullc20d, method = "vi")

rm(of2perturb20d,of2p20d,of2p20dm,nullperturb20d,nullp20dm,nullp20d,of2c20d,nullc20d)

#FIFTH TRIAL
of2perturb20e <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20em <- as.matrix(of2perturb20e)
of2p20e <- graph_from_adjacency_matrix(of2p20em, mode = "undirected", weighted = TRUE)
nullperturb20e <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20em <- as.matrix(nullperturb20e)
nullp20e <- graph_from_adjacency_matrix(nullp20em, mode = "undirected", weighted = TRUE)

of2c20e <- cluster_louvain(of2p20e, weights = E(of2p20e)$weight)
nullc20e <- cluster_louvain(nullp20e, weights = E(nullp20e)$weight)

VIC20[6] <- compare(of2c,of2c20e,method = "vi")
VICrandom20[6] <- compare(nullc, nullc20e, method = "vi")

rm(of2perturb20e,of2p20e,of2p20em,nullperturb20e,nullp20em,nullp20e,of2c20e,nullc20e)

#SIXTH TRIAL
of2perturb20f <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20fm <- as.matrix(of2perturb20f)
of2p20f <- graph_from_adjacency_matrix(of2p20fm, mode = "undirected", weighted = TRUE)
nullperturb20f <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20fm <- as.matrix(nullperturb20f)
nullp20f <- graph_from_adjacency_matrix(nullp20fm, mode = "undirected", weighted = TRUE)

of2c20f <- cluster_louvain(of2p20f, weights = E(of2p20f)$weight)
nullc20f <- cluster_louvain(nullp20f, weights = E(nullp20f)$weight)

VIC20[7] <- compare(of2c,of2c20f,method = "vi")
VICrandom20[7] <- compare(nullc, nullc20f, method = "vi")

rm(of2perturb20f,of2p20f,of2p20fm,nullperturb20f,nullp20fm,nullp20f,of2c20f,nullc20f)

#SEVENTH TRIAL
of2perturb20g <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20gm <- as.matrix(of2perturb20g)
of2p20g <- graph_from_adjacency_matrix(of2p20gm, mode = "undirected", weighted = TRUE)
nullperturb20g <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20gm <- as.matrix(nullperturb20g)
nullp20g <- graph_from_adjacency_matrix(nullp20gm, mode = "undirected", weighted = TRUE)

of2c20g <- cluster_louvain(of2p20g, weights = E(of2p20g)$weight)
nullc20g <- cluster_louvain(nullp20g, weights = E(nullp20g)$weight)

VIC20[8] <- compare(of2c,of2c20g,method = "vi")
VICrandom20[8] <- compare(nullc, nullc20g, method = "vi")

rm(of2perturb20g,of2p20g,of2p20gm,nullperturb20g,nullp20gm,nullp20g,of2c20g,nullc20g)

#EIGHTH TRIAL
of2perturb20h <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20hm <- as.matrix(of2perturb20h)
of2p20h <- graph_from_adjacency_matrix(of2p20hm, mode = "undirected", weighted = TRUE)
nullperturb20h <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20hm <- as.matrix(nullperturb20h)
nullp20h <- graph_from_adjacency_matrix(nullp20hm, mode = "undirected", weighted = TRUE)

of2c20h <- cluster_louvain(of2p20h, weights = E(of2p20h)$weight)
nullc20h <- cluster_louvain(nullp20h, weights = E(nullp20h)$weight)

VIC20[9] <- compare(of2c,of2c20h,method = "vi")
VICrandom20[9] <- compare(nullc, nullc20h, method = "vi")

rm(of2perturb20h,of2p20h,of2p20hm,nullperturb20h,nullp20hm,nullp20h,of2c20h,nullc20h)

#NINTH TRIAL
of2perturb20i <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20im <- as.matrix(of2perturb20i)
of2p20i <- graph_from_adjacency_matrix(of2p20im, mode = "undirected", weighted = TRUE)
nullperturb20i <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20im <- as.matrix(nullperturb20i)
nullp20i <- graph_from_adjacency_matrix(nullp20im, mode = "undirected", weighted = TRUE)

of2c20i <- cluster_louvain(of2p20i, weights = E(of2p20i)$weight)
nullc20i <- cluster_louvain(nullp20i, weights = E(nullp20i)$weight)

VIC20[10] <- compare(of2c,of2c20i,method = "vi")
VICrandom20[10] <- compare(nullc, nullc20i, method = "vi")

rm(of2perturb20i,of2p20i,of2p20im,nullperturb20i,nullp20im,nullp20i,of2c20i,nullc20i)

#TENTH TRIAL
of2perturb20j <- rewireR(of2adj, 863478, dist = "NegBinom")
of2p20jm <- as.matrix(of2perturb20j)
of2p20j <- graph_from_adjacency_matrix(of2p20jm, mode = "undirected", weighted = TRUE)
nullperturb20j <- rewireR(nullmatrix, 863478, dist = "NegBinom")
nullp20jm <- as.matrix(nullperturb20j)
nullp20j <- graph_from_adjacency_matrix(nullp20jm, mode = "undirected", weighted = TRUE)

of2c20j <- cluster_louvain(of2p20j, weights = E(of2p20j)$weight)
nullc20j <- cluster_louvain(nullp20j, weights = E(nullp20j)$weight)

VIC20[11] <- compare(of2c,of2c20j,method = "vi")
VICrandom20[11] <- compare(nullc, nullc20j, method = "vi")

rm(of2perturb20j,of2p20j,of2p20jm,nullperturb20j,nullp20jm,nullp20j,of2c20j,nullc20j)

#############25% Perturbation Level######################

#FIRST TRIAL
of2perturb25a <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25am <- as.matrix(of2perturb25a)
of2p25a <- graph_from_adjacency_matrix(of2p25am, mode = "undirected", weighted = TRUE)
nullperturb25a <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25am <- as.matrix(nullperturb25a)
nullp25a <- graph_from_adjacency_matrix(nullp25am, mode = "undirected", weighted = TRUE)

of2c25a <- cluster_louvain(of2p25a, weights = E(of2p25a)$weight)
nullc25a <- cluster_louvain(nullp25a, weights = E(nullp25a)$weight)

VIC25[2] <- compare(of2c,of2c25a,method = "vi")
VICrandom25[2] <- compare(nullc, nullc25a, method = "vi")

rm(of2perturb25a,of2p25a,of2p25am,nullperturb25a,nullp25am,nullp25a,of2c25a,nullc25a)

#SECOND TRIAL
of2perturb25b <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25bm <- as.matrix(of2perturb25b)
of2p25b <- graph_from_adjacency_matrix(of2p25bm, mode = "undirected", weighted = TRUE)
nullperturb25b <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25bm <- as.matrix(nullperturb25b)
nullp25b <- graph_from_adjacency_matrix(nullp25bm, mode = "undirected", weighted = TRUE)

of2c25b <- cluster_louvain(of2p25b, weights = E(of2p25b)$weight)
nullc25b <- cluster_louvain(nullp25b, weights = E(nullp25b)$weight)

VIC25[3] <- compare(of2c,of2c25b,method = "vi")
VICrandom25[3] <- compare(nullc, nullc25b, method = "vi")

rm(of2perturb25b,of2p25b,of2p25bm,nullperturb25b,nullp25bm,nullp25b,of2c25b,nullc25b)

#THIRD TRIAL
of2perturb25c <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25cm <- as.matrix(of2perturb25c)
of2p25c <- graph_from_adjacency_matrix(of2p25cm, mode = "undirected", weighted = TRUE)
nullperturb25c <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25cm <- as.matrix(nullperturb25c)
nullp25c <- graph_from_adjacency_matrix(nullp25cm, mode = "undirected", weighted = TRUE)

of2c25c <- cluster_louvain(of2p25c, weights = E(of2p25c)$weight)
nullc25c <- cluster_louvain(nullp25c, weights = E(nullp25c)$weight)

VIC25[4] <- compare(of2c,of2c25c,method = "vi")
VICrandom25[4] <- compare(nullc, nullc25c, method = "vi")

rm(of2perturb25c,of2p25c,of2p25cm,nullperturb25c,nullp25cm,nullp25c,of2c25c,nullc25c)

#FOURTH TRIAL
of2perturb25d <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25dm <- as.matrix(of2perturb25d)
of2p25d <- graph_from_adjacency_matrix(of2p25dm, mode = "undirected", weighted = TRUE)
nullperturb25d <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25dm <- as.matrix(nullperturb25d)
nullp25d <- graph_from_adjacency_matrix(nullp25dm, mode = "undirected", weighted = TRUE)

of2c25d <- cluster_louvain(of2p25d, weights = E(of2p25d)$weight)
nullc25d <- cluster_louvain(nullp25d, weights = E(nullp25d)$weight)

VIC25[5] <- compare(of2c,of2c25d,method = "vi")
VICrandom25[5] <- compare(nullc, nullc25d, method = "vi")

rm(of2perturb25d,of2p25d,of2p25dm,nullperturb25d,nullp25dm,nullp25d,of2c25d,nullc25d)

#FIFTH TRIAL
of2perturb25e <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25em <- as.matrix(of2perturb25e)
of2p25e <- graph_from_adjacency_matrix(of2p25em, mode = "undirected", weighted = TRUE)
nullperturb25e <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25em <- as.matrix(nullperturb25e)
nullp25e <- graph_from_adjacency_matrix(nullp25em, mode = "undirected", weighted = TRUE)

of2c25e <- cluster_louvain(of2p25e, weights = E(of2p25e)$weight)
nullc25e <- cluster_louvain(nullp25e, weights = E(nullp25e)$weight)

VIC25[6] <- compare(of2c,of2c25e,method = "vi")
VICrandom25[6] <- compare(nullc, nullc25e, method = "vi")

rm(of2perturb25e,of2p25e,of2p25em,nullperturb25e,nullp25em,nullp25e,of2c25e,nullc25e)

#SIXTH TRIAL
of2perturb25f <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25fm <- as.matrix(of2perturb25f)
of2p25f <- graph_from_adjacency_matrix(of2p25fm, mode = "undirected", weighted = TRUE)
nullperturb25f <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25fm <- as.matrix(nullperturb25f)
nullp25f <- graph_from_adjacency_matrix(nullp25fm, mode = "undirected", weighted = TRUE)

of2c25f <- cluster_louvain(of2p25f, weights = E(of2p25f)$weight)
nullc25f <- cluster_louvain(nullp25f, weights = E(nullp25f)$weight)

VIC25[7] <- compare(of2c,of2c25f,method = "vi")
VICrandom25[7] <- compare(nullc, nullc25f, method = "vi")

rm(of2perturb25f,of2p25f,of2p25fm,nullperturb25f,nullp25fm,nullp25f,of2c25f,nullc25f)

#SEVENTH TRIAL
of2perturb25g <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25gm <- as.matrix(of2perturb25g)
of2p25g <- graph_from_adjacency_matrix(of2p25gm, mode = "undirected", weighted = TRUE)
nullperturb25g <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25gm <- as.matrix(nullperturb25g)
nullp25g <- graph_from_adjacency_matrix(nullp25gm, mode = "undirected", weighted = TRUE)

of2c25g <- cluster_louvain(of2p25g, weights = E(of2p25g)$weight)
nullc25g <- cluster_louvain(nullp25g, weights = E(nullp25g)$weight)

VIC25[8] <- compare(of2c,of2c25g,method = "vi")
VICrandom25[8] <- compare(nullc, nullc25g, method = "vi")

rm(of2perturb25g,of2p25g,of2p25gm,nullperturb25g,nullp25gm,nullp25g,of2c25g,nullc25g)

#EIGHTH TRIAL
of2perturb25h <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25hm <- as.matrix(of2perturb25h)
of2p25h <- graph_from_adjacency_matrix(of2p25hm, mode = "undirected", weighted = TRUE)
nullperturb25h <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25hm <- as.matrix(nullperturb25h)
nullp25h <- graph_from_adjacency_matrix(nullp25hm, mode = "undirected", weighted = TRUE)

of2c25h <- cluster_louvain(of2p25h, weights = E(of2p25h)$weight)
nullc25h <- cluster_louvain(nullp25h, weights = E(nullp25h)$weight)

VIC25[9] <- compare(of2c,of2c25h,method = "vi")
VICrandom25[9] <- compare(nullc, nullc25h, method = "vi")

rm(of2perturb25h,of2p25h,of2p25hm,nullperturb25h,nullp25hm,nullp25h,of2c25h,nullc25h)

#NINTH TRIAL
of2perturb25i <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25im <- as.matrix(of2perturb25i)
of2p25i <- graph_from_adjacency_matrix(of2p25im, mode = "undirected", weighted = TRUE)
nullperturb25i <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25im <- as.matrix(nullperturb25i)
nullp25i <- graph_from_adjacency_matrix(nullp25im, mode = "undirected", weighted = TRUE)

of2c25i <- cluster_louvain(of2p25i, weights = E(of2p25i)$weight)
nullc25i <- cluster_louvain(nullp25i, weights = E(nullp25i)$weight)

VIC25[10] <- compare(of2c,of2c25i,method = "vi")
VICrandom25[10] <- compare(nullc, nullc25i, method = "vi")

rm(of2perturb25i,of2p25i,of2p25im,nullperturb25i,nullp25im,nullp25i,of2c25i,nullc25i)

#TENTH TRIAL
of2perturb25j <- rewireR(of2adj, 1079346, dist = "NegBinom")
of2p25jm <- as.matrix(of2perturb25j)
of2p25j <- graph_from_adjacency_matrix(of2p25jm, mode = "undirected", weighted = TRUE)
nullperturb25j <- rewireR(nullmatrix, 1079346, dist = "NegBinom")
nullp25jm <- as.matrix(nullperturb25j)
nullp25j <- graph_from_adjacency_matrix(nullp25jm, mode = "undirected", weighted = TRUE)

of2c25j <- cluster_louvain(of2p25j, weights = E(of2p25j)$weight)
nullc25j <- cluster_louvain(nullp25j, weights = E(nullp25j)$weight)

VIC25[11] <- compare(of2c,of2c25j,method = "vi")
VICrandom25[11] <- compare(nullc, nullc25j, method = "vi")

rm(of2perturb25j,of2p25j,of2p25jm,nullperturb25j,nullp25jm,nullp25j,of2c25j,nullc25j)

############30% Perturbation Level#######################

#FIRST TRIAL
of2perturb30a <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30am <- as.matrix(of2perturb30a)
of2p30a <- graph_from_adjacency_matrix(of2p30am, mode = "undirected", weighted = TRUE)
nullperturb30a <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30am <- as.matrix(nullperturb30a)
nullp30a <- graph_from_adjacency_matrix(nullp30am, mode = "undirected", weighted = TRUE)

of2c30a <- cluster_louvain(of2p30a, weights = E(of2p30a)$weight)
nullc30a <- cluster_louvain(nullp30a, weights = E(nullp30a)$weight)

VIC30[2] <- compare(of2c,of2c30a,method = "vi")
VICrandom30[2] <- compare(nullc, nullc30a, method = "vi")

rm(of2perturb30a,of2p30a,of2p30am,nullperturb30a,nullp30am,nullp30a,of2c30a,nullc30a)

#SECOND TRIAL
of2perturb30b <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30bm <- as.matrix(of2perturb30b)
of2p30b <- graph_from_adjacency_matrix(of2p30bm, mode = "undirected", weighted = TRUE)
nullperturb30b <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30bm <- as.matrix(nullperturb30b)
nullp30b <- graph_from_adjacency_matrix(nullp30bm, mode = "undirected", weighted = TRUE)

of2c30b <- cluster_louvain(of2p30b, weights = E(of2p30b)$weight)
nullc30b <- cluster_louvain(nullp30b, weights = E(nullp30b)$weight)

VIC30[3] <- compare(of2c,of2c30b,method = "vi")
VICrandom30[3] <- compare(nullc, nullc30b, method = "vi")

rm(of2perturb30b,of2p30b,of2p30bm,nullperturb30b,nullp30bm,nullp30b,of2c30b,nullc30b)

#THIRD TRIAL
of2perturb30c <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30cm <- as.matrix(of2perturb30c)
of2p30c <- graph_from_adjacency_matrix(of2p30cm, mode = "undirected", weighted = TRUE)
nullperturb30c <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30cm <- as.matrix(nullperturb30c)
nullp30c <- graph_from_adjacency_matrix(nullp30cm, mode = "undirected", weighted = TRUE)

of2c30c <- cluster_louvain(of2p30c, weights = E(of2p30c)$weight)
nullc30c <- cluster_louvain(nullp30c, weights = E(nullp30c)$weight)

VIC30[4] <- compare(of2c,of2c30c,method = "vi")
VICrandom30[4] <- compare(nullc, nullc30c, method = "vi")

rm(of2perturb30c,of2p30c,of2p30cm,nullperturb30c,nullp30cm,nullp30c,of2c30c,nullc30c)

#FOURTH TRIAL
of2perturb30d <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30dm <- as.matrix(of2perturb30d)
of2p30d <- graph_from_adjacency_matrix(of2p30dm, mode = "undirected", weighted = TRUE)
nullperturb30d <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30dm <- as.matrix(nullperturb30d)
nullp30d <- graph_from_adjacency_matrix(nullp30dm, mode = "undirected", weighted = TRUE)

of2c30d <- cluster_louvain(of2p30d, weights = E(of2p30d)$weight)
nullc30d <- cluster_louvain(nullp30d, weights = E(nullp30d)$weight)

VIC30[5] <- compare(of2c,of2c30d,method = "vi")
VICrandom30[5] <- compare(nullc, nullc30d, method = "vi")

rm(of2perturb30d,of2p30d,of2p30dm,nullperturb30d,nullp30dm,nullp30d,of2c30d,nullc30d)

#FIFTH TRIAL
of2perturb30e <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30em <- as.matrix(of2perturb30e)
of2p30e <- graph_from_adjacency_matrix(of2p30em, mode = "undirected", weighted = TRUE)
nullperturb30e <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30em <- as.matrix(nullperturb30e)
nullp30e <- graph_from_adjacency_matrix(nullp30em, mode = "undirected", weighted = TRUE)

of2c30e <- cluster_louvain(of2p30e, weights = E(of2p30e)$weight)
nullc30e <- cluster_louvain(nullp30e, weights = E(nullp30e)$weight)

VIC30[6] <- compare(of2c,of2c30e,method = "vi")
VICrandom30[6] <- compare(nullc, nullc30e, method = "vi")

rm(of2perturb30e,of2p30e,of2p30em,nullperturb30e,nullp30em,nullp30e,of2c30e,nullc30e)

#SIXTH TRIAL
of2perturb30f <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30fm <- as.matrix(of2perturb30f)
of2p30f <- graph_from_adjacency_matrix(of2p30fm, mode = "undirected", weighted = TRUE)
nullperturb30f <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30fm <- as.matrix(nullperturb30f)
nullp30f <- graph_from_adjacency_matrix(nullp30fm, mode = "undirected", weighted = TRUE)

of2c30f <- cluster_louvain(of2p30f, weights = E(of2p30f)$weight)
nullc30f <- cluster_louvain(nullp30f, weights = E(nullp30f)$weight)

VIC30[7] <- compare(of2c,of2c30f,method = "vi")
VICrandom30[7] <- compare(nullc, nullc30f, method = "vi")

rm(of2perturb30f,of2p30f,of2p30fm,nullperturb30f,nullp30fm,nullp30f,of2c30f,nullc30f)

#SEVENTH TRIAL
of2perturb30g <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30gm <- as.matrix(of2perturb30g)
of2p30g <- graph_from_adjacency_matrix(of2p30gm, mode = "undirected", weighted = TRUE)
nullperturb30g <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30gm <- as.matrix(nullperturb30g)
nullp30g <- graph_from_adjacency_matrix(nullp30gm, mode = "undirected", weighted = TRUE)

of2c30g <- cluster_louvain(of2p30g, weights = E(of2p30g)$weight)
nullc30g <- cluster_louvain(nullp30g, weights = E(nullp30g)$weight)

VIC30[8] <- compare(of2c,of2c30g,method = "vi")
VICrandom30[8] <- compare(nullc, nullc30g, method = "vi")

rm(of2perturb30g,of2p30g,of2p30gm,nullperturb30g,nullp30gm,nullp30g,of2c30g,nullc30g)

#EIGHTH TRIAL
of2perturb30h <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30hm <- as.matrix(of2perturb30h)
of2p30h <- graph_from_adjacency_matrix(of2p30hm, mode = "undirected", weighted = TRUE)
nullperturb30h <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30hm <- as.matrix(nullperturb30h)
nullp30h <- graph_from_adjacency_matrix(nullp30hm, mode = "undirected", weighted = TRUE)

of2c30h <- cluster_louvain(of2p30h, weights = E(of2p30h)$weight)
nullc30h <- cluster_louvain(nullp30h, weights = E(nullp30h)$weight)

VIC30[9] <- compare(of2c,of2c30h,method = "vi")
VICrandom30[9] <- compare(nullc, nullc30h, method = "vi")

rm(of2perturb30h,of2p30h,of2p30hm,nullperturb30h,nullp30hm,nullp30h,of2c30h,nullc30h)

#NINTH TRIAL
of2perturb30i <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30im <- as.matrix(of2perturb30i)
of2p30i <- graph_from_adjacency_matrix(of2p30im, mode = "undirected", weighted = TRUE)
nullperturb30i <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30im <- as.matrix(nullperturb30i)
nullp30i <- graph_from_adjacency_matrix(nullp30im, mode = "undirected", weighted = TRUE)

of2c30i <- cluster_louvain(of2p30i, weights = E(of2p30i)$weight)
nullc30i <- cluster_louvain(nullp30i, weights = E(nullp30i)$weight)

VIC30[10] <- compare(of2c,of2c30i,method = "vi")
VICrandom30[10] <- compare(nullc, nullc30i, method = "vi")

rm(of2perturb30i,of2p30i,of2p30im,nullperturb30i,nullp30im,nullp30i,of2c30i,nullc30i)

#TENTH TRIAL
of2perturb30j <- rewireR(of2adj, 1295217, dist = "NegBinom")
of2p30jm <- as.matrix(of2perturb30j)
of2p30j <- graph_from_adjacency_matrix(of2p30jm, mode = "undirected", weighted = TRUE)
nullperturb30j <- rewireR(nullmatrix, 1295217, dist = "NegBinom")
nullp30jm <- as.matrix(nullperturb30j)
nullp30j <- graph_from_adjacency_matrix(nullp30jm, mode = "undirected", weighted = TRUE)

of2c30j <- cluster_louvain(of2p30j, weights = E(of2p30j)$weight)
nullc30j <- cluster_louvain(nullp30j, weights = E(nullp30j)$weight)

VIC30[11] <- compare(of2c,of2c30j,method = "vi")
VICrandom30[11] <- compare(nullc, nullc30j, method = "vi")

rm(of2perturb30j,of2p30j,of2p30jm,nullperturb30j,nullp30jm,nullp30j,of2c30j,nullc30j)

##############35% Perturbation Level#####################

#FIRST TRIAL
of2perturb35a <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35am <- as.matrix(of2perturb35a)
of2p35a <- graph_from_adjacency_matrix(of2p35am, mode = "undirected", weighted = TRUE)
nullperturb35a <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35am <- as.matrix(nullperturb35a)
nullp35a <- graph_from_adjacency_matrix(nullp35am, mode = "undirected", weighted = TRUE)

of2c35a <- cluster_louvain(of2p35a, weights = E(of2p35a)$weight)
nullc35a <- cluster_louvain(nullp35a, weights = E(nullp35a)$weight)

VIC35[2] <- compare(of2c,of2c35a,method = "vi")
VICrandom35[2] <- compare(nullc, nullc35a, method = "vi")

rm(of2perturb35a,of2p35a,of2p35am,nullperturb35a,nullp35am,nullp35a,of2c35a,nullc35a)

#SECOND TRIAL
of2perturb35b <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35bm <- as.matrix(of2perturb35b)
of2p35b <- graph_from_adjacency_matrix(of2p35bm, mode = "undirected", weighted = TRUE)
nullperturb35b <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35bm <- as.matrix(nullperturb35b)
nullp35b <- graph_from_adjacency_matrix(nullp35bm, mode = "undirected", weighted = TRUE)

of2c35b <- cluster_louvain(of2p35b, weights = E(of2p35b)$weight)
nullc35b <- cluster_louvain(nullp35b, weights = E(nullp35b)$weight)

VIC35[3] <- compare(of2c,of2c35b,method = "vi")
VICrandom35[3] <- compare(nullc, nullc35b, method = "vi")

rm(of2perturb35b,of2p35b,of2p35bm,nullperturb35b,nullp35bm,nullp35b,of2c35b,nullc35b)

#THIRD TRIAL
of2perturb35c <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35cm <- as.matrix(of2perturb35c)
of2p35c <- graph_from_adjacency_matrix(of2p35cm, mode = "undirected", weighted = TRUE)
nullperturb35c <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35cm <- as.matrix(nullperturb35c)
nullp35c <- graph_from_adjacency_matrix(nullp35cm, mode = "undirected", weighted = TRUE)

of2c35c <- cluster_louvain(of2p35c, weights = E(of2p35c)$weight)
nullc35c <- cluster_louvain(nullp35c, weights = E(nullp35c)$weight)

VIC35[4] <- compare(of2c,of2c35c,method = "vi")
VICrandom35[4] <- compare(nullc, nullc35c, method = "vi")

rm(of2perturb35c,of2p35c,of2p35cm,nullperturb35c,nullp35cm,nullp35c,of2c35c,nullc35c)

#FOURTH TRIAL
of2perturb35d <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35dm <- as.matrix(of2perturb35d)
of2p35d <- graph_from_adjacency_matrix(of2p35dm, mode = "undirected", weighted = TRUE)
nullperturb35d <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35dm <- as.matrix(nullperturb35d)
nullp35d <- graph_from_adjacency_matrix(nullp35dm, mode = "undirected", weighted = TRUE)

of2c35d <- cluster_louvain(of2p35d, weights = E(of2p35d)$weight)
nullc35d <- cluster_louvain(nullp35d, weights = E(nullp35d)$weight)

VIC35[5] <- compare(of2c,of2c35d,method = "vi")
VICrandom35[5] <- compare(nullc, nullc35d, method = "vi")

rm(of2perturb35d,of2p35d,of2p35dm,nullperturb35d,nullp35dm,nullp35d,of2c35d,nullc35d)

#FIFTH TRIAL
of2perturb35e <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35em <- as.matrix(of2perturb35e)
of2p35e <- graph_from_adjacency_matrix(of2p35em, mode = "undirected", weighted = TRUE)
nullperturb35e <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35em <- as.matrix(nullperturb35e)
nullp35e <- graph_from_adjacency_matrix(nullp35em, mode = "undirected", weighted = TRUE)

of2c35e <- cluster_louvain(of2p35e, weights = E(of2p35e)$weight)
nullc35e <- cluster_louvain(nullp35e, weights = E(nullp35e)$weight)

VIC35[6] <- compare(of2c,of2c35e,method = "vi")
VICrandom35[6] <- compare(nullc, nullc35e, method = "vi")

rm(of2perturb35e,of2p35e,of2p35em,nullperturb35e,nullp35em,nullp35e,of2c35e,nullc35e)

#SIXTH TRIAL
of2perturb35f <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35fm <- as.matrix(of2perturb35f)
of2p35f <- graph_from_adjacency_matrix(of2p35fm, mode = "undirected", weighted = TRUE)
nullperturb35f <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35fm <- as.matrix(nullperturb35f)
nullp35f <- graph_from_adjacency_matrix(nullp35fm, mode = "undirected", weighted = TRUE)

of2c35f <- cluster_louvain(of2p35f, weights = E(of2p35f)$weight)
nullc35f <- cluster_louvain(nullp35f, weights = E(nullp35f)$weight)

VIC35[7] <- compare(of2c,of2c35f,method = "vi")
VICrandom35[7] <- compare(nullc, nullc35f, method = "vi")

rm(of2perturb35f,of2p35f,of2p35fm,nullperturb35f,nullp35fm,nullp35f,of2c35f,nullc35f)

#SEVENTH TRIAL
of2perturb35g <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35gm <- as.matrix(of2perturb35g)
of2p35g <- graph_from_adjacency_matrix(of2p35gm, mode = "undirected", weighted = TRUE)
nullperturb35g <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35gm <- as.matrix(nullperturb35g)
nullp35g <- graph_from_adjacency_matrix(nullp35gm, mode = "undirected", weighted = TRUE)

of2c35g <- cluster_louvain(of2p35g, weights = E(of2p35g)$weight)
nullc35g <- cluster_louvain(nullp35g, weights = E(nullp35g)$weight)

VIC35[8] <- compare(of2c,of2c35g,method = "vi")
VICrandom35[8] <- compare(nullc, nullc35g, method = "vi")

rm(of2perturb35g,of2p35g,of2p35gm,nullperturb35g,nullp35gm,nullp35g,of2c35g,nullc35g)

#EIGHTH TRIAL
of2perturb35h <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35hm <- as.matrix(of2perturb35h)
of2p35h <- graph_from_adjacency_matrix(of2p35hm, mode = "undirected", weighted = TRUE)
nullperturb35h <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35hm <- as.matrix(nullperturb35h)
nullp35h <- graph_from_adjacency_matrix(nullp35hm, mode = "undirected", weighted = TRUE)

of2c35h <- cluster_louvain(of2p35h, weights = E(of2p35h)$weight)
nullc35h <- cluster_louvain(nullp35h, weights = E(nullp35h)$weight)

VIC35[9] <- compare(of2c,of2c35h,method = "vi")
VICrandom35[9] <- compare(nullc, nullc35h, method = "vi")

rm(of2perturb35h,of2p35h,of2p35hm,nullperturb35h,nullp35hm,nullp35h,of2c35h,nullc35h)

#NINTH TRIAL
of2perturb35i <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35im <- as.matrix(of2perturb35i)
of2p35i <- graph_from_adjacency_matrix(of2p35im, mode = "undirected", weighted = TRUE)
nullperturb35i <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35im <- as.matrix(nullperturb35i)
nullp35i <- graph_from_adjacency_matrix(nullp35im, mode = "undirected", weighted = TRUE)

of2c35i <- cluster_louvain(of2p35i, weights = E(of2p35i)$weight)
nullc35i <- cluster_louvain(nullp35i, weights = E(nullp35i)$weight)

VIC35[10] <- compare(of2c,of2c35i,method = "vi")
VICrandom35[10] <- compare(nullc, nullc35i, method = "vi")

rm(of2perturb35i,of2p35i,of2p35im,nullperturb35i,nullp35im,nullp35i,of2c35i,nullc35i)

#TENTH TRIAL
of2perturb35j <- rewireR(of2adj, 1511083, dist = "NegBinom")
of2p35jm <- as.matrix(of2perturb35j)
of2p35j <- graph_from_adjacency_matrix(of2p35jm, mode = "undirected", weighted = TRUE)
nullperturb35j <- rewireR(nullmatrix, 1511083, dist = "NegBinom")
nullp35jm <- as.matrix(nullperturb35j)
nullp35j <- graph_from_adjacency_matrix(nullp35jm, mode = "undirected", weighted = TRUE)

of2c35j <- cluster_louvain(of2p35j, weights = E(of2p35j)$weight)
nullc35j <- cluster_louvain(nullp35j, weights = E(nullp35j)$weight)

VIC35[11] <- compare(of2c,of2c35j,method = "vi")
VICrandom35[11] <- compare(nullc, nullc35j, method = "vi")

rm(of2perturb35j,of2p35j,of2p35jm,nullperturb35j,nullp35jm,nullp35j,of2c35j,nullc35j)

#############40% Perturbation Level#####################

#FIRST TRIAL
of2perturb40a <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40am <- as.matrix(of2perturb40a)
of2p40a <- graph_from_adjacency_matrix(of2p40am, mode = "undirected", weighted = TRUE)
nullperturb40a <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40am <- as.matrix(nullperturb40a)
nullp40a <- graph_from_adjacency_matrix(nullp40am, mode = "undirected", weighted = TRUE)

of2c40a <- cluster_louvain(of2p40a, weights = E(of2p40a)$weight)
nullc40a <- cluster_louvain(nullp40a, weights = E(nullp40a)$weight)

VIC40[2] <- compare(of2c,of2c40a,method = "vi")
VICrandom40[2] <- compare(nullc, nullc40a, method = "vi")

rm(of2perturb40a,of2p40a,of2p40am,nullperturb40a,nullp40am,nullp40a,of2c40a,nullc40a)

#SECOND TRIAL
of2perturb40b <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40bm <- as.matrix(of2perturb40b)
of2p40b <- graph_from_adjacency_matrix(of2p40bm, mode = "undirected", weighted = TRUE)
nullperturb40b <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40bm <- as.matrix(nullperturb40b)
nullp40b <- graph_from_adjacency_matrix(nullp40bm, mode = "undirected", weighted = TRUE)

of2c40b <- cluster_louvain(of2p40b, weights = E(of2p40b)$weight)
nullc40b <- cluster_louvain(nullp40b, weights = E(nullp40b)$weight)

VIC40[3] <- compare(of2c,of2c40b,method = "vi")
VICrandom40[3] <- compare(nullc, nullc40b, method = "vi")

rm(of2perturb40b,of2p40b,of2p40bm,nullperturb40b,nullp40bm,nullp40b,of2c40b,nullc40b)

#THIRD TRIAL
of2perturb40c <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40cm <- as.matrix(of2perturb40c)
of2p40c <- graph_from_adjacency_matrix(of2p40cm, mode = "undirected", weighted = TRUE)
nullperturb40c <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40cm <- as.matrix(nullperturb40c)
nullp40c <- graph_from_adjacency_matrix(nullp40cm, mode = "undirected", weighted = TRUE)

of2c40c <- cluster_louvain(of2p40c, weights = E(of2p40c)$weight)
nullc40c <- cluster_louvain(nullp40c, weights = E(nullp40c)$weight)

VIC40[4] <- compare(of2c,of2c40c,method = "vi")
VICrandom40[4] <- compare(nullc, nullc40c, method = "vi")

rm(of2perturb40c,of2p40c,of2p40cm,nullperturb40c,nullp40cm,nullp40c,of2c40c,nullc40c)

#FOURTH TRIAL
of2perturb40d <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40dm <- as.matrix(of2perturb40d)
of2p40d <- graph_from_adjacency_matrix(of2p40dm, mode = "undirected", weighted = TRUE)
nullperturb40d <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40dm <- as.matrix(nullperturb40d)
nullp40d <- graph_from_adjacency_matrix(nullp40dm, mode = "undirected", weighted = TRUE)

of2c40d <- cluster_louvain(of2p40d, weights = E(of2p40d)$weight)
nullc40d <- cluster_louvain(nullp40d, weights = E(nullp40d)$weight)

VIC40[5] <- compare(of2c,of2c40d,method = "vi")
VICrandom40[5] <- compare(nullc, nullc40d, method = "vi")

rm(of2perturb40d,of2p40d,of2p40dm,nullperturb40d,nullp40dm,nullp40d,of2c40d,nullc40d)

#FIFTH TRIAL
of2perturb40e <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40em <- as.matrix(of2perturb40e)
of2p40e <- graph_from_adjacency_matrix(of2p40em, mode = "undirected", weighted = TRUE)
nullperturb40e <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40em <- as.matrix(nullperturb40e)
nullp40e <- graph_from_adjacency_matrix(nullp40em, mode = "undirected", weighted = TRUE)

of2c40e <- cluster_louvain(of2p40e, weights = E(of2p40e)$weight)
nullc40e <- cluster_louvain(nullp40e, weights = E(nullp40e)$weight)

VIC40[6] <- compare(of2c,of2c40e,method = "vi")
VICrandom40[6] <- compare(nullc, nullc40e, method = "vi")

rm(of2perturb40e,of2p40e,of2p40em,nullperturb40e,nullp40em,nullp40e,of2c40e,nullc40e)

#SIXTH TRIAL
of2perturb40f <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40fm <- as.matrix(of2perturb40f)
of2p40f <- graph_from_adjacency_matrix(of2p40fm, mode = "undirected", weighted = TRUE)
nullperturb40f <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40fm <- as.matrix(nullperturb40f)
nullp40f <- graph_from_adjacency_matrix(nullp40fm, mode = "undirected", weighted = TRUE)

of2c40f <- cluster_louvain(of2p40f, weights = E(of2p40f)$weight)
nullc40f <- cluster_louvain(nullp40f, weights = E(nullp40f)$weight)

VIC40[7] <- compare(of2c,of2c40f,method = "vi")
VICrandom40[7] <- compare(nullc, nullc40f, method = "vi")

rm(of2perturb40f,of2p40f,of2p40fm,nullperturb40f,nullp40fm,nullp40f,of2c40f,nullc40f)

#SEVENTH TRIAL
of2perturb40g <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40gm <- as.matrix(of2perturb40g)
of2p40g <- graph_from_adjacency_matrix(of2p40gm, mode = "undirected", weighted = TRUE)
nullperturb40g <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40gm <- as.matrix(nullperturb40g)
nullp40g <- graph_from_adjacency_matrix(nullp40gm, mode = "undirected", weighted = TRUE)

of2c40g <- cluster_louvain(of2p40g, weights = E(of2p40g)$weight)
nullc40g <- cluster_louvain(nullp40g, weights = E(nullp40g)$weight)

VIC40[8] <- compare(of2c,of2c40g,method = "vi")
VICrandom40[8] <- compare(nullc, nullc40g, method = "vi")

rm(of2perturb40g,of2p40g,of2p40gm,nullperturb40g,nullp40gm,nullp40g,of2c40g,nullc40g)

#EIGHTH TRIAL
of2perturb40h <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40hm <- as.matrix(of2perturb40h)
of2p40h <- graph_from_adjacency_matrix(of2p40hm, mode = "undirected", weighted = TRUE)
nullperturb40h <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40hm <- as.matrix(nullperturb40h)
nullp40h <- graph_from_adjacency_matrix(nullp40hm, mode = "undirected", weighted = TRUE)

of2c40h <- cluster_louvain(of2p40h, weights = E(of2p40h)$weight)
nullc40h <- cluster_louvain(nullp40h, weights = E(nullp40h)$weight)

VIC40[9] <- compare(of2c,of2c40h,method = "vi")
VICrandom40[9] <- compare(nullc, nullc40h, method = "vi")

rm(of2perturb40h,of2p40h,of2p40hm,nullperturb40h,nullp40hm,nullp40h,of2c40h,nullc40h)

#NINTH TRIAL
of2perturb40i <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40im <- as.matrix(of2perturb40i)
of2p40i <- graph_from_adjacency_matrix(of2p40im, mode = "undirected", weighted = TRUE)
nullperturb40i <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40im <- as.matrix(nullperturb40i)
nullp40i <- graph_from_adjacency_matrix(nullp40im, mode = "undirected", weighted = TRUE)

of2c40i <- cluster_louvain(of2p40i, weights = E(of2p40i)$weight)
nullc40i <- cluster_louvain(nullp40i, weights = E(nullp40i)$weight)

VIC40[10] <- compare(of2c,of2c40i,method = "vi")
VICrandom40[10] <- compare(nullc, nullc40i, method = "vi")

rm(of2perturb40i,of2p40i,of2p40im,nullperturb40i,nullp40im,nullp40i,of2c40i,nullc40i)

#TENTH TRIAL
of2perturb40j <- rewireR(of2adj, 1726956, dist = "NegBinom")
of2p40jm <- as.matrix(of2perturb40j)
of2p40j <- graph_from_adjacency_matrix(of2p40jm, mode = "undirected", weighted = TRUE)
nullperturb40j <- rewireR(nullmatrix, 1726956, dist = "NegBinom")
nullp40jm <- as.matrix(nullperturb40j)
nullp40j <- graph_from_adjacency_matrix(nullp40jm, mode = "undirected", weighted = TRUE)

of2c40j <- cluster_louvain(of2p40j, weights = E(of2p40j)$weight)
nullc40j <- cluster_louvain(nullp40j, weights = E(nullp40j)$weight)

VIC40[11] <- compare(of2c,of2c40j,method = "vi")
VICrandom40[11] <- compare(nullc, nullc40j, method = "vi")

rm(of2perturb40j,of2p40j,of2p40jm,nullperturb40j,nullp40jm,nullp40j,of2c40j,nullc40j)

############45% Perturbation Level######################

#FIRST TRIAL
of2perturb45a <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45am <- as.matrix(of2perturb45a)
of2p45a <- graph_from_adjacency_matrix(of2p45am, mode = "undirected", weighted = TRUE)
nullperturb45a <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45am <- as.matrix(nullperturb45a)
nullp45a <- graph_from_adjacency_matrix(nullp45am, mode = "undirected", weighted = TRUE)

of2c45a <- cluster_louvain(of2p45a, weights = E(of2p45a)$weight)
nullc45a <- cluster_louvain(nullp45a, weights = E(nullp45a)$weight)

VIC45[2] <- compare(of2c,of2c45a,method = "vi")
VICrandom45[2] <- compare(nullc, nullc45a, method = "vi")

rm(of2perturb45a,of2p45a,of2p45am,nullperturb45a,nullp45am,nullp45a,of2c45a,nullc45a)

#SECOND TRIAL
of2perturb45b <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45bm <- as.matrix(of2perturb45b)
of2p45b <- graph_from_adjacency_matrix(of2p45bm, mode = "undirected", weighted = TRUE)
nullperturb45b <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45bm <- as.matrix(nullperturb45b)
nullp45b <- graph_from_adjacency_matrix(nullp45bm, mode = "undirected", weighted = TRUE)

of2c45b <- cluster_louvain(of2p45b, weights = E(of2p45b)$weight)
nullc45b <- cluster_louvain(nullp45b, weights = E(nullp45b)$weight)

VIC45[3] <- compare(of2c,of2c45b,method = "vi")
VICrandom45[3] <- compare(nullc, nullc45b, method = "vi")

rm(of2perturb45b,of2p45b,of2p45bm,nullperturb45b,nullp45bm,nullp45b,of2c45b,nullc45b)

#THIRD TRIAL
of2perturb45c <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45cm <- as.matrix(of2perturb45c)
of2p45c <- graph_from_adjacency_matrix(of2p45cm, mode = "undirected", weighted = TRUE)
nullperturb45c <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45cm <- as.matrix(nullperturb45c)
nullp45c <- graph_from_adjacency_matrix(nullp45cm, mode = "undirected", weighted = TRUE)

of2c45c <- cluster_louvain(of2p45c, weights = E(of2p45c)$weight)
nullc45c <- cluster_louvain(nullp45c, weights = E(nullp45c)$weight)

VIC45[4] <- compare(of2c,of2c45c,method = "vi")
VICrandom45[4] <- compare(nullc, nullc45c, method = "vi")

rm(of2perturb45c,of2p45c,of2p45cm,nullperturb45c,nullp45cm,nullp45c,of2c45c,nullc45c)

#FOURTH TRIAL
of2perturb45d <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45dm <- as.matrix(of2perturb45d)
of2p45d <- graph_from_adjacency_matrix(of2p45dm, mode = "undirected", weighted = TRUE)
nullperturb45d <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45dm <- as.matrix(nullperturb45d)
nullp45d <- graph_from_adjacency_matrix(nullp45dm, mode = "undirected", weighted = TRUE)

of2c45d <- cluster_louvain(of2p45d, weights = E(of2p45d)$weight)
nullc45d <- cluster_louvain(nullp45d, weights = E(nullp45d)$weight)

VIC45[5] <- compare(of2c,of2c45d,method = "vi")
VICrandom45[5] <- compare(nullc, nullc45d, method = "vi")

rm(of2perturb45d,of2p45d,of2p45dm,nullperturb45d,nullp45dm,nullp45d,of2c45d,nullc45d)

#FIFTH TRIAL
of2perturb45e <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45em <- as.matrix(of2perturb45e)
of2p45e <- graph_from_adjacency_matrix(of2p45em, mode = "undirected", weighted = TRUE)
nullperturb45e <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45em <- as.matrix(nullperturb45e)
nullp45e <- graph_from_adjacency_matrix(nullp45em, mode = "undirected", weighted = TRUE)

of2c45e <- cluster_louvain(of2p45e, weights = E(of2p45e)$weight)
nullc45e <- cluster_louvain(nullp45e, weights = E(nullp45e)$weight)

VIC45[6] <- compare(of2c,of2c45e,method = "vi")
VICrandom45[6] <- compare(nullc, nullc45e, method = "vi")

rm(of2perturb45e,of2p45e,of2p45em,nullperturb45e,nullp45em,nullp45e,of2c45e,nullc45e)

#SIXTH TRIAL
of2perturb45f <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45fm <- as.matrix(of2perturb45f)
of2p45f <- graph_from_adjacency_matrix(of2p45fm, mode = "undirected", weighted = TRUE)
nullperturb45f <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45fm <- as.matrix(nullperturb45f)
nullp45f <- graph_from_adjacency_matrix(nullp45fm, mode = "undirected", weighted = TRUE)

of2c45f <- cluster_louvain(of2p45f, weights = E(of2p45f)$weight)
nullc45f <- cluster_louvain(nullp45f, weights = E(nullp45f)$weight)

VIC45[7] <- compare(of2c,of2c45f,method = "vi")
VICrandom45[7] <- compare(nullc, nullc45f, method = "vi")

rm(of2perturb45f,of2p45f,of2p45fm,nullperturb45f,nullp45fm,nullp45f,of2c45f,nullc45f)

#SEVENTH TRIAL
of2perturb45g <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45gm <- as.matrix(of2perturb45g)
of2p45g <- graph_from_adjacency_matrix(of2p45gm, mode = "undirected", weighted = TRUE)
nullperturb45g <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45gm <- as.matrix(nullperturb45g)
nullp45g <- graph_from_adjacency_matrix(nullp45gm, mode = "undirected", weighted = TRUE)

of2c45g <- cluster_louvain(of2p45g, weights = E(of2p45g)$weight)
nullc45g <- cluster_louvain(nullp45g, weights = E(nullp45g)$weight)

VIC45[8] <- compare(of2c,of2c45g,method = "vi")
VICrandom45[8] <- compare(nullc, nullc45g, method = "vi")

rm(of2perturb45g,of2p45g,of2p45gm,nullperturb45g,nullp45gm,nullp45g,of2c45g,nullc45g)

#EIGHTH TRIAL
of2perturb45h <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45hm <- as.matrix(of2perturb45h)
of2p45h <- graph_from_adjacency_matrix(of2p45hm, mode = "undirected", weighted = TRUE)
nullperturb45h <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45hm <- as.matrix(nullperturb45h)
nullp45h <- graph_from_adjacency_matrix(nullp45hm, mode = "undirected", weighted = TRUE)

of2c45h <- cluster_louvain(of2p45h, weights = E(of2p45h)$weight)
nullc45h <- cluster_louvain(nullp45h, weights = E(nullp45h)$weight)

VIC45[9] <- compare(of2c,of2c45h,method = "vi")
VICrandom45[9] <- compare(nullc, nullc45h, method = "vi")

rm(of2perturb45h,of2p45h,of2p45hm,nullperturb45h,nullp45hm,nullp45h,of2c45h,nullc45h)

#NINTH TRIAL
of2perturb45i <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45im <- as.matrix(of2perturb45i)
of2p45i <- graph_from_adjacency_matrix(of2p45im, mode = "undirected", weighted = TRUE)
nullperturb45i <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45im <- as.matrix(nullperturb45i)
nullp45i <- graph_from_adjacency_matrix(nullp45im, mode = "undirected", weighted = TRUE)

of2c45i <- cluster_louvain(of2p45i, weights = E(of2p45i)$weight)
nullc45i <- cluster_louvain(nullp45i, weights = E(nullp45i)$weight)

VIC45[10] <- compare(of2c,of2c45i,method = "vi")
VICrandom45[10] <- compare(nullc, nullc45i, method = "vi")

rm(of2perturb45i,of2p45i,of2p45im,nullperturb45i,nullp45im,nullp45i,of2c45i,nullc45i)

#TENTH TRIAL
of2perturb45j <- rewireR(of2adj, 1942821, dist = "NegBinom")
of2p45jm <- as.matrix(of2perturb45j)
of2p45j <- graph_from_adjacency_matrix(of2p45jm, mode = "undirected", weighted = TRUE)
nullperturb45j <- rewireR(nullmatrix, 1942821, dist = "NegBinom")
nullp45jm <- as.matrix(nullperturb45j)
nullp45j <- graph_from_adjacency_matrix(nullp45jm, mode = "undirected", weighted = TRUE)

of2c45j <- cluster_louvain(of2p45j, weights = E(of2p45j)$weight)
nullc45j <- cluster_louvain(nullp45j, weights = E(nullp45j)$weight)

VIC45[11] <- compare(of2c,of2c45j,method = "vi")
VICrandom45[11] <- compare(nullc, nullc45j, method = "vi")

rm(of2perturb45j,of2p45j,of2p45jm,nullperturb45j,nullp45jm,nullp45j,of2c45j,nullc45j)

##############50% Perturbation Level#####################

#FIRST TRIAL
of2perturb50a <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50am <- as.matrix(of2perturb50a)
of2p50a <- graph_from_adjacency_matrix(of2p50am, mode = "undirected", weighted = TRUE)
nullperturb50a <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50am <- as.matrix(nullperturb50a)
nullp50a <- graph_from_adjacency_matrix(nullp50am, mode = "undirected", weighted = TRUE)

of2c50a <- cluster_louvain(of2p50a, weights = E(of2p50a)$weight)
nullc50a <- cluster_louvain(nullp50a, weights = E(nullp50a)$weight)

VIC50[2] <- compare(of2c,of2c50a,method = "vi")
VICrandom50[2] <- compare(nullc, nullc50a, method = "vi")

rm(of2perturb50a,of2p50a,of2p50am,nullperturb50a,nullp50am,nullp50a,of2c50a,nullc50a)

#SECOND TRIAL
of2perturb50b <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50bm <- as.matrix(of2perturb50b)
of2p50b <- graph_from_adjacency_matrix(of2p50bm, mode = "undirected", weighted = TRUE)
nullperturb50b <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50bm <- as.matrix(nullperturb50b)
nullp50b <- graph_from_adjacency_matrix(nullp50bm, mode = "undirected", weighted = TRUE)

of2c50b <- cluster_louvain(of2p50b, weights = E(of2p50b)$weight)
nullc50b <- cluster_louvain(nullp50b, weights = E(nullp50b)$weight)

VIC50[3] <- compare(of2c,of2c50b,method = "vi")
VICrandom50[3] <- compare(nullc, nullc50b, method = "vi")

rm(of2perturb50b,of2p50b,of2p50bm,nullperturb50b,nullp50bm,nullp50b,of2c50b,nullc50b)

#THIRD TRIAL
of2perturb50c <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50cm <- as.matrix(of2perturb50c)
of2p50c <- graph_from_adjacency_matrix(of2p50cm, mode = "undirected", weighted = TRUE)
nullperturb50c <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50cm <- as.matrix(nullperturb50c)
nullp50c <- graph_from_adjacency_matrix(nullp50cm, mode = "undirected", weighted = TRUE)

of2c50c <- cluster_louvain(of2p50c, weights = E(of2p50c)$weight)
nullc50c <- cluster_louvain(nullp50c, weights = E(nullp50c)$weight)

VIC50[4] <- compare(of2c,of2c50c,method = "vi")
VICrandom50[4] <- compare(nullc, nullc50c, method = "vi")

rm(of2perturb50c,of2p50c,of2p50cm,nullperturb50c,nullp50cm,nullp50c,of2c50c,nullc50c)

#FOURTH TRIAL
of2perturb50d <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50dm <- as.matrix(of2perturb50d)
of2p50d <- graph_from_adjacency_matrix(of2p50dm, mode = "undirected", weighted = TRUE)
nullperturb50d <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50dm <- as.matrix(nullperturb50d)
nullp50d <- graph_from_adjacency_matrix(nullp50dm, mode = "undirected", weighted = TRUE)

of2c50d <- cluster_louvain(of2p50d, weights = E(of2p50d)$weight)
nullc50d <- cluster_louvain(nullp50d, weights = E(nullp50d)$weight)

VIC50[5] <- compare(of2c,of2c50d,method = "vi")
VICrandom50[5] <- compare(nullc, nullc50d, method = "vi")

rm(of2perturb50d,of2p50d,of2p50dm,nullperturb50d,nullp50dm,nullp50d,of2c50d,nullc50d)

#FIFTH TRIAL
of2perturb50e <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50em <- as.matrix(of2perturb50e)
of2p50e <- graph_from_adjacency_matrix(of2p50em, mode = "undirected", weighted = TRUE)
nullperturb50e <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50em <- as.matrix(nullperturb50e)
nullp50e <- graph_from_adjacency_matrix(nullp50em, mode = "undirected", weighted = TRUE)

of2c50e <- cluster_louvain(of2p50e, weights = E(of2p50e)$weight)
nullc50e <- cluster_louvain(nullp50e, weights = E(nullp50e)$weight)

VIC50[6] <- compare(of2c,of2c50e,method = "vi")
VICrandom50[6] <- compare(nullc, nullc50e, method = "vi")

rm(of2perturb50e,of2p50e,of2p50em,nullperturb50e,nullp50em,nullp50e,of2c50e,nullc50e)

#SIXTH TRIAL
of2perturb50f <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50fm <- as.matrix(of2perturb50f)
of2p50f <- graph_from_adjacency_matrix(of2p50fm, mode = "undirected", weighted = TRUE)
nullperturb50f <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50fm <- as.matrix(nullperturb50f)
nullp50f <- graph_from_adjacency_matrix(nullp50fm, mode = "undirected", weighted = TRUE)

of2c50f <- cluster_louvain(of2p50f, weights = E(of2p50f)$weight)
nullc50f <- cluster_louvain(nullp50f, weights = E(nullp50f)$weight)

VIC50[7] <- compare(of2c,of2c50f,method = "vi")
VICrandom50[7] <- compare(nullc, nullc50f, method = "vi")

rm(of2perturb50f,of2p50f,of2p50fm,nullperturb50f,nullp50fm,nullp50f,of2c50f,nullc50f)

#SEVENTH TRIAL
of2perturb50g <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50gm <- as.matrix(of2perturb50g)
of2p50g <- graph_from_adjacency_matrix(of2p50gm, mode = "undirected", weighted = TRUE)
nullperturb50g <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50gm <- as.matrix(nullperturb50g)
nullp50g <- graph_from_adjacency_matrix(nullp50gm, mode = "undirected", weighted = TRUE)

of2c50g <- cluster_louvain(of2p50g, weights = E(of2p50g)$weight)
nullc50g <- cluster_louvain(nullp50g, weights = E(nullp50g)$weight)

VIC50[8] <- compare(of2c,of2c50g,method = "vi")
VICrandom50[8] <- compare(nullc, nullc50g, method = "vi")

rm(of2perturb50g,of2p50g,of2p50gm,nullperturb50g,nullp50gm,nullp50g,of2c50g,nullc50g)

#EIGHTH TRIAL
of2perturb50h <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50hm <- as.matrix(of2perturb50h)
of2p50h <- graph_from_adjacency_matrix(of2p50hm, mode = "undirected", weighted = TRUE)
nullperturb50h <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50hm <- as.matrix(nullperturb50h)
nullp50h <- graph_from_adjacency_matrix(nullp50hm, mode = "undirected", weighted = TRUE)

of2c50h <- cluster_louvain(of2p50h, weights = E(of2p50h)$weight)
nullc50h <- cluster_louvain(nullp50h, weights = E(nullp50h)$weight)

VIC50[9] <- compare(of2c,of2c50h,method = "vi")
VICrandom50[9] <- compare(nullc, nullc50h, method = "vi")

rm(of2perturb50h,of2p50h,of2p50hm,nullperturb50h,nullp50hm,nullp50h,of2c50h,nullc50h)

#NINTH TRIAL
of2perturb50i <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50im <- as.matrix(of2perturb50i)
of2p50i <- graph_from_adjacency_matrix(of2p50im, mode = "undirected", weighted = TRUE)
nullperturb50i <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50im <- as.matrix(nullperturb50i)
nullp50i <- graph_from_adjacency_matrix(nullp50im, mode = "undirected", weighted = TRUE)

of2c50i <- cluster_louvain(of2p50i, weights = E(of2p50i)$weight)
nullc50i <- cluster_louvain(nullp50i, weights = E(nullp50i)$weight)

VIC50[10] <- compare(of2c,of2c50i,method = "vi")
VICrandom50[10] <- compare(nullc, nullc50i, method = "vi")

rm(of2perturb50i,of2p50i,of2p50im,nullperturb50i,nullp50im,nullp50i,of2c50i,nullc50i)

#TENTH TRIAL
of2perturb50j <- rewireR(of2adj, 2158695, dist = "NegBinom")
of2p50jm <- as.matrix(of2perturb50j)
of2p50j <- graph_from_adjacency_matrix(of2p50jm, mode = "undirected", weighted = TRUE)
nullperturb50j <- rewireR(nullmatrix, 2158695, dist = "NegBinom")
nullp50jm <- as.matrix(nullperturb50j)
nullp50j <- graph_from_adjacency_matrix(nullp50jm, mode = "undirected", weighted = TRUE)

of2c50j <- cluster_louvain(of2p50j, weights = E(of2p50j)$weight)
nullc50j <- cluster_louvain(nullp50j, weights = E(nullp50j)$weight)

VIC50[11] <- compare(of2c,of2c50j,method = "vi")
VICrandom50[11] <- compare(nullc, nullc50j, method = "vi")

rm(of2perturb50j,of2p50j,of2p50jm,nullperturb50j,nullp50jm,nullp50j,of2c50j,nullc50j)

###############55% Perturbation Level################

#FIRST TRIAL
of2perturb55a <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55am <- as.matrix(of2perturb55a)
of2p55a <- graph_from_adjacency_matrix(of2p55am, mode = "undirected", weighted = TRUE)
nullperturb55a <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55am <- as.matrix(nullperturb55a)
nullp55a <- graph_from_adjacency_matrix(nullp55am, mode = "undirected", weighted = TRUE)

of2c55a <- cluster_louvain(of2p55a, weights = E(of2p55a)$weight)
nullc55a <- cluster_louvain(nullp55a, weights = E(nullp55a)$weight)

VIC55[2] <- compare(of2c,of2c55a,method = "vi")
VICrandom55[2] <- compare(nullc, nullc55a, method = "vi")

rm(of2perturb55a,of2p55a,of2p55am,nullperturb55a,nullp55am,nullp55a,of2c55a,nullc55a)

#SECOND TRIAL
of2perturb55b <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55bm <- as.matrix(of2perturb55b)
of2p55b <- graph_from_adjacency_matrix(of2p55bm, mode = "undirected", weighted = TRUE)
nullperturb55b <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55bm <- as.matrix(nullperturb55b)
nullp55b <- graph_from_adjacency_matrix(nullp55bm, mode = "undirected", weighted = TRUE)

of2c55b <- cluster_louvain(of2p55b, weights = E(of2p55b)$weight)
nullc55b <- cluster_louvain(nullp55b, weights = E(nullp55b)$weight)

VIC55[3] <- compare(of2c,of2c55b,method = "vi")
VICrandom55[3] <- compare(nullc, nullc55b, method = "vi")

rm(of2perturb55b,of2p55b,of2p55bm,nullperturb55b,nullp55bm,nullp55b,of2c55b,nullc55b)

#THIRD TRIAL
of2perturb55c <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55cm <- as.matrix(of2perturb55c)
of2p55c <- graph_from_adjacency_matrix(of2p55cm, mode = "undirected", weighted = TRUE)
nullperturb55c <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55cm <- as.matrix(nullperturb55c)
nullp55c <- graph_from_adjacency_matrix(nullp55cm, mode = "undirected", weighted = TRUE)

of2c55c <- cluster_louvain(of2p55c, weights = E(of2p55c)$weight)
nullc55c <- cluster_louvain(nullp55c, weights = E(nullp55c)$weight)

VIC55[4] <- compare(of2c,of2c55c,method = "vi")
VICrandom55[4] <- compare(nullc, nullc55c, method = "vi")

rm(of2perturb55c,of2p55c,of2p55cm,nullperturb55c,nullp55cm,nullp55c,of2c55c,nullc55c)

#FOURTH TRIAL
of2perturb55d <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55dm <- as.matrix(of2perturb55d)
of2p55d <- graph_from_adjacency_matrix(of2p55dm, mode = "undirected", weighted = TRUE)
nullperturb55d <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55dm <- as.matrix(nullperturb55d)
nullp55d <- graph_from_adjacency_matrix(nullp55dm, mode = "undirected", weighted = TRUE)

of2c55d <- cluster_louvain(of2p55d, weights = E(of2p55d)$weight)
nullc55d <- cluster_louvain(nullp55d, weights = E(nullp55d)$weight)

VIC55[5] <- compare(of2c,of2c55d,method = "vi")
VICrandom55[5] <- compare(nullc, nullc55d, method = "vi")

rm(of2perturb55d,of2p55d,of2p55dm,nullperturb55d,nullp55dm,nullp55d,of2c55d,nullc55d)

#FIFTH TRIAL
of2perturb55e <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55em <- as.matrix(of2perturb55e)
of2p55e <- graph_from_adjacency_matrix(of2p55em, mode = "undirected", weighted = TRUE)
nullperturb55e <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55em <- as.matrix(nullperturb55e)
nullp55e <- graph_from_adjacency_matrix(nullp55em, mode = "undirected", weighted = TRUE)

of2c55e <- cluster_louvain(of2p55e, weights = E(of2p55e)$weight)
nullc55e <- cluster_louvain(nullp55e, weights = E(nullp55e)$weight)

VIC55[6] <- compare(of2c,of2c55e,method = "vi")
VICrandom55[6] <- compare(nullc, nullc55e, method = "vi")

rm(of2perturb55e,of2p55e,of2p55em,nullperturb55e,nullp55em,nullp55e,of2c55e,nullc55e)

#SIXTH TRIAL
of2perturb55f <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55fm <- as.matrix(of2perturb55f)
of2p55f <- graph_from_adjacency_matrix(of2p55fm, mode = "undirected", weighted = TRUE)
nullperturb55f <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55fm <- as.matrix(nullperturb55f)
nullp55f <- graph_from_adjacency_matrix(nullp55fm, mode = "undirected", weighted = TRUE)

of2c55f <- cluster_louvain(of2p55f, weights = E(of2p55f)$weight)
nullc55f <- cluster_louvain(nullp55f, weights = E(nullp55f)$weight)

VIC55[7] <- compare(of2c,of2c55f,method = "vi")
VICrandom55[7] <- compare(nullc, nullc55f, method = "vi")

rm(of2perturb55f,of2p55f,of2p55fm,nullperturb55f,nullp55fm,nullp55f,of2c55f,nullc55f)

#SEVENTH TRIAL
of2perturb55g <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55gm <- as.matrix(of2perturb55g)
of2p55g <- graph_from_adjacency_matrix(of2p55gm, mode = "undirected", weighted = TRUE)
nullperturb55g <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55gm <- as.matrix(nullperturb55g)
nullp55g <- graph_from_adjacency_matrix(nullp55gm, mode = "undirected", weighted = TRUE)

of2c55g <- cluster_louvain(of2p55g, weights = E(of2p55g)$weight)
nullc55g <- cluster_louvain(nullp55g, weights = E(nullp55g)$weight)

VIC55[8] <- compare(of2c,of2c55g,method = "vi")
VICrandom55[8] <- compare(nullc, nullc55g, method = "vi")

rm(of2perturb55g,of2p55g,of2p55gm,nullperturb55g,nullp55gm,nullp55g,of2c55g,nullc55g)

#EIGHTH TRIAL
of2perturb55h <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55hm <- as.matrix(of2perturb55h)
of2p55h <- graph_from_adjacency_matrix(of2p55hm, mode = "undirected", weighted = TRUE)
nullperturb55h <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55hm <- as.matrix(nullperturb55h)
nullp55h <- graph_from_adjacency_matrix(nullp55hm, mode = "undirected", weighted = TRUE)

of2c55h <- cluster_louvain(of2p55h, weights = E(of2p55h)$weight)
nullc55h <- cluster_louvain(nullp55h, weights = E(nullp55h)$weight)

VIC55[9] <- compare(of2c,of2c55h,method = "vi")
VICrandom55[9] <- compare(nullc, nullc55h, method = "vi")

rm(of2perturb55h,of2p55h,of2p55hm,nullperturb55h,nullp55hm,nullp55h,of2c55h,nullc55h)

#NINTH TRIAL
of2perturb55i <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55im <- as.matrix(of2perturb55i)
of2p55i <- graph_from_adjacency_matrix(of2p55im, mode = "undirected", weighted = TRUE)
nullperturb55i <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55im <- as.matrix(nullperturb55i)
nullp55i <- graph_from_adjacency_matrix(nullp55im, mode = "undirected", weighted = TRUE)

of2c55i <- cluster_louvain(of2p55i, weights = E(of2p55i)$weight)
nullc55i <- cluster_louvain(nullp55i, weights = E(nullp55i)$weight)

VIC55[10] <- compare(of2c,of2c55i,method = "vi")
VICrandom55[10] <- compare(nullc, nullc55i, method = "vi")

rm(of2perturb55i,of2p55i,of2p55im,nullperturb55i,nullp55im,nullp55i,of2c55i,nullc55i)

#TENTH TRIAL
of2perturb55j <- rewireR(of2adj, 2374559, dist = "NegBinom")
of2p55jm <- as.matrix(of2perturb55j)
of2p55j <- graph_from_adjacency_matrix(of2p55jm, mode = "undirected", weighted = TRUE)
nullperturb55j <- rewireR(nullmatrix, 2374559, dist = "NegBinom")
nullp55jm <- as.matrix(nullperturb55j)
nullp55j <- graph_from_adjacency_matrix(nullp55jm, mode = "undirected", weighted = TRUE)

of2c55j <- cluster_louvain(of2p55j, weights = E(of2p55j)$weight)
nullc55j <- cluster_louvain(nullp55j, weights = E(nullp55j)$weight)

VIC55[11] <- compare(of2c,of2c55j,method = "vi")
VICrandom55[11] <- compare(nullc, nullc55j, method = "vi")

rm(of2perturb55j,of2p55j,of2p55jm,nullperturb55j,nullp55jm,nullp55j,of2c55j,nullc55j)

##############60% Perturbation Level#######################

#FIRST TRIAL
of2perturb60a <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60am <- as.matrix(of2perturb60a)
of2p60a <- graph_from_adjacency_matrix(of2p60am, mode = "undirected", weighted = TRUE)
nullperturb60a <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60am <- as.matrix(nullperturb60a)
nullp60a <- graph_from_adjacency_matrix(nullp60am, mode = "undirected", weighted = TRUE)

of2c60a <- cluster_louvain(of2p60a, weights = E(of2p60a)$weight)
nullc60a <- cluster_louvain(nullp60a, weights = E(nullp60a)$weight)

VIC60[2] <- compare(of2c,of2c60a,method = "vi")
VICrandom60[2] <- compare(nullc, nullc60a, method = "vi")

rm(of2perturb60a,of2p60a,of2p60am,nullperturb60a,nullp60am,nullp60a,of2c60a,nullc60a)

#SECOND TRIAL
of2perturb60b <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60bm <- as.matrix(of2perturb60b)
of2p60b <- graph_from_adjacency_matrix(of2p60bm, mode = "undirected", weighted = TRUE)
nullperturb60b <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60bm <- as.matrix(nullperturb60b)
nullp60b <- graph_from_adjacency_matrix(nullp60bm, mode = "undirected", weighted = TRUE)

of2c60b <- cluster_louvain(of2p60b, weights = E(of2p60b)$weight)
nullc60b <- cluster_louvain(nullp60b, weights = E(nullp60b)$weight)

VIC60[3] <- compare(of2c,of2c60b,method = "vi")
VICrandom60[3] <- compare(nullc, nullc60b, method = "vi")

rm(of2perturb60b,of2p60b,of2p60bm,nullperturb60b,nullp60bm,nullp60b,of2c60b,nullc60b)

#THIRD TRIAL
of2perturb60c <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60cm <- as.matrix(of2perturb60c)
of2p60c <- graph_from_adjacency_matrix(of2p60cm, mode = "undirected", weighted = TRUE)
nullperturb60c <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60cm <- as.matrix(nullperturb60c)
nullp60c <- graph_from_adjacency_matrix(nullp60cm, mode = "undirected", weighted = TRUE)

of2c60c <- cluster_louvain(of2p60c, weights = E(of2p60c)$weight)
nullc60c <- cluster_louvain(nullp60c, weights = E(nullp60c)$weight)

VIC60[4] <- compare(of2c,of2c60c,method = "vi")
VICrandom60[4] <- compare(nullc, nullc60c, method = "vi")

rm(of2perturb60c,of2p60c,of2p60cm,nullperturb60c,nullp60cm,nullp60c,of2c60c,nullc60c)

#FOURTH TRIAL
of2perturb60d <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60dm <- as.matrix(of2perturb60d)
of2p60d <- graph_from_adjacency_matrix(of2p60dm, mode = "undirected", weighted = TRUE)
nullperturb60d <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60dm <- as.matrix(nullperturb60d)
nullp60d <- graph_from_adjacency_matrix(nullp60dm, mode = "undirected", weighted = TRUE)

of2c60d <- cluster_louvain(of2p60d, weights = E(of2p60d)$weight)
nullc60d <- cluster_louvain(nullp60d, weights = E(nullp60d)$weight)

VIC60[5] <- compare(of2c,of2c60d,method = "vi")
VICrandom60[5] <- compare(nullc, nullc60d, method = "vi")

rm(of2perturb60d,of2p60d,of2p60dm,nullperturb60d,nullp60dm,nullp60d,of2c60d,nullc60d)

#FIFTH TRIAL
of2perturb60e <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60em <- as.matrix(of2perturb60e)
of2p60e <- graph_from_adjacency_matrix(of2p60em, mode = "undirected", weighted = TRUE)
nullperturb60e <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60em <- as.matrix(nullperturb60e)
nullp60e <- graph_from_adjacency_matrix(nullp60em, mode = "undirected", weighted = TRUE)

of2c60e <- cluster_louvain(of2p60e, weights = E(of2p60e)$weight)
nullc60e <- cluster_louvain(nullp60e, weights = E(nullp60e)$weight)

VIC60[6] <- compare(of2c,of2c60e,method = "vi")
VICrandom60[6] <- compare(nullc, nullc60e, method = "vi")

rm(of2perturb60e,of2p60e,of2p60em,nullperturb60e,nullp60em,nullp60e,of2c60e,nullc60e)

#SIXTH TRIAL
of2perturb60f <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60fm <- as.matrix(of2perturb60f)
of2p60f <- graph_from_adjacency_matrix(of2p60fm, mode = "undirected", weighted = TRUE)
nullperturb60f <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60fm <- as.matrix(nullperturb60f)
nullp60f <- graph_from_adjacency_matrix(nullp60fm, mode = "undirected", weighted = TRUE)

of2c60f <- cluster_louvain(of2p60f, weights = E(of2p60f)$weight)
nullc60f <- cluster_louvain(nullp60f, weights = E(nullp60f)$weight)

VIC60[7] <- compare(of2c,of2c60f,method = "vi")
VICrandom60[7] <- compare(nullc, nullc60f, method = "vi")

rm(of2perturb60f,of2p60f,of2p60fm,nullperturb60f,nullp60fm,nullp60f,of2c60f,nullc60f)

#SEVENTH TRIAL
of2perturb60g <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60gm <- as.matrix(of2perturb60g)
of2p60g <- graph_from_adjacency_matrix(of2p60gm, mode = "undirected", weighted = TRUE)
nullperturb60g <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60gm <- as.matrix(nullperturb60g)
nullp60g <- graph_from_adjacency_matrix(nullp60gm, mode = "undirected", weighted = TRUE)

of2c60g <- cluster_louvain(of2p60g, weights = E(of2p60g)$weight)
nullc60g <- cluster_louvain(nullp60g, weights = E(nullp60g)$weight)

VIC60[8] <- compare(of2c,of2c60g,method = "vi")
VICrandom60[8] <- compare(nullc, nullc60g, method = "vi")

rm(of2perturb60g,of2p60g,of2p60gm,nullperturb60g,nullp60gm,nullp60g,of2c60g,nullc60g)

#EIGHTH TRIAL
of2perturb60h <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60hm <- as.matrix(of2perturb60h)
of2p60h <- graph_from_adjacency_matrix(of2p60hm, mode = "undirected", weighted = TRUE)
nullperturb60h <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60hm <- as.matrix(nullperturb60h)
nullp60h <- graph_from_adjacency_matrix(nullp60hm, mode = "undirected", weighted = TRUE)

of2c60h <- cluster_louvain(of2p60h, weights = E(of2p60h)$weight)
nullc60h <- cluster_louvain(nullp60h, weights = E(nullp60h)$weight)

VIC60[9] <- compare(of2c,of2c60h,method = "vi")
VICrandom60[9] <- compare(nullc, nullc60h, method = "vi")

rm(of2perturb60h,of2p60h,of2p60hm,nullperturb60h,nullp60hm,nullp60h,of2c60h,nullc60h)

#NINTH TRIAL
of2perturb60i <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60im <- as.matrix(of2perturb60i)
of2p60i <- graph_from_adjacency_matrix(of2p60im, mode = "undirected", weighted = TRUE)
nullperturb60i <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60im <- as.matrix(nullperturb60i)
nullp60i <- graph_from_adjacency_matrix(nullp60im, mode = "undirected", weighted = TRUE)

of2c60i <- cluster_louvain(of2p60i, weights = E(of2p60i)$weight)
nullc60i <- cluster_louvain(nullp60i, weights = E(nullp60i)$weight)

VIC60[10] <- compare(of2c,of2c60i,method = "vi")
VICrandom60[10] <- compare(nullc, nullc60i, method = "vi")

rm(of2perturb60i,of2p60i,of2p60im,nullperturb60i,nullp60im,nullp60i,of2c60i,nullc60i)

#TENTH TRIAL
of2perturb60j <- rewireR(of2adj, 2590434, dist = "NegBinom")
of2p60jm <- as.matrix(of2perturb60j)
of2p60j <- graph_from_adjacency_matrix(of2p60jm, mode = "undirected", weighted = TRUE)
nullperturb60j <- rewireR(nullmatrix, 2590434, dist = "NegBinom")
nullp60jm <- as.matrix(nullperturb60j)
nullp60j <- graph_from_adjacency_matrix(nullp60jm, mode = "undirected", weighted = TRUE)

of2c60j <- cluster_louvain(of2p60j, weights = E(of2p60j)$weight)
nullc60j <- cluster_louvain(nullp60j, weights = E(nullp60j)$weight)

VIC60[11] <- compare(of2c,of2c60j,method = "vi")
VICrandom60[11] <- compare(nullc, nullc60j, method = "vi")

rm(of2perturb60j,of2p60j,of2p60jm,nullperturb60j,nullp60jm,nullp60j,of2c60j,nullc60j)

###############PLOTTING VI CURVES AND TESTING######################

#Adding the perturbation level as the first entry of each vector
VIC0[1] <- 0
VIC5[1] <- 5
VIC10[1] <- 10
VIC15[1] <- 15
VIC20[1] <- 20
VIC25[1] <- 25
VIC30[1] <- 30
VIC35[1] <- 35
VIC40[1] <- 40
VIC45[1] <- 45
VIC50[1] <- 50
VIC55[1] <- 55
VIC60[1] <- 60

VICrandom0[1] <- 0
VICrandom5[1] <- 5
VICrandom10[1] <- 10
VICrandom15[1] <- 15
VICrandom20[1] <- 20
VICrandom25[1] <- 25
VICrandom30[1] <- 30
VICrandom35[1] <- 35
VICrandom40[1] <- 40
VICrandom45[1] <- 45
VICrandom50[1] <- 50
VICrandom55[1] <- 55
VICrandom60[1] <- 60

#Creating matrices containing all the VI values
VIC <- cbind(VIC0,VIC5,VIC10,VIC15,VIC20,VIC25,VIC30,VIC35,VIC40,VIC45,VIC50,VIC55,VIC60)
VICrandom <- cbind(VICrandom0,VICrandom5,VICrandom10,VICrandom15,VICrandom20,VICrandom25,VICrandom30,VICrandom35,VICrandom40,VICrandom45,VICrandom50,VICrandom55,VICrandom60)

#Computing the area under the VI curves for the of2 network and the random model
robinAUC(of2, VIC[2:11,1:13], VICrandom[2:11,1:13], measure = "vi", verbose = FALSE)
#Calculating the Bayes Factor using GP regression
robinGPTest(VIC[2:11,1:13], VICrandom[2:11,1:13])
#Performing the ITP test
robinFDATest(of2, VIC[2:11,1:13], VICrandom[2:11,1:13], measure = "vi", legend = c("VIC", "VICrandom"), verbose = FALSE)

#Setting up the variables to plot the VI curves
VICplot <- c(colMeans(VIC[2:11,1:13]))
VICrandomplot <- c(colMeans(VICrandom[2:11,1:13]))
values <- c(VICplot, VICrandomplot)
perturb <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6)
group <- c("VIC","VIC","VIC","VIC","VIC","VIC","VIC","VIC","VIC","VIC","VIC","VIC","VIC","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom","VICrandom")
curves <- data.frame(group,perturb,values)

#Plotting the VI curves using ggplot2
ggplot(curves, aes(x = perturb, y= values, colour = group, fill = group)) +
  geom_line() +
  geom_point(size = 2, shape = 21) +
  ggtitle("VI curves for network of2") +
  xlab("Perturbation level") +
  ylab("VI") +
  labs(colour = "Model", fill = "Model")

########Testing the robustness of the compressed network######

##############Creating super node network####################

#Importing the functions from the SuperNode.R code from Stanley (2018)
source('SuperNode.R')
#Compressing of2 into a super node representation with 500 super nodes
Out=SuperNode(of2,500)
#Storing the node to super node conversion key
snodes<- Out$SNAssn

#Finding the seeds of the super node network and sorting them into numerical order
seeds = SeedsFromCore(of2, 500)
seedsort = sort(seeds, decreasing = FALSE)
extras = setdiff(Out$SNAssn,seeds)

#Cluster the super node network and store the list of community memberships
snlv <- cluster_louvain(Out$SNNet)
memb <- snlv$membership

#Find which community in the clustering on the super node network each node in the original network belongs to
clust = c()
for(i in c(1:2939)) {
  sn = Out$SNAssn[i]
  if(! sn %in% seeds){
    clust[i] = i
  } else{
    clust[i] = memb[match(sn,seedsort)]
  }
}

#Calculating the VI between the clustering on the original network and the clustering on the super node network
compare(clustlv,clust,method = "vi")

#Plotting the super node network
plot.igraph(Out$SNNet, vertex.size = 7, vertex.label = NA)

#Plotting the original network colour coded by super node membership
#Plotting the super node network colour coded by clustering membership
my_collsn <- as.numeric(as.factor(snlv$membership))
my_colSuper <- as.numeric(as.factor(Out$SNAssn))
plot.igraph(Out$SNNet, vertex.color = my_collsn, vertex.size = 7, vertex.label = NA)
plot.igraph(of2, vertex.color = my_collsn, vertex.size = 10, vertex.label = NA)

##################Testing the robustness#######################

#Finding the adjacency matrix of the super node network
SNadj <- as_adjacency_matrix(Out$SNNet, type = "both", attr = "weight", edges = FALSE, names = FALSE, sparse = FALSE)
#Finding the null random model corresponding to the super node network
SNrewire <- rewireR(SNadj,124750,dist="NegBinom")
#Creating a graph of the null random model
SNnullmatrix <- as.matrix(SNrewire)
SNnull <- graph_from_adjacency_matrix(SNnullmatrix, mode = "undirected", weighted = TRUE)

#Storing the edge weights of the super node network and null random model
SNw<- E(Out$SNNet)$weight
SNnullw<- E(SNnull)$weight

#Clustering the super node network and null random model
SNc <- cluster_louvain(Out$SNNet, weights = SNw)
SNnullc <- cluster_louvain(SNnull, weights = SNnullw)

#Creating vectors to store the VI between the original and perturbed networks
VICSN0 <- c(0,0,0,0,0,0,0,0,0,0,0)
VICSN5 <- c()
VICSN10 <- c()
VICSN15 <- c()
VICSN20 <- c()
VICSN25 <- c()
VICSN30 <- c()
VICSN35 <- c()
VICSN40 <- c()
VICSN45 <- c()
VICSN50 <- c()
VICSN55 <- c()
VICSN60 <- c()

VICSNrandom0 <- c(0,0,0,0,0,0,0,0,0,0,0)
VICSNrandom5 <- c()
VICSNrandom10 <- c()
VICSNrandom15 <- c()
VICSNrandom20 <- c()
VICSNrandom25 <- c()
VICSNrandom30 <- c()
VICSNrandom35 <- c()
VICSNrandom40 <- c()
VICSNrandom45 <- c()
VICSNrandom50 <- c()
VICSNrandom55 <- c()
VICSNrandom60 <- c()

#################5% Perturbation level#####################

#First trial at the 5% perturbation level

#Perturbing 5% of the edges in the super node network and the null model 
SNperturb5a <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5am <- as.matrix(SNperturb5a)
SNp5a <- graph_from_adjacency_matrix(SNp5am, mode = "undirected", weighted = TRUE)
SNnullperturb5a <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5am <- as.matrix(SNnullperturb5a)
SNnullp5a <- graph_from_adjacency_matrix(SNnullp5am, mode = "undirected", weighted = TRUE)

#Clustering the perturbed networks
SNc5a <- cluster_louvain(SNp5a, weights = E(SNp5a)$weight)
SNnullc5a <- cluster_louvain(SNnullp5a, weights = E(SNnullp5a)$weight)

#Calculating the VI between the original networks and the perturbed networks
VICSN5[2] <- compare(SNc,SNc5a,method = "vi")
VICSNrandom5[2] <- compare(SNnullc, SNnullc5a, method = "vi")

#Removing the variables used for this section of the analysis from the global environment to save space
rm(SNperturb5a,SNp5a,SNp5am,SNnullperturb5a,SNnullp5am,SNnullp5a,SNc5a,SNnullc5a)

#I now repeat the same process ten times at each perturbation level

#SECOND TRIAL
SNperturb5b <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5bm <- as.matrix(SNperturb5b)
SNp5b <- graph_from_adjacency_matrix(SNp5bm, mode = "undirected", weighted = TRUE)
SNnullperturb5b <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5bm <- as.matrix(SNnullperturb5b)
SNnullp5b <- graph_from_adjacency_matrix(SNnullp5bm, mode = "undirected", weighted = TRUE)

SNc5b <- cluster_louvain(SNp5b, weights = E(SNp5b)$weight)
SNnullc5b <- cluster_louvain(SNnullp5b, weights = E(SNnullp5b)$weight)

VICSN5[3] <- compare(SNc,SNc5b,method = "vi")
VICSNrandom5[3] <- compare(SNnullc, SNnullc5b, method = "vi")

rm(SNperturb5b,SNp5b,SNp5bm,SNnullperturb5b,SNnullp5bm,SNnullp5b,SNc5b,SNnullc5b)

#THIRD TRIAL
SNperturb5c <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5cm <- as.matrix(SNperturb5c)
SNp5c <- graph_from_adjacency_matrix(SNp5cm, mode = "undirected", weighted = TRUE)
SNnullperturb5c <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5cm <- as.matrix(SNnullperturb5c)
SNnullp5c <- graph_from_adjacency_matrix(SNnullp5cm, mode = "undirected", weighted = TRUE)

SNc5c <- cluster_louvain(SNp5c, weights = E(SNp5c)$weight)
SNnullc5c <- cluster_louvain(SNnullp5c, weights = E(SNnullp5c)$weight)

VICSN5[4] <- compare(SNc,SNc5c,method = "vi")
VICSNrandom5[4] <- compare(SNnullc, SNnullc5c, method = "vi")

rm(SNperturb5c,SNp5c,SNp5cm,SNnullperturb5c,SNnullp5cm,SNnullp5c,SNc5c,SNnullc5c)

#FOURTH TRIAL
SNperturb5d <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5dm <- as.matrix(SNperturb5d)
SNp5d <- graph_from_adjacency_matrix(SNp5dm, mode = "undirected", weighted = TRUE)
SNnullperturb5d <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5dm <- as.matrix(SNnullperturb5d)
SNnullp5d <- graph_from_adjacency_matrix(SNnullp5dm, mode = "undirected", weighted = TRUE)

SNc5d <- cluster_louvain(SNp5d, weights = E(SNp5d)$weight)
SNnullc5d <- cluster_louvain(SNnullp5d, weights = E(SNnullp5d)$weight)

VICSN5[5] <- compare(SNc,SNc5d,method = "vi")
VICSNrandom5[5] <- compare(SNnullc, SNnullc5d, method = "vi")

rm(SNperturb5d,SNp5d,SNp5dm,SNnullperturb5d,SNnullp5dm,SNnullp5d,SNc5d,SNnullc5d)

#FIFTH TRIAL
SNperturb5e <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5em <- as.matrix(SNperturb5e)
SNp5e <- graph_from_adjacency_matrix(SNp5em, mode = "undirected", weighted = TRUE)
SNnullperturb5e <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5em <- as.matrix(SNnullperturb5e)
SNnullp5e <- graph_from_adjacency_matrix(SNnullp5em, mode = "undirected", weighted = TRUE)

SNc5e <- cluster_louvain(SNp5e, weights = E(SNp5e)$weight)
SNnullc5e <- cluster_louvain(SNnullp5e, weights = E(SNnullp5e)$weight)

VICSN5[6] <- compare(SNc,SNc5e,method = "vi")
VICSNrandom5[6] <- compare(SNnullc, SNnullc5e, method = "vi")

rm(SNperturb5e,SNp5e,SNp5em,SNnullperturb5e,SNnullp5em,SNnullp5e,SNc5e,SNnullc5e)

#SIXTH TRIAL
SNperturb5f <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5fm <- as.matrix(SNperturb5f)
SNp5f <- graph_from_adjacency_matrix(SNp5fm, mode = "undirected", weighted = TRUE)
SNnullperturb5f <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5fm <- as.matrix(SNnullperturb5f)
SNnullp5f <- graph_from_adjacency_matrix(SNnullp5fm, mode = "undirected", weighted = TRUE)

SNc5f <- cluster_louvain(SNp5f, weights = E(SNp5f)$weight)
SNnullc5f <- cluster_louvain(SNnullp5f, weights = E(SNnullp5f)$weight)

VICSN5[7] <- compare(SNc,SNc5f,method = "vi")
VICSNrandom5[7] <- compare(SNnullc, SNnullc5f, method = "vi")

rm(SNperturb5f,SNp5f,SNp5fm,SNnullperturb5f,SNnullp5fm,SNnullp5f,SNc5f,SNnullc5f)

#SEVENTH TRIAL
SNperturb5g <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5gm <- as.matrix(SNperturb5g)
SNp5g <- graph_from_adjacency_matrix(SNp5gm, mode = "undirected", weighted = TRUE)
SNnullperturb5g <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5gm <- as.matrix(SNnullperturb5g)
SNnullp5g <- graph_from_adjacency_matrix(SNnullp5gm, mode = "undirected", weighted = TRUE)

SNc5g <- cluster_louvain(SNp5g, weights = E(SNp5g)$weight)
SNnullc5g <- cluster_louvain(SNnullp5g, weights = E(SNnullp5g)$weight)

VICSN5[8] <- compare(SNc,SNc5g,method = "vi")
VICSNrandom5[8] <- compare(SNnullc, SNnullc5g, method = "vi")

rm(SNperturb5g,SNp5g,SNp5gm,SNnullperturb5g,SNnullp5gm,SNnullp5g,SNc5g,SNnullc5g)

#EIGHTH TRIAL
SNperturb5h <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5hm <- as.matrix(SNperturb5h)
SNp5h <- graph_from_adjacency_matrix(SNp5hm, mode = "undirected", weighted = TRUE)
SNnullperturb5h <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5hm <- as.matrix(SNnullperturb5h)
SNnullp5h <- graph_from_adjacency_matrix(SNnullp5hm, mode = "undirected", weighted = TRUE)

SNc5h <- cluster_louvain(SNp5h, weights = E(SNp5h)$weight)
SNnullc5h <- cluster_louvain(SNnullp5h, weights = E(SNnullp5h)$weight)

VICSN5[9] <- compare(SNc,SNc5h,method = "vi")
VICSNrandom5[9] <- compare(SNnullc, SNnullc5h, method = "vi")

rm(SNperturb5h,SNp5h,SNp5hm,SNnullperturb5h,SNnullp5hm,SNnullp5h,SNc5h,SNnullc5h)

#NINTH TRIAL
SNperturb5i <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5im <- as.matrix(SNperturb5i)
SNp5i <- graph_from_adjacency_matrix(SNp5im, mode = "undirected", weighted = TRUE)
SNnullperturb5i <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5im <- as.matrix(SNnullperturb5i)
SNnullp5i <- graph_from_adjacency_matrix(SNnullp5im, mode = "undirected", weighted = TRUE)

SNc5i <- cluster_louvain(SNp5i, weights = E(SNp5i)$weight)
SNnullc5i <- cluster_louvain(SNnullp5i, weights = E(SNnullp5i)$weight)

VICSN5[10] <- compare(SNc,SNc5i,method = "vi")
VICSNrandom5[10] <- compare(SNnullc, SNnullc5i, method = "vi")

rm(SNperturb5i,SNp5i,SNp5im,SNnullperturb5i,SNnullp5im,SNnullp5i,SNc5i,SNnullc5i)

#TENTH TRIAL
SNperturb5j <- rewireR(SNadj, 6237, dist = "NegBinom")
SNp5jm <- as.matrix(SNperturb5j)
SNp5j <- graph_from_adjacency_matrix(SNp5jm, mode = "undirected", weighted = TRUE)
SNnullperturb5j <- rewireR(SNnullmatrix, 6237, dist = "NegBinom")
SNnullp5jm <- as.matrix(SNnullperturb5j)
SNnullp5j <- graph_from_adjacency_matrix(SNnullp5jm, mode = "undirected", weighted = TRUE)

SNc5j <- cluster_louvain(SNp5j, weights = E(SNp5j)$weight)
SNnullc5j <- cluster_louvain(SNnullp5j, weights = E(SNnullp5j)$weight)

VICSN5[11] <- compare(SNc,SNc5j,method = "vi")
VICSNrandom5[11] <- compare(SNnullc, SNnullc5j, method = "vi")

rm(SNperturb5j,SNp5j,SNp5jm,SNnullperturb5j,SNnullp5jm,SNnullp5j,SNc5j,SNnullc5j)

###############10% Perturbation level#####################

#FIRST TRIAL
SNperturb10a <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10am <- as.matrix(SNperturb10a)
SNp10a <- graph_from_adjacency_matrix(SNp10am, mode = "undirected", weighted = TRUE)
SNnullperturb10a <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10am <- as.matrix(SNnullperturb10a)
SNnullp10a <- graph_from_adjacency_matrix(SNnullp10am, mode = "undirected", weighted = TRUE)

SNc10a <- cluster_louvain(SNp10a, weights = E(SNp10a)$weight)
SNnullc10a <- cluster_louvain(SNnullp10a, weights = E(SNnullp10a)$weight)

VICSN10[2] <- compare(SNc,SNc10a,method = "vi")
VICSNrandom10[2] <- compare(SNnullc, SNnullc10a, method = "vi")

rm(SNperturb10a,SNp10a,SNp10am,SNnullperturb10a,SNnullp10am,SNnullp10a,SNc10a,SNnullc10a)

#SECOND TRIAL
SNperturb10b <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10bm <- as.matrix(SNperturb10b)
SNp10b <- graph_from_adjacency_matrix(SNp10bm, mode = "undirected", weighted = TRUE)
SNnullperturb10b <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10bm <- as.matrix(SNnullperturb10b)
SNnullp10b <- graph_from_adjacency_matrix(SNnullp10bm, mode = "undirected", weighted = TRUE)

SNc10b <- cluster_louvain(SNp10b, weights = E(SNp10b)$weight)
SNnullc10b <- cluster_louvain(SNnullp10b, weights = E(SNnullp10b)$weight)

VICSN10[3] <- compare(SNc,SNc10b,method = "vi")
VICSNrandom10[3] <- compare(SNnullc, SNnullc10b, method = "vi")

rm(SNperturb10b,SNp10b,SNp10bm,SNnullperturb10b,SNnullp10bm,SNnullp10b,SNc10b,SNnullc10b)

#THIRD TRIAL
SNperturb10c <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10cm <- as.matrix(SNperturb10c)
SNp10c <- graph_from_adjacency_matrix(SNp10cm, mode = "undirected", weighted = TRUE)
SNnullperturb10c <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10cm <- as.matrix(SNnullperturb10c)
SNnullp10c <- graph_from_adjacency_matrix(SNnullp10cm, mode = "undirected", weighted = TRUE)

SNc10c <- cluster_louvain(SNp10c, weights = E(SNp10c)$weight)
SNnullc10c <- cluster_louvain(SNnullp10c, weights = E(SNnullp10c)$weight)

VICSN10[4] <- compare(SNc,SNc10c,method = "vi")
VICSNrandom10[4] <- compare(SNnullc, SNnullc10c, method = "vi")

rm(SNperturb10c,SNp10c,SNp10cm,SNnullperturb10c,SNnullp10cm,SNnullp10c,SNc10c,SNnullc10c)

#FOURTH TRIAL
SNperturb10d <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10dm <- as.matrix(SNperturb10d)
SNp10d <- graph_from_adjacency_matrix(SNp10dm, mode = "undirected", weighted = TRUE)
SNnullperturb10d <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10dm <- as.matrix(SNnullperturb10d)
SNnullp10d <- graph_from_adjacency_matrix(SNnullp10dm, mode = "undirected", weighted = TRUE)

SNc10d <- cluster_louvain(SNp10d, weights = E(SNp10d)$weight)
SNnullc10d <- cluster_louvain(SNnullp10d, weights = E(SNnullp10d)$weight)

VICSN10[5] <- compare(SNc,SNc10d,method = "vi")
VICSNrandom10[5] <- compare(SNnullc, SNnullc10d, method = "vi")

rm(SNperturb10d,SNp10d,SNp10dm,SNnullperturb10d,SNnullp10dm,SNnullp10d,SNc10d,SNnullc10d)

#FIFTH TRIAL
SNperturb10e <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10em <- as.matrix(SNperturb10e)
SNp10e <- graph_from_adjacency_matrix(SNp10em, mode = "undirected", weighted = TRUE)
SNnullperturb10e <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10em <- as.matrix(SNnullperturb10e)
SNnullp10e <- graph_from_adjacency_matrix(SNnullp10em, mode = "undirected", weighted = TRUE)

SNc10e <- cluster_louvain(SNp10e, weights = E(SNp10e)$weight)
SNnullc10e <- cluster_louvain(SNnullp10e, weights = E(SNnullp10e)$weight)

VICSN10[6] <- compare(SNc,SNc10e,method = "vi")
VICSNrandom10[6] <- compare(SNnullc, SNnullc10e, method = "vi")

rm(SNperturb10e,SNp10e,SNp10em,SNnullperturb10e,SNnullp10em,SNnullp10e,SNc10e,SNnullc10e)

#SIXTH TRIAL
SNperturb10f <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10fm <- as.matrix(SNperturb10f)
SNp10f <- graph_from_adjacency_matrix(SNp10fm, mode = "undirected", weighted = TRUE)
SNnullperturb10f <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10fm <- as.matrix(SNnullperturb10f)
SNnullp10f <- graph_from_adjacency_matrix(SNnullp10fm, mode = "undirected", weighted = TRUE)

SNc10f <- cluster_louvain(SNp10f, weights = E(SNp10f)$weight)
SNnullc10f <- cluster_louvain(SNnullp10f, weights = E(SNnullp10f)$weight)

VICSN10[7] <- compare(SNc,SNc10f,method = "vi")
VICSNrandom10[7] <- compare(SNnullc, SNnullc10f, method = "vi")

rm(SNperturb10f,SNp10f,SNp10fm,SNnullperturb10f,SNnullp10fm,SNnullp10f,SNc10f,SNnullc10f)

#SEVENTH TRIAL
SNperturb10g <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10gm <- as.matrix(SNperturb10g)
SNp10g <- graph_from_adjacency_matrix(SNp10gm, mode = "undirected", weighted = TRUE)
SNnullperturb10g <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10gm <- as.matrix(SNnullperturb10g)
SNnullp10g <- graph_from_adjacency_matrix(SNnullp10gm, mode = "undirected", weighted = TRUE)

SNc10g <- cluster_louvain(SNp10g, weights = E(SNp10g)$weight)
SNnullc10g <- cluster_louvain(SNnullp10g, weights = E(SNnullp10g)$weight)

VICSN10[8] <- compare(SNc,SNc10g,method = "vi")
VICSNrandom10[8] <- compare(SNnullc, SNnullc10g, method = "vi")

rm(SNperturb10g,SNp10g,SNp10gm,SNnullperturb10g,SNnullp10gm,SNnullp10g,SNc10g,SNnullc10g)

#EIGHTH TRIAL
SNperturb10h <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10hm <- as.matrix(SNperturb10h)
SNp10h <- graph_from_adjacency_matrix(SNp10hm, mode = "undirected", weighted = TRUE)
SNnullperturb10h <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10hm <- as.matrix(SNnullperturb10h)
SNnullp10h <- graph_from_adjacency_matrix(SNnullp10hm, mode = "undirected", weighted = TRUE)

SNc10h <- cluster_louvain(SNp10h, weights = E(SNp10h)$weight)
SNnullc10h <- cluster_louvain(SNnullp10h, weights = E(SNnullp10h)$weight)

VICSN10[9] <- compare(SNc,SNc10h,method = "vi")
VICSNrandom10[9] <- compare(SNnullc, SNnullc10h, method = "vi")

rm(SNperturb10h,SNp10h,SNp10hm,SNnullperturb10h,SNnullp10hm,SNnullp10h,SNc10h,SNnullc10h)

#NINTH TRIAL
SNperturb10i <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10im <- as.matrix(SNperturb10i)
SNp10i <- graph_from_adjacency_matrix(SNp10im, mode = "undirected", weighted = TRUE)
SNnullperturb10i <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10im <- as.matrix(SNnullperturb10i)
SNnullp10i <- graph_from_adjacency_matrix(SNnullp10im, mode = "undirected", weighted = TRUE)

SNc10i <- cluster_louvain(SNp10i, weights = E(SNp10i)$weight)
SNnullc10i <- cluster_louvain(SNnullp10i, weights = E(SNnullp10i)$weight)

VICSN10[10] <- compare(SNc,SNc10i,method = "vi")
VICSNrandom10[10] <- compare(SNnullc, SNnullc10i, method = "vi")

rm(SNperturb10i,SNp10i,SNp10im,SNnullperturb10i,SNnullp10im,SNnullp10i,SNc10i,SNnullc10i)

#TENTH TRIAL
SNperturb10j <- rewireR(SNadj, 12475, dist = "NegBinom")
SNp10jm <- as.matrix(SNperturb10j)
SNp10j <- graph_from_adjacency_matrix(SNp10jm, mode = "undirected", weighted = TRUE)
SNnullperturb10j <- rewireR(SNnullmatrix, 12475, dist = "NegBinom")
SNnullp10jm <- as.matrix(SNnullperturb10j)
SNnullp10j <- graph_from_adjacency_matrix(SNnullp10jm, mode = "undirected", weighted = TRUE)

SNc10j <- cluster_louvain(SNp10j, weights = E(SNp10j)$weight)
SNnullc10j <- cluster_louvain(SNnullp10j, weights = E(SNnullp10j)$weight)

VICSN10[11] <- compare(SNc,SNc10j,method = "vi")
VICSNrandom10[11] <- compare(SNnullc, SNnullc10j, method = "vi")

rm(SNperturb10j,SNp10j,SNp10jm,SNnullperturb10j,SNnullp10jm,SNnullp10j,SNc10j,SNnullc10j)

################15% Perturbation Level#####################

#FIRST TRIAL
SNperturb15a <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15am <- as.matrix(SNperturb15a)
SNp15a <- graph_from_adjacency_matrix(SNp15am, mode = "undirected", weighted = TRUE)
SNnullperturb15a <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15am <- as.matrix(SNnullperturb15a)
SNnullp15a <- graph_from_adjacency_matrix(SNnullp15am, mode = "undirected", weighted = TRUE)

SNc15a <- cluster_louvain(SNp15a, weights = E(SNp15a)$weight)
SNnullc15a <- cluster_louvain(SNnullp15a, weights = E(SNnullp15a)$weight)

VICSN15[2] <- compare(SNc,SNc15a,method = "vi")
VICSNrandom15[2] <- compare(SNnullc, SNnullc15a, method = "vi")

rm(SNperturb15a,SNp15a,SNp15am,SNnullperturb15a,SNnullp15am,SNnullp15a,SNc15a,SNnullc15a)

#SECOND TRIAL
SNperturb15b <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15bm <- as.matrix(SNperturb15b)
SNp15b <- graph_from_adjacency_matrix(SNp15bm, mode = "undirected", weighted = TRUE)
SNnullperturb15b <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15bm <- as.matrix(SNnullperturb15b)
SNnullp15b <- graph_from_adjacency_matrix(SNnullp15bm, mode = "undirected", weighted = TRUE)

SNc15b <- cluster_louvain(SNp15b, weights = E(SNp15b)$weight)
SNnullc15b <- cluster_louvain(SNnullp15b, weights = E(SNnullp15b)$weight)

VICSN15[3] <- compare(SNc,SNc15b,method = "vi")
VICSNrandom15[3] <- compare(SNnullc, SNnullc15b, method = "vi")

rm(SNperturb15b,SNp15b,SNp15bm,SNnullperturb15b,SNnullp15bm,SNnullp15b,SNc15b,SNnullc15b)

#THIRD TRIAL
SNperturb15c <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15cm <- as.matrix(SNperturb15c)
SNp15c <- graph_from_adjacency_matrix(SNp15cm, mode = "undirected", weighted = TRUE)
SNnullperturb15c <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15cm <- as.matrix(SNnullperturb15c)
SNnullp15c <- graph_from_adjacency_matrix(SNnullp15cm, mode = "undirected", weighted = TRUE)

SNc15c <- cluster_louvain(SNp15c, weights = E(SNp15c)$weight)
SNnullc15c <- cluster_louvain(SNnullp15c, weights = E(SNnullp15c)$weight)

VICSN15[4] <- compare(SNc,SNc15c,method = "vi")
VICSNrandom15[4] <- compare(SNnullc, SNnullc15c, method = "vi")

rm(SNperturb15c,SNp15c,SNp15cm,SNnullperturb15c,SNnullp15cm,SNnullp15c,SNc15c,SNnullc15c)

#FOURTH TRIAL
SNperturb15d <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15dm <- as.matrix(SNperturb15d)
SNp15d <- graph_from_adjacency_matrix(SNp15dm, mode = "undirected", weighted = TRUE)
SNnullperturb15d <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15dm <- as.matrix(SNnullperturb15d)
SNnullp15d <- graph_from_adjacency_matrix(SNnullp15dm, mode = "undirected", weighted = TRUE)

SNc15d <- cluster_louvain(SNp15d, weights = E(SNp15d)$weight)
SNnullc15d <- cluster_louvain(SNnullp15d, weights = E(SNnullp15d)$weight)

VICSN15[5] <- compare(SNc,SNc15d,method = "vi")
VICSNrandom15[5] <- compare(SNnullc, SNnullc15d, method = "vi")

rm(SNperturb15d,SNp15d,SNp15dm,SNnullperturb15d,SNnullp15dm,SNnullp15d,SNc15d,SNnullc15d)

#FIFTH TRIAL
SNperturb15e <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15em <- as.matrix(SNperturb15e)
SNp15e <- graph_from_adjacency_matrix(SNp15em, mode = "undirected", weighted = TRUE)
SNnullperturb15e <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15em <- as.matrix(SNnullperturb15e)
SNnullp15e <- graph_from_adjacency_matrix(SNnullp15em, mode = "undirected", weighted = TRUE)

SNc15e <- cluster_louvain(SNp15e, weights = E(SNp15e)$weight)
SNnullc15e <- cluster_louvain(SNnullp15e, weights = E(SNnullp15e)$weight)

VICSN15[6] <- compare(SNc,SNc15e,method = "vi")
VICSNrandom15[6] <- compare(SNnullc, SNnullc15e, method = "vi")

rm(SNperturb15e,SNp15e,SNp15em,SNnullperturb15e,SNnullp15em,SNnullp15e,SNc15e,SNnullc15e)

#SIXTH TRIAL
SNperturb15f <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15fm <- as.matrix(SNperturb15f)
SNp15f <- graph_from_adjacency_matrix(SNp15fm, mode = "undirected", weighted = TRUE)
SNnullperturb15f <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15fm <- as.matrix(SNnullperturb15f)
SNnullp15f <- graph_from_adjacency_matrix(SNnullp15fm, mode = "undirected", weighted = TRUE)

SNc15f <- cluster_louvain(SNp15f, weights = E(SNp15f)$weight)
SNnullc15f <- cluster_louvain(SNnullp15f, weights = E(SNnullp15f)$weight)

VICSN15[7] <- compare(SNc,SNc15f,method = "vi")
VICSNrandom15[7] <- compare(SNnullc, SNnullc15f, method = "vi")

rm(SNperturb15f,SNp15f,SNp15fm,SNnullperturb15f,SNnullp15fm,SNnullp15f,SNc15f,SNnullc15f)

#SEVENTH TRIAL
SNperturb15g <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15gm <- as.matrix(SNperturb15g)
SNp15g <- graph_from_adjacency_matrix(SNp15gm, mode = "undirected", weighted = TRUE)
SNnullperturb15g <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15gm <- as.matrix(SNnullperturb15g)
SNnullp15g <- graph_from_adjacency_matrix(SNnullp15gm, mode = "undirected", weighted = TRUE)

SNc15g <- cluster_louvain(SNp15g, weights = E(SNp15g)$weight)
SNnullc15g <- cluster_louvain(SNnullp15g, weights = E(SNnullp15g)$weight)

VICSN15[8] <- compare(SNc,SNc15g,method = "vi")
VICSNrandom15[8] <- compare(SNnullc, SNnullc15g, method = "vi")

rm(SNperturb15g,SNp15g,SNp15gm,SNnullperturb15g,SNnullp15gm,SNnullp15g,SNc15g,SNnullc15g)

#EIGHTH TRIAL
SNperturb15h <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15hm <- as.matrix(SNperturb15h)
SNp15h <- graph_from_adjacency_matrix(SNp15hm, mode = "undirected", weighted = TRUE)
SNnullperturb15h <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15hm <- as.matrix(SNnullperturb15h)
SNnullp15h <- graph_from_adjacency_matrix(SNnullp15hm, mode = "undirected", weighted = TRUE)

SNc15h <- cluster_louvain(SNp15h, weights = E(SNp15h)$weight)
SNnullc15h <- cluster_louvain(SNnullp15h, weights = E(SNnullp15h)$weight)

VICSN15[9] <- compare(SNc,SNc15h,method = "vi")
VICSNrandom15[9] <- compare(SNnullc, SNnullc15h, method = "vi")

rm(SNperturb15h,SNp15h,SNp15hm,SNnullperturb15h,SNnullp15hm,SNnullp15h,SNc15h,SNnullc15h)

#NINTH TRIAL
SNperturb15i <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15im <- as.matrix(SNperturb15i)
SNp15i <- graph_from_adjacency_matrix(SNp15im, mode = "undirected", weighted = TRUE)
SNnullperturb15i <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15im <- as.matrix(SNnullperturb15i)
SNnullp15i <- graph_from_adjacency_matrix(SNnullp15im, mode = "undirected", weighted = TRUE)

SNc15i <- cluster_louvain(SNp15i, weights = E(SNp15i)$weight)
SNnullc15i <- cluster_louvain(SNnullp15i, weights = E(SNnullp15i)$weight)

VICSN15[10] <- compare(SNc,SNc15i,method = "vi")
VICSNrandom15[10] <- compare(SNnullc, SNnullc15i, method = "vi")

rm(SNperturb15i,SNp15i,SNp15im,SNnullperturb15i,SNnullp15im,SNnullp15i,SNc15i,SNnullc15i)

#TENTH TRIAL
SNperturb15j <- rewireR(SNadj, 18712, dist = "NegBinom")
SNp15jm <- as.matrix(SNperturb15j)
SNp15j <- graph_from_adjacency_matrix(SNp15jm, mode = "undirected", weighted = TRUE)
SNnullperturb15j <- rewireR(SNnullmatrix, 18712, dist = "NegBinom")
SNnullp15jm <- as.matrix(SNnullperturb15j)
SNnullp15j <- graph_from_adjacency_matrix(SNnullp15jm, mode = "undirected", weighted = TRUE)

SNc15j <- cluster_louvain(SNp15j, weights = E(SNp15j)$weight)
SNnullc15j <- cluster_louvain(SNnullp15j, weights = E(SNnullp15j)$weight)

VICSN15[11] <- compare(SNc,SNc15j,method = "vi")
VICSNrandom15[11] <- compare(SNnullc, SNnullc15j, method = "vi")

rm(SNperturb15j,SNp15j,SNp15jm,SNnullperturb15j,SNnullp15jm,SNnullp15j,SNc15j,SNnullc15j)

################20% Perturbation Level#######################

#FIRST TRIAL
SNperturb20a <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20am <- as.matrix(SNperturb20a)
SNp20a <- graph_from_adjacency_matrix(SNp20am, mode = "undirected", weighted = TRUE)
SNnullperturb20a <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20am <- as.matrix(SNnullperturb20a)
SNnullp20a <- graph_from_adjacency_matrix(SNnullp20am, mode = "undirected", weighted = TRUE)

SNc20a <- cluster_louvain(SNp20a, weights = E(SNp20a)$weight)
SNnullc20a <- cluster_louvain(SNnullp20a, weights = E(SNnullp20a)$weight)

VICSN20[2] <- compare(SNc,SNc20a,method = "vi")
VICSNrandom20[2] <- compare(SNnullc, SNnullc20a, method = "vi")

rm(SNperturb20a,SNp20a,SNp20am,SNnullperturb20a,SNnullp20am,SNnullp20a,SNc20a,SNnullc20a)

#SECOND TRIAL
SNperturb20b <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20bm <- as.matrix(SNperturb20b)
SNp20b <- graph_from_adjacency_matrix(SNp20bm, mode = "undirected", weighted = TRUE)
SNnullperturb20b <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20bm <- as.matrix(SNnullperturb20b)
SNnullp20b <- graph_from_adjacency_matrix(SNnullp20bm, mode = "undirected", weighted = TRUE)

SNc20b <- cluster_louvain(SNp20b, weights = E(SNp20b)$weight)
SNnullc20b <- cluster_louvain(SNnullp20b, weights = E(SNnullp20b)$weight)

VICSN20[3] <- compare(SNc,SNc20b,method = "vi")
VICSNrandom20[3] <- compare(SNnullc, SNnullc20b, method = "vi")

rm(SNperturb20b,SNp20b,SNp20bm,SNnullperturb20b,SNnullp20bm,SNnullp20b,SNc20b,SNnullc20b)

#THIRD TRIAL
SNperturb20c <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20cm <- as.matrix(SNperturb20c)
SNp20c <- graph_from_adjacency_matrix(SNp20cm, mode = "undirected", weighted = TRUE)
SNnullperturb20c <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20cm <- as.matrix(SNnullperturb20c)
SNnullp20c <- graph_from_adjacency_matrix(SNnullp20cm, mode = "undirected", weighted = TRUE)

SNc20c <- cluster_louvain(SNp20c, weights = E(SNp20c)$weight)
SNnullc20c <- cluster_louvain(SNnullp20c, weights = E(SNnullp20c)$weight)

VICSN20[4] <- compare(SNc,SNc20c,method = "vi")
VICSNrandom20[4] <- compare(SNnullc, SNnullc20c, method = "vi")

rm(SNperturb20c,SNp20c,SNp20cm,SNnullperturb20c,SNnullp20cm,SNnullp20c,SNc20c,SNnullc20c)

#FOURTH TRIAL
SNperturb20d <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20dm <- as.matrix(SNperturb20d)
SNp20d <- graph_from_adjacency_matrix(SNp20dm, mode = "undirected", weighted = TRUE)
SNnullperturb20d <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20dm <- as.matrix(SNnullperturb20d)
SNnullp20d <- graph_from_adjacency_matrix(SNnullp20dm, mode = "undirected", weighted = TRUE)

SNc20d <- cluster_louvain(SNp20d, weights = E(SNp20d)$weight)
SNnullc20d <- cluster_louvain(SNnullp20d, weights = E(SNnullp20d)$weight)

VICSN20[5] <- compare(SNc,SNc20d,method = "vi")
VICSNrandom20[5] <- compare(SNnullc, SNnullc20d, method = "vi")

rm(SNperturb20d,SNp20d,SNp20dm,SNnullperturb20d,SNnullp20dm,SNnullp20d,SNc20d,SNnullc20d)

#FIFTH TRIAL
SNperturb20e <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20em <- as.matrix(SNperturb20e)
SNp20e <- graph_from_adjacency_matrix(SNp20em, mode = "undirected", weighted = TRUE)
SNnullperturb20e <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20em <- as.matrix(SNnullperturb20e)
SNnullp20e <- graph_from_adjacency_matrix(SNnullp20em, mode = "undirected", weighted = TRUE)

SNc20e <- cluster_louvain(SNp20e, weights = E(SNp20e)$weight)
SNnullc20e <- cluster_louvain(SNnullp20e, weights = E(SNnullp20e)$weight)

VICSN20[6] <- compare(SNc,SNc20e,method = "vi")
VICSNrandom20[6] <- compare(SNnullc, SNnullc20e, method = "vi")

rm(SNperturb20e,SNp20e,SNp20em,SNnullperturb20e,SNnullp20em,SNnullp20e,SNc20e,SNnullc20e)

#SIXTH TRIAL
SNperturb20f <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20fm <- as.matrix(SNperturb20f)
SNp20f <- graph_from_adjacency_matrix(SNp20fm, mode = "undirected", weighted = TRUE)
SNnullperturb20f <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20fm <- as.matrix(SNnullperturb20f)
SNnullp20f <- graph_from_adjacency_matrix(SNnullp20fm, mode = "undirected", weighted = TRUE)

SNc20f <- cluster_louvain(SNp20f, weights = E(SNp20f)$weight)
SNnullc20f <- cluster_louvain(SNnullp20f, weights = E(SNnullp20f)$weight)

VICSN20[7] <- compare(SNc,SNc20f,method = "vi")
VICSNrandom20[7] <- compare(SNnullc, SNnullc20f, method = "vi")

rm(SNperturb20f,SNp20f,SNp20fm,SNnullperturb20f,SNnullp20fm,SNnullp20f,SNc20f,SNnullc20f)

#SEVENTH TRIAL
SNperturb20g <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20gm <- as.matrix(SNperturb20g)
SNp20g <- graph_from_adjacency_matrix(SNp20gm, mode = "undirected", weighted = TRUE)
SNnullperturb20g <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20gm <- as.matrix(SNnullperturb20g)
SNnullp20g <- graph_from_adjacency_matrix(SNnullp20gm, mode = "undirected", weighted = TRUE)

SNc20g <- cluster_louvain(SNp20g, weights = E(SNp20g)$weight)
SNnullc20g <- cluster_louvain(SNnullp20g, weights = E(SNnullp20g)$weight)

VICSN20[8] <- compare(SNc,SNc20g,method = "vi")
VICSNrandom20[8] <- compare(SNnullc, SNnullc20g, method = "vi")

rm(SNperturb20g,SNp20g,SNp20gm,SNnullperturb20g,SNnullp20gm,SNnullp20g,SNc20g,SNnullc20g)

#EIGHTH TRIAL
SNperturb20h <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20hm <- as.matrix(SNperturb20h)
SNp20h <- graph_from_adjacency_matrix(SNp20hm, mode = "undirected", weighted = TRUE)
SNnullperturb20h <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20hm <- as.matrix(SNnullperturb20h)
SNnullp20h <- graph_from_adjacency_matrix(SNnullp20hm, mode = "undirected", weighted = TRUE)

SNc20h <- cluster_louvain(SNp20h, weights = E(SNp20h)$weight)
SNnullc20h <- cluster_louvain(SNnullp20h, weights = E(SNnullp20h)$weight)

VICSN20[9] <- compare(SNc,SNc20h,method = "vi")
VICSNrandom20[9] <- compare(SNnullc, SNnullc20h, method = "vi")

rm(SNperturb20h,SNp20h,SNp20hm,SNnullperturb20h,SNnullp20hm,SNnullp20h,SNc20h,SNnullc20h)

#NINTH TRIAL
SNperturb20i <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20im <- as.matrix(SNperturb20i)
SNp20i <- graph_from_adjacency_matrix(SNp20im, mode = "undirected", weighted = TRUE)
SNnullperturb20i <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20im <- as.matrix(SNnullperturb20i)
SNnullp20i <- graph_from_adjacency_matrix(SNnullp20im, mode = "undirected", weighted = TRUE)

SNc20i <- cluster_louvain(SNp20i, weights = E(SNp20i)$weight)
SNnullc20i <- cluster_louvain(SNnullp20i, weights = E(SNnullp20i)$weight)

VICSN20[10] <- compare(SNc,SNc20i,method = "vi")
VICSNrandom20[10] <- compare(SNnullc, SNnullc20i, method = "vi")

rm(SNperturb20i,SNp20i,SNp20im,SNnullperturb20i,SNnullp20im,SNnullp20i,SNc20i,SNnullc20i)

#TENTH TRIAL
SNperturb20j <- rewireR(SNadj, 24950, dist = "NegBinom")
SNp20jm <- as.matrix(SNperturb20j)
SNp20j <- graph_from_adjacency_matrix(SNp20jm, mode = "undirected", weighted = TRUE)
SNnullperturb20j <- rewireR(SNnullmatrix, 24950, dist = "NegBinom")
SNnullp20jm <- as.matrix(SNnullperturb20j)
SNnullp20j <- graph_from_adjacency_matrix(SNnullp20jm, mode = "undirected", weighted = TRUE)

SNc20j <- cluster_louvain(SNp20j, weights = E(SNp20j)$weight)
SNnullc20j <- cluster_louvain(SNnullp20j, weights = E(SNnullp20j)$weight)

VICSN20[11] <- compare(SNc,SNc20j,method = "vi")
VICSNrandom20[11] <- compare(SNnullc, SNnullc20j, method = "vi")

rm(SNperturb20j,SNp20j,SNp20jm,SNnullperturb20j,SNnullp20jm,SNnullp20j,SNc20j,SNnullc20j)

#############25% Perturbation Level######################

#FIRST TRIAL
SNperturb25a <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25am <- as.matrix(SNperturb25a)
SNp25a <- graph_from_adjacency_matrix(SNp25am, mode = "undirected", weighted = TRUE)
SNnullperturb25a <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25am <- as.matrix(SNnullperturb25a)
SNnullp25a <- graph_from_adjacency_matrix(SNnullp25am, mode = "undirected", weighted = TRUE)

SNc25a <- cluster_louvain(SNp25a, weights = E(SNp25a)$weight)
SNnullc25a <- cluster_louvain(SNnullp25a, weights = E(SNnullp25a)$weight)

VICSN25[2] <- compare(SNc,SNc25a,method = "vi")
VICSNrandom25[2] <- compare(SNnullc, SNnullc25a, method = "vi")

rm(SNperturb25a,SNp25a,SNp25am,SNnullperturb25a,SNnullp25am,SNnullp25a,SNc25a,SNnullc25a)

#SECOND TRIAL
SNperturb25b <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25bm <- as.matrix(SNperturb25b)
SNp25b <- graph_from_adjacency_matrix(SNp25bm, mode = "undirected", weighted = TRUE)
SNnullperturb25b <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25bm <- as.matrix(SNnullperturb25b)
SNnullp25b <- graph_from_adjacency_matrix(SNnullp25bm, mode = "undirected", weighted = TRUE)

SNc25b <- cluster_louvain(SNp25b, weights = E(SNp25b)$weight)
SNnullc25b <- cluster_louvain(SNnullp25b, weights = E(SNnullp25b)$weight)

VICSN25[3] <- compare(SNc,SNc25b,method = "vi")
VICSNrandom25[3] <- compare(SNnullc, SNnullc25b, method = "vi")

rm(SNperturb25b,SNp25b,SNp25bm,SNnullperturb25b,SNnullp25bm,SNnullp25b,SNc25b,SNnullc25b)

#THIRD TRIAL
SNperturb25c <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25cm <- as.matrix(SNperturb25c)
SNp25c <- graph_from_adjacency_matrix(SNp25cm, mode = "undirected", weighted = TRUE)
SNnullperturb25c <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25cm <- as.matrix(SNnullperturb25c)
SNnullp25c <- graph_from_adjacency_matrix(SNnullp25cm, mode = "undirected", weighted = TRUE)

SNc25c <- cluster_louvain(SNp25c, weights = E(SNp25c)$weight)
SNnullc25c <- cluster_louvain(SNnullp25c, weights = E(SNnullp25c)$weight)

VICSN25[4] <- compare(SNc,SNc25c,method = "vi")
VICSNrandom25[4] <- compare(SNnullc, SNnullc25c, method = "vi")

rm(SNperturb25c,SNp25c,SNp25cm,SNnullperturb25c,SNnullp25cm,SNnullp25c,SNc25c,SNnullc25c)

#FOURTH TRIAL
SNperturb25d <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25dm <- as.matrix(SNperturb25d)
SNp25d <- graph_from_adjacency_matrix(SNp25dm, mode = "undirected", weighted = TRUE)
SNnullperturb25d <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25dm <- as.matrix(SNnullperturb25d)
SNnullp25d <- graph_from_adjacency_matrix(SNnullp25dm, mode = "undirected", weighted = TRUE)

SNc25d <- cluster_louvain(SNp25d, weights = E(SNp25d)$weight)
SNnullc25d <- cluster_louvain(SNnullp25d, weights = E(SNnullp25d)$weight)

VICSN25[5] <- compare(SNc,SNc25d,method = "vi")
VICSNrandom25[5] <- compare(SNnullc, SNnullc25d, method = "vi")

rm(SNperturb25d,SNp25d,SNp25dm,SNnullperturb25d,SNnullp25dm,SNnullp25d,SNc25d,SNnullc25d)

#FIFTH TRIAL
SNperturb25e <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25em <- as.matrix(SNperturb25e)
SNp25e <- graph_from_adjacency_matrix(SNp25em, mode = "undirected", weighted = TRUE)
SNnullperturb25e <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25em <- as.matrix(SNnullperturb25e)
SNnullp25e <- graph_from_adjacency_matrix(SNnullp25em, mode = "undirected", weighted = TRUE)

SNc25e <- cluster_louvain(SNp25e, weights = E(SNp25e)$weight)
SNnullc25e <- cluster_louvain(SNnullp25e, weights = E(SNnullp25e)$weight)

VICSN25[6] <- compare(SNc,SNc25e,method = "vi")
VICSNrandom25[6] <- compare(SNnullc, SNnullc25e, method = "vi")

rm(SNperturb25e,SNp25e,SNp25em,SNnullperturb25e,SNnullp25em,SNnullp25e,SNc25e,SNnullc25e)

#SIXTH TRIAL
SNperturb25f <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25fm <- as.matrix(SNperturb25f)
SNp25f <- graph_from_adjacency_matrix(SNp25fm, mode = "undirected", weighted = TRUE)
SNnullperturb25f <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25fm <- as.matrix(SNnullperturb25f)
SNnullp25f <- graph_from_adjacency_matrix(SNnullp25fm, mode = "undirected", weighted = TRUE)

SNc25f <- cluster_louvain(SNp25f, weights = E(SNp25f)$weight)
SNnullc25f <- cluster_louvain(SNnullp25f, weights = E(SNnullp25f)$weight)

VICSN25[7] <- compare(SNc,SNc25f,method = "vi")
VICSNrandom25[7] <- compare(SNnullc, SNnullc25f, method = "vi")

rm(SNperturb25f,SNp25f,SNp25fm,SNnullperturb25f,SNnullp25fm,SNnullp25f,SNc25f,SNnullc25f)

#SEVENTH TRIAL
SNperturb25g <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25gm <- as.matrix(SNperturb25g)
SNp25g <- graph_from_adjacency_matrix(SNp25gm, mode = "undirected", weighted = TRUE)
SNnullperturb25g <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25gm <- as.matrix(SNnullperturb25g)
SNnullp25g <- graph_from_adjacency_matrix(SNnullp25gm, mode = "undirected", weighted = TRUE)

SNc25g <- cluster_louvain(SNp25g, weights = E(SNp25g)$weight)
SNnullc25g <- cluster_louvain(SNnullp25g, weights = E(SNnullp25g)$weight)

VICSN25[8] <- compare(SNc,SNc25g,method = "vi")
VICSNrandom25[8] <- compare(SNnullc, SNnullc25g, method = "vi")

rm(SNperturb25g,SNp25g,SNp25gm,SNnullperturb25g,SNnullp25gm,SNnullp25g,SNc25g,SNnullc25g)

#EIGHTH TRIAL
SNperturb25h <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25hm <- as.matrix(SNperturb25h)
SNp25h <- graph_from_adjacency_matrix(SNp25hm, mode = "undirected", weighted = TRUE)
SNnullperturb25h <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25hm <- as.matrix(SNnullperturb25h)
SNnullp25h <- graph_from_adjacency_matrix(SNnullp25hm, mode = "undirected", weighted = TRUE)

SNc25h <- cluster_louvain(SNp25h, weights = E(SNp25h)$weight)
SNnullc25h <- cluster_louvain(SNnullp25h, weights = E(SNnullp25h)$weight)

VICSN25[9] <- compare(SNc,SNc25h,method = "vi")
VICSNrandom25[9] <- compare(SNnullc, SNnullc25h, method = "vi")

rm(SNperturb25h,SNp25h,SNp25hm,SNnullperturb25h,SNnullp25hm,SNnullp25h,SNc25h,SNnullc25h)

#NINTH TRIAL
SNperturb25i <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25im <- as.matrix(SNperturb25i)
SNp25i <- graph_from_adjacency_matrix(SNp25im, mode = "undirected", weighted = TRUE)
SNnullperturb25i <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25im <- as.matrix(SNnullperturb25i)
SNnullp25i <- graph_from_adjacency_matrix(SNnullp25im, mode = "undirected", weighted = TRUE)

SNc25i <- cluster_louvain(SNp25i, weights = E(SNp25i)$weight)
SNnullc25i <- cluster_louvain(SNnullp25i, weights = E(SNnullp25i)$weight)

VICSN25[10] <- compare(SNc,SNc25i,method = "vi")
VICSNrandom25[10] <- compare(SNnullc, SNnullc25i, method = "vi")

rm(SNperturb25i,SNp25i,SNp25im,SNnullperturb25i,SNnullp25im,SNnullp25i,SNc25i,SNnullc25i)

#TENTH TRIAL
SNperturb25j <- rewireR(SNadj, 31187, dist = "NegBinom")
SNp25jm <- as.matrix(SNperturb25j)
SNp25j <- graph_from_adjacency_matrix(SNp25jm, mode = "undirected", weighted = TRUE)
SNnullperturb25j <- rewireR(SNnullmatrix, 31187, dist = "NegBinom")
SNnullp25jm <- as.matrix(SNnullperturb25j)
SNnullp25j <- graph_from_adjacency_matrix(SNnullp25jm, mode = "undirected", weighted = TRUE)

SNc25j <- cluster_louvain(SNp25j, weights = E(SNp25j)$weight)
SNnullc25j <- cluster_louvain(SNnullp25j, weights = E(SNnullp25j)$weight)

VICSN25[11] <- compare(SNc,SNc25j,method = "vi")
VICSNrandom25[11] <- compare(SNnullc, SNnullc25j, method = "vi")

rm(SNperturb25j,SNp25j,SNp25jm,SNnullperturb25j,SNnullp25jm,SNnullp25j,SNc25j,SNnullc25j)

############30% Perturbation Level#######################

#FIRST TRIAL
SNperturb30a <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30am <- as.matrix(SNperturb30a)
SNp30a <- graph_from_adjacency_matrix(SNp30am, mode = "undirected", weighted = TRUE)
SNnullperturb30a <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30am <- as.matrix(SNnullperturb30a)
SNnullp30a <- graph_from_adjacency_matrix(SNnullp30am, mode = "undirected", weighted = TRUE)

SNc30a <- cluster_louvain(SNp30a, weights = E(SNp30a)$weight)
SNnullc30a <- cluster_louvain(SNnullp30a, weights = E(SNnullp30a)$weight)

VICSN30[2] <- compare(SNc,SNc30a,method = "vi")
VICSNrandom30[2] <- compare(SNnullc, SNnullc30a, method = "vi")

rm(SNperturb30a,SNp30a,SNp30am,SNnullperturb30a,SNnullp30am,SNnullp30a,SNc30a,SNnullc30a)

#SECOND TRIAL
SNperturb30b <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30bm <- as.matrix(SNperturb30b)
SNp30b <- graph_from_adjacency_matrix(SNp30bm, mode = "undirected", weighted = TRUE)
SNnullperturb30b <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30bm <- as.matrix(SNnullperturb30b)
SNnullp30b <- graph_from_adjacency_matrix(SNnullp30bm, mode = "undirected", weighted = TRUE)

SNc30b <- cluster_louvain(SNp30b, weights = E(SNp30b)$weight)
SNnullc30b <- cluster_louvain(SNnullp30b, weights = E(SNnullp30b)$weight)

VICSN30[3] <- compare(SNc,SNc30b,method = "vi")
VICSNrandom30[3] <- compare(SNnullc, SNnullc30b, method = "vi")

rm(SNperturb30b,SNp30b,SNp30bm,SNnullperturb30b,SNnullp30bm,SNnullp30b,SNc30b,SNnullc30b)

#THIRD TRIAL
SNperturb30c <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30cm <- as.matrix(SNperturb30c)
SNp30c <- graph_from_adjacency_matrix(SNp30cm, mode = "undirected", weighted = TRUE)
SNnullperturb30c <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30cm <- as.matrix(SNnullperturb30c)
SNnullp30c <- graph_from_adjacency_matrix(SNnullp30cm, mode = "undirected", weighted = TRUE)

SNc30c <- cluster_louvain(SNp30c, weights = E(SNp30c)$weight)
SNnullc30c <- cluster_louvain(SNnullp30c, weights = E(SNnullp30c)$weight)

VICSN30[4] <- compare(SNc,SNc30c,method = "vi")
VICSNrandom30[4] <- compare(SNnullc, SNnullc30c, method = "vi")

rm(SNperturb30c,SNp30c,SNp30cm,SNnullperturb30c,SNnullp30cm,SNnullp30c,SNc30c,SNnullc30c)

#FOURTH TRIAL
SNperturb30d <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30dm <- as.matrix(SNperturb30d)
SNp30d <- graph_from_adjacency_matrix(SNp30dm, mode = "undirected", weighted = TRUE)
SNnullperturb30d <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30dm <- as.matrix(SNnullperturb30d)
SNnullp30d <- graph_from_adjacency_matrix(SNnullp30dm, mode = "undirected", weighted = TRUE)

SNc30d <- cluster_louvain(SNp30d, weights = E(SNp30d)$weight)
SNnullc30d <- cluster_louvain(SNnullp30d, weights = E(SNnullp30d)$weight)

VICSN30[5] <- compare(SNc,SNc30d,method = "vi")
VICSNrandom30[5] <- compare(SNnullc, SNnullc30d, method = "vi")

rm(SNperturb30d,SNp30d,SNp30dm,SNnullperturb30d,SNnullp30dm,SNnullp30d,SNc30d,SNnullc30d)

#FIFTH TRIAL
SNperturb30e <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30em <- as.matrix(SNperturb30e)
SNp30e <- graph_from_adjacency_matrix(SNp30em, mode = "undirected", weighted = TRUE)
SNnullperturb30e <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30em <- as.matrix(SNnullperturb30e)
SNnullp30e <- graph_from_adjacency_matrix(SNnullp30em, mode = "undirected", weighted = TRUE)

SNc30e <- cluster_louvain(SNp30e, weights = E(SNp30e)$weight)
SNnullc30e <- cluster_louvain(SNnullp30e, weights = E(SNnullp30e)$weight)

VICSN30[6] <- compare(SNc,SNc30e,method = "vi")
VICSNrandom30[6] <- compare(SNnullc, SNnullc30e, method = "vi")

rm(SNperturb30e,SNp30e,SNp30em,SNnullperturb30e,SNnullp30em,SNnullp30e,SNc30e,SNnullc30e)

#SIXTH TRIAL
SNperturb30f <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30fm <- as.matrix(SNperturb30f)
SNp30f <- graph_from_adjacency_matrix(SNp30fm, mode = "undirected", weighted = TRUE)
SNnullperturb30f <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30fm <- as.matrix(SNnullperturb30f)
SNnullp30f <- graph_from_adjacency_matrix(SNnullp30fm, mode = "undirected", weighted = TRUE)

SNc30f <- cluster_louvain(SNp30f, weights = E(SNp30f)$weight)
SNnullc30f <- cluster_louvain(SNnullp30f, weights = E(SNnullp30f)$weight)

VICSN30[7] <- compare(SNc,SNc30f,method = "vi")
VICSNrandom30[7] <- compare(SNnullc, SNnullc30f, method = "vi")

rm(SNperturb30f,SNp30f,SNp30fm,SNnullperturb30f,SNnullp30fm,SNnullp30f,SNc30f,SNnullc30f)

#SEVENTH TRIAL
SNperturb30g <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30gm <- as.matrix(SNperturb30g)
SNp30g <- graph_from_adjacency_matrix(SNp30gm, mode = "undirected", weighted = TRUE)
SNnullperturb30g <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30gm <- as.matrix(SNnullperturb30g)
SNnullp30g <- graph_from_adjacency_matrix(SNnullp30gm, mode = "undirected", weighted = TRUE)

SNc30g <- cluster_louvain(SNp30g, weights = E(SNp30g)$weight)
SNnullc30g <- cluster_louvain(SNnullp30g, weights = E(SNnullp30g)$weight)

VICSN30[8] <- compare(SNc,SNc30g,method = "vi")
VICSNrandom30[8] <- compare(SNnullc, SNnullc30g, method = "vi")

rm(SNperturb30g,SNp30g,SNp30gm,SNnullperturb30g,SNnullp30gm,SNnullp30g,SNc30g,SNnullc30g)

#EIGHTH TRIAL
SNperturb30h <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30hm <- as.matrix(SNperturb30h)
SNp30h <- graph_from_adjacency_matrix(SNp30hm, mode = "undirected", weighted = TRUE)
SNnullperturb30h <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30hm <- as.matrix(SNnullperturb30h)
SNnullp30h <- graph_from_adjacency_matrix(SNnullp30hm, mode = "undirected", weighted = TRUE)

SNc30h <- cluster_louvain(SNp30h, weights = E(SNp30h)$weight)
SNnullc30h <- cluster_louvain(SNnullp30h, weights = E(SNnullp30h)$weight)

VICSN30[9] <- compare(SNc,SNc30h,method = "vi")
VICSNrandom30[9] <- compare(SNnullc, SNnullc30h, method = "vi")

rm(SNperturb30h,SNp30h,SNp30hm,SNnullperturb30h,SNnullp30hm,SNnullp30h,SNc30h,SNnullc30h)

#NINTH TRIAL
SNperturb30i <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30im <- as.matrix(SNperturb30i)
SNp30i <- graph_from_adjacency_matrix(SNp30im, mode = "undirected", weighted = TRUE)
SNnullperturb30i <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30im <- as.matrix(SNnullperturb30i)
SNnullp30i <- graph_from_adjacency_matrix(SNnullp30im, mode = "undirected", weighted = TRUE)

SNc30i <- cluster_louvain(SNp30i, weights = E(SNp30i)$weight)
SNnullc30i <- cluster_louvain(SNnullp30i, weights = E(SNnullp30i)$weight)

VICSN30[10] <- compare(SNc,SNc30i,method = "vi")
VICSNrandom30[10] <- compare(SNnullc, SNnullc30i, method = "vi")

rm(SNperturb30i,SNp30i,SNp30im,SNnullperturb30i,SNnullp30im,SNnullp30i,SNc30i,SNnullc30i)

#TENTH TRIAL
SNperturb30j <- rewireR(SNadj, 37425, dist = "NegBinom")
SNp30jm <- as.matrix(SNperturb30j)
SNp30j <- graph_from_adjacency_matrix(SNp30jm, mode = "undirected", weighted = TRUE)
SNnullperturb30j <- rewireR(SNnullmatrix, 37425, dist = "NegBinom")
SNnullp30jm <- as.matrix(SNnullperturb30j)
SNnullp30j <- graph_from_adjacency_matrix(SNnullp30jm, mode = "undirected", weighted = TRUE)

SNc30j <- cluster_louvain(SNp30j, weights = E(SNp30j)$weight)
SNnullc30j <- cluster_louvain(SNnullp30j, weights = E(SNnullp30j)$weight)

VICSN30[11] <- compare(SNc,SNc30j,method = "vi")
VICSNrandom30[11] <- compare(SNnullc, SNnullc30j, method = "vi")

rm(SNperturb30j,SNp30j,SNp30jm,SNnullperturb30j,SNnullp30jm,SNnullp30j,SNc30j,SNnullc30j)

##############35% Perturbation Level#####################

#FIRST TRIAL
SNperturb35a <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35am <- as.matrix(SNperturb35a)
SNp35a <- graph_from_adjacency_matrix(SNp35am, mode = "undirected", weighted = TRUE)
SNnullperturb35a <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35am <- as.matrix(SNnullperturb35a)
SNnullp35a <- graph_from_adjacency_matrix(SNnullp35am, mode = "undirected", weighted = TRUE)

SNc35a <- cluster_louvain(SNp35a, weights = E(SNp35a)$weight)
SNnullc35a <- cluster_louvain(SNnullp35a, weights = E(SNnullp35a)$weight)

VICSN35[2] <- compare(SNc,SNc35a,method = "vi")
VICSNrandom35[2] <- compare(SNnullc, SNnullc35a, method = "vi")

rm(SNperturb35a,SNp35a,SNp35am,SNnullperturb35a,SNnullp35am,SNnullp35a,SNc35a,SNnullc35a)

#SECOND TRIAL
SNperturb35b <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35bm <- as.matrix(SNperturb35b)
SNp35b <- graph_from_adjacency_matrix(SNp35bm, mode = "undirected", weighted = TRUE)
SNnullperturb35b <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35bm <- as.matrix(SNnullperturb35b)
SNnullp35b <- graph_from_adjacency_matrix(SNnullp35bm, mode = "undirected", weighted = TRUE)

SNc35b <- cluster_louvain(SNp35b, weights = E(SNp35b)$weight)
SNnullc35b <- cluster_louvain(SNnullp35b, weights = E(SNnullp35b)$weight)

VICSN35[3] <- compare(SNc,SNc35b,method = "vi")
VICSNrandom35[3] <- compare(SNnullc, SNnullc35b, method = "vi")

rm(SNperturb35b,SNp35b,SNp35bm,SNnullperturb35b,SNnullp35bm,SNnullp35b,SNc35b,SNnullc35b)

#THIRD TRIAL
SNperturb35c <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35cm <- as.matrix(SNperturb35c)
SNp35c <- graph_from_adjacency_matrix(SNp35cm, mode = "undirected", weighted = TRUE)
SNnullperturb35c <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35cm <- as.matrix(SNnullperturb35c)
SNnullp35c <- graph_from_adjacency_matrix(SNnullp35cm, mode = "undirected", weighted = TRUE)

SNc35c <- cluster_louvain(SNp35c, weights = E(SNp35c)$weight)
SNnullc35c <- cluster_louvain(SNnullp35c, weights = E(SNnullp35c)$weight)

VICSN35[4] <- compare(SNc,SNc35c,method = "vi")
VICSNrandom35[4] <- compare(SNnullc, SNnullc35c, method = "vi")

rm(SNperturb35c,SNp35c,SNp35cm,SNnullperturb35c,SNnullp35cm,SNnullp35c,SNc35c,SNnullc35c)

#FOURTH TRIAL
SNperturb35d <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35dm <- as.matrix(SNperturb35d)
SNp35d <- graph_from_adjacency_matrix(SNp35dm, mode = "undirected", weighted = TRUE)
SNnullperturb35d <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35dm <- as.matrix(SNnullperturb35d)
SNnullp35d <- graph_from_adjacency_matrix(SNnullp35dm, mode = "undirected", weighted = TRUE)

SNc35d <- cluster_louvain(SNp35d, weights = E(SNp35d)$weight)
SNnullc35d <- cluster_louvain(SNnullp35d, weights = E(SNnullp35d)$weight)

VICSN35[5] <- compare(SNc,SNc35d,method = "vi")
VICSNrandom35[5] <- compare(SNnullc, SNnullc35d, method = "vi")

rm(SNperturb35d,SNp35d,SNp35dm,SNnullperturb35d,SNnullp35dm,SNnullp35d,SNc35d,SNnullc35d)

#FIFTH TRIAL
SNperturb35e <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35em <- as.matrix(SNperturb35e)
SNp35e <- graph_from_adjacency_matrix(SNp35em, mode = "undirected", weighted = TRUE)
SNnullperturb35e <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35em <- as.matrix(SNnullperturb35e)
SNnullp35e <- graph_from_adjacency_matrix(SNnullp35em, mode = "undirected", weighted = TRUE)

SNc35e <- cluster_louvain(SNp35e, weights = E(SNp35e)$weight)
SNnullc35e <- cluster_louvain(SNnullp35e, weights = E(SNnullp35e)$weight)

VICSN35[6] <- compare(SNc,SNc35e,method = "vi")
VICSNrandom35[6] <- compare(SNnullc, SNnullc35e, method = "vi")

rm(SNperturb35e,SNp35e,SNp35em,SNnullperturb35e,SNnullp35em,SNnullp35e,SNc35e,SNnullc35e)

#SIXTH TRIAL
SNperturb35f <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35fm <- as.matrix(SNperturb35f)
SNp35f <- graph_from_adjacency_matrix(SNp35fm, mode = "undirected", weighted = TRUE)
SNnullperturb35f <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35fm <- as.matrix(SNnullperturb35f)
SNnullp35f <- graph_from_adjacency_matrix(SNnullp35fm, mode = "undirected", weighted = TRUE)

SNc35f <- cluster_louvain(SNp35f, weights = E(SNp35f)$weight)
SNnullc35f <- cluster_louvain(SNnullp35f, weights = E(SNnullp35f)$weight)

VICSN35[7] <- compare(SNc,SNc35f,method = "vi")
VICSNrandom35[7] <- compare(SNnullc, SNnullc35f, method = "vi")

rm(SNperturb35f,SNp35f,SNp35fm,SNnullperturb35f,SNnullp35fm,SNnullp35f,SNc35f,SNnullc35f)

#SEVENTH TRIAL
SNperturb35g <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35gm <- as.matrix(SNperturb35g)
SNp35g <- graph_from_adjacency_matrix(SNp35gm, mode = "undirected", weighted = TRUE)
SNnullperturb35g <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35gm <- as.matrix(SNnullperturb35g)
SNnullp35g <- graph_from_adjacency_matrix(SNnullp35gm, mode = "undirected", weighted = TRUE)

SNc35g <- cluster_louvain(SNp35g, weights = E(SNp35g)$weight)
SNnullc35g <- cluster_louvain(SNnullp35g, weights = E(SNnullp35g)$weight)

VICSN35[8] <- compare(SNc,SNc35g,method = "vi")
VICSNrandom35[8] <- compare(SNnullc, SNnullc35g, method = "vi")

rm(SNperturb35g,SNp35g,SNp35gm,SNnullperturb35g,SNnullp35gm,SNnullp35g,SNc35g,SNnullc35g)

#EIGHTH TRIAL
SNperturb35h <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35hm <- as.matrix(SNperturb35h)
SNp35h <- graph_from_adjacency_matrix(SNp35hm, mode = "undirected", weighted = TRUE)
SNnullperturb35h <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35hm <- as.matrix(SNnullperturb35h)
SNnullp35h <- graph_from_adjacency_matrix(SNnullp35hm, mode = "undirected", weighted = TRUE)

SNc35h <- cluster_louvain(SNp35h, weights = E(SNp35h)$weight)
SNnullc35h <- cluster_louvain(SNnullp35h, weights = E(SNnullp35h)$weight)

VICSN35[9] <- compare(SNc,SNc35h,method = "vi")
VICSNrandom35[9] <- compare(SNnullc, SNnullc35h, method = "vi")

rm(SNperturb35h,SNp35h,SNp35hm,SNnullperturb35h,SNnullp35hm,SNnullp35h,SNc35h,SNnullc35h)

#NINTH TRIAL
SNperturb35i <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35im <- as.matrix(SNperturb35i)
SNp35i <- graph_from_adjacency_matrix(SNp35im, mode = "undirected", weighted = TRUE)
SNnullperturb35i <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35im <- as.matrix(SNnullperturb35i)
SNnullp35i <- graph_from_adjacency_matrix(SNnullp35im, mode = "undirected", weighted = TRUE)

SNc35i <- cluster_louvain(SNp35i, weights = E(SNp35i)$weight)
SNnullc35i <- cluster_louvain(SNnullp35i, weights = E(SNnullp35i)$weight)

VICSN35[10] <- compare(SNc,SNc35i,method = "vi")
VICSNrandom35[10] <- compare(SNnullc, SNnullc35i, method = "vi")

rm(SNperturb35i,SNp35i,SNp35im,SNnullperturb35i,SNnullp35im,SNnullp35i,SNc35i,SNnullc35i)

#TENTH TRIAL
SNperturb35j <- rewireR(SNadj, 43662, dist = "NegBinom")
SNp35jm <- as.matrix(SNperturb35j)
SNp35j <- graph_from_adjacency_matrix(SNp35jm, mode = "undirected", weighted = TRUE)
SNnullperturb35j <- rewireR(SNnullmatrix, 43662, dist = "NegBinom")
SNnullp35jm <- as.matrix(SNnullperturb35j)
SNnullp35j <- graph_from_adjacency_matrix(SNnullp35jm, mode = "undirected", weighted = TRUE)

SNc35j <- cluster_louvain(SNp35j, weights = E(SNp35j)$weight)
SNnullc35j <- cluster_louvain(SNnullp35j, weights = E(SNnullp35j)$weight)

VICSN35[11] <- compare(SNc,SNc35j,method = "vi")
VICSNrandom35[11] <- compare(SNnullc, SNnullc35j, method = "vi")

rm(SNperturb35j,SNp35j,SNp35jm,SNnullperturb35j,SNnullp35jm,SNnullp35j,SNc35j,SNnullc35j)

#############40% Perturbation Level#####################

#FIRST TRIAL
SNperturb40a <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40am <- as.matrix(SNperturb40a)
SNp40a <- graph_from_adjacency_matrix(SNp40am, mode = "undirected", weighted = TRUE)
SNnullperturb40a <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40am <- as.matrix(SNnullperturb40a)
SNnullp40a <- graph_from_adjacency_matrix(SNnullp40am, mode = "undirected", weighted = TRUE)

SNc40a <- cluster_louvain(SNp40a, weights = E(SNp40a)$weight)
SNnullc40a <- cluster_louvain(SNnullp40a, weights = E(SNnullp40a)$weight)

VICSN40[2] <- compare(SNc,SNc40a,method = "vi")
VICSNrandom40[2] <- compare(SNnullc, SNnullc40a, method = "vi")

rm(SNperturb40a,SNp40a,SNp40am,SNnullperturb40a,SNnullp40am,SNnullp40a,SNc40a,SNnullc40a)

#SECOND TRIAL
SNperturb40b <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40bm <- as.matrix(SNperturb40b)
SNp40b <- graph_from_adjacency_matrix(SNp40bm, mode = "undirected", weighted = TRUE)
SNnullperturb40b <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40bm <- as.matrix(SNnullperturb40b)
SNnullp40b <- graph_from_adjacency_matrix(SNnullp40bm, mode = "undirected", weighted = TRUE)

SNc40b <- cluster_louvain(SNp40b, weights = E(SNp40b)$weight)
SNnullc40b <- cluster_louvain(SNnullp40b, weights = E(SNnullp40b)$weight)

VICSN40[3] <- compare(SNc,SNc40b,method = "vi")
VICSNrandom40[3] <- compare(SNnullc, SNnullc40b, method = "vi")

rm(SNperturb40b,SNp40b,SNp40bm,SNnullperturb40b,SNnullp40bm,SNnullp40b,SNc40b,SNnullc40b)

#THIRD TRIAL
SNperturb40c <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40cm <- as.matrix(SNperturb40c)
SNp40c <- graph_from_adjacency_matrix(SNp40cm, mode = "undirected", weighted = TRUE)
SNnullperturb40c <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40cm <- as.matrix(SNnullperturb40c)
SNnullp40c <- graph_from_adjacency_matrix(SNnullp40cm, mode = "undirected", weighted = TRUE)

SNc40c <- cluster_louvain(SNp40c, weights = E(SNp40c)$weight)
SNnullc40c <- cluster_louvain(SNnullp40c, weights = E(SNnullp40c)$weight)

VICSN40[4] <- compare(SNc,SNc40c,method = "vi")
VICSNrandom40[4] <- compare(SNnullc, SNnullc40c, method = "vi")

rm(SNperturb40c,SNp40c,SNp40cm,SNnullperturb40c,SNnullp40cm,SNnullp40c,SNc40c,SNnullc40c)

#FOURTH TRIAL
SNperturb40d <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40dm <- as.matrix(SNperturb40d)
SNp40d <- graph_from_adjacency_matrix(SNp40dm, mode = "undirected", weighted = TRUE)
SNnullperturb40d <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40dm <- as.matrix(SNnullperturb40d)
SNnullp40d <- graph_from_adjacency_matrix(SNnullp40dm, mode = "undirected", weighted = TRUE)

SNc40d <- cluster_louvain(SNp40d, weights = E(SNp40d)$weight)
SNnullc40d <- cluster_louvain(SNnullp40d, weights = E(SNnullp40d)$weight)

VICSN40[5] <- compare(SNc,SNc40d,method = "vi")
VICSNrandom40[5] <- compare(SNnullc, SNnullc40d, method = "vi")

rm(SNperturb40d,SNp40d,SNp40dm,SNnullperturb40d,SNnullp40dm,SNnullp40d,SNc40d,SNnullc40d)

#FIFTH TRIAL
SNperturb40e <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40em <- as.matrix(SNperturb40e)
SNp40e <- graph_from_adjacency_matrix(SNp40em, mode = "undirected", weighted = TRUE)
SNnullperturb40e <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40em <- as.matrix(SNnullperturb40e)
SNnullp40e <- graph_from_adjacency_matrix(SNnullp40em, mode = "undirected", weighted = TRUE)

SNc40e <- cluster_louvain(SNp40e, weights = E(SNp40e)$weight)
SNnullc40e <- cluster_louvain(SNnullp40e, weights = E(SNnullp40e)$weight)

VICSN40[6] <- compare(SNc,SNc40e,method = "vi")
VICSNrandom40[6] <- compare(SNnullc, SNnullc40e, method = "vi")

rm(SNperturb40e,SNp40e,SNp40em,SNnullperturb40e,SNnullp40em,SNnullp40e,SNc40e,SNnullc40e)

#SIXTH TRIAL
SNperturb40f <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40fm <- as.matrix(SNperturb40f)
SNp40f <- graph_from_adjacency_matrix(SNp40fm, mode = "undirected", weighted = TRUE)
SNnullperturb40f <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40fm <- as.matrix(SNnullperturb40f)
SNnullp40f <- graph_from_adjacency_matrix(SNnullp40fm, mode = "undirected", weighted = TRUE)

SNc40f <- cluster_louvain(SNp40f, weights = E(SNp40f)$weight)
SNnullc40f <- cluster_louvain(SNnullp40f, weights = E(SNnullp40f)$weight)

VICSN40[7] <- compare(SNc,SNc40f,method = "vi")
VICSNrandom40[7] <- compare(SNnullc, SNnullc40f, method = "vi")

rm(SNperturb40f,SNp40f,SNp40fm,SNnullperturb40f,SNnullp40fm,SNnullp40f,SNc40f,SNnullc40f)

#SEVENTH TRIAL
SNperturb40g <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40gm <- as.matrix(SNperturb40g)
SNp40g <- graph_from_adjacency_matrix(SNp40gm, mode = "undirected", weighted = TRUE)
SNnullperturb40g <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40gm <- as.matrix(SNnullperturb40g)
SNnullp40g <- graph_from_adjacency_matrix(SNnullp40gm, mode = "undirected", weighted = TRUE)

SNc40g <- cluster_louvain(SNp40g, weights = E(SNp40g)$weight)
SNnullc40g <- cluster_louvain(SNnullp40g, weights = E(SNnullp40g)$weight)

VICSN40[8] <- compare(SNc,SNc40g,method = "vi")
VICSNrandom40[8] <- compare(SNnullc, SNnullc40g, method = "vi")

rm(SNperturb40g,SNp40g,SNp40gm,SNnullperturb40g,SNnullp40gm,SNnullp40g,SNc40g,SNnullc40g)

#EIGHTH TRIAL
SNperturb40h <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40hm <- as.matrix(SNperturb40h)
SNp40h <- graph_from_adjacency_matrix(SNp40hm, mode = "undirected", weighted = TRUE)
SNnullperturb40h <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40hm <- as.matrix(SNnullperturb40h)
SNnullp40h <- graph_from_adjacency_matrix(SNnullp40hm, mode = "undirected", weighted = TRUE)

SNc40h <- cluster_louvain(SNp40h, weights = E(SNp40h)$weight)
SNnullc40h <- cluster_louvain(SNnullp40h, weights = E(SNnullp40h)$weight)

VICSN40[9] <- compare(SNc,SNc40h,method = "vi")
VICSNrandom40[9] <- compare(SNnullc, SNnullc40h, method = "vi")

rm(SNperturb40h,SNp40h,SNp40hm,SNnullperturb40h,SNnullp40hm,SNnullp40h,SNc40h,SNnullc40h)

#NINTH TRIAL
SNperturb40i <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40im <- as.matrix(SNperturb40i)
SNp40i <- graph_from_adjacency_matrix(SNp40im, mode = "undirected", weighted = TRUE)
SNnullperturb40i <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40im <- as.matrix(SNnullperturb40i)
SNnullp40i <- graph_from_adjacency_matrix(SNnullp40im, mode = "undirected", weighted = TRUE)

SNc40i <- cluster_louvain(SNp40i, weights = E(SNp40i)$weight)
SNnullc40i <- cluster_louvain(SNnullp40i, weights = E(SNnullp40i)$weight)

VICSN40[10] <- compare(SNc,SNc40i,method = "vi")
VICSNrandom40[10] <- compare(SNnullc, SNnullc40i, method = "vi")

rm(SNperturb40i,SNp40i,SNp40im,SNnullperturb40i,SNnullp40im,SNnullp40i,SNc40i,SNnullc40i)

#TENTH TRIAL
SNperturb40j <- rewireR(SNadj, 49900, dist = "NegBinom")
SNp40jm <- as.matrix(SNperturb40j)
SNp40j <- graph_from_adjacency_matrix(SNp40jm, mode = "undirected", weighted = TRUE)
SNnullperturb40j <- rewireR(SNnullmatrix, 49900, dist = "NegBinom")
SNnullp40jm <- as.matrix(SNnullperturb40j)
SNnullp40j <- graph_from_adjacency_matrix(SNnullp40jm, mode = "undirected", weighted = TRUE)

SNc40j <- cluster_louvain(SNp40j, weights = E(SNp40j)$weight)
SNnullc40j <- cluster_louvain(SNnullp40j, weights = E(SNnullp40j)$weight)

VICSN40[11] <- compare(SNc,SNc40j,method = "vi")
VICSNrandom40[11] <- compare(SNnullc, SNnullc40j, method = "vi")

rm(SNperturb40j,SNp40j,SNp40jm,SNnullperturb40j,SNnullp40jm,SNnullp40j,SNc40j,SNnullc40j)

############45% Perturbation Level######################

#FIRST TRIAL
SNperturb45a <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45am <- as.matrix(SNperturb45a)
SNp45a <- graph_from_adjacency_matrix(SNp45am, mode = "undirected", weighted = TRUE)
SNnullperturb45a <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45am <- as.matrix(SNnullperturb45a)
SNnullp45a <- graph_from_adjacency_matrix(SNnullp45am, mode = "undirected", weighted = TRUE)

SNc45a <- cluster_louvain(SNp45a, weights = E(SNp45a)$weight)
SNnullc45a <- cluster_louvain(SNnullp45a, weights = E(SNnullp45a)$weight)

VICSN45[2] <- compare(SNc,SNc45a,method = "vi")
VICSNrandom45[2] <- compare(SNnullc, SNnullc45a, method = "vi")

rm(SNperturb45a,SNp45a,SNp45am,SNnullperturb45a,SNnullp45am,SNnullp45a,SNc45a,SNnullc45a)

#SECOND TRIAL
SNperturb45b <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45bm <- as.matrix(SNperturb45b)
SNp45b <- graph_from_adjacency_matrix(SNp45bm, mode = "undirected", weighted = TRUE)
SNnullperturb45b <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45bm <- as.matrix(SNnullperturb45b)
SNnullp45b <- graph_from_adjacency_matrix(SNnullp45bm, mode = "undirected", weighted = TRUE)

SNc45b <- cluster_louvain(SNp45b, weights = E(SNp45b)$weight)
SNnullc45b <- cluster_louvain(SNnullp45b, weights = E(SNnullp45b)$weight)

VICSN45[3] <- compare(SNc,SNc45b,method = "vi")
VICSNrandom45[3] <- compare(SNnullc, SNnullc45b, method = "vi")

rm(SNperturb45b,SNp45b,SNp45bm,SNnullperturb45b,SNnullp45bm,SNnullp45b,SNc45b,SNnullc45b)

#THIRD TRIAL
SNperturb45c <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45cm <- as.matrix(SNperturb45c)
SNp45c <- graph_from_adjacency_matrix(SNp45cm, mode = "undirected", weighted = TRUE)
SNnullperturb45c <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45cm <- as.matrix(SNnullperturb45c)
SNnullp45c <- graph_from_adjacency_matrix(SNnullp45cm, mode = "undirected", weighted = TRUE)

SNc45c <- cluster_louvain(SNp45c, weights = E(SNp45c)$weight)
SNnullc45c <- cluster_louvain(SNnullp45c, weights = E(SNnullp45c)$weight)

VICSN45[4] <- compare(SNc,SNc45c,method = "vi")
VICSNrandom45[4] <- compare(SNnullc, SNnullc45c, method = "vi")

rm(SNperturb45c,SNp45c,SNp45cm,SNnullperturb45c,SNnullp45cm,SNnullp45c,SNc45c,SNnullc45c)

#FOURTH TRIAL
SNperturb45d <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45dm <- as.matrix(SNperturb45d)
SNp45d <- graph_from_adjacency_matrix(SNp45dm, mode = "undirected", weighted = TRUE)
SNnullperturb45d <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45dm <- as.matrix(SNnullperturb45d)
SNnullp45d <- graph_from_adjacency_matrix(SNnullp45dm, mode = "undirected", weighted = TRUE)

SNc45d <- cluster_louvain(SNp45d, weights = E(SNp45d)$weight)
SNnullc45d <- cluster_louvain(SNnullp45d, weights = E(SNnullp45d)$weight)

VICSN45[5] <- compare(SNc,SNc45d,method = "vi")
VICSNrandom45[5] <- compare(SNnullc, SNnullc45d, method = "vi")

rm(SNperturb45d,SNp45d,SNp45dm,SNnullperturb45d,SNnullp45dm,SNnullp45d,SNc45d,SNnullc45d)

#FIFTH TRIAL
SNperturb45e <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45em <- as.matrix(SNperturb45e)
SNp45e <- graph_from_adjacency_matrix(SNp45em, mode = "undirected", weighted = TRUE)
SNnullperturb45e <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45em <- as.matrix(SNnullperturb45e)
SNnullp45e <- graph_from_adjacency_matrix(SNnullp45em, mode = "undirected", weighted = TRUE)

SNc45e <- cluster_louvain(SNp45e, weights = E(SNp45e)$weight)
SNnullc45e <- cluster_louvain(SNnullp45e, weights = E(SNnullp45e)$weight)

VICSN45[6] <- compare(SNc,SNc45e,method = "vi")
VICSNrandom45[6] <- compare(SNnullc, SNnullc45e, method = "vi")

rm(SNperturb45e,SNp45e,SNp45em,SNnullperturb45e,SNnullp45em,SNnullp45e,SNc45e,SNnullc45e)

#SIXTH TRIAL
SNperturb45f <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45fm <- as.matrix(SNperturb45f)
SNp45f <- graph_from_adjacency_matrix(SNp45fm, mode = "undirected", weighted = TRUE)
SNnullperturb45f <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45fm <- as.matrix(SNnullperturb45f)
SNnullp45f <- graph_from_adjacency_matrix(SNnullp45fm, mode = "undirected", weighted = TRUE)

SNc45f <- cluster_louvain(SNp45f, weights = E(SNp45f)$weight)
SNnullc45f <- cluster_louvain(SNnullp45f, weights = E(SNnullp45f)$weight)

VICSN45[7] <- compare(SNc,SNc45f,method = "vi")
VICSNrandom45[7] <- compare(SNnullc, SNnullc45f, method = "vi")

rm(SNperturb45f,SNp45f,SNp45fm,SNnullperturb45f,SNnullp45fm,SNnullp45f,SNc45f,SNnullc45f)

#SEVENTH TRIAL
SNperturb45g <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45gm <- as.matrix(SNperturb45g)
SNp45g <- graph_from_adjacency_matrix(SNp45gm, mode = "undirected", weighted = TRUE)
SNnullperturb45g <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45gm <- as.matrix(SNnullperturb45g)
SNnullp45g <- graph_from_adjacency_matrix(SNnullp45gm, mode = "undirected", weighted = TRUE)

SNc45g <- cluster_louvain(SNp45g, weights = E(SNp45g)$weight)
SNnullc45g <- cluster_louvain(SNnullp45g, weights = E(SNnullp45g)$weight)

VICSN45[8] <- compare(SNc,SNc45g,method = "vi")
VICSNrandom45[8] <- compare(SNnullc, SNnullc45g, method = "vi")

rm(SNperturb45g,SNp45g,SNp45gm,SNnullperturb45g,SNnullp45gm,SNnullp45g,SNc45g,SNnullc45g)

#EIGHTH TRIAL
SNperturb45h <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45hm <- as.matrix(SNperturb45h)
SNp45h <- graph_from_adjacency_matrix(SNp45hm, mode = "undirected", weighted = TRUE)
SNnullperturb45h <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45hm <- as.matrix(SNnullperturb45h)
SNnullp45h <- graph_from_adjacency_matrix(SNnullp45hm, mode = "undirected", weighted = TRUE)

SNc45h <- cluster_louvain(SNp45h, weights = E(SNp45h)$weight)
SNnullc45h <- cluster_louvain(SNnullp45h, weights = E(SNnullp45h)$weight)

VICSN45[9] <- compare(SNc,SNc45h,method = "vi")
VICSNrandom45[9] <- compare(SNnullc, SNnullc45h, method = "vi")

rm(SNperturb45h,SNp45h,SNp45hm,SNnullperturb45h,SNnullp45hm,SNnullp45h,SNc45h,SNnullc45h)

#NINTH TRIAL
SNperturb45i <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45im <- as.matrix(SNperturb45i)
SNp45i <- graph_from_adjacency_matrix(SNp45im, mode = "undirected", weighted = TRUE)
SNnullperturb45i <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45im <- as.matrix(SNnullperturb45i)
SNnullp45i <- graph_from_adjacency_matrix(SNnullp45im, mode = "undirected", weighted = TRUE)

SNc45i <- cluster_louvain(SNp45i, weights = E(SNp45i)$weight)
SNnullc45i <- cluster_louvain(SNnullp45i, weights = E(SNnullp45i)$weight)

VICSN45[10] <- compare(SNc,SNc45i,method = "vi")
VICSNrandom45[10] <- compare(SNnullc, SNnullc45i, method = "vi")

rm(SNperturb45i,SNp45i,SNp45im,SNnullperturb45i,SNnullp45im,SNnullp45i,SNc45i,SNnullc45i)

#TENTH TRIAL
SNperturb45j <- rewireR(SNadj, 56137, dist = "NegBinom")
SNp45jm <- as.matrix(SNperturb45j)
SNp45j <- graph_from_adjacency_matrix(SNp45jm, mode = "undirected", weighted = TRUE)
SNnullperturb45j <- rewireR(SNnullmatrix, 56137, dist = "NegBinom")
SNnullp45jm <- as.matrix(SNnullperturb45j)
SNnullp45j <- graph_from_adjacency_matrix(SNnullp45jm, mode = "undirected", weighted = TRUE)

SNc45j <- cluster_louvain(SNp45j, weights = E(SNp45j)$weight)
SNnullc45j <- cluster_louvain(SNnullp45j, weights = E(SNnullp45j)$weight)

VICSN45[11] <- compare(SNc,SNc45j,method = "vi")
VICSNrandom45[11] <- compare(SNnullc, SNnullc45j, method = "vi")

rm(SNperturb45j,SNp45j,SNp45jm,SNnullperturb45j,SNnullp45jm,SNnullp45j,SNc45j,SNnullc45j)

##############50% Perturbation Level#####################

#FIRST TRIAL
SNperturb50a <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50am <- as.matrix(SNperturb50a)
SNp50a <- graph_from_adjacency_matrix(SNp50am, mode = "undirected", weighted = TRUE)
SNnullperturb50a <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50am <- as.matrix(SNnullperturb50a)
SNnullp50a <- graph_from_adjacency_matrix(SNnullp50am, mode = "undirected", weighted = TRUE)

SNc50a <- cluster_louvain(SNp50a, weights = E(SNp50a)$weight)
SNnullc50a <- cluster_louvain(SNnullp50a, weights = E(SNnullp50a)$weight)

VICSN50[2] <- compare(SNc,SNc50a,method = "vi")
VICSNrandom50[2] <- compare(SNnullc, SNnullc50a, method = "vi")

rm(SNperturb50a,SNp50a,SNp50am,SNnullperturb50a,SNnullp50am,SNnullp50a,SNc50a,SNnullc50a)

#SECOND TRIAL
SNperturb50b <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50bm <- as.matrix(SNperturb50b)
SNp50b <- graph_from_adjacency_matrix(SNp50bm, mode = "undirected", weighted = TRUE)
SNnullperturb50b <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50bm <- as.matrix(SNnullperturb50b)
SNnullp50b <- graph_from_adjacency_matrix(SNnullp50bm, mode = "undirected", weighted = TRUE)

SNc50b <- cluster_louvain(SNp50b, weights = E(SNp50b)$weight)
SNnullc50b <- cluster_louvain(SNnullp50b, weights = E(SNnullp50b)$weight)

VICSN50[3] <- compare(SNc,SNc50b,method = "vi")
VICSNrandom50[3] <- compare(SNnullc, SNnullc50b, method = "vi")

rm(SNperturb50b,SNp50b,SNp50bm,SNnullperturb50b,SNnullp50bm,SNnullp50b,SNc50b,SNnullc50b)

#THIRD TRIAL
SNperturb50c <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50cm <- as.matrix(SNperturb50c)
SNp50c <- graph_from_adjacency_matrix(SNp50cm, mode = "undirected", weighted = TRUE)
SNnullperturb50c <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50cm <- as.matrix(SNnullperturb50c)
SNnullp50c <- graph_from_adjacency_matrix(SNnullp50cm, mode = "undirected", weighted = TRUE)

SNc50c <- cluster_louvain(SNp50c, weights = E(SNp50c)$weight)
SNnullc50c <- cluster_louvain(SNnullp50c, weights = E(SNnullp50c)$weight)

VICSN50[4] <- compare(SNc,SNc50c,method = "vi")
VICSNrandom50[4] <- compare(SNnullc, SNnullc50c, method = "vi")

rm(SNperturb50c,SNp50c,SNp50cm,SNnullperturb50c,SNnullp50cm,SNnullp50c,SNc50c,SNnullc50c)

#FOURTH TRIAL
SNperturb50d <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50dm <- as.matrix(SNperturb50d)
SNp50d <- graph_from_adjacency_matrix(SNp50dm, mode = "undirected", weighted = TRUE)
SNnullperturb50d <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50dm <- as.matrix(SNnullperturb50d)
SNnullp50d <- graph_from_adjacency_matrix(SNnullp50dm, mode = "undirected", weighted = TRUE)

SNc50d <- cluster_louvain(SNp50d, weights = E(SNp50d)$weight)
SNnullc50d <- cluster_louvain(SNnullp50d, weights = E(SNnullp50d)$weight)

VICSN50[5] <- compare(SNc,SNc50d,method = "vi")
VICSNrandom50[5] <- compare(SNnullc, SNnullc50d, method = "vi")

rm(SNperturb50d,SNp50d,SNp50dm,SNnullperturb50d,SNnullp50dm,SNnullp50d,SNc50d,SNnullc50d)

#FIFTH TRIAL
SNperturb50e <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50em <- as.matrix(SNperturb50e)
SNp50e <- graph_from_adjacency_matrix(SNp50em, mode = "undirected", weighted = TRUE)
SNnullperturb50e <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50em <- as.matrix(SNnullperturb50e)
SNnullp50e <- graph_from_adjacency_matrix(SNnullp50em, mode = "undirected", weighted = TRUE)

SNc50e <- cluster_louvain(SNp50e, weights = E(SNp50e)$weight)
SNnullc50e <- cluster_louvain(SNnullp50e, weights = E(SNnullp50e)$weight)

VICSN50[6] <- compare(SNc,SNc50e,method = "vi")
VICSNrandom50[6] <- compare(SNnullc, SNnullc50e, method = "vi")

rm(SNperturb50e,SNp50e,SNp50em,SNnullperturb50e,SNnullp50em,SNnullp50e,SNc50e,SNnullc50e)

#SIXTH TRIAL
SNperturb50f <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50fm <- as.matrix(SNperturb50f)
SNp50f <- graph_from_adjacency_matrix(SNp50fm, mode = "undirected", weighted = TRUE)
SNnullperturb50f <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50fm <- as.matrix(SNnullperturb50f)
SNnullp50f <- graph_from_adjacency_matrix(SNnullp50fm, mode = "undirected", weighted = TRUE)

SNc50f <- cluster_louvain(SNp50f, weights = E(SNp50f)$weight)
SNnullc50f <- cluster_louvain(SNnullp50f, weights = E(SNnullp50f)$weight)

VICSN50[7] <- compare(SNc,SNc50f,method = "vi")
VICSNrandom50[7] <- compare(SNnullc, SNnullc50f, method = "vi")

rm(SNperturb50f,SNp50f,SNp50fm,SNnullperturb50f,SNnullp50fm,SNnullp50f,SNc50f,SNnullc50f)

#SEVENTH TRIAL
SNperturb50g <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50gm <- as.matrix(SNperturb50g)
SNp50g <- graph_from_adjacency_matrix(SNp50gm, mode = "undirected", weighted = TRUE)
SNnullperturb50g <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50gm <- as.matrix(SNnullperturb50g)
SNnullp50g <- graph_from_adjacency_matrix(SNnullp50gm, mode = "undirected", weighted = TRUE)

SNc50g <- cluster_louvain(SNp50g, weights = E(SNp50g)$weight)
SNnullc50g <- cluster_louvain(SNnullp50g, weights = E(SNnullp50g)$weight)

VICSN50[8] <- compare(SNc,SNc50g,method = "vi")
VICSNrandom50[8] <- compare(SNnullc, SNnullc50g, method = "vi")

rm(SNperturb50g,SNp50g,SNp50gm,SNnullperturb50g,SNnullp50gm,SNnullp50g,SNc50g,SNnullc50g)

#EIGHTH TRIAL
SNperturb50h <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50hm <- as.matrix(SNperturb50h)
SNp50h <- graph_from_adjacency_matrix(SNp50hm, mode = "undirected", weighted = TRUE)
SNnullperturb50h <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50hm <- as.matrix(SNnullperturb50h)
SNnullp50h <- graph_from_adjacency_matrix(SNnullp50hm, mode = "undirected", weighted = TRUE)

SNc50h <- cluster_louvain(SNp50h, weights = E(SNp50h)$weight)
SNnullc50h <- cluster_louvain(SNnullp50h, weights = E(SNnullp50h)$weight)

VICSN50[9] <- compare(SNc,SNc50h,method = "vi")
VICSNrandom50[9] <- compare(SNnullc, SNnullc50h, method = "vi")

rm(SNperturb50h,SNp50h,SNp50hm,SNnullperturb50h,SNnullp50hm,SNnullp50h,SNc50h,SNnullc50h)

#NINTH TRIAL
SNperturb50i <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50im <- as.matrix(SNperturb50i)
SNp50i <- graph_from_adjacency_matrix(SNp50im, mode = "undirected", weighted = TRUE)
SNnullperturb50i <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50im <- as.matrix(SNnullperturb50i)
SNnullp50i <- graph_from_adjacency_matrix(SNnullp50im, mode = "undirected", weighted = TRUE)

SNc50i <- cluster_louvain(SNp50i, weights = E(SNp50i)$weight)
SNnullc50i <- cluster_louvain(SNnullp50i, weights = E(SNnullp50i)$weight)

VICSN50[10] <- compare(SNc,SNc50i,method = "vi")
VICSNrandom50[10] <- compare(SNnullc, SNnullc50i, method = "vi")

rm(SNperturb50i,SNp50i,SNp50im,SNnullperturb50i,SNnullp50im,SNnullp50i,SNc50i,SNnullc50i)

#TENTH TRIAL
SNperturb50j <- rewireR(SNadj, 62375, dist = "NegBinom")
SNp50jm <- as.matrix(SNperturb50j)
SNp50j <- graph_from_adjacency_matrix(SNp50jm, mode = "undirected", weighted = TRUE)
SNnullperturb50j <- rewireR(SNnullmatrix, 62375, dist = "NegBinom")
SNnullp50jm <- as.matrix(SNnullperturb50j)
SNnullp50j <- graph_from_adjacency_matrix(SNnullp50jm, mode = "undirected", weighted = TRUE)

SNc50j <- cluster_louvain(SNp50j, weights = E(SNp50j)$weight)
SNnullc50j <- cluster_louvain(SNnullp50j, weights = E(SNnullp50j)$weight)

VICSN50[11] <- compare(SNc,SNc50j,method = "vi")
VICSNrandom50[11] <- compare(SNnullc, SNnullc50j, method = "vi")

rm(SNperturb50j,SNp50j,SNp50jm,SNnullperturb50j,SNnullp50jm,SNnullp50j,SNc50j,SNnullc50j)

###############55% Perturbation Level################

#FIRST TRIAL
SNperturb55a <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55am <- as.matrix(SNperturb55a)
SNp55a <- graph_from_adjacency_matrix(SNp55am, mode = "undirected", weighted = TRUE)
SNnullperturb55a <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55am <- as.matrix(SNnullperturb55a)
SNnullp55a <- graph_from_adjacency_matrix(SNnullp55am, mode = "undirected", weighted = TRUE)

SNc55a <- cluster_louvain(SNp55a, weights = E(SNp55a)$weight)
SNnullc55a <- cluster_louvain(SNnullp55a, weights = E(SNnullp55a)$weight)

VICSN55[2] <- compare(SNc,SNc55a,method = "vi")
VICSNrandom55[2] <- compare(SNnullc, SNnullc55a, method = "vi")

rm(SNperturb55a,SNp55a,SNp55am,SNnullperturb55a,SNnullp55am,SNnullp55a,SNc55a,SNnullc55a)

#SECOND TRIAL
SNperturb55b <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55bm <- as.matrix(SNperturb55b)
SNp55b <- graph_from_adjacency_matrix(SNp55bm, mode = "undirected", weighted = TRUE)
SNnullperturb55b <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55bm <- as.matrix(SNnullperturb55b)
SNnullp55b <- graph_from_adjacency_matrix(SNnullp55bm, mode = "undirected", weighted = TRUE)

SNc55b <- cluster_louvain(SNp55b, weights = E(SNp55b)$weight)
SNnullc55b <- cluster_louvain(SNnullp55b, weights = E(SNnullp55b)$weight)

VICSN55[3] <- compare(SNc,SNc55b,method = "vi")
VICSNrandom55[3] <- compare(SNnullc, SNnullc55b, method = "vi")

rm(SNperturb55b,SNp55b,SNp55bm,SNnullperturb55b,SNnullp55bm,SNnullp55b,SNc55b,SNnullc55b)

#THIRD TRIAL
SNperturb55c <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55cm <- as.matrix(SNperturb55c)
SNp55c <- graph_from_adjacency_matrix(SNp55cm, mode = "undirected", weighted = TRUE)
SNnullperturb55c <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55cm <- as.matrix(SNnullperturb55c)
SNnullp55c <- graph_from_adjacency_matrix(SNnullp55cm, mode = "undirected", weighted = TRUE)

SNc55c <- cluster_louvain(SNp55c, weights = E(SNp55c)$weight)
SNnullc55c <- cluster_louvain(SNnullp55c, weights = E(SNnullp55c)$weight)

VICSN55[4] <- compare(SNc,SNc55c,method = "vi")
VICSNrandom55[4] <- compare(SNnullc, SNnullc55c, method = "vi")

rm(SNperturb55c,SNp55c,SNp55cm,SNnullperturb55c,SNnullp55cm,SNnullp55c,SNc55c,SNnullc55c)

#FOURTH TRIAL
SNperturb55d <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55dm <- as.matrix(SNperturb55d)
SNp55d <- graph_from_adjacency_matrix(SNp55dm, mode = "undirected", weighted = TRUE)
SNnullperturb55d <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55dm <- as.matrix(SNnullperturb55d)
SNnullp55d <- graph_from_adjacency_matrix(SNnullp55dm, mode = "undirected", weighted = TRUE)

SNc55d <- cluster_louvain(SNp55d, weights = E(SNp55d)$weight)
SNnullc55d <- cluster_louvain(SNnullp55d, weights = E(SNnullp55d)$weight)

VICSN55[5] <- compare(SNc,SNc55d,method = "vi")
VICSNrandom55[5] <- compare(SNnullc, SNnullc55d, method = "vi")

rm(SNperturb55d,SNp55d,SNp55dm,SNnullperturb55d,SNnullp55dm,SNnullp55d,SNc55d,SNnullc55d)

#FIFTH TRIAL
SNperturb55e <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55em <- as.matrix(SNperturb55e)
SNp55e <- graph_from_adjacency_matrix(SNp55em, mode = "undirected", weighted = TRUE)
SNnullperturb55e <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55em <- as.matrix(SNnullperturb55e)
SNnullp55e <- graph_from_adjacency_matrix(SNnullp55em, mode = "undirected", weighted = TRUE)

SNc55e <- cluster_louvain(SNp55e, weights = E(SNp55e)$weight)
SNnullc55e <- cluster_louvain(SNnullp55e, weights = E(SNnullp55e)$weight)

VICSN55[6] <- compare(SNc,SNc55e,method = "vi")
VICSNrandom55[6] <- compare(SNnullc, SNnullc55e, method = "vi")

rm(SNperturb55e,SNp55e,SNp55em,SNnullperturb55e,SNnullp55em,SNnullp55e,SNc55e,SNnullc55e)

#SIXTH TRIAL
SNperturb55f <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55fm <- as.matrix(SNperturb55f)
SNp55f <- graph_from_adjacency_matrix(SNp55fm, mode = "undirected", weighted = TRUE)
SNnullperturb55f <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55fm <- as.matrix(SNnullperturb55f)
SNnullp55f <- graph_from_adjacency_matrix(SNnullp55fm, mode = "undirected", weighted = TRUE)

SNc55f <- cluster_louvain(SNp55f, weights = E(SNp55f)$weight)
SNnullc55f <- cluster_louvain(SNnullp55f, weights = E(SNnullp55f)$weight)

VICSN55[7] <- compare(SNc,SNc55f,method = "vi")
VICSNrandom55[7] <- compare(SNnullc, SNnullc55f, method = "vi")

rm(SNperturb55f,SNp55f,SNp55fm,SNnullperturb55f,SNnullp55fm,SNnullp55f,SNc55f,SNnullc55f)

#SEVENTH TRIAL
SNperturb55g <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55gm <- as.matrix(SNperturb55g)
SNp55g <- graph_from_adjacency_matrix(SNp55gm, mode = "undirected", weighted = TRUE)
SNnullperturb55g <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55gm <- as.matrix(SNnullperturb55g)
SNnullp55g <- graph_from_adjacency_matrix(SNnullp55gm, mode = "undirected", weighted = TRUE)

SNc55g <- cluster_louvain(SNp55g, weights = E(SNp55g)$weight)
SNnullc55g <- cluster_louvain(SNnullp55g, weights = E(SNnullp55g)$weight)

VICSN55[8] <- compare(SNc,SNc55g,method = "vi")
VICSNrandom55[8] <- compare(SNnullc, SNnullc55g, method = "vi")

rm(SNperturb55g,SNp55g,SNp55gm,SNnullperturb55g,SNnullp55gm,SNnullp55g,SNc55g,SNnullc55g)

#EIGHTH TRIAL
SNperturb55h <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55hm <- as.matrix(SNperturb55h)
SNp55h <- graph_from_adjacency_matrix(SNp55hm, mode = "undirected", weighted = TRUE)
SNnullperturb55h <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55hm <- as.matrix(SNnullperturb55h)
SNnullp55h <- graph_from_adjacency_matrix(SNnullp55hm, mode = "undirected", weighted = TRUE)

SNc55h <- cluster_louvain(SNp55h, weights = E(SNp55h)$weight)
SNnullc55h <- cluster_louvain(SNnullp55h, weights = E(SNnullp55h)$weight)

VICSN55[9] <- compare(SNc,SNc55h,method = "vi")
VICSNrandom55[9] <- compare(SNnullc, SNnullc55h, method = "vi")

rm(SNperturb55h,SNp55h,SNp55hm,SNnullperturb55h,SNnullp55hm,SNnullp55h,SNc55h,SNnullc55h)

#NINTH TRIAL
SNperturb55i <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55im <- as.matrix(SNperturb55i)
SNp55i <- graph_from_adjacency_matrix(SNp55im, mode = "undirected", weighted = TRUE)
SNnullperturb55i <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55im <- as.matrix(SNnullperturb55i)
SNnullp55i <- graph_from_adjacency_matrix(SNnullp55im, mode = "undirected", weighted = TRUE)

SNc55i <- cluster_louvain(SNp55i, weights = E(SNp55i)$weight)
SNnullc55i <- cluster_louvain(SNnullp55i, weights = E(SNnullp55i)$weight)

VICSN55[10] <- compare(SNc,SNc55i,method = "vi")
VICSNrandom55[10] <- compare(SNnullc, SNnullc55i, method = "vi")

rm(SNperturb55i,SNp55i,SNp55im,SNnullperturb55i,SNnullp55im,SNnullp55i,SNc55i,SNnullc55i)

#TENTH TRIAL
SNperturb55j <- rewireR(SNadj, 68612, dist = "NegBinom")
SNp55jm <- as.matrix(SNperturb55j)
SNp55j <- graph_from_adjacency_matrix(SNp55jm, mode = "undirected", weighted = TRUE)
SNnullperturb55j <- rewireR(SNnullmatrix, 68612, dist = "NegBinom")
SNnullp55jm <- as.matrix(SNnullperturb55j)
SNnullp55j <- graph_from_adjacency_matrix(SNnullp55jm, mode = "undirected", weighted = TRUE)

SNc55j <- cluster_louvain(SNp55j, weights = E(SNp55j)$weight)
SNnullc55j <- cluster_louvain(SNnullp55j, weights = E(SNnullp55j)$weight)

VICSN55[11] <- compare(SNc,SNc55j,method = "vi")
VICSNrandom55[11] <- compare(SNnullc, SNnullc55j, method = "vi")

rm(SNperturb55j,SNp55j,SNp55jm,SNnullperturb55j,SNnullp55jm,SNnullp55j,SNc55j,SNnullc55j)

##############60% Perturbation Level#######################

#FIRST TRIAL
SNperturb60a <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60am <- as.matrix(SNperturb60a)
SNp60a <- graph_from_adjacency_matrix(SNp60am, mode = "undirected", weighted = TRUE)
SNnullperturb60a <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60am <- as.matrix(SNnullperturb60a)
SNnullp60a <- graph_from_adjacency_matrix(SNnullp60am, mode = "undirected", weighted = TRUE)

SNc60a <- cluster_louvain(SNp60a, weights = E(SNp60a)$weight)
SNnullc60a <- cluster_louvain(SNnullp60a, weights = E(SNnullp60a)$weight)

VICSN60[2] <- compare(SNc,SNc60a,method = "vi")
VICSNrandom60[2] <- compare(SNnullc, SNnullc60a, method = "vi")

rm(SNperturb60a,SNp60a,SNp60am,SNnullperturb60a,SNnullp60am,SNnullp60a,SNc60a,SNnullc60a)

#SECOND TRIAL
SNperturb60b <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60bm <- as.matrix(SNperturb60b)
SNp60b <- graph_from_adjacency_matrix(SNp60bm, mode = "undirected", weighted = TRUE)
SNnullperturb60b <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60bm <- as.matrix(SNnullperturb60b)
SNnullp60b <- graph_from_adjacency_matrix(SNnullp60bm, mode = "undirected", weighted = TRUE)

SNc60b <- cluster_louvain(SNp60b, weights = E(SNp60b)$weight)
SNnullc60b <- cluster_louvain(SNnullp60b, weights = E(SNnullp60b)$weight)

VICSN60[3] <- compare(SNc,SNc60b,method = "vi")
VICSNrandom60[3] <- compare(SNnullc, SNnullc60b, method = "vi")

rm(SNperturb60b,SNp60b,SNp60bm,SNnullperturb60b,SNnullp60bm,SNnullp60b,SNc60b,SNnullc60b)

#THIRD TRIAL
SNperturb60c <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60cm <- as.matrix(SNperturb60c)
SNp60c <- graph_from_adjacency_matrix(SNp60cm, mode = "undirected", weighted = TRUE)
SNnullperturb60c <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60cm <- as.matrix(SNnullperturb60c)
SNnullp60c <- graph_from_adjacency_matrix(SNnullp60cm, mode = "undirected", weighted = TRUE)

SNc60c <- cluster_louvain(SNp60c, weights = E(SNp60c)$weight)
SNnullc60c <- cluster_louvain(SNnullp60c, weights = E(SNnullp60c)$weight)

VICSN60[4] <- compare(SNc,SNc60c,method = "vi")
VICSNrandom60[4] <- compare(SNnullc, SNnullc60c, method = "vi")

rm(SNperturb60c,SNp60c,SNp60cm,SNnullperturb60c,SNnullp60cm,SNnullp60c,SNc60c,SNnullc60c)

#FOURTH TRIAL
SNperturb60d <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60dm <- as.matrix(SNperturb60d)
SNp60d <- graph_from_adjacency_matrix(SNp60dm, mode = "undirected", weighted = TRUE)
SNnullperturb60d <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60dm <- as.matrix(SNnullperturb60d)
SNnullp60d <- graph_from_adjacency_matrix(SNnullp60dm, mode = "undirected", weighted = TRUE)

SNc60d <- cluster_louvain(SNp60d, weights = E(SNp60d)$weight)
SNnullc60d <- cluster_louvain(SNnullp60d, weights = E(SNnullp60d)$weight)

VICSN60[5] <- compare(SNc,SNc60d,method = "vi")
VICSNrandom60[5] <- compare(SNnullc, SNnullc60d, method = "vi")

rm(SNperturb60d,SNp60d,SNp60dm,SNnullperturb60d,SNnullp60dm,SNnullp60d,SNc60d,SNnullc60d)

#FIFTH TRIAL
SNperturb60e <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60em <- as.matrix(SNperturb60e)
SNp60e <- graph_from_adjacency_matrix(SNp60em, mode = "undirected", weighted = TRUE)
SNnullperturb60e <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60em <- as.matrix(SNnullperturb60e)
SNnullp60e <- graph_from_adjacency_matrix(SNnullp60em, mode = "undirected", weighted = TRUE)

SNc60e <- cluster_louvain(SNp60e, weights = E(SNp60e)$weight)
SNnullc60e <- cluster_louvain(SNnullp60e, weights = E(SNnullp60e)$weight)

VICSN60[6] <- compare(SNc,SNc60e,method = "vi")
VICSNrandom60[6] <- compare(SNnullc, SNnullc60e, method = "vi")

rm(SNperturb60e,SNp60e,SNp60em,SNnullperturb60e,SNnullp60em,SNnullp60e,SNc60e,SNnullc60e)

#SIXTH TRIAL
SNperturb60f <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60fm <- as.matrix(SNperturb60f)
SNp60f <- graph_from_adjacency_matrix(SNp60fm, mode = "undirected", weighted = TRUE)
SNnullperturb60f <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60fm <- as.matrix(SNnullperturb60f)
SNnullp60f <- graph_from_adjacency_matrix(SNnullp60fm, mode = "undirected", weighted = TRUE)

SNc60f <- cluster_louvain(SNp60f, weights = E(SNp60f)$weight)
SNnullc60f <- cluster_louvain(SNnullp60f, weights = E(SNnullp60f)$weight)

VICSN60[7] <- compare(SNc,SNc60f,method = "vi")
VICSNrandom60[7] <- compare(SNnullc, SNnullc60f, method = "vi")

rm(SNperturb60f,SNp60f,SNp60fm,SNnullperturb60f,SNnullp60fm,SNnullp60f,SNc60f,SNnullc60f)

#SEVENTH TRIAL
SNperturb60g <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60gm <- as.matrix(SNperturb60g)
SNp60g <- graph_from_adjacency_matrix(SNp60gm, mode = "undirected", weighted = TRUE)
SNnullperturb60g <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60gm <- as.matrix(SNnullperturb60g)
SNnullp60g <- graph_from_adjacency_matrix(SNnullp60gm, mode = "undirected", weighted = TRUE)

SNc60g <- cluster_louvain(SNp60g, weights = E(SNp60g)$weight)
SNnullc60g <- cluster_louvain(SNnullp60g, weights = E(SNnullp60g)$weight)

VICSN60[8] <- compare(SNc,SNc60g,method = "vi")
VICSNrandom60[8] <- compare(SNnullc, SNnullc60g, method = "vi")

rm(SNperturb60g,SNp60g,SNp60gm,SNnullperturb60g,SNnullp60gm,SNnullp60g,SNc60g,SNnullc60g)

#EIGHTH TRIAL
SNperturb60h <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60hm <- as.matrix(SNperturb60h)
SNp60h <- graph_from_adjacency_matrix(SNp60hm, mode = "undirected", weighted = TRUE)
SNnullperturb60h <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60hm <- as.matrix(SNnullperturb60h)
SNnullp60h <- graph_from_adjacency_matrix(SNnullp60hm, mode = "undirected", weighted = TRUE)

SNc60h <- cluster_louvain(SNp60h, weights = E(SNp60h)$weight)
SNnullc60h <- cluster_louvain(SNnullp60h, weights = E(SNnullp60h)$weight)

VICSN60[9] <- compare(SNc,SNc60h,method = "vi")
VICSNrandom60[9] <- compare(SNnullc, SNnullc60h, method = "vi")

rm(SNperturb60h,SNp60h,SNp60hm,SNnullperturb60h,SNnullp60hm,SNnullp60h,SNc60h,SNnullc60h)

#NINTH TRIAL
SNperturb60i <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60im <- as.matrix(SNperturb60i)
SNp60i <- graph_from_adjacency_matrix(SNp60im, mode = "undirected", weighted = TRUE)
SNnullperturb60i <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60im <- as.matrix(SNnullperturb60i)
SNnullp60i <- graph_from_adjacency_matrix(SNnullp60im, mode = "undirected", weighted = TRUE)

SNc60i <- cluster_louvain(SNp60i, weights = E(SNp60i)$weight)
SNnullc60i <- cluster_louvain(SNnullp60i, weights = E(SNnullp60i)$weight)

VICSN60[10] <- compare(SNc,SNc60i,method = "vi")
VICSNrandom60[10] <- compare(SNnullc, SNnullc60i, method = "vi")

rm(SNperturb60i,SNp60i,SNp60im,SNnullperturb60i,SNnullp60im,SNnullp60i,SNc60i,SNnullc60i)

#TENTH TRIAL
SNperturb60j <- rewireR(SNadj, 74850, dist = "NegBinom")
SNp60jm <- as.matrix(SNperturb60j)
SNp60j <- graph_from_adjacency_matrix(SNp60jm, mode = "undirected", weighted = TRUE)
SNnullperturb60j <- rewireR(SNnullmatrix, 74850, dist = "NegBinom")
SNnullp60jm <- as.matrix(SNnullperturb60j)
SNnullp60j <- graph_from_adjacency_matrix(SNnullp60jm, mode = "undirected", weighted = TRUE)

SNc60j <- cluster_louvain(SNp60j, weights = E(SNp60j)$weight)
SNnullc60j <- cluster_louvain(SNnullp60j, weights = E(SNnullp60j)$weight)

VICSN60[11] <- compare(SNc,SNc60j,method = "vi")
VICSNrandom60[11] <- compare(SNnullc, SNnullc60j, method = "vi")

rm(SNperturb60j,SNp60j,SNp60jm,SNnullperturb60j,SNnullp60jm,SNnullp60j,SNc60j,SNnullc60j)

###############PLOTTING VI CURVES AND TESTING######################

#Adding the perturbation level as the first value of each vector 
VICSN0[1] <- 0
VICSN5[1] <- 5
VICSN10[1] <- 10
VICSN15[1] <- 15
VICSN20[1] <- 20
VICSN25[1] <- 25
VICSN30[1] <- 30
VICSN35[1] <- 35
VICSN40[1] <- 40
VICSN45[1] <- 45
VICSN50[1] <- 50
VICSN55[1] <- 55
VICSN60[1] <- 60

VICSNrandom0[1] <- 0
VICSNrandom5[1] <- 5
VICSNrandom10[1] <- 10
VICSNrandom15[1] <- 15
VICSNrandom20[1] <- 20
VICSNrandom25[1] <- 25
VICSNrandom30[1] <- 30
VICSNrandom35[1] <- 35
VICSNrandom40[1] <- 40
VICSNrandom45[1] <- 45
VICSNrandom50[1] <- 50
VICSNrandom55[1] <- 55
VICSNrandom60[1] <- 60

#Creating matrices containing all the VI values
VICSN <- cbind(VICSN0,VICSN5,VICSN10,VICSN15,VICSN20,VICSN25,VICSN30,VICSN35,VICSN40,VICSN45,VICSN50,VICSN55,VICSN60)
VICSNrandom <- cbind(VICSNrandom0,VICSNrandom5,VICSNrandom10,VICSNrandom15,VICSNrandom20,VICSNrandom25,VICSNrandom30,VICSNrandom35,VICSNrandom40,VICSNrandom45,VICSNrandom50,VICSNrandom55,VICSNrandom60)

#Computing the area under the VI curves for the of2 network and the random model
robinAUC(Out$SNNet, VICSN[2:11,1:13], VICSNrandom[2:11,1:13], measure = "vi", verbose = FALSE)
#Calculating the Bayes Factor using GP regression
robinGPTest(VICSN[2:11,1:13], VICSNrandom[2:11,1:13])
#Performing the ITP test
robinFDATest(Out$SNNet, VICSN[2:11,1:13], VICSNrandom[2:11,1:13], measure = "vi", legend = c("VICSN", "VICSNrandom"), verbose = FALSE)

#Setting up the variables to plot the VI curves
VICSNplot <- c(colMeans(VICSN[2:11,1:13]))
VICSNrandomplot <- c(colMeans(VICSNrandom[2:11,1:13]))
valuesSN <- c(VICSNplot, VICSNrandomplot)
perturbSN <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6)
groupSN <- c("VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSN","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom","VICSNrandom")
curvesSN <- data.frame(groupSN,perturbSN,valuesSN)

#Plotting the VI curves using ggplot2
ggplot(curvesSN, aes(x = perturbSN, y= valuesSN, colour = groupSN, fill = groupSN)) +
  geom_line() +
  geom_point(size = 2, shape = 21) +
  ggtitle("VI curves for super node representation") +
  xlab("Perturbation level") +
  ylab("VI") +
  labs(colour = "Model", fill = "Model")
