#install.packages("igraph")
#install.packages("network")
#install.packages("intergraph")
#install.packages("RColorBrewer")

library(igraph)
library(network)
library(intergraph)
library(RColorBrewer)

#EXERCISE 1 Build and analyse a small network from facebook 
#getwd()
setwd("~/ADS/Neural Network/Part2_social")

nodes <- read.csv("Facebook_att.csv", header = TRUE)
links <- read.csv("Facebook_edge.csv", header = TRUE)
facebook <- graph_from_data_frame(d = links, vertices = nodes$NodeID, directed = FALSE)
facebook_net <- asNetwork(facebook)

# add attributes of vertexes (nodes) using igraph 
# Load in other attributes.
facebook_net %v% "vertex.names"        <- nodes$vertex.names
facebook_net %v% "sex"                 <- nodes$sex
facebook_net %v% "friend_count"        <- nodes$friend_count
facebook_net %v% "group"               <- nodes$group
facebook_net %v% "mutual_friend_count" <- nodes$mutual_friend_count
facebook_net %v% "na"                  <- nodes$na
facebook_net %v% "relationship_status" <- nodes$relationship_status

summary(facebook_net)
plot(facebook_net, vertex.cex = 1.2, main = "Basic plot of Douglas's Facebook friends")
#question 1

#Measure node-level metrics 
degree(facebook, v = 1, mode = "total") #friend nr 1
closeness(facebook, v = 1, normalized = TRUE) 
  #you get a warning message, due to the fact that the graph is disconnected
betweenness(facebook, v = 1, directed = FALSE, normalized = TRUE)
eigen_centrality(facebook)$vector[1]
transitivity(facebook, type = "localundirected", vids = 1)

deg <- degree(facebook, mode = "all")
cls <- closeness(facebook, normalized = TRUE)
btw <- betweenness(facebook, directed = FALSE, normalized = TRUE)
eig <- eigen_centrality(facebook)
lcl <- transitivity(facebook, type = "local")

centrality_v1<- (cbind(mean(deg),mean(cls),mean(btw)))
centrality_v1
sqrt(var(btw))
#question 2 vertex 1 in the network

head(sort(deg, decreasing = TRUE),6)
head(sort(btw, decreasing = TRUE),5)
head(sort(cls, decreasing = TRUE),5)
head(sort(eig$vector, decreasing = TRUE),5)

#question 3
par(mfrow = c(2,2))
plot(deg, cls, main = "Degree versus closeness", 
     xlab = "Degree", ylab = "Closeness") 
plot(deg, btw, main = "Degree versus betweenness",
     xlab = "Degree", ylab = "Betweenness")
plot(deg, eig$vector, main = "Degree versus eigenvector",
     xlab = "Degree", ylab = "Eigenvector")
plot(deg, lcl, main = "Degree versus local clustering",
     xlab = "Degree", ylab = "Local clustering")


#measure group-level metrics
deg_dist <- degree_distribution(facebook)
barplot(deg_dist)
edge_density(facebook)
diameter(facebook, directed = FALSE)
transitivity(facebook, type = "global")
centr_degree(facebook)$centralization

### Detect the component and community 
components(facebook)
components_fb <-components(facebook)
barplot(components_fb$csize, ylab = " csize frequency", xlab = "no number of cluster")

# node attributes identifying the different groups 
group <- as.factor(get.vertex.attribute(facebook_net, attrname = 'group'))

pal <- brewer.pal(nlevels(group), "Set1")

plot(facebook_net, vertex.col = pal[group], vertex.cex = 1.5, 
     main = "Plot of Facebook Data colored by friend type")
legend(x = "bottomleft", legend = levels(group), col = pal, 
       pch = 19, pt.cex = 1.2, bty = "n", title = "Friend type")

#sex
group_sex <- as.factor(get.vertex.attribute(facebook_net, attrname = 'sex'))

pal <- brewer.pal(nlevels(group_sex), "Set1")

plot(facebook_net, vertex.col = pal[group_sex], vertex.cex = 1.5, 
     main = "Plot of Facebook Data colored by friend type")
legend(x = "bottomleft", legend = levels(group_sex), col = pal, 
       pch = 19, pt.cex = 1.2, bty = "n", title = "Friend type")
#relationship status 
group_rela <- as.factor(get.vertex.attribute(facebook_net, attrname = 'relationship_status'))

pal <- brewer.pal(nlevels(group_rela), "Set1")

plot(facebook_net, vertex.col = pal[group_rela], vertex.cex = 1.5, 
     main = "Plot of Facebook Data colored by friend type")
legend(x = "bottomleft", legend = levels(group_rela), col = pal, 
       pch = 19, pt.cex = 1.2, bty = "n", title = "Friend type")

#denisty of each group of Douglas friends
sapply(levels(group), function(x) {
  y <- get.inducedSubgraph(facebook_net, 
                           which(facebook_net %v% "group" == x))
  paste0("Density for ", x, " friends is ", 
         edge_density(asIgraph(y)))
})


############################
#exercise 2 

cw <- cluster_walktrap(facebook)
plot(cw, facebook, vertex.label = V(facebook)$group, 
     main = "Walktrap")
ceb <- cluster_edge_betweenness(facebook)
plot(ceb, facebook, vertex.label = V(facebook)$group, 
     main = "Edge Betweenness")
cfg <- cluster_fast_greedy(facebook)
plot(cfg, facebook, vertex.label = V(facebook)$group,
     main = "Fast Greedy")
clp <- cluster_label_prop(facebook)
plot(clp, facebook, vertex.label = V(facebook)$group,
     main = "Label Prop")

ER <- sample_gnp(100, 1/100)
plot(ER, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "ER Random Network: G(N,p) model")
ER_10 <- sample_gnp(100, 10/100)
plot(ER_10, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "ER Random Network: G(N,p) model")
ER_5 <- sample_gnp(100, 5/100)
plot(ER_5, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "ER Random Network: G(N,p) model")

ER_100 <- sample_gnp(100, 10/100)
plot(ER_100, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "ER Random Network: G(N,p) model 100")
ER_1000 <- sample_gnp(1000, 10/100)
plot(ER_1000, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "ER Random Network: G(N,p) model 1000")
ER_500 <- sample_gnp(500, 10/100)
plot(ER_500, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "ER Random Network: G(N,p) model ")
ER_200 <- sample_gnp(200, 10/100)
ER_2000<- sample_gnp(2000, 10/100)

transitivity(ER_100)
transitivity(ER_200)
transitivity(ER_500)
transitivity(ER_1000)
transitivity(ER_2000)


### Small world model 
Regular<-watts.strogatz.game(dim=1,size=300,nei=6, p=0)
plot(Regular, layout=layout.circle, vertex.label=NA, vertex.size=5, main= "Network with zero rewiring probability ")
SW1<-watts.strogatz.game(dim=1,size=300,nei=6, p=0.001)
plot(SW1, layout=layout.circle, vertex.label=NA, vertex.size=5, main= "Network with 0.001 rewiring probability ")
SW2<-watts.strogatz.game(dim=1,size=300,nei=6, p=0.01)
plot(SW2, layout=layout.circle, vertex.label=NA, vertex.size=5, main= "Network with 0.01 rewiring probability ")
SW3<-watts.strogatz.game(dim=1,size=300,nei=6, p=0.1)
plot(SW3, layout=layout.circle, vertex.label=NA, vertex.size=5, main= "Network with 0.1 rewiring probability ")

#Q 17
dia <- c()
for (x in 20:500){
  result <- sample_gnp(x, 0.001)
  dia <- c(dia,diameter(result))
}

results <- data.frame(x=seq(20, 500))
results$diam <- dia
results$logdia = log(results$diam)
plot(x=results$x, y=results$logdia, xlab = "N", ylab= "Log(diameter)")

transitivity(Regular)
transitivity(SW1)
transitivity(SW2)
transitivity(SW3)
SW2_5<-watts.strogatz.game(dim=1,size=300,nei=6, p=0.02)
transitivity(SW2_5)

set.seed(1)
avg.stat <- function(nei, p) {
  result <- replicate(1000, {
    wsg <- watts.strogatz.game(1, 300, nei, p)
    c(average.path.length(wsg),
      transitivity(wsg))
  })
  apply(result, 1, quantile, probs = c(0.5, 0.05, 0.95))
}
nei <- 6
p <- 2^-seq(0, 10, len = 21)

result <- sapply(p, avg.stat, nei = nei)
result <- t(result / rep(avg.stat(nei, 0)[1,], each = 3))
par(mar=c(3.2, 2, 0.2, 0.2), mgp=c(2, 1, 0))
matplot(p, result, type = "l", log = "x", xaxt = "n", ylab = "",lty = rep(c(1,2,2),2), col=rep(c(1,2), each=3))
axis(1, at = 2^- (0:10), labels =  c(1, parse(text = paste(2, 1:10, sep ="^-",collapse = ";"))))
legend("bottomleft", c("average path length", "clustering coefficient"),lty = 1, col = c(1, 2))


g0 <- barabasi.game(100, power = 0.5, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE, zero.appeal = 1, directed = FALSE,algorithm ="psumtree", start.graph = NULL)
plot(g0, vertex.label= NA, edge.arrow.size=0.02,vertex.size =5, main = "Scale-free network model, power=0.5")

g1 <- barabasi.game(100, power = 0.5, m = NULL, out.dist = NULL, out.seq = NULL, out.pref = FALSE, zero.appeal = 1, directed = FALSE,algorithm ="psumtree", start.graph = NULL)
g1Net<-asNetwork(g1)
VS = 3+ 0.5*degree(g1)
plot(g1, vertex.label= NA, edge.arrow.size=0.02,vertex.size =VS, main = "Scale-free network model, power=0.5")

# Q 23
setwd("~/ADS/Neural Network/Part2_social")
getwd()
snap_data <- readRDS('C:\\Users\\mirth\\OneDrive\\Documenten\\ADS\\Neural Network\\Part2_social\\Attributes.rds')
Amaz_data<- readRDS('C:/Users/mirth/OneDrive/Documenten/ADS/Neural Network/Part2_social/amazon_edge.rds')
snap <- graph_from_data_frame(d = Amaz_data, vertices = snap_data$NodeID, directed = FALSE)
snap_net<-asNetwork(snap)
summary(snap_net)

#plot(snap_net, vertex.cex = 1.2, main = "plot of snap net")

transitivity(snap)
#diameter(snap, directed = FALSE)
edge_density(snap)
centr_degree(snap)$centralization
centr_degree(snap)
eigen_centrality(snap)$vector[1]

deg_dist_a <- degree_distribution(snap)
barplot(deg_dist_a)

#########################################################################################
########################################################################################
#########################################################################################
#exercise 3
setwd("~/ADS/Neural Network/Part2_social")
#install.packages("readxl")
library(readxl)
class_net <- read_excel("C:\\Users\\mirth\\OneDrive\\Documenten\\ADS\\Neural Network\\Part2_social\\Class_network_survey.xlsx")
class_nodes <- readRDS("C:\\Users\\mirth\\OneDrive\\Documenten\\ADS\\Neural Network\\Part2_social\\Class_attributes.rds")
class_edge <- readRDS('C:\\Users\\mirth\\OneDrive\\Documenten\\ADS\\Neural Network\\Part2_social\\class_nei.rds')

class(class_nodes)

class_surv<-graph_from_data_frame(d = class_edge, vertices = class_nodes$ID, directed = FALSE)
class_s<-asNetwork(class_surv)

summary(class_s)

plot(class_s, vertex.cex = 1.2, main = "Class survey network")

deg<-degree(class_surv, mode = "all")
mean_deg<-mean(deg)

share_paper<- mean(class_nodes$Share_Paper)
share_veg<- mean(class_nodes$Share_Veg)
share_youtube<-mean(class_nodes$Share_Youtube)
share<-cbind(share_paper, share_veg, share_youtube)
barplot(share, deg)
plot(deg, class_nodes$Share_Paper)
plot(class_nodes$Share_Paper, deg)
plot(deg, share)

cShare_Paper<-cor(deg, class_nodes$Share_Paper)
cs_veg<-cor(deg, class_nodes$Share_Veg)
CShare_Yout<-cor(deg, class_nodes$Share_Youtube)

cTh_paper<-cor(deg, class_nodes$Th_Paper)
cTh_veg<-cor(deg, class_nodes$Th_Veg)
cTh_Yout<-cor(deg, class_nodes$Th_Youtube)
cor_deg<-cbind(cShare_Paper,cs_veg, CShare_Yout, cTh_paper,cTh_veg,cTh_Yout)
cor_deg

#install.packages("tidyr")
library(tidyr)
#install.packages("dplyr")
library(dplyr)
g <- class_edge %>%
  rename("ID" = matrix..) %>%
  pivot_longer(cols = starts_with("V")) %>%
  filter(!is.na(value)) %>%
  dplyr::select(ID, value) %>%
  graph_from_data_frame(vertices = class_nodes)

deg<-degree(g, mode = "all")

deg_dist <- degree_distribution(g)
deg <- degree(g)
sort(deg, decreasing=TRUE)
edge_density(g)
diameter(g, directed = FALSE)
transitivity(g, type = "global")
centr_degree(g)$centralization
cen<-centralization.degree(g)

cls <- closeness(g, normalized = TRUE)
btw <- betweenness(g, directed = FALSE, normalized = TRUE)
sort(cls, decreasing=FALSE)
sort(btw, decreasing=FALSE)

nodes_zonder <- class_nodes[class_nodes$ID != 68, ]
edges_zonder <- class_edge[class_edge$matrix.. != 68, ]
g_zonder <- edges_zonder %>%
  rename("ID" = matrix..) %>%
  pivot_longer(cols = starts_with("V")) %>%
  filter(!is.na(value)) %>%
  dplyr::select(ID, value) %>%
  graph_from_data_frame(vertices = nodes_zonder)

degree(g_zonder)
edge_density(g_zonder)
diameter(g_zonder, directed = FALSE)
transitivity(g_zonder, type = "global")
centr_degree(g_zonder)$centralization


g_net<-asNetwork(g)

as.matrix.network(g_net)
delete.edges(g_net,68)              #Remove an edge
as.matrix.network(g_net)
delete.vertices(g_net,68)  #Remove a vertex
as.matrix.network(g_net)

g_graph<- asIgraph(g_net)
class(g_graph)

deg_g<-degree(g_graph)
sort(deg_g, decreasing=FALSE)

edge_density(g_graph)
diameter(g_graph, directed = FALSE)
transitivity(g_graph, type = "global")
centr_degree(g_graph)$centralization



