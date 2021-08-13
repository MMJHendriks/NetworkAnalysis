
###############################################################################
## Required packages
## Install them first if you have not done so already
###############################################################################

library(tidyverse)
library(igraph)
library(magrittr)

###############################################################################
## The simulated class data
###############################################################################

nodes <- readRDS("class_nei.rds")
vertices <- readRDS("Attributes.rds")

g <- nodes %>%
  rename("ID" = matrix..) %>%
  pivot_longer(cols = starts_with("V")) %>%
  filter(!is.na(value)) %>%
  dplyr::select(ID, value) %>%
  graph_from_data_frame(vertices = vertices)
###############################################################################
## Linear threshold model
###############################################################################

## Don't forget to set the thresholds, see below how to do it.

#LT <- function(g, activated, plot = TRUE,q1=q1,q2=q2,xi=xi) {  #q2=1, b=1
  
  ## Transform the graph to an undirected graph if it is undirected
  if (is_directed(g)) {
    g <- as.undirected(g)
    cat(paste0("The inserted (directed) graph is for simplicity transformed to a undirected graph.\n"))
  }
  
  if (is.null(E(g)$weight)) {
    stop("Specify the individual activation thresholds as V(g)$threshold")
  }
  
  ## If you want a plot, fix the layout of the plot.
  if (plot) {
    l <- layout.fruchterman.reingold(g)
  }
  
  ## Empty list for the output
  out <- list()
  
  ## Set all nodes to inactive  
  V(g)$active <- FALSE
  ## Except the activated nodes
  V(g)$active[activated] <- TRUE
  ## If changes occur, continue the while loop (set to true initially)
  changes <- TRUE
  ## Number of active nodes
  n_active <- sum(V(g)$active)
  ## Timepoint
  t <- 1
  
  while(changes) {
    ## Set changes to false
    changes <- FALSE
    ## And set the time parameter to two
    t <- t + 1
    
    if (plot) { ## If you want a plot
      ## Set the colors of the nodes to red if they are active
      V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
      ## And plot the graph
      plot(g, layout = l)
    }
    
    ## For all nodes that are not yet active
    for (v in V(g)[!V(g)$active]) {
      ## If the proportion of active neighbors is greater than your threshold, the node
      ## gets activated
      if(xi>q1){ 
        V(g)$active[v]<-TRUE
      }
      if (q2<xi &x1<q1){ 
        if (q2/q1> x1){ 
          V(g)$active[v]<-FALSE
        }
        else{ 
        V(g)$active[v]<-TRUE
        }
      } 
      if(q2<xi){ 
        V(g)$active[v]<-FALSE
      }
    } 
  #    if (sum(neighbors(g, v)$active) / max(length(neighbors(g, v)), 1) > V(g)$threshold[v]) {
  #      V(g)$active[v] <- TRUE
  #    }
  #  }
    ## The number of active nodes at timepoint t
    n_active[t] <- sum(V(g)$active)
    
    ## If the number of active nodes is still increasing, continue the while loop
    if (n_active[t] > n_active[t-1]) {
      changes <- TRUE
    }
  }
  
  ## The final version of the plot
  if (plot) {
    V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
    plot(g, layout = l)
  }
  # Return the output
  c(active = n_active)
}




###############################################################################
## Linear threshold model - greedy algorithm - only slightly different from
## the IC greedy algorithm
###############################################################################

LT_greedy <- function(g, k = 3, q2=q2,b=b){  #q2=0.1, b=1) 
  
  ## Make the graph undirected if it isn't already
  if (is_directed(g)) {
    g <- as.undirected(g)
    cat("To work with the LT model, the graph is transformed to an undirected graph.")
  }
  ## Empty object for the number of spreads
  spreads  <- NA
  ## No one is in the solution yet (just a placeholder)
  solution <- rep(FALSE, vcount(g))
  node_solution <- NA
  ## Ignore the nodes without any ties
  ignore <- degree(g, mode = "all") == 0
  
  ## Set a progress bar, so that you have any idea how long it will take
  pb <- txtProgressBar(min = 1, max = k * (vcount(g) - sum(ignore)) - (k*(k-1)/2), style = 3)
  counter <- 0
  
  ## For the number of nodes that you want in the final solution
  for (i in 1:k) {
    ## No best node yet
    best_node <- -1
    ## And also no spread
    best_spread <- -Inf
    ## Select the nodes that you want to loop over (ignore those that are)
    ## already chosen as best nodes, and the nodes without ties
    nodes <- V(g)[!(solution | ignore)]
    
    ## For all nodes, check how much other nodes will get infected
    ## Store the node with most influence in best_node, and its spread
    ## in the best_spread
    for (node in nodes) {
      LT_solution <- LT(g, c(which(solution), node), plot = F)
      spread <- LT_solution[length(LT_solution)]
      if (spread > best_spread) {
        best_spread <- spread
        best_node   <- node
      }
      ## Update the progress bar
      setTxtProgressBar(pb, value = (counter <- counter + 1))
    }
    ## Store the best node to true
    solution[best_node] <- TRUE
    node_solution[i] <- best_node
    ## And record the spread
    spreads[i] <- best_spread
  }
  cat("\n")
  ## Return the output
  c(nodes = node_solution, 
    spread = spreads)
}

## Youtube adoption threshold
#V(g)$threshold <- vertices$Th_Youtube
#LT_greedy(g, k = 3)

## Vegetarian recipe adoption threshold
#V(g)$threshold <- vertices$Th_Veg
#LT_greedy(g, k = 3)

## Paper adoption threshold
#V(g)$threshold <- vertices$Th_Paper
#LT_greedy(g, k = 3)




###########################################################

library(igraph)
library(network)
library(intergraph)
library(RColorBrewer)
library(tidyr)
library(dplyr)

nodes <- readRDS("class_nei.rds")
vertices <- readRDS("Attributes.rds")

g <- nodes %>%
  rename("ID" = matrix..) %>%
  pivot_longer(cols = starts_with("V")) %>%
  filter(!is.na(value)) %>%
  dplyr::select(ID, value) %>%
  graph_from_data_frame(vertices = vertices)


#vertices <- readRDS('C:/Users/evast/MasterADS/Network Analysis/Attributes.rds')
#edges <- readRDS('C:/Users/evast/MasterADS/Network Analysis/class_nei.rds')

#g <- edges %>%
#  rename("ID" = matrix..) %>%
#  pivot_longer(cols = starts_with("V")) %>%
#  filter(!is.na(value)) %>%
#  dplyr::select(ID, value) %>%
#  graph_from_data_frame(vertices = vertices)

## Vegetarian recipe adoption threshold
V(g)$threshold <- vertices$Th_Veg +0.2
V(g)$threshold
V(g)$threshold2 <- vertices$Th_Veg / 4
V(g)$threshold2

E(g)$weight <- rep(vertices$Share_Veg, times = degree(g, mode = "out"))
E(g)$weight

###############################################################################
## Linear threshold model
###############################################################################

## Don't forget to set the thresholds, see below how to do it.

LT <- function(g, activated, plot = TRUE) {
  ## Transform the graph to an undirected graph if it is undirected
  if (is_directed(g)) {
    g <- as.undirected(g)
    cat(paste0("The inserted (directed) graph is for simplicity transformed to a undirected graph.\n"))
  }
  if (is.null(E(g)$weight)) {
    stop("Specify the individual activation thresholds as V(g)$threshold")
  }
  if (plot) {
    l <- layout.fruchterman.reingold(g)
  }
  out <- list()  ## Empty list for the output
  ## Set all nodes to inactive  
  V(g)$active <- FALSE
  ## Except the activated nodes
  V(g)$active[activated] <- TRUE
  ## If changes occur, continue the while loop (set to true initially)
  changes <- TRUE

  n_active <- sum(V(g)$active)  ## Number of active nodes
 ## Timepoint
  t <- 1
  while(changes) {
    changes <- FALSE ## Set changes to false
    ## And set the time parameter 
    t <- t + 1
    if (plot) { 
      V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
      plot(g, layout = l)
    }

    ## For all nodes that are not yet active
    for (v in V(g)[!V(g)$active]) {
      xi <- (sum(neighbors(g, v)$active) / max(length(neighbors(g, v)), 1))
      q1 <- V(g)$threshold[v]
      q2 <- V(g)$threshold2[v]
      ## If the proportion of active neighbors is greater than your threshold, the node
      ## gets activated, also allow for a small probability of below threshold adoption
      if (xi > q1) { #xi>q1  --> 1
        V(g)$active[v] <- TRUE}
      if (q2 < xi & xi < q1) { # q2 < xi < q1 --> pi
        if (q2+((q1-q2)/2) > xi) {
          V(g)$active[v] <- FALSE
        } 
        else {
          V(g)$active[v] <- TRUE
        } 
      } 
      else { ## xi<q2 --> 0
        V(g)$active[v] <- FALSE
      }
    }
    ## The number of active nodes at timepoint t
    n_active[t] <- sum(V(g)$active)
    ## If the number of active nodes is still increasing, continue the while loop
    if (n_active[t] > n_active[t-1]) {
      changes <- TRUE
    }
  }
  if (plot) {
    V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
    plot(g, layout = l)
  }
  c(active = n_active)# Return the output
}


###
## Vegetarian recipe adoption threshold
V(g)$threshold <- vertices$Th_Veg +0.6
V(g)$threshold
V(g)$threshold2 <- vertices$Th_Veg / 2
V(g)$threshold2

E(g)$weight <- rep(vertices$Share_Veg, times = degree(g, mode = "out"))
E(g)$weight


LT <- function(g, activated, plot = TRUE) {
  
  ## Transform the graph to an undirected graph if it is undirected
  if (is_directed(g)) {
    g <- as.undirected(g)
    cat(paste0("The inserted (directed) graph is for simplicity transformed to a undirected graph.\n"))
  }
  
  if (is.null(E(g)$weight)) {
    stop("Specify the individual activation thresholds as V(g)$threshold")
  }
  
  ## If you want a plot, fix the layout of the plot.
  if (plot) {
    l <- layout.fruchterman.reingold(g)
  }
  
  ## Empty list for the output
  out <- list()
  
  ## Set all nodes to inactive  
  V(g)$active <- FALSE
  ## Except the activated nodes
  V(g)$active[activated] <- TRUE
  ## If changes occur, continue the while loop (set to true initially)
  changes <- TRUE
  ## Number of active nodes
  n_active <- sum(V(g)$active)
  ## Timepoint
  t <- 1
  
  while(changes) {
    ## Set changes to false
    changes <- FALSE
    ## And set the time parameter to two
    t <- t + 1
    
    if (plot) { ## If you want a plot
      ## Set the colors of the nodes to red if they are active
      V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
      ## And plot the graph
      plot(g, layout = l)
    }
    
    ## For all nodes that are not yet active
    for (v in V(g)[!V(g)$active]) {
      xi <- (sum(neighbors(g, v)$active) / max(length(neighbors(g, v)), 1))
      q1 <- V(g)$threshold[v]
      q2 <- V(g)$threshold2[v]
      
      ## If the proportion of active neighbors is greater than your threshold, the node
      ## gets activated
      if (xi > q1) { #xi>q1  --> 1
        V(g)$active[v] <- TRUE}
      if (q2 < xi & xi < q1) { # q2 < xi < q1 --> pi
          V(g)$active[v] <- FALSE
      } 
      else { ## xi<q2 --> 0
        V(g)$active[v] <- TRUE
      }
    }
    ## The number of active nodes at timepoint t
    n_active[t] <- sum(V(g)$active)
    
    ## If the number of active nodes is still increasing, continue the while loop
    if (n_active[t] > n_active[t-1]) {
      changes <- TRUE
    }
  }
  
  ## The final version of the plot
  if (plot) {
    V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
    plot(g, layout = l)
  }
  # Return the output
  c(active = n_active)
}


LT_greedy(g, k=3)
LT(g, 21)


x=c(0,1,2,3,4,5,6)
#y1 = c(1,28,52,70,75,76,76)  #threshold +0.2 threshold2 /4
#y1 = c(1,35,75,86,87,87)#threshold +0.4 threshold2 /2
#y2 = c(1,42,67,68,68,68)  #threshold +0.2 threshold2 /2
y1=c(1,42,85,98,98,98,98) #threshold +0.01 threshold2 /2
y2= c(1,30,61,87,91,93,93) #threshold +0.2 threshold2 /2
y3=c(1,14,36,58,74,85,89) #threshold +0.4 threshold2 /2
y4=c(1,5,21,48,61,64,67)#threshold +0.6 threshold2 /2
plot(x=x,y=y1, type='l', xlab="time", ylab="Number of active nodes", main="Diffusion curve Vegetarian recipe")
lines(x, y1, col='pink')
lines(x, y2, col='red')
lines(x, y3, col='purple')
lines(x, y4, col='blue')
legend(x = 0, y= 73, legend=c("q1 lowest", "q1 low", "q1 high", "q1 highest"),
       col=c("pink", "red", "purple", "blue"), lty=1:1, cex=0.8)



for (i in c(-1,0,0.1,0.5,1,1.5,2)){ 
  k=5
  V(g)$threshold<-vertices$Th_Veg
  x<-LT_greedy(g,k=k,q2=i,b=0.1)
  print(x)
  
  
  title=c('q=',i)
  plot(x=c(1:5), y=x[6:10], xlab="Number of K", ylab="Number of affected nodes", main=title)}








#######################

V(g)$threshold <- vertices$Th_Veg
E(g)$weight <- rep(vertices$Share_Veg, times = degree(g, mode = "out"))

LT_S <- function(g, activated, plot = TRUE) {
  
  ## Transform the graph to an undirected graph if it is undirected
  if (is_directed(g)) {
    g <- as.undirected(g)
    cat(paste0("The inserted (directed) graph is for simplicity transformed to a undirected graph.\n"))
  }
  
  if (is.null(E(g)$weight)) {
    stop("Specify the individual activation thresholds as V(g)$threshold")
  }
  
  ## If you want a plot, fix the layout of the plot.
  if (plot) {
    l <- layout.fruchterman.reingold(g)
  }
  
  ## Empty list for the output
  out <- list()
  
  ## Set all nodes to inactive  
  V(g)$active <- FALSE
  ## Except the activated nodes
  V(g)$active[activated] <- TRUE
  ## If changes occur, continue the while loop (set to true initially)
  changes <- TRUE
  ## Number of active nodes
  n_active <- sum(V(g)$active)
  ## Timepoint
  t <- 1
  
  while(changes) {
    ## Set changes to false
    changes <- FALSE
    ## And set the time parameter to two
    t <- t + 1
    
    if (plot) { ## If you want a plot
      ## Set the colors of the nodes to red if they are active
      V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
      ## And plot the graph
      plot(g, layout = l)
    }
    
    ## For all nodes that are not yet active
    for (v in V(g)[!V(g)$active]) {
      xi <- (sum(neighbors(g, v)$active) / max(length(neighbors(g, v)), 1))
      q1 <- V(g)$threshold[v]
      
      ## If the proportion of active neighbors is greater than your threshold, the node
      ## gets activated
     S_log<-(1/(1+exp(-(xi-q1)/b)))
     if (xi>q1){ 
       V(g)$active[v]<-TRUE
     }
     else if(S_log>q2){ 
       V(g)$active[v]<-S_log
     }
     else{ 
       V(g)$active[v]<-FALSE}
    }
    ## The number of active nodes at timepoint t
    n_active[t] <- sum(V(g)$active)
    
    ## If the number of active nodes is still increasing, continue the while loop
    if (n_active[t] > n_active[t-1]) {
      changes <- TRUE
    }
  }
  
  ## The final version of the plot
  if (plot) {
    V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
    plot(g, layout = l, type="l")
  }
  # Return the output
  c(active = n_active)
}



for (i in c(-1,0,0.1,0.5,1,1.5,2)){ 
  k=5
  V(g)$threshold<-vertices$Th_Veg
  x<-LT_greedy(g,k=k,q2=0.1,b=i)
  print(x)
  
  title=c('b=',i)
  plot(x=c(1:5), y=x[6:10], xlab="Number of K", ylab="Number of affected nodes", main=title)}




## Vegetarian recipe adoption threshold
V(g)$threshold <- vertices$Th_Veg -0.1
V(g)$threshold
V(g)$threshold2 <- vertices$Th_Veg / 2
V(g)$threshold2
LT_greedy(g, k=3)
LT(g, 68)


###############################


LT_s <- function(g, activated, plot = TRUE) {
  
  ## Transform the graph to an undirected graph if it is undirected
  if (is_directed(g)) {
    g <- as.undirected(g)
    cat(paste0("The inserted (directed) graph is for simplicity transformed to a undirected graph.\n"))
  }
  
  if (is.null(E(g)$weight)) {
    stop("Specify the individual activation thresholds as V(g)$threshold")
  }
  
  ## If you want a plot, fix the layout of the plot.
  if (plot) {
    l <- layout.fruchterman.reingold(g)
  }
  
  ## Empty list for the output
  out <- list()
  
  ## Set all nodes to inactive  
  V(g)$active <- FALSE
  ## Except the activated nodes
  V(g)$active[activated] <- TRUE
  ## If changes occur, continue the while loop (set to true initially)
  changes <- TRUE
  ## Number of active nodes
  n_active <- sum(V(g)$active)
  ## Timepoint
  t <- 1
  
  while(changes) {
    ## Set changes to false
    changes <- FALSE
    ## And set the time parameter to two
    t <- t + 1
    
    if (plot) { ## If you want a plot
      ## Set the colors of the nodes to red if they are active
      V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
      ## And plot the graph
      plot(g, layout = l)
    }
    
    ## For all nodes that are not yet active
    for (v in V(g)[!V(g)$active]) {
      xi <- (sum(neighbors(g, v)$active) / max(length(neighbors(g, v)), 1))
      q1 <- V(g)$threshold[v] -0.01
      b <- 0.1
      
      S_log <- 1 / (1 + exp(- (xi-q1)/b))
      
      ## If the proportion of active neighbors is greater than your threshold, the node
      ## gets activated
      if (S_log > q1) { 
        V(g)$active[v] <- TRUE 
      } 
      else {
        V(g)$active[v] <- FALSE
      } 
    }
    ## The number of active nodes at timepoint t
    n_active[t] <- sum(V(g)$active)
    
    ## If the number of active nodes is still increasing, continue the while loop
    if (n_active[t] > n_active[t-1]) {
      changes <- TRUE
    }
  }
  
  ## The final version of the plot
  if (plot) {
    V(g)$color <- ifelse(V(g)$active, "red", "lightblue")
    plot(g, layout = l)
  }
  # Return the output
  c(active = n_active)
}

LT_sgreedy <- function(g, k = 3) {
  
  ## Make the graph undirected if it isn't already
  if (is_directed(g)) {
    g <- as.undirected(g)
    cat("To work with the LT model, the graph is transformed to an undirected graph.")
  }
  ## Empty object for the number of spreads
  spreads  <- NA
  ## No one is in the solution yet (just a placeholder)
  solution <- rep(FALSE, vcount(g))
  node_solution <- NA
  ## Ignore the nodes without any ties
  ignore <- degree(g, mode = "all") == 0
  
  ## Set a progress bar, so that you have any idea how long it will take
  pb <- txtProgressBar(min = 1, max = k * (vcount(g) - sum(ignore)) - (k*(k-1)/2), style = 3)
  counter <- 0
  
  ## For the number of nodes that you want in the final solution
  for (i in 1:k) {
    ## No best node yet
    best_node <- -1
    ## And also no spread
    best_spread <- -Inf
    ## Select the nodes that you want to loop over (ignore those that are)
    ## already chosen as best nodes, and the nodes without ties
    nodes <- V(g)[!(solution | ignore)]
    
    ## For all nodes, check how much other nodes will get infected
    ## Store the node with most influence in best_node, and its spread
    ## in the best_spread
    for (node in nodes) {
      LT_solution <- LT_s(g, c(which(solution), node), plot = F) ## !!
      spread <- LT_solution[length(LT_solution)]
      if (spread > best_spread) {
        best_spread <- spread
        best_node   <- node
      }
      ## Update the progress bar
      setTxtProgressBar(pb, value = (counter <- counter + 1))
    }
    ## Store the best node to true
    solution[best_node] <- TRUE
    node_solution[i] <- best_node
    ## And record the spread
    spreads[i] <- best_spread
  }
  cat("\n")
  ## Return the output
  c(nodes = node_solution, 
    spread = spreads)
}

V(g)$threshold <- vertices$Th_Veg
V(g)$threshold

E(g)$weight <- rep(vertices$Share_Veg, times = degree(g, mode = "out"))
E(g)$weight



LT_sgreedy(g, k=2)
LT_s(g, c(8))
LT_s(g, c(15))


