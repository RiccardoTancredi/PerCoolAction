library(igraph) |> suppressPackageStartupMessages()
library(ggplot2) |> suppressPackageStartupMessages()

# Main to generate a network with a multiplicative process and preferential attachment

# Model 1

N_it <- 3  # Number of iterations
m <- 2    # Multiplicative factor 


# Initialize the star graph
g <- make_graph(~ 1--2, 1--3, 1--4, 1--5)
# plot(g, 
#      vertex.label=NA, 
#      vertex.color="green",
#      vertex.size=2,
#      edge.color="red",
#      edge.width=1)

#

# Generate the network with a multiplicative process and preferential attachment
for (t in 1:N_it) {
  for (node in seq(1:vcount(g))) {
    vcount(g)
    node_degree <- degree(g, v = node)
    mk = m*node_degree
    for( j in seq(1:mk)){
      g <- add.vertices(g, nv = 1)
      new_node <- vcount(g)
      g <- add_edges(g, edges = c(node, new_node ))
    }
  }
}

plot(g, 
     vertex.label=NA, 
     vertex.color="green",
     vertex.size=2,
     edge.color="red",
     edge.width=1)




# Model 2

N_it <- 2 # Number of iterations
m <- 2    # Multiplicative factor 
nn_list <- list()
nn_list_new <- list()
# Initialize the star graph
g <- make_graph(~ 1--2, 1--3, 1--4, 1--5)


# Generate the network with a multiplicative process and preferential attachment

for (t in 1:N_it) {
  states <- logical(vcount(g))
  states <- TRUE
  V(g)$state <- states
  degree_sequence <- degree(g)
  for( s in (1:vcount(g) ) ) {
    nn_list[[length(nn_list) + 1]] <- as_ids(neighbors(g, s, mode = c( "all")))
  }
  g <- delete_edges(g, E(g))
 
  
  # loop su tutti i nodi (mi serve per sapere il degree)
  for (node in (1:vcount(g))) {
    mk = m*degree_sequence[node]
    # loop sui figli generati 
    for( j in (1:mk)){
      g <- add.vertices(g, nv = 1, state=FALSE, name = vcount(g) + 1, busy=FALSE)
      new_node <- vcount(g)
      g <- add_edges(g, edges = c(node, new_node )) 
    }
  }
   degree_sequence_new <- degree(g)
# genero una lista di vicini aggiornata
  for( d in (1:vcount(g) ) ) {
    nn_list_new[[length(nn_list_new) + 1]] <- as_ids(neighbors(g, d, mode = c( "all")))
  }
  # loop su tutti i nodi per vedere quali sono quelli già esistenti
  for( l in seq(1:vcount(g))){
    if(V(g)$state[l] == TRUE) {
      # scorre sui quelli che erano i vicini di l allo step precedente
      for( p in nn_list[[l]]) { 
        #Qua sotto dovrei collegare tanti figli di l quanto era prima il degree di l
        shortest_path <- shortest_paths(g, from = l, to = p)
        if (length(shortest_path$vpath[[1]]) != 4) {
          p <- as.integer(p)
          random_son_l <- sample((1:degree_sequence_new[l]), 1) 
          while (V(g)$busy[as.integer(nn_list_new[[l]][random_son_l])] == TRUE) {
          random_son_l <- sample((1:degree_sequence_new[l]), 1) 
          }
          son_l <- nn_list_new[[l]][random_son_l]
          V(g)$busy[as.integer(nn_list_new[[l]][random_son_l])] <- TRUE
          random_son_p <- sample((1:degree_sequence_new[p]), 1) 
          while (V(g)$busy[as.integer(nn_list_new[[p]][random_son_p])] == TRUE) {
           random_son_p <- sample((1:degree_sequence_new[p]), 1) 
          }
          V(g)$busy[as.integer(nn_list_new[[p]][random_son_p])] <- TRUE
          son_p <- nn_list_new[[p]][random_son_p]
          cat(son_p, "\n")
          g <- add_edges(g, edges = c(son_l, son_p))
        }
      }
    }
  }
  nn_list <- NULL
  nn_list_new <- NULL
}

plot(g, 
     vertex.label=NA, 
     vertex.color="green",
     vertex.size=2,
     edge.color="red",
     edge.width=1)




# Auxiliary graph 

lb_max <- 5
N_box <- c()
for (lb in 2:lb_max){
  # auxiliary graph 
  h <- delete_edges(g, E(g))
  for (v1 in (1:vcount(h))){
    nn <- neighbors(graph = h, v = v1)
    for (v2 in (1:vcount(h))){
      shortest_path <- shortest_paths(g, from = v1, to = v2)$vpath[[1]]
      dist_ij <- length(shortest_path)
      
      if((dist_ij > lb) && !(v2 %in% nn )) {
        h <- add_edges(h, edges = c(v1, v2))
      }
    }
  }

  # plot(h)


  # Box covering: greedy-color 
  color <- rep(-1, vcount(h))
  B <- lapply(1:vcount(h), function(i) neighborhood(h, order = 1, nodes = i))
  color[1] <- 0
  palette <- 0:(vcount(h)-1)

  for (v in V(g)){
    v <- as.numeric(v)
    if (v == 1){
      next
    }
    else{
      neig <- B[[v]][[1]][-1]
      col_neig <- color[neig] # (-1,-1,-1)
      minimum <- min(palette[!(palette %in% col_neig)])  
      color[v] <- minimum
    }
  }
  N_box <- append(N_box, max(color) + 1)
}
cat("Number of boxes =", N_box, "\n")

# linear fit:
# N_b(l_B)/N ~ l_b^(-d_b)
# log(N_b(l_B)/N) ~ -d_b*log(l_b)
y <- log(N_box/vcount(g))
x <- log(2:lb_max)
model <- lm(y ~ x)

summary(model)

pred <- predict(model, data.frame(x=x))
plot(x, y)
abline(model, col="red", lty=2, lw=4)
