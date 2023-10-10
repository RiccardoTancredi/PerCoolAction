library(igraph)


get_graph <- function(model, only_bcc = T, add_spins = T, param1 = NA, layout_design = T) {
  name <- NULL
  # Erdos-Renyi
  if (model == "ER" || model == "er" || model == "erdos") {
    P <- param1 / N
    graph <- graph <- erdos.renyi.game(N, P, type = "gnp")
    name <- "Erdos_Renyi"
    nickname <- "er"
  } else if (model == "BA" || model == "ba" || model == "scalefree") {
    # Barabasi-Albert model (scale free) with a power that is x^-power
    graph <- barabasi.game(n = N, directed = F, power = param1)
    name <- "Barabasi-Albert"
    nickname <- "ba"
    graph$m0 <- 1
  } else if (model == "fully" || model == "full") {
    graph <- make_full_graph(N, directed = FALSE, loops = FALSE)
    name <- "Fully Connected"
    nickname <- "fc"
  } else {
    cat("Wrong model name\n")
    return()
  }

  # Fully connected graph
  # graph <- make_full_graph(N, directed = FALSE, loops = FALSE)

  # Generate a small-world network using the Watts-Strogatz model
  # k <- 4       # Number of neighbors to connect to initially
  # p <- 0.5/N     # Probability of rewiring edges
  # graph <- watts.strogatz.game(1, N, k, p)

  graph$name <- name
  graph$nickname <- nickname
  graph$param <- param1
  graph$N <- N
  if (only_bcc) {
    graph <- keep_only_bcc(graph)
  }

  if (add_spins) {
    graph$spins <- generate_spins(N)
  }
  # common things
  if (layout_design) {
    layout_ <- layout_nicely(graph)
    graph$layout_design <- layout_
  }
  return(graph)
}

generate_spins <- function(n) {
  spin_values <- sample(c(1, -1), n, replace = TRUE)
  return(spin_values)
}

generate_colors <- function(spins, pos_color = "green", neg_color = "red") {
  node_colors <- ifelse(spins == 1, pos_color, neg_color)
  return(node_colors)
}

get_bcc_size <- function(graph) {
  size <- max(components(graph)$csize)
  return(size)
}

keep_only_bcc <- function(graph) {
  clts <- clusters(graph)
  best_clst <- which.max(clts$csize)
  nodes <- clusters(graph)$membership != best_clst
  graph <- delete_vertices(graph, nodes)
  return(graph)
}

get_Tc <- function(graph) {
  k1 <- get_avg_degree(graph)
  k2 <- get_moment_2(graph)
  return(1 / (-0.5 * log(1 - 2 * k1 / k2)))
}

get_magnetization <- function(graph) {
  # clts <- clusters(graph)


  # magnets <- rep(0, clts$no)
  # for (i in 1:clts$no) {
  #   # # Se la dimnesione del gruppo è 1 (nodo singolo), allora non conta per la magnetizzazione
  #   # if (clts$csize[i] == 1) {
  #   #   next
  #   # }

  #   # indice di T o F in base se l'indice del cluster a cui appartengono è uguale ad i
  #   nodes <- clusters(graph)$membership == i

  #   # Sommiamo gli spins dei nodi che soddisfano la condizione sopra
  #   magnets[i] <- sum(graph$spins[nodes])
  #   magnets[i] <- abs(magnets[i])

  #   # cat(i, " ", magnets[i], "\n")
  # }
  # # cat("magnet: ", sum(magnets), "\n")
  # return(sum(magnets))

  clts <- clusters(graph)
  best_clst <- which.max(clts$csize)
  # indice di T o F in base se l'indice del cluster a cui appartengono è uguale ad best_clst
  nodes <- clusters(graph)$membership == best_clst
  magnets <- sum(graph$spins[nodes])
  # cat(magnets, graph$spins[nodes], "\n")
  return(magnets)
}

get_energy <- function(graph) {
  n <- vcount(graph)
  M <- get_avg_degree(graph) * n / 2
  J <- 1
  # cat("M:", M, "avg_deg:", avg_degree, "\n")

  allsum <- 0
  for (edge in E(graph)) {
    nodes <- ends(graph, edge)
    node1 <- nodes[[1]]
    node2 <- nodes[[2]]

    spin1 <- graph$spins[node1]
    spin2 <- graph$spins[node2]
    # cat("Edge:", edge, "Nodes:", node1, " spin1:", spin1, "and", node2, "spin2:", spin2, "\n")

    # cat("A allsum (", allsum, ") aggiungo ", spin1 * spin2 * J)
    allsum <- allsum + spin1 * spin2 * J
    # cat(" per un totale di ", allsum, "\n")
  }
  return(M - allsum)
  # return(- allsum) # per semplicità per ora è questa
}

get_avg_degree <- function(graph) {
  # v2:
  # dd <- degree_distribution(graph)
  # return(sum(dd*seq(0, length(dd) - 1)))
  # v1:
  return(mean(degree(graph)))
  # degrees <- degree(graph)
  # degree_probs <- table(degrees) / length(degrees)
  # average_degree <- sum(unname(degree_probs) * as.integer(names(degree_probs)))
  # average_degree
}

get_moment_2 <- function(graph) {
  node_degrees <- degree(graph)
  return(sum(node_degrees^2) / length(node_degrees))
}

get_moment_4 <- function(graph) {
  node_degrees <- degree(graph)
  return(sum(node_degrees^4) / length(node_degrees))
}

extract_index <- function(graph) {
  index <- sample.int(vcount(graph), 1) # we select which node will have the spin changed
  return(index)
}

# calculate the difference between the configuration with the node index flipped (the old configuration)
# and the actual configuration (so we are assuming graph$spins[index] is already flipped)
calculate_delta <- function(graph, spins, index) {
  nb <- neighbors(graph, index)
  J <- 1
  # the two it's because we have to subtract the old energy and then add the new one that is its opposite (so - (A) + (-A))
  delta <- 2 * sum(spins[index] * spins[nb] * J * -1)
  return(delta)
}

calculate_delta_m <- function(graph, spins, index) {
  delta <- 2 * spins[index]
  return(delta)
}

plot_graph <- function(graph, set_size = T, name = "graph") {
  if (set_size) {
    options(repr.plot.width = 12, repr.plot.height = 12)
  }
  if (is.null(graph$layout_design)) {
    colors <- "cornflowerblue"
  } else {
    colors <- generate_colors(graph$spins)
  }


  if (is.null(graph$layout_design)) {
    layout_ <- layout_nicely(graph)
  } else {
    layout_ <- graph$layout_design
  }

  plot(graph, main = name, edge.width = 10, vertex.color = colors, vertex.size = 10, layout = layout_)
}


simulate_ising <- function(graph, kbt, nstep, plot_every = 0) {
  energies <- rep(-Inf, nstep)
  energies[1] <- get_energy(graph)

  magnets <- rep(-Inf, nstep)
  magnets[1] <- get_magnetization(graph)

  suscep <- rep(-Inf, nstep)
  heat <- rep(-Inf, nstep)
  suscep[1] <- var(magnets) / kbt
  heat[1] <- var(heat) / (kbt^2)

  last_m_avg <- mean(graph$spins)

  spins <- graph$spins

  randoms <- runif(nstep)

  for (i in 2:nstep) {
    index <- extract_index(graph)
    # cat("Ho estratto l'indice", index, "che ha spin:", graph$spins[index], "\n")

    # let's change the index
    spins[index] <- spins[index] * -1

    # let's calculate the delta
    delta <- calculate_delta(graph, spins, index)

    # the prob of happening
    prob <- exp(-delta / kbt)
    random_number <- randoms[i]

    change <- F
    if (prob > 1 || prob > random_number) {
      change <- T
    }
    # cat("Index:", index, " Delta:", delta, "prob:", prob, "random:", random_number, "change:", change, "\n")

    graph$spins <- spins


    if (change) {
      energies[i] <- energies[i - 1] + delta
      magnets[i] <- magnets[i - 1] + calculate_delta_m(graph, spins, index)

      # it should be (<m^2> - <m>^2)/kbT but <m^2> is 1 always (because m is only -1 or 1)
      foo <- (last_m_avg + (spins[index] - spins[index] * (-1)) / N)
      suscep[i] <- (1 - foo^2) / kbt
      last_m_avg <- foo

      bar <- energies[i] / N
      heat[i] <- (1 - bar^2) / (kbt^2)
    } else {
      energies[i] <- energies[i - 1]
      magnets[i] <- magnets[i - 1]
      spins[index] <- spins[index] * -1

      suscep[i] <- suscep[i - 1]
      heat[i] <- heat[i - 1]
    }

    # if ((plot_every > 0) & ((i %% plot_every) == 0)) {
    #   plot_graph(graph)
    # }

    # cat(spins[index], "\n")


    # heat[i] <- var(heat) / (kbt^2)
  }

  graph$spins <- spins
  graph$energies <- energies
  graph$magnetizations <- magnets
  graph$susceptibility <- suscep
  graph$heat <- heat

  return(graph)
}


multi_T_ising <- function(nsteps, kbts, graph, perc = 0.02) {
  tsteps <- length(kbts)
  avgs_en <- rep(-Inf, tsteps)
  avgs_m <- rep(-Inf, tsteps)
  avg_suscep <- rep(-Inf, tsteps)
  avg_heat <- rep(-Inf, tsteps)
  # kbts <- rep(0, tsteps)
  # all_energies <- matrix(0, tsteps, nsteps)
  # all_magnets <- matrix(0, tsteps, nsteps)

  LAST_N <- as.integer(nsteps * perc)
  # PARAM1 <- 1.8
  # graph <- get_graph(model = model, only_bcc = T, add_spins = T, param1 = PARAM1)

  # graph <- delete_vertices(graph, which(degrees(graph) == 1))

  TC <- get_Tc(graph)
  cat("Avg degree:", get_avg_degree(graph), "\n")
  cat("Second moment:", get_moment_2(graph), "\n")
  cat("TC:", TC, "\n")

  # KBT <- TC * 2
  # kbt0 <- 3

  for (i in 1:tsteps) {
    kbt <- kbts[i]
    # print(kbt)
    # cat("Inizio Ising\n")
    graph <- simulate_ising(graph, kbt = kbt, nstep = nsteps)
    # cat("Fine Ising\n")
    energies <- graph$energies
    magnetizations <- graph$magnetizations
    suscep <- graph$susceptibility
    heat <- graph$heat

    avgs_en[i] <- mean(energies[(nsteps - LAST_N):nsteps])
    avgs_m[i] <- mean(magnetizations[(nsteps - LAST_N):nsteps])
    avg_suscep[i] <- mean(suscep[(nsteps - LAST_N):nsteps])
    avg_heat[i] <- mean(heat[(nsteps - LAST_N):nsteps])

    message(i, " - Ho calcolato la simulazione con kbt: ", kbt, appendLF = F)

    graph$energies <- NULL
    graph$magnetizations <- NULL
    graph$avg_suscep <- NULL
    graph$avg_heat <- NULL

    gc()
  }

  graph$avgs_en <- avgs_en
  graph$avgs_m <- avgs_m
  graph$avg_sus <- avg_suscep
  graph$avg_heat <- avg_heat
  # graph$all_energies <- all_energies
  # graph$all_magnets <- all_magnets
  graph$kbts <- kbts

  return(graph)
}


get_scalefree_with_constraint <- function(n, power, min_degree) {
  i <- 0
  while (TRUE) {
    tryCatch(
      {
        sequence <- sample(min_degree:n, n, replace = TRUE, prob = (min_degree:n)^-power)
        graph <- degree.sequence.game(sequence, method = "vl")
        graph$name <- "Scale Free"
        graph$nickname <- "sf"
        graph$m0 <- min_degree
        graph$param <- power
        graph$N <- N
        message("Dopo ", i, " tentativi sono riuscito a generare il grafo", appendLF = F)
        return(graph)
      },
      warning = function(w) {
        print("Warning! Something went wrong but I can continue")
      },
      error = function(e) {
        # cat("Errore nella generazione, riprovo\n")

      },
      finally = {
      }
    )
    i <- i + 1
  }
  # real_n <- n
  # while (TRUE) {
  #   n <- n + 1
  #   tryCatch(
  #     {
  #       x <- min_degree:n
  #       y <- (x^-power) * n / sum((x^-power))
  #       ry <- round(y)
  #       xi_cut <- which(ry != 0)
  #       ry <- ry[xi_cut]
  #       x <- x[xi_cut]

  #       delta <- real_n - sum(ry)
  #       cat("Devo togliere", delta, "nodi\n")
  #       ry[1] <- ry[1] + delta
  #       # cat("sumt", sum(ry), "\n")

  #       degree_sequence <- rep(x, ry)
  #       message(degree_sequence)
  #       graph <- degree.sequence.game(degree_sequence, method = "vl")

  #       graph$name <- "Barabasi-Albert"
  #       graph$nickname <- "ba"
  #       graph$m0 <- min_degree
  #       graph$param <- power
  #       graph$N <- N

  #       return(graph)
  #     },
  #     warning = function(w) {
  #       print("Warning! Something went wrong but I can continue")
  #     },
  #     error = function(e) {
  #       cat("Errore nella generazione, riprovo\n")
  #     },
  #     finally = {
  #     }
  #   )
  # }
}
