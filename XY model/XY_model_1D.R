library("igraph") |> suppressPackageStartupMessages()
source("const.R")

# Set up your bot
bot_token <- API
bot_chat_id <- Ric_id # Riccardo
bot_chat_id_2 <- Anto_id # Antonio

dim <- 1 # 1D
nei <- 3
p <- 0.2

H_matrix_calc <- function(N, mat, angles.samples){
    H <- rep(0, N)
    for (i in 1:N){
        for (j in (mat[[i]])){
            # These are all the n.n.
            H[i] <- H[i] - cos(angles.samples[i] - angles.samples[j]) 
        }
    }
    co <- mean(cos(angles.samples))
    si <- mean(sin(angles.samples))
    return(list(H = H, co = co, si = si))
}

compute_delta_H <- function(N, mat, angles.samples, T, H){
    change_pos <- sample(N, 1)
    new_val <- runif(1, -pi/6, pi/6) + angles.samples[change_pos]
    old_val <- angles.samples[change_pos]

    neighbors <- mat[change_pos][[1]]

    old <- - sum(cos(angles.samples[change_pos] - angles.samples[neighbors]))
    new <- - sum(cos(new_val - angles.samples[neighbors]))
    
    delta_H <- new - old

    # Metropolis passage
    if (delta_H < 0 || runif(n = 1, min = 0, max = 1) < exp(-delta_H/T)){
        angles.samples[change_pos] <- new_val
        H <- H + delta_H
    }
    return(list(angles.samples = angles.samples, H = H))
}


get_specific_heat <- function(all_H, T, nei, N){
    K <- 1/T
    mu <- mean(all_H*2/(N*nei))
    return(K^2 * (1-mu/K - mu^2))
}

# Step 1: Generate a Regular Lattice
generate_regular_lattice <- function(dim, size, nei) {
    g <- watts.strogatz.game(dim = dim, size = size, nei = nei, p = 0, loops = FALSE, multiple = FALSE)
    return(g)
}

# Step 2: Rewire Edges
rewire_edges <- function(mat, p, size, nei) {
    l <- (1:size)
    for (k in 1:nei){
        for (i in 1:size){
            # neig <- which(mat[i, l] == 1, arr.ind = TRUE)
            j <- (i + k) %% size
            if (j == 0) {
                j <- 1
            }
            new_pos <- sample(l, 1)
            if (mat[i, new_pos] == 1 || new_pos == i){
                next
            }
            else{
                if (runif(1) < p){
                    mat[i, j] <- 0; mat[j, i] <- 0
                    mat[i, new_pos] <- 1; mat[new_pos, i] <- 1
                }
            }
        } 
    }
    return(mat)
}

# Step 3: Create a Watts-Strogatz Network
create_WS_network <- function(dim, size, nei, p) {
    # Generate a regular lattice
    g <- generate_regular_lattice(dim, size, nei)
    mat <- rewire_edges(get.adjacency(g), p, size, nei)
    # Rewire edges
    adjacency_list <- lapply(1:nrow(mat), function(i) which(mat[i,] == 1))    
    return(adjacency_list)
}

NN <- c(1600)#, 200, 400)
Temperature <- seq(2, 2.4, 0.05)
# probs <- seq(0, 1, 0.1)
p <- 0.2
runs <- 1:1000
iter_per_step <- 1:5000
network_realization <- 100

for (N in NN){
    magn_N_2 <- c(); magn_N_4 <- c(); magn_N <- c(); heats <- c()
    for (n in 1:network_realization){
        m <- c(); m.2 <- c(); m.4 <- c()
    specific_heat <- c()

    mat <- create_WS_network(dim = dim, size = N, nei = nei, p = p)
    angles.samples <- runif(n = N, min = -pi, max = pi) 

    result <- H_matrix_calc(N = N, mat = mat, angles.samples = angles.samples)
    H_matrix <- result$H; co <- result$co; si <- result$si
    H <- sum(H_matrix)/2     
    
    for (T in Temperature){
        magn <- c(); all_H <- c(); m_squared <- c(); m_fourth <- c()
        for (run in runs){
            for (iter in iter_per_step){
                res <- compute_delta_H(N, mat, angles.samples, T, H)
                angles.samples <- res$angles.samples; H <- res$H
            }
            co <- mean(cos(angles.samples))
            si <- mean(sin(angles.samples))
            all_H <- append(all_H, H)
            
            # Compute magnetization and higer orders of m
            m_temp <- c(co, si)
            m_2 <- sum(m_temp^2)
            m_4 <- m_2^2

            # Save results for each run
            magn <- append(magn, norm(m_temp, type="2"))
            m_squared <- append(m_squared, m_2)
            m_fourth <- append(m_fourth, m_4) 
        }

        specific_heat <- append(specific_heat, get_specific_heat(all_H, T, nei, N))
        # Take the average for each run, so to have only one value per temperature
        m <- append(m, mean(magn))
        m.2 <- append(m.2, mean(m_squared))
        m.4 <- append(m.4, mean(m_fourth))
    }
    magn_N <- append(magn_N, m)
    magn_N_2 <- append(magn_N_2, m.2)
    magn_N_4 <- append(magn_N_4, m.4)
    heats <- append(heats, specific_heat)
}
    cat("Saving results for N =", N, "\n")
    res_m_2 <- rep(0, length(Temperature))
    res_m_4 <- rep(0, length(Temperature))
    res_m <- rep(0, length(Temperature))
    res_c <- rep(0, length(Temperature))
    for (i in 1:length(Temperature)){
        res_m_2[i] <- mean(magn_N_2[seq(i, length(magn_N_2), length(Temperature))])
        res_m_4[i] <- mean(magn_N_4[seq(i, length(magn_N_4), length(Temperature))])
        res_m[i] <- mean(magn_N[seq(i, length(magn_N), length(Temperature))])
        res_c[i] <- mean(heats[seq(i, length(heats), length(Temperature))])
    }
    U.N.T <- 1 - res_m_4/(3*res_m_2^2)
    # res_m <- res_m*(N^(1/4))
    folder <- 'res/'
    file_path_U <- paste(folder, 'U_', N, '.txt', sep='')
    file_path_c <- paste(folder, 'c_', N, '.txt', sep='')
    file_path_mN <- paste(folder, 'mN_', N, '.txt', sep='')
    write.table(U.N.T, file = file_path_U, col.names = FALSE, row.names = FALSE)
    write.table(res_c, file = file_path_c, col.names = FALSE, row.names = FALSE)
    write.table(res_m, file = file_path_mN, col.names = FALSE, row.names = FALSE)
}
header = paste('#  netw iters =', network_realization,' p =',p, 'newCv_newU', "Delta_T =", diff(Temperature)[1])
writeLines(header, 'res/metadata.txt')
cat("Done!\n")

# Send a message to the bot
message <- "Your R script has finished executing."
command <- sprintf('curl -X POST "https://api.telegram.org/bot%s/sendMessage" -d "chat_id=%s&text=%s"', bot_token, bot_chat_id, message)

# system(command)

folder_path <- "res/"  # Replace with the actual path to your folder
# Get a list of all files in the folder
files <- list.files(path = folder_path, full.names = TRUE)
# Iterate through the files and send them
for (file in files) {
    command <- sprintf('curl -F chat_id=%s -F document=@%s "https://api.telegram.org/bot%s/sendDocument"', bot_chat_id, file, bot_token)
    system(command)
    command <- sprintf('curl -F chat_id=%s -F document=@%s "https://api.telegram.org/bot%s/sendDocument"', bot_chat_id_2, file, bot_token)
    system(command)
    print("File sent.")
}
