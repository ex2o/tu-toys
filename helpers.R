cartesian_product <- function(M){
  purrr::invoke(tidyr::expand_grid, M, .name_repair="minimal")}

powerset_coalitions <- function(N) {
  list_01 <- rep(list(c(F,T)),N)
  coalitions <- cartesian_product(list_01)
  colnames(coalitions) <- paste0('x',1:N)
  coalitions
}

powerset_children <- function(coalitions) {
  coalitions <- dplyr::select(coalitions, starts_with('x'))
  C <- nrow(coalitions)
  player_names <- colnames(coalitions)
  result <- list()
  for (k in 1:C) {
    r <- coalitions[k,] # current row
    R <- coalitions[1:(k-1),] # previous rows
    #downset <- which(apply(R, MAR=1, FUN=function(rr){all(rr <= r)}))
    children <- which(apply(R, MAR=1, FUN=function(rr){
      r <- unlist(r)
      sum(rr[r]) == sum(r)-1 & sum(!rr[!r]) == sum(!r)
    }))
    result[[k]] <- children
  }
  return(result)
}

v_equal_value <- function(game) {
  coalitions <- dplyr::select(game, starts_with('x'))
  rowSums(coalitions)
}

marginal_contributions <- function(game) {
  C <- nrow(game)
  result <- list()
  for (k in 1:C) {
    coalition <- game[k,]
    children <- coalition$children[[1]]
    n_children <- length(children) 
    if (n_children == 0) {
      m <- numeric(0)
    } else {
      m <- rep(NA,n_children)
      for (j in seq_along(children)) {
        child_v <- game[children[j],]$v
        m[j] <- coalition$v - child_v
      }
    }
    result[[k]] <- m
  }
  result
}

shapley_weights <- function(game) {
  coalitions <- dplyr::select(game, starts_with('x'))
  sizes <- rowSums(coalitions)
  N <- ncol(coalitions)
  1/N/choose(N-1,sizes)
}

# We can write any mu as mu = sum(b*v) for coefficients b
# This function computes the weighting distribution for marginal contributions,
# given a choice of b
b_to_weights <- function(b) {
  n <- length(b)
  w <- rep(NA, n)
  w[n] <- b[n]/n
  for (t in (n-1):1) {
    w[t] <- (b[t] + (n-t)*w[t+1])/t 
  }
  return(w)
}

# This function goes in the opposite direction from b_to_weights
weights_to_b <- function(w) {
  n <- length(w)
  b <- rep(NA, n)
  b[n] <- n*w[n]
  for (t in (n-1):1) {
    b[t] <- t*w[t] - (n-t)*w[t+1]
  }
  return(b)
}

check_symmetric_weights <- function(w,N) {
  sum(w*choose(N-1,0:(N-1))) == 1
}

b_to_mu <- function(game,b) {
  teams <- game %>% dplyr::select(starts_with("x"))
  v <- game$v[-1]
  size_s <- rowSums(teams)
  return(sum(v*b[size_s]))
}

# find subsets with a given player (or set of indices)
find_subsets <- function(player, game) {
  n <- length(player)
  game[rowSums(game[,player]) == n,]
}

# given a game with binary subset matrix S, sort it by subset size, and keep subset sizes
sort_by_size <- function(game) {
  game %>% rowwise() %>%  
    mutate(size = sum(c_across(starts_with("x")))) %>% 
    arrange(size)
}

is_child <- function(x,y,size_x,size_y) {
  size_x == size_y-1 && all(x<=y) 
}

powerset_lattice_weights <- function(game, w) {
  S <- game %>% select(starts_with("x"))
  size <- game$size
  n <- 2^N
  M <- matrix(0, nrow=n, ncol=n)
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      if (is_child(S[j,], S[i,], size[j], size[i])) {
        # assign weight of the subset,
        # since weight vector w has w[|S|+1] = weight(|S|) for |S| = 0,...,N-1
        M[i,j] <- w[size[i]]
      }
    }
  }
  return(M)
}

powerset_lattice_binary <- function(game) {
  S <- game %>% select(starts_with("x"))
  size <- game$size
  n <- 2^N
  M <- matrix(0, nrow=n, ncol=n)
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      if (is_child(S[j,], S[i,], size[j], size[i])) {
        M[i,j] <- 1
      }
    }
  }
  return(M)
}
