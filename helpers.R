cartesian_product <- function(M){
  purrr::invoke(tidyr::expand_grid, M, .name_repair="minimal")}

powerset_coalitions <- function(N) {
  list_01 <- rep(list(c(F,T)),N)
  coalitions <- cartesian_product(list_01)
  colnames(coalitions) <- paste0('x',1:N)
  coalitions
}

powerset_parents <- function(coalitions) {
  C <- nrow(coalitions)
  result <- list()
  for (k in 1:C) {
    r <- coalitions[k,] # current row
    R <- coalitions[1:(k-1),] # previous rows
    #downset <- which(apply(R, MAR=1, FUN=function(rr){all(rr <= r)}))
    parent <- which(apply(R, MAR=1, FUN=function(rr){
      r <- unlist(r)
      sum(rr[r]) == sum(r)-1 & sum(!rr[!r]) == sum(!r)
    }))
    result[[k]] <- parent
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
    parents <- coalition$parents[[1]]
    n_parents <- length(parents) 
    if (n_parents == 0) {
      m <- numeric(0)
    } else {
      m <- rep(NA,n_parents)
      for (j in seq_along(parents)) {
        parent_v <- game[parents[j],]$v
        m[j] <- coalition$v - parent_v
      }
    }
    result[[k]] <- m
  }
  result
}
