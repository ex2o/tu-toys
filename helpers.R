cartesian_product <- function(M){
  purrr::invoke(tidyr::expand_grid, M, .name_repair="minimal")}

powerset_coalitions <- function(N) {
  list_01 <- rep(list(c(0,1)),N)
  coalitions <- cartesian_product(list_01)
  colnames(coalitions) <- paste0('p',1:N)
  coalitions
}

powerset_order <- function(coalitions) {
  C <- nrow(coalitions)
  result <- list()
  for (k in 1:C) {
    r <- coalitions[k,] # current row
    R <- coalitions[1:(k-1),] # previous rows
    downset <- which(apply(R, MAR=1, FUN=function(rr){all(rr <= r)}))
    result[[k]] <- downset
  }
  return(result)
}

v_equal_value_game <- function(game) {
  coalitions <- dplyr::select(game, starts_with('p'))
  rowSums(coalitions)
}
