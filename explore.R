library(purrr)
library(tidyr)
library(dplyr)
source("helpers.R")

# Simpler approach ---------------------------------------------------------

N <- 3
game <- powerset_coalitions(N)
game$v <- v_equal_value(game)

# weight vector w has w[|S|+1] = weight(|S|) for |S| = 0,...,N-1
w <- b_to_weights(c(rep(0,N-1),1)) # Banzhaf = rep(1,N)/2^(N-1) # or Shapley = b_to_weights(c(rep(0,N-1),1))
b <- weights_to_b(w)
check_symmetric_weights(w,N)
mu <- b_to_mu(game, b)
mu > 0 # check mu is positive

# we represent the lattice using a lower-triangular weight matrix (size 2^(2N))
# the matrix generating function assumes coalitions are sorted by size
game <- sort_by_size(game)
M <- powerset_lattice_weights(game,w)
V <- matrix(game$v)
sol <- rep(NA, N)
for (i in 1:N) {
  S_with_i <- game[[paste0("x",i)]] == 1
  M_i <- M[S_with_i, !S_with_i]
  sol[i] <- sum(M_i%*%V[S_with_i] - M_i%*%V[!S_with_i])
}
sol; sum(sol); sum(sol) == mu

# an alternative that avoids non-binary lattice matrix (for semivalues)
game <- sort_by_size(game)
M <- powerset_lattice_binary(game)
V <- matrix(game$v)
sol <- rep(NA, N)
for (i in 1:N) {
  S_with_i <- game[[paste0("x",i)]] == 1
  w_i <- w[game$size[S_with_i]]
  M_i <- M[S_with_i, !S_with_i]
  delta <- M_i%*%V[S_with_i] - M_i%*%V[!S_with_i]
  sol[i] <- sum(w_i*delta)
}
sol; sum(sol); sum(sol) == mu

# Experimenting with mu functions -----------------------------------------

N <- 3
game <- powerset_coalitions(N)
game$parents <- powerset_parents(game)
game$v <- v_equal_value(game)
coalitions <- dplyr::select(game, starts_with('x'))
banz <- (2*rowSums(coalitions) - N)/(2^(N-1))
# Due to symmetry, only team size matters
sort(banz);sum(banz)

banz*game$v;sum(banz*game$v)
v <- c(0,1,1,1,2,2,1,2)
v[1] <- 0
sort(banz)*v;sum(sort(banz)*v)
plot(sort(banz), type='l')

# Sum to sum of Banzhaf values (creates Banzhaf weights)
b_to_weights(c(-1/4,1/4,3/4))
weights_to_b(c(1/4,1/4,1/4))

# TODO: AWESOME: insight into Banzhaf value - sample from binomial distribution!
plot(weights_to_b(rep(1/2^99,100))*choose(99,0:99)) # TODO: SO YEAH! can tend to sum to 0
plot(choose(99,0:99)/(2^99)) # weights are constant for Banzhaf so middle coalitions have more influence
# so yeah, just get the binomial parameters right...
hist(rbinom(100000,102,1/2), freq=F, breaks=50) # try other p
lines(choose(99,0:99)/(2^99), col="red") # Banzhaf
# 1. Method 1: binomial sample coalition sizes S, then compute all for size S (for all features for fairness)
# 2. Method 2: uniform sample coalitions, then binomial emerges
lines(40:60,choose(99,39:59)/(2^99), col="green", lwd=3)

# Avg weight of unit coalitions
b_to_weights(c(1/3,0,0))

# Avg final contribution weights
b_to_weights(c(0,-1/3, 1))

# Sum to average of first few coalitions only
b_to_weights(c(rep(1/10,10),rep(0,90)))

# Sum to average of first and last coalitions
b_to_weights(c(1/5,0,0,0,1))

# Shapley style weights come from b[n] = 1 and rest 0
b_to_weights(c(rep(0,2),1))*choose(2,0:2) # multiple by num coalitions
b_to_weights(c(rep(0,99),1))*choose(99,0:99)

# TODO: AWESOME. if b gives equal weight to all coalitions, then the weights drop dramatically mid way
b_to_weights(rep(1,5))
w <- b_to_weights(rep(1/100,100))
plot(w/w[1]*choose(99,0:99))
abline(v=50)
plot(w/w[1])
sum(w/w[1]>0.00000001)

# If b gives binomial weights to coalitions
plot(b_to_weights(choose(99,0:99))*choose(99,0:99))

# TODO: OOOH. If b gives 1/binomial weights to coalitions
# In other words, b gives Shapley weights to coalitions
plot(b_to_weights(1/choose(99,0:99))*choose(99,0:99))
plot(b_to_weights(b_to_weights(c(rep(0,99),1)))*choose(99,0:99))
# go one deeper
b <- b_to_weights(b_to_weights(b_to_weights(c(rep(0,99),1))))
plot(b*choose(99,0:99))
plot(b_to_weights(b)*choose(99,0:99)) # go one deeper 
# but if we apply the rescaling before going one deeper... WTH
w <- b_to_weights(b_to_weights(b_to_weights(c(rep(0,99),1)))*choose(99,0:99))*choose(99,0:99)
plot(w/w[1], type = 'l', col='black')
w2 <- b_to_weights(rep(1/100,100))*choose(99,0:99)
w3 <- b_to_weights(choose(99,0:99))*choose(99,0:99)
lines(w2/w2[1], col = "blue")
lines(w3/w3[1], col = "red")

# if mu gives weight equal to coalition size, then decrease is slower
b_to_weights((1:5))

# TODO: AWESOME!! effect of decreasing b by coalition size
b_to_weights(1/(1:10)^(1/4))
b_to_weights(1/(1:10)^(1/2))
b_to_weights(1/(1:10))
b_to_weights(1/((1:10)^2))
b_to_weights(1/((1:10)^3))
sum(b_to_weights(1/(1:100)) > 0.01)
b_to_weights(10:1)
plot(b_to_weights(1/(1:10)), type='l')
lines(b_to_weights(1/((1:10)^2)))
lines(b_to_weights(1/((1:10)^3)))

# effect of high b for last t coalitions only
b_to_weights(c(rep(0,7),rep(1,3)))
w <- b_to_weights(c(rep(0,70),rep(1,30)))
round(w/w[1],3)

# test that it works
weights_to_b(b_to_weights(c(rep(0,3),1)))
weights_to_b(b_to_weights(rep(1,5)))

# effect of the banzhaf weights on mu 
# (3 player game -> 4 subsets containing i, 3 subset sizes, all equal weight)
weights_to_b(rep(1/4,3))

# effect of shapley weights on mu
weights_to_b(c(1/3,1/6,1/3))
weights_to_b(c(0.333, 0.167, 0.333))

# effect of increasing weights by coalition size
weights_to_b(1:5)

# effect of decreasing weights by coalition size
weights_to_b(1/(1:10)^(1/4))
weights_to_b(1/(1:10)^(1/2))
weights_to_b(1/(1:10))
weights_to_b(1/((1:10)^2))
weights_to_b(1/((1:10)^3))
weights_to_b(1/((1:10)^4))
w <- weights_to_b(1/((1:10)^5))
round(w/w[1],3)

# Note the weights technically need to be positive
# effect of weighting first t coalitions only, few features
weights_to_b(c(rep(1,3),rep(0,7)))

# effect of weighting first t coalitions only, many features
weights_to_b(c(rep(1/100,5),rep(0,95)))

# TODO: COOL effect of weighting last t coalitions only, many features
weights_to_b(c(rep(0,95),rep(1/5,5)))*choose(99,0:99)

# TODO: GOOD effect of weighting last t coalitions only, few features
weights_to_b(c(rep(0,7),rep(1/3,3)))

# effect of weighting last t coalitions only, increasing weight
weights_to_b(c(rep(0,100),1:3))


# Old approaches ----------------------------------------------------------

# 1. provide a set of sub-coalitions (list of binary vectors) 
# 2. provide parents in the partial order (list of parent indices)
# 3. provide a characteristic function on 1
# 4. compute the marginal contributions using 1 and 2
# 5. provide a weighting distribution or vector of weights on 3
# 6. compute the weighted average using 3 and 4
# 7. compute a transformation (e.g., identity or share function form) of 5

# 1. sub-coalitions -------------------------------------------------------

N <- 3
game <- powerset_coalitions(N)

# 2. partial order children ------------------------------------------------

game$children <- powerset_children(game)

for (i in 1:nrow(game)) {
  
  current <- game[i,]
  print(current)
  print(game[current$children[[1]],])
  cat("\n\n")
  
}

# 3. characteristic function ----------------------------------------------

game$v <- v_equal_value(game)

# 4. compute marginal contributions ---------------------------------------

game$marginals <- marginal_contributions(game)
#View(game)

# 5. Provide weighting distribution ---------------------------------------

game$shap_w <- shapley_weights(game)

allocation <- function(game, player="x1", weights="shap_w") {
  S <- game %>% dplyr::filter(.data[[player]] == T)
  
}

allocation(game, "x1")$parents


# 0. Specify sets as binary strings with value column (B,v) 
#    -- matrix
# 1. Specify function for the weight w (for contributions) for each subset size 
#    OR specify the b coefficient (for mu) for each subset size (one gives the other).
#    (more generally, w can be a function of the set) 
#    -- weight(B)
# 2. For each index (player) i 
#    -- that is, column of B (or, more generally, i can be a set of columns):
#      1) Find all subsets S containing i
#         -- subsets(B,player)
#      2) For each S:
#         i) Find all "i-children" C
#            -- children(S,i)
#         ii) For each (i,C):
#             (1) compute marginal contribution
#             (2) weight?

# # 2.??
# # we need a map from player index to columns?? (usually just identity)
# player_map <- list("x1"=1,"x2"=2,"x3"=3) # e.g., list("x1"=c(1,2),"x2"=3)
# S <- lapply(player_map, find_subsets, game=game)

# 100 1 
# 010 1
# 001 1
# 110 2
# 101 2
# 011 2
# 111 3

# 1000 1
# 0100 2
# 0010 3
# 0001 4
# 1100 12
# 1010 13
# 1001 14
# 0110 23
# 0101 24
# 0011 34
# 1110 123
# 1101 124
# 1011 134
# 0111 234
# 1111 1234


#{}
#{3}
#{2}
#{2,3}
#{1}
#{1,3}
#{1,2}
#{1,2,3}


