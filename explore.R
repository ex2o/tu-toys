# 1. provide a set of sub-coalitions (list of binary vectors) 
# 2. provide parents in the partial order (list of parent indices)
# 3. provide a characteristic function on 1
# 4. compute the marginal contributions using 1 and 2
# 5. provide a weighting distribution or vector of weights on 3
# 6. compute the weighted average using 3 and 4
# 7. compute a transformation (e.g., identity or share function form) of 5

library(purrr)
library(tidyr)
source("helpers.R")

# 1. sub-coalitions -------------------------------------------------------

N <- 3
game <- powerset_coalitions(N)

# 2. partial order parents ------------------------------------------------

game$parents <- powerset_parents(game)

# 3. characteristic function ----------------------------------------------

game$v <- v_equal_value(game)

# 4. compute marginal contributions ---------------------------------------

game$marginals <- marginal_contributions(game)
#View(game)

# 5. Provide weighting distribution ---------------------------------------



