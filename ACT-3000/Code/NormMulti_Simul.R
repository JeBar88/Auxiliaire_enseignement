### ACT-3000
## Chapitre 3 partie 3 : Loi normal multivariée
## Jérémie Barde

##### Simulation : Exemple 10 #####
m <- 1e6
d <- 3
mu <- c(15, 10, 16)
s <- matrix(c(4, 2.4, -2, 
              2.4, 9, -4.5,
              -2, -4.5, 25), ncol = 3, byrow = TRUE) 

mB <- t(chol(s))

Z <- replicate(d, rnorm(m))
X <- t(mu + (mB %*% t(Z)))

apply(X, 2, mean)
cov(X)

S <- rowSums(X)
var(S)
