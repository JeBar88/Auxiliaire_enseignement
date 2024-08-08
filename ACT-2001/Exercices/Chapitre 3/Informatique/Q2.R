### ACT-2001
## Chapitre 3 informatique : Question 2
## Jérémie Barde

m <- 1E5
a <- 2
lam <- 0.2
mu <- log(10) - 0.32
s <- 0.8

set.seed(2019)
U <- matrix(runif(2 *m), ncol = 2, byrow = T)
head(U)
tail(U)

x1 <- qgamma(U[,1], a, lam)
head(x1)
tail(x1)

x2 <- qlnorm(U[,2], mu, s)
head(x2, 4)
tail(x2,1)
