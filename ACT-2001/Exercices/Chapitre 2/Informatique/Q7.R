### ACT-2001
## Chapitre 2 informatique : Question 7
## Jérémie Barde

a <- 0.8
mu <- c(0.1, -0.3)
s <- c(0.2, 0.1)
k <- c(0.0001, 0.01, 0.5, 0.99, 0.9999)

### C)
Fx <- function(x) a*pnorm(x, mu[1], s[1]) + (1-a)*pnorm(x, mu[2], s[2])
Fx(-2)
Fx(2)

VaR <- function(k) optimize(function(x) abs(Fx(x) - k), c(-2, 2))$min
sapply(k, VaR)
