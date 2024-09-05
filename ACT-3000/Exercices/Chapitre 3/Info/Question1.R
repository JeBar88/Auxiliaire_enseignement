### ACT-3000
## Chapitre 3 info : Question 1
## Jérémie Barde

### Variables
mu <- log(100) - 0.5
s <- 1
a <- 1.5
lam <- 50
b <- 0.01
u <- 0.95

### (b)
VaRS <- function(u) qexp(u, b) + qpareto(u, a, lam) + qlnorm(u, mu, s)
VaRS(0.99)

### (c)
Fs <- function(x) optimize(function(u) abs(VaRS(u) - x), c(0, 1))$min
Fs(2000)
