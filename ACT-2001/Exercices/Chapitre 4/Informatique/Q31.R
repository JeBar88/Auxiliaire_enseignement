### ACT-2001
## Chapitre 4 informatique : Question 31
## Jérémie Barde

library(actuar)
### Chapitre 4 : Question 31
mu <- log(1000) - 0.405
s <- sqrt(0.405 * 2)
q <-  0.2

## a)
Esp <- q * mlnorm(1, mu, s)
Var <- q * (mlnorm(2, mu, s) - mlnorm(1, mu, s)^2) + q * (1-q) * mlnorm(1, mu, s)^2

qnorm(0.95)

Fx <- function(x) dbinom(0, 1, q) + dbinom(1, 1, q) * plnorm(x, mu, s)
Fx(7.9831)

VaR <- optimise(function(x) abs(Fx(x) - 0.99), c(0, 10000))$min

0.2 * 1000/(1-0.99) * exp(pnorm((log(VaR) - mu - s^2)/0.9))
