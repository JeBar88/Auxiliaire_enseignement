### ACT-2001
## Chapitre 2 informatique : Question 27
## Jérémie Barde

## a)
b <- 1e4
q <- 0.0012
k <- 0.995
x <- c(0, b)
fx <- c(1-q, q)
Fx <- cumsum(fx)

VaRX <- b*qbinom(k, 1, q)
TVaRX <- 1/(1-k) * sum(pmax(x - VaRX, 0) * fx)
cbind(VaRX, TVaRX)

## b)
VaRK <- qbinom(k, 200, q)
VaRS <- b*VaRK
TVaRS <- b * (1/(1-k) * sum(pmax(0:200 - VaRK, 0) * dbinom(0:200, 200, q)) + VaRK)
cbind(VaRS, TVaRS)