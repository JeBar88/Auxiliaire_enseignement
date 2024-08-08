### ACT-2001 ###
## Examen final informatique 2017 : Question 3
n <- 5
q <- 0.3

M <- function(x) sum(dbinom(0:5, n, q) * ppois(x, 0.1 * ((0:5) + 1)))
m <- function(x) sum(dbinom(0:5, n, q) * dpois(x, 0.1 * ((0:5) + 1)))

Fm <- sapply(0:100, M)
fm <- sapply(0:100, m)

VaRM <- min(which(Fm >= 0.9)) - 1
TVaRM <- VaRM + sum(pmax(0:100 - VaRM, 0) * fm)/(1-0.9) 
cbind(VaRM, TVaRM)

VaRX <- 1000 * VaRM
TVaRX <- 1000 * TVaRM
cbind(VaRX, TVaRX)
    
S <- function(x) sum(dbinom(0:5, n, q) * ppois(x, 10 * ((0:5) + 1)))
s <- function(x) sum(dbinom(0:5, n, q) * dpois(x, 10 * ((0:5) + 1)))

Fs <- sapply(0:100, S)
fs <- sapply(0:100, s)

VaRS <- min(which(Fs >= 0.9)) - 1

VaRWn <- 10 * (min(which(Fs >= 0.9)) - 1)
TVaRWn <- 10 * VaRS + sum(pmax(0:100 - VaRS, 0) * fs)/(1-0.9) 
cbind(VaRWn, TVaRWn)
