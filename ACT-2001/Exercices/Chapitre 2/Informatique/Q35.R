### ACT-2001
## Chapitre 2 informatique : Question 35
## Jérémie Barde

x <- c(0, 10)
fX <- c(0.9, 0.1)
FX <- cumsum(fX)

s <- c(0, 10, 20)
fS <- c(0.81, 0.18, 0.01)
FS <- cumsum(fS)

# (a)
VaRX <- function(k) x[min(which(FX >= k))]
sapply(c(0.5, 0.85, 0.99), VaRX)

VaRS <- function(k) s[min(which(FS >= k))]
sapply(c(0.5, 0.85, 0.99), VaRS)

# (b)
TVaRX <- function(k) sum(pmax(x - VaRX(k), 0) * fX)/(1-k) + VaRX(k)
sapply(c(0.5, 0.85, 0.99), TVaRX)

TVaRS <- function(k) sum(pmax(s - VaRS(k), 0) * fS)/(1-k) + VaRS(k)
sapply(c(0.5, 0.85, 0.99), TVaRS)

# Utiliser la formule pour passer à la LTVaR