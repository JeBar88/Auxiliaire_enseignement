##########################
### Question 1 section 1.6
##########################
b <- c(0.1, 0.2)
q <- b[1]/b[2]
vk <- 0:10000
### a)
fj <- dgeom(vk, q)

Fs <- function(x) sum(fj * pgamma(x, 2 + vk, b[2]))

VaRS <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 2e2))$min
VaRS(0.99)

TVaRS <- function(k) sum(fj * (vk+2)/b[2] * pgamma(VaRS(k),vk + 3, b[2], low = F))/(1 - k)
TVaRS(0.99)

CTVaRX1 <- (1/b[1] - (1/b[1] * pexp(VaRS(0.99), b[1]) - VaRS(0.99) * exp(-b[1] * VaRS(0.99))))/(1 - 0.99)

