###############################
### Chapitre 3 exercice trad 24
###############################
library(Runuran)
library(actuar)
nfft <- 2^20
b <- c(0.1, 0.2, 1/(1/0.1 + 1/0.2))
a <- c(2.5, 5)
g <- 1
aS <- sum(a) - g
q <- (b/max(b))[-which.max((b/max(b)))]

j1 <- dnbinom(0:(nfft - 1), a[1] - g, q[1])
j2 <- dnbinom(0:(nfft - 1), g, q[2])

kt <- fft(j1) * fft(j2) 
k <- Re(fft(kt, T))/nfft

Fs <- function(x) sum(k[1:1001] * pgamma(x, aS + 0:1000, max(b)))

### b)
E <- sum(k[1:1001] * (aS + 0:1000)/max(b))
Var <- sum(k[1:1001] * ((aS + 0:1000) + (aS + 0:1000)^2)/max(b)^2) - E^2
cbind(E, Var)

### c)
# Avec optimize
VaR <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 1000))$min
VaR(0.95)

# Par interpolation
gen <- pinv.new(cdf=Fs, lb=0, ub=Inf)
VaR2 <- function(u) uq(gen, u)

TVaR <- function(u) 1/(1 - u) * sum(k[1:1001] * (aS + 0:1000)/max(b) * pgamma(VaR2(u), aS + 0:1000 + 1, max(b), low = F))
TVaR(0.95)








