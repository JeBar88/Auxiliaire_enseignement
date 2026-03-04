### ACT-3000
## Chapitre 3 info : Question 7
## Jérémie Barde

### Foncition utiles
dPoTei <- function(m1, m2, lam1, lam2, a){
  f <- matrix(numeric(length(m1) * length(m2)), ncol = length(m2))
  
  for (i in m1){
    for (j in m2){
      
      k <- 0:min(i, j)
      
      f[i + 1, j + 1] <-  sum(dpois(k, a) * dpois(i - k, lam1 - a) * dpois(j - k, lam2 - a))
    }
  }
  f
}
repart <- function(x, y, fxy){
  f <- matrix(numeric(length(x)*length(y)), ncol = length(y))
  for (i in x) {
    for (j in y) {
      f[i + 1, j + 1] <- sum(fxy[0:i + 1, 0:j + 1])
    }
  }
  f
}

### Variables
a <- 2:1
b <- 1/100
lam <- c(10, 7)
a0 <- 2
k <- 0:50

### a)
fm12 <- dPoTei(k, k, lam[1], lam[2], a0)
fm12[8, 6]

### b)
Fm12 <- repart(k, k, fm12)
Fm12[6, 11]
Fm12[6, 11] - Fm12[4, 11] - Fm12[6, 5] + Fm12[4, 5]

### c)
EM12 <- sum(outer(k, k)*fm12)
sum(outer(pmin(k, 5), k*I(k > 4))*fm12)

### d)
CovM12 <- a0
CovX12 <- prod(a/b)*CovM12

### e)
perm <- expand.grid(k, k)
am <- a[1]*perm[, 1] + a[2]*perm[, 2]
fm <- as.vector(fm12)

Fs <- function(x) fm[1] + sum(fm[-1]*pgamma(x, am[-1], b))
1 - Fs(2e3)

# ou trouer la loi de S

### f)
VaRS <- function(u) optimize(function(x) abs(Fs(x) - u), c(0, 1e4))$min
VaRS(0.995)

TVaRS <- function(u){
  VaR <- VaRS(u)
  sum(fm*am/b*pgamma(VaR, am + 1, b, low=FALSE))/(1 - u)
}
TVaRS(0.995)
  


















