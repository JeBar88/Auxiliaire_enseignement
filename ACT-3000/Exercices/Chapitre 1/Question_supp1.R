### ACT-3000
## Chapitre 3 : Q25 À changer
## Jérémie Barde

### Fonction utils
dPoTei <- function(m1, m2, lam1, lam2, a){
  f <- matrix(0, nrow=length(m1), ncol = length(m2))
  
  for (i in m1){
    for (j in m2){
      
      k <- 0:min(i, j)
      
      f[i + 1, j + 1] <-  sum(dpois(k, a) * dpois(i - k, lam1 - a) * dpois(j - k, lam2 - a))
    }
  }
  f
}
repart <- function(x, y, fxy){
  f <- matrix(0, nrow=length(x), ncol = length(y))
  for (i in x) {
    for (j in y) {
      I <- matrix(numeric(length(x)*length(y)), ncol = length(y))
      I[0:i + 1, 0:j + 1] <- 1
      f[i + 1, j + 1] <- sum(I * fxy)
    }
  }
  f
}

### Variable
lam <- c(10, 7)
a <- c(2, 1)
b <- 1/100
g <- c(0, 0.5, 2)
m1 <- 0:100
m2 <- 0:100

### a)
fm12 <- dPoTei(m1, m2, lam[1], lam[2], g[3])
sum(fm12)
fm12[8, 6]

### b)
sum(fm12[0:5 + 1, 0:10 + 1])
## ou
Fm12 <- repart(m1, m2, fm12)
Fm12[6, 11]

a1 <- 3 + 1
b1 <- 5 + 1
a2 <- 4 + 1
b2 <- 10 + 1

Fm12[6, 11] - Fm12[4, 11] - Fm12[6, 5] + Fm12[4, 5]
Fm12[b1, b2] - Fm12[a1, b2] - Fm12[b1, a2] + Fm12[a1, a2]

### c)
Em12 <- sum(outer(m1, m2, '*')*fm12)
Em12S <- sum(outer(pmin(m1, 5), m2*I(m2 > 4), '*')*fm12)

### d)
covM12 <- Em12 - prod(lam)
EB <- mgamma(1, a, b)
VarB <- mgamma(2, a, b) - EX^2
VarX <- lam*VarB + lam*EB^2

covX12 <- prod(EB)*covM12
rhoX12 <- covX12/sqrt(prod(VarX))

### e)
## Fonction de masse sous forme de vecteur
fm <- as.vector(fm12)

## Ajustement des paramètre
comb <- expand.grid(m1, m2)
am <- comb[, 1]*a[1] + comb[, 2]*a[2]

## Calculer Fs
Fs <- function(x) fm[1] + sum(fm[-1]*pgamma(x, am[-1], b))
1 - Fs(2000)

### e)
u <- 0.995
VaR <- optimize(function(x) abs(Fs(x) - u), c(0, 1e4))$min
TVaR <- sum(fm[-1] * am[-1]/b * pgamma(VaR, am[-1] + 1, b, low=FALSE))/(1 - u)
cbind(VaR, TVaR)







