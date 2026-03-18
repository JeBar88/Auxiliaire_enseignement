### ACT-3000
## Loi poisson Teicher composée Expo/Gamma, bêta indentique
## Jérémie Barde

#### Loi discrète composée -- B~Bi~Ga(a, b) ####
k <- 0:3
a <- c(2.1, 3.4)
b <- 0.1

fm12 <- matrix(c(0.30, 0.10, 0.02, 0.03,
                 0.20, 0.10, 0.01, 0.01,
                 0.04, 0.02, 0.05, 0.03,
                 0.01, 0.02, 0.01, 0.05), ncol = 4, byrow =TRUE)

fm <- as.vector(fm12)

## Marignales pour M
fm1 <- apply(fm12, 1, sum)
fm2 <- apply(fm12, 2, sum)

Emi <- c(sum(k*fm1), sum(k*fm2))
Ex <- a/b*Emi 
Es <- sum(Ex)

## Vecteur (X1, X2)
comb <- expand.grid(0:3, 0:3)
am1 <- a[1]* comb[, 1]
am2 <- a[2]* comb[, 2]
am <- am1 + am2

Fx12 <- function(x1, x2) fm[1] + sum(fm[-1] * pgamma(x1, am1[-1], b)*pgamma(x2, am2[-1], b))
Fx12(2, 5)

## V.a. S
Fs <- function(x) fm[1] + sum(fm[-1] * pgamma(x, am[-1], b))
Fs(400)

Es_test <- sum(fm*am/b)

SL <- function(d) sum(fm*(am/b * pgamma(d, am + 1, b, low=FALSE) - d*pgamma(d, am, b, low=FALSE)))
VaRS <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 400))$min
TVaRS <- function(k){
  VaR <- VaRS(k)
  1/(1 - k)*sum(fm*am/b * pgamma(VaR, am + 1, b, low = F))
}

SL(100) # 3.630182
VaRS(0.95)
TVaRS(0.95)

##### Poisson Teicher composée #####
### Fonction utile
dPoTei <- function(k1, k2, lam1, lam2, a0){
  f <- function(i, j){
    k <- 0:min(i, j)
    sum(dpois(k, a0) * dpois(i - k, lam1 - a0) * dpois(j - k, lam2 - a0))
  }
  outer(k1, k2, Vectorize(f))
}

### Variable
m1 <- 0:20
m2 <- 0:20
lam <- c(3, 2)
g0 <- 1
a <- c(2, 3)
b <- 0.2

### Fonction de densité de la Teicher
fm12 <- dPoTei(m1, m2, lam[1], lam[2], g0)

### Vecteur de densité pij
fm <- as.vector(fm12)

### Vecteur de paramèetre a*
comb <- expand.grid(m1, m2)
am <- colSums(t(as.matrix(comb)) * a)
# ou
am <- as.vector(as.matrix(comb) %*% matrix(a))
# ou
am1 <- a[1]*comb[, 1]
am2 <- a[2]*comb[, 2]
am <- am1 + am2

### Repartition de FS
Fs <- function(x) fm[1] + sum(fm[-1]*pgamma(x, am[-1], b))
Fs(100)



