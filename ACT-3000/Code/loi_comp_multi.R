### ACT-3000
## Loi poisson Teicher composée Expo/Gamma, bêta indentique
## Jérémie Barde

##### Loi discrète arbitraire #####
a <- c(2.1, 3.4)
b <- 0.1

fm12 <- matrix(c(0.30, 0.10, 0.02, 0.03,
                 0.20, 0.10, 0.01, 0.01,
                 0.04, 0.02, 0.05, 0.03,
                 0.01, 0.02, 0.01, 0.05), ncol = 4, byrow =TRUE)

fm <- as.vector(fm12)

cas <- expand.grid(0:3, 0:3)[-1, ]
am1 <- a[1]* cas[, 1]
am2 <- a[2]* cas[, 2]
am <- am1 + am2

Fx12 <- function(x1, x2) fm[1] + sum(fm[-1] * pgamma(x1, am1, b)*pgamma(x2, am2, b))
Fs <- function(x) fm[1] + sum(fm[-1] * pgamma(x, am, b))

SL <- function(d) sum(fm[-1]*(am/b * pgamma(d, am + 1, b, low=F) - d*pgamma(d, am, b, low=F)))
VaRS <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 500))$min
TVaRS <- function(k){
  VaR <- VaRS(k)
  1/(1-k)*sum(fm[-1]*am/b * pgamma(VaR, am + 1, b, low = F))
}

SL(100)
VaRS(0.95)
TVaRS(0.95)

##### Avec poisson Teicher #####
### Fonction utile
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
cas <- expand.grid(m1, m2)[-1, ]
am <- colSums(t(as.matrix(cas)) * a)
# ou
am <- as.vector(as.matrix(cas) %*% matrix(a))
# ou
am1 <- a[1]*cas[, 1]
am2 <- a[2]*cas[, 2]
am <- am1 + am2

### Repartition de FS
FS <- function(x) fm[1] + sum(fm[-1]*pgamma(x, am, b))
FS(100)

max(3/4*5 - 2, 0)
max(5 - 4/3*2, 0)



