### ACT-3000
## Exemple de base v.a discrètes multivariée
## # Loi discrète quelconque
## # Loi Poisson de Teicher
## # Loi Bernoulli CGMR
## # Loi Binomial bivariée de Marshall-Olking
## Jérémie Barde

#### Package utile ####
library(actuar)

#### Fonction utile ####
repart_2d <- function(x, y, fxy){
  f <- function(i, j) sum(fxy[0:i + 1, 0:j + 1])
  outer(x, y, Vectorize(f))
}
repart_3d <- function(k, d, f) {
  res <- array(0, dim(f))
  
  for (i in k) {
    for (j in k) {
      for (l in k) {
        res[i + 1, j + 1, l + 1] <- sum(f[0:i + 1, 0:j + 1, 0:l + 1])
      }
    }
  }
  res
}

# Explication: On calule les sommmes cumulatives en fixant les dimensions, comme 
# apply mélange les dimensions on utilise aperm pour les remettres dans le bonne autre.
repart_nd <- function(f) {
  d <- length(dim(f))
  for (i in 1:d) {
    f <- apply(f, (1:d)[-i], cumsum) |> aperm(order(c(i, (1:d)[-i])))
  }
  f
}

ds_2d <- function(k, f){
  fs <- function(s) sum(sapply(0:s, function(i) f[i + 1, s - i + 1]))
  sapply(k, fs)
}
ds_3d <- function(k, f){
  fs <- function(s){
    res <- 0
    for (k1 in 0:s) {
      for (k2 in 0:(s - k1)) {
        res <- res + f[k1 + 1, k2 + 1, (s - k1 - k2) + 1]
      }
    }
    res
  }
  sapply(k, fs)
}
ds <- function(f) {
  index <- arrayInd(seq_along(f), dim(f)) - 1
  sums <- rowSums(index)
  sapply(0:max(sums), function(s) sum(f[sums == s]))
}

dPoTeibiv <- function(k1, k2, lam, a0){
  f <- function(i, j){
    k <- 0:min(i, j)
    sum(dpois(k, a0) * dpois(i - k, lam[1] - a0) * dpois(j - k, lam[2] - a0))
  }
  outer(k1, k2, Vectorize(f))
}
dPoTeibiv_old <- function(m1, m2, lam1, lam2, a0){
  f <- matrix(numeric(length(m1) * length(m2)), ncol = length(m2))

  for (i in m1){
    for (j in m2){

      k <- 0:min(i, j)

      f[i + 1, j + 1] <-  sum(dpois(k, a0) * dpois(i - k, lam1 - a0) * dpois(j - k, lam2 - a0))
    }
  }
  f
}
dBernCGMR <- function(qt, q0) {
  f <- (1 - q0) * Reduce(outer, lapply(qt, function(p) c(1 - p, p)))
  f[length(f)] <- f[length(f)] + q0
  f
}
dBernCGMR_old <- function(qt, q0){
  d <- length(qt)
  perm <- expand.grid(rep(list(0:1), d))
  pr <- apply(perm, 1, function(k) prod(dbinom(k, 1, qt)))
  f_vec <- (1 - q0)*pr + q0*I(apply(perm, 1, prod) == 1)
  f <- array(f_vec, dim = rep(2, d))
  f
}
dBinBivMO <- function(k1, k2, n, q, p11) {
  p10 <- q[1] - p11
  p01 <- q[2] - p11
  p00 <- 1 - p01 - p10 - p11
  f <- function(i, j) {
    l <- max(i + j - n, 0):min(i, j)
    sum((factorial(n)*p11^l*p10^(i - l)*p01^(j - l)*p00^(n - i - j + l))/
          (factorial(l)*factorial(i - l)*factorial(j - l)*factorial(n - i - j + l)))
  }
  outer(k1, k2, Vectorize(f))
}

#### Loi discrète quelconque ####
### Exemple 1 -- Loi discrete 2d ###
fm12 <- matrix(c(0.5, 0.05, 0.05,
                 0.05, 0.3, 0.01,
                 0.01, 0.02, 0.01), ncol = 3, byrow = TRUE)
k <- 0:2

## Calculer les marginales
fm1 <- rowSums(fm12)
fm2 <- colSums(fm12)

## Calcule EMi et VarMi
Em1 <- sum(fm1 * k)
Em2 <- sum(fm2 * k)
Em12 <- sum(outer(k, k, "*")*fm12)

Vm1 <- sum(k^2*fm12) - Em1^2
Vm2 <- sum(k^2*fm12) - Em2^2

## Coveriance et correlation de Pearson
cov <- Em12 - Em1*Em2
pp <- cov/sqrt(Vm1*Vm2)

## Fonction de repartition
Fm12 <- repart_2d(k, k, fm12)
# Fm12 <- repart_nd(fm12)

## V.a. S = M1 + M2
kmax <- max(k)*2
fm12_0 <- matrix(0, ncol=kmax + 1, nrow=kmax + 1)
fm12_0[k + 1, k + 1] <- fm12

fs <- ds_2d(0:kmax, fm12_0)
# ou
fs <- ds(fm12)
c("ES_test"=sum(0:kmax*fs), "ES"=Em1 + Em2)

### Exemple 2 - nd ###
# pmf d'une Bernouilli CGMR avec q0 <- 0.1 et qt <- c(0.5, 0.25, 0.15, 0.25, 0.20)
# fm <- matrix(c(0.5, 0.05, 0.05,
#                0.05, 0.3, 0.01,
#                0.01, 0.02, 0.01), ncol = 3, byrow = TRUE)

# fm <- array(c(0.083, 0.031, 0.054, 0.054, 0.050, 0.005, 0.035,
#                  0.047, 0.051, 0.117, 0.044, 0.011, 0.114, 0.032,
#                  0.025, 0.088, 0.028, 0.013, 0.005, 0.003, 0.009,
#                  0.014, 0.001, 0.019, 0.024, 0.012, 0.031), dim=rep(3, 3))

fm <- array(c(0.172125, 0.172125, 0.057375, 0.057375,
               0.030375, 0.030375, 0.010125, 0.010125,
               0.057375, 0.057375, 0.019125, 0.019125,
               0.010125, 0.010125, 0.003375, 0.003375,
               0.04303125, 0.04303125, 0.01434375, 0.01434375,
               0.00759375, 0.00759375, 0.00253125, 0.00253125,
               0.01434375, 0.01434375, 0.00478125, 0.00478125,
               0.00253125, 0.00253125, 0.00084375, 0.10084375), dim = rep(2, 5))

d <- length(dim(fm))
k <- 0:(dim(fm)[1] - 1)

## Calculer les marginales
# Marginale univariée
fmi <- sapply(1:d, function(i) apply(fm, i, sum))

# Marginale bivariée
comb <- combn(d, 2)
fmij <- mapply(function(i, j) apply(fm, c(i, j), sum), comb[1 ,], comb[2, ], SIMPLIFY = FALSE)
names(fmij) <- paste0('fm', comb[1, ], comb[2, ])
sapply(fmij, sum)

# Marginale trivariée
# comb3 <- combn(d, 3)
# fmijl <- mapply(function(i, j, l) apply(fm, c(i, j, l), sum), comb3[1 ,], comb3[2, ], comb3[3, ], SIMPLIFY = FALSE)
# names(fmijl) <- paste0('fm', comb3[1, ], comb3[2, ], comb3[3, ])
# sapply(fmijl, sum)

# Marginale quadrivariée
# comb4 <- combn(d, 4)
# fmijlh <- mapply(function(i, j, l, h) apply(fm, c(i, j, l), sum), comb4[1 ,], comb4[2, ], comb4[3, ], comb4[4, ], SIMPLIFY = FALSE)
# names(fmijlh) <- paste0('fm', comb4[1, ], comb4[2, ], comb4[3, ], comb4[4, ])
# sapply(fmijlh, sum)


## Calcule EMi et VarMi
Emi <- apply(fmi, 2, function(f) sum(k*f))
Em_prod <- sum(Reduce("%o%", rep(list(k), d))*fm)

Vmi <- sapply(1:d, function(i) sum(k^2*fmi[, i]) - Emi[i]^2)

## Coveriance et correlation de Pearson
covmij <- sapply(1:ncol(comb), function(i) sum(outer(k, k)*fmij[[i]]) - prod(Emi[comb[, i]]))
names(covmij) <- paste0('covm', comb[1, ], comb[2, ])

ppmij <- sapply(1:ncol(comb), function(i) covmij[i]/sqrt(prod(Vmi[comb[, i]])))
names(ppmij) <- paste0('pp', comb[1, ], comb[2, ])

## Fonction de repartition
Fm <- repart_nd(fm)

## V.a. S
kmax <- max(k)*d
k2 <- 0:kmax
fs <- ds(fm)
sum(fs)
c("Es_test"=sum(k2*fs), "Es"=sum(Emi))

#### Loi Poisson bivariée de Teicher ####
nfft <- 2^13
k <- 0:(nfft - 1)
m1 <- 0:50
m2 <- 0:50
a0 <- 1.5
lam <- c(2, 3)

### pmf et cdf
fm12 <- dPoTeibiv(m1, m2, lam, a0)
Fm12 <- repart_nd(fm12)

### Marginales
fm1 <- rowSums(fm12)
fm2 <- colSums(fm12)

### Covariance et correlation pearson
cov <- a0
pp <- a0/sqrt(prod(lam))

### Espérance conditionnelle
lam[1] - a0 + a0/lam[2]*0:10
lam[2] - a0 + a0/lam[1]*0:10

### Espérance de M1M2
Em_prod <- sum(outer(m1, m2, '*') * fm12)
cov <- Em_prod - prod(lam)
## E[min(M1, 7) * M2 * 1[M2 < = 6]]
sum(outer(pmin(m1, 7), m2 * I(m2 <= 6), '*')*fm12)

### V.a. N = M1 + M2
Es <- sum(lam)

## Méthode 1 : N = K1 + K2 + 2K0
fk12 <- dpois(k, sum(lam - a0))
f2k0 <- dpois(k/2, a0)
fnt <- fft(fk12) * fft(f2k0)
fn <- Re(fft(fnt, TRUE))/nfft
cbind("Es_test"=sum(k*fn), Es)

## Méthode 2 : utiliser directement fgp (Recommandé de Jé)
fi <- numeric(nfft)
fi[2] <- 1 
fnt <- exp((lam[1] - a0)*(fft(fi) - 1)) * exp((lam[2] - a0)*(fft(fi) - 1)) * exp(a0*(fft(fi)^2 - 1))
fn2 <- Re(fft(fnt, TRUE))/nfft
cbind("Es_test"=sum(k*fn2), Es)

## Méthode 3 : Connaitre la loi de N (Preuve)
lamN <- sum(lam) - a0
fb <- numeric(nfft)
fb[2:3] <- c((sum(lam) - 2*a0)/lamN, a0/lamN)
fnt <- exp(lamN*(fft(fb) - 1))
fn3 <- Re(fft(fnt, T))/nfft
cbind("Es_test"=sum(k*fn3), Es)

## Méthode 4 : 'a la main'
kmax <- max(m1)*2
k2 <- 0:kmax
fn4 <- ds(fm12)
cbind("Es_test"=sum(k2*fn4), Es)

Fn <- cumsum(fn)
u <- 0.99
VaRN <- k[min(which(Fn > u))]
TVaRN <- VaRN + 1/(1 - u)*sum(pmax(k - VaRN, 0)*fn)
cbind(VaRN, TVaRN) 

##### Loi Bernoulli CGMR ####
nfft <- 2^3
k <- 0:(nfft - 1) 
d <- 3
q0 <- 0.1
q <- c(0.2, 0.4, 0.7)

### Calculer les q tildes
qt <- 1 - (1 - q)/(1 - q0)

### pmf et cdf
fi <- dBernCGMR(qt, q0)
Fi <- repart_nd(fi)

### Calculer les marginales
# Marginale univariée
fmi <- sapply(1:d, function(i) apply(fi, i, sum))

### Marginale bivariée
comb <- combn(d, 2)
fij <- mapply(function(i, j) apply(fi, c(i, j), sum), comb[1 ,], comb[2, ], SIMPLIFY = FALSE)
names(fij) <- paste0('fm', comb[1, ], comb[2, ])
sapply(fij, sum)

### Espérance
Eii <- q

### Coveriance et correlation de Pearson
Covij <- sapply(1:d, function(i) sum(fij[[i]][2, 2]) - prod(Eii[comb[, i]]))
names(Covij) <- paste0('Cov', comb[1, ], comb[2, ])

## Verif
Covi <- mapply(function(i, j) (1 - q0)*qt[i]*qt[j] + q0 - q[i]*q[j], comb[1, ], comb[2, ])
names(Covi) <- paste0("Cov", comb[1, ], comb[2, ])

### V.a. S
Es <- sum(q)
fgpI <- function(t, q) (1 - q) + q*t

## Méthode   : Utiliser la fgp
f <- numeric(nfft) ; f[2] <- 1 
fst <- (1 - q0)*apply(outer(fft(f), qt, fgpI), 1, prod) + q0*fft(f)^d
fs <- Re(fft(fst, TRUE))/nfft
cbind("Es_test"=sum(k*fs), Es)

## Méthode 2 : utiliser la méthode brute
kmax <- d
k2 <- 0:kmax
fs <- ds(fi)
cbind("Es_test"=sum(k2*fs), Es)

#### Loi Binomial bivariée de Marshall-Olking ####
nfft <- 2^10
n <- 3
q <- c(0.2, 0.4)
p11 <- 0.1
k <- 0:n

### pmf et cdf
fm <- dBinBivMO(k, k, n, q, p11)
Fm <- repart_nd(fm)

### Marginales
fm1 <- rowSums(fm)
fm2 <- colSums(fm)

### Espérance
Emi <- n*q

### Coveriance et correlation de Pearson
Covm <- sum(outer(k, k, '*')*fm) - prod(Emi)

## Verif
Covm_test <- n*(p11 - prod(q))
cbind(Covm_test, Covm)

### V.a. S
Es <- sum(Emi)

## Méthode 1: Utiliser la fgp
p10 <- q[1] - p11
p01 <- q[2] - p11
p00 <- 1 - p01 - p10 - p11

fi <- numeric(nfft) ; fi[2] <- 1
fit <- fft(fi)
fst <- (p00 + p10*fit + p01*fit + p11*fit^2)^n
fs <- Re(fft(fst, TRUE))/nfft
cbind("Es_test"=sum(0:(nfft - 1)*fs), Es)

## Méthode 2 : utiliser la méthode brute
k2 <- 0:6
fs <- ds(fm)
cbind("Es_test"=sum(k2*fs), Es)


