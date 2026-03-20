### ACT-3000
## Exemple de base pour une v.a discrète bivariée
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
repart <- function(f) {
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

ds <- function(k, d, f) {
  index <- arrayInd(seq_along(f), dim(f)) - 1
  sums <- rowSums(index)
  sapply(k, function(s) sum(f[sums == s]))
}

#### Exemple 1 -- Loi discrete ####
fm12 <- matrix(c(0.5, 0.05, 0.05,
                 0.05, 0.3, 0.01,
                 0.01, 0.02, 0.01), ncol = 3, byrow = T)
k <- 0:2

## Calculer les marginales
fm1 <- apply(fm12, 1, sum)
fm2 <- apply(fm12, 2, sum)

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
fs <- ds(0:kmax, d, fm12)
c("ES_test"=sum(0:kmax*fs), "ES"=Em1 + Em2)

#### Exemple 2 -- Loi discrete dimension d ####
fm123 <- array(c(0.083, 0.031, 0.054, 0.054, 0.050, 0.005, 0.035, 
                 0.047, 0.051, 0.117, 0.044, 0.011, 0.114, 0.032,
                 0.025, 0.088, 0.028, 0.013, 0.005, 0.003, 0.009,
                 0.014, 0.001, 0.019, 0.024, 0.012, 0.031), dim=rep(3, 3))
d <- 3
k <- 0:2

## Calculer les marginales
fmi <- sapply(1:d, function(i) apply(fm123, i, sum))

comb <- combn(d, 2)
fmij <- mapply(function(i, j) apply(fm123, c(i, j), sum), comb[1 ,], comb[2, ], SIMPLIFY = FALSE)
names(fmij) <- paste0('fm', comb[1, ], comb[2, ])


## Calcule EMi et VarMi
Emi <- apply(fmi, 2, function(f) sum(k*f))
Em123 <- sum(Reduce("%o%", rep(list(k), d))*fm123)

Vmi <- sapply(1:d, function(i) sum(k^2*fmi[, i]) - Emi[i]^2)

## Coveriance et correlation de Pearson
covmij <- sapply(1:ncol(comb), function(i) sum(outer(k, k)*fmij[[i]]) - prod(Emi[comb[, i]]))
names(covmij) <- paste0('covm', comb[1, ], comb[2, ])

ppmij <- sapply(1:d, function(i) covmij[i]/sqrt(prod(Vmi[comb[, i]])))
names(ppmij) <- paste0('pp', comb[1, ], comb[2, ])


## Fonction de repartition
Fm123 <- repart_3d(k, d, fm123)
# Fm123 <- repart_nd(fm123)

## V.a. S = M1 + M2 + M3
kmax <- max(k)*3
fm123_0 <- array(0, dim = rep(kmax + 1, d))
fm123_0[k + 1, k + 1, k + 1] <- fm123

k2 <- 0:kmax
fs <- ds_3d(k2, fm123_0)
# ou
fs <- ds(k2, d, fm123)
c("Es_test"=sum(k2*fs), "Es"=sum(Emi))

#### Exemple 3 - dimension 5 ####
# pmf d'une Bernouilli CGMR avec q0 <- 0.1 et qt <- c(0.5, 0.25, 0.15, 0.25, 0.20)
fm5 <- array(c(0.172125, 0.172125, 0.057375, 0.057375,
              0.030375, 0.030375, 0.010125, 0.010125,
              0.057375, 0.057375, 0.019125, 0.019125,
              0.010125, 0.010125, 0.003375, 0.003375,
              0.04303125, 0.04303125, 0.01434375, 0.01434375,
              0.00759375, 0.00759375, 0.00253125, 0.00253125,
              0.01434375, 0.01434375, 0.00478125, 0.00478125,
              0.00253125, 0.00253125, 0.00084375, 0.10084375), dim = rep(2, 5))
d <- 5
k <- 0:1

## Calculer les marginales
# Marginale univariée
fmi <- sapply(1:d, function(i) apply(fm5, i, sum))

# Marginale bivariée
comb <- combn(d, 2)
fmij <- mapply(function(i, j) apply(fm5, c(i, j), sum), comb[1 ,], comb[2, ], SIMPLIFY = FALSE)
names(fmij) <- paste0('fm', comb[1, ], comb[2, ])
sapply(fmij, sum)

# Marginale trivariée
comb3 <- combn(d, 3)
fmijl <- mapply(function(i, j, l) apply(fm5, c(i, j, l), sum), comb3[1 ,], comb3[2, ], comb3[3, ], SIMPLIFY = FALSE)
names(fmijl) <- paste0('fm', comb3[1, ], comb3[2, ], comb3[3, ])
sapply(fmijl, sum)

# Marginale quadrivariée
comb4 <- combn(d, 4)
fmijlh <- mapply(function(i, j, l, h) apply(fm5, c(i, j, l), sum), comb4[1 ,], comb4[2, ], comb4[3, ], comb4[4, ], SIMPLIFY = FALSE)
names(fmijlh) <- paste0('fm', comb4[1, ], comb4[2, ], comb4[3, ], comb4[4, ])
sapply(fmijlh, sum)


## Calcule EMi et VarMi
Emi <- apply(fmi, 2, function(f) sum(k*f))
Em12345 <- sum(Reduce("%o%", rep(list(k), d))*fm5)

Vmi <- sapply(1:d, function(i) sum(k^2*fmi[, i]) - Emi[i]^2)

## Coveriance et correlation de Pearson
covmij <- sapply(1:ncol(comb), function(i) sum(outer(k, k)*fmij[[i]]) - prod(Emi[comb[, i]]))
names(covmij) <- paste0('covm', comb[1, ], comb[2, ])



ppmij <- sapply(1:ncol(comb), function(i) covmij[i]/sqrt(prod(Vmi[comb[, i]])))
names(ppmij) <- paste0('pp', comb[1, ], comb[2, ])


## Fonction de repartition
Fm5 <- repart(fm5)

## V.a. S = M1 + M2 + M3
kmax <- max(k)*d
fm5_0 <- array(0, dim = rep(kmax + 1, d))
fm5_0[k + 1, k + 1, k + 1, k + 1, k + 1] <- fm5

k2 <- 0:kmax
fs <- ds(k2, d, fm5)
sum(fs)
c("Es_test"=sum(k2*fs), "Es"=sum(Emi))



