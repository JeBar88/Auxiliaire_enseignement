### ACT-3000
## Fonctions utiles
## Jérémie Barde

##### Erlang généralisée ####
derlgen <- function(x, b){
  P <- sapply(seq_along(b), function(i) prod(b[-i]/(b[-i] - b[i])))
  sum(P*dexp(x, b))
}
perlgen <- function(x, b){
  P <- sapply(seq_along(b), function(i) prod(b[-i]/(b[-i] - b[i])))
  sum(P*pexp(x, b))
}
terlgen <- function(VaR, b){
  P <- sapply(seq_along(b), function(i) prod(b[-i]/(b[-i] - b[i])))
  A <- exp(-b*VaR)*(VaR + 1/b)
  1/(1 - u)*sum(P*A)
}

#### Répartition multivariée lois disrètes ####
#### Répartition multivariée lois disrètes ####
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
repart_nd <- function(f){
  d <- length(dim(f))
  for (i in 1:d) {
    f <- apply(f, (1:d)[-i], cumsum) |> aperm(order(c(i, (1:d)[-i])))
  }
  f
}

#### Aggrégation ####
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
#### Fonction de masse de probabilité multivariée discrète ####
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

#### Comonotonne et antimonotonne ####
dcomo <- function(x1, x2, Fx1, Fx2){
  Fx1 <- c(0, Fx1)
  Fx2 <- c(0, Fx2)
  
  pcomo <- function(x1, x2) pmin(Fx1[x1 + 2], Fx2[x2 + 2])
  f <- function(x1, x2) pcomo(x1, x2) - pcomo(x1-1, x2) - pcomo(x1, x2-1) + pcomo(x1-1, x2-1)
  
  outer(x1, x2, Vectorize(f))
}

danti <- function(x1, x2, Fx1, Fx2){
  Fx1 <- c(0, Fx1)
  Fx2 <- c(0, Fx2)
  
  panti <- function(x1, x2) pmax(Fx1[x1 + 2] + Fx2[x2 + 2] - 1, 0)
  f <- function(x1, x2) panti(x1, x2) - panti(x1-1, x2) - panti(x1, x2-1) + panti(x1-1, x2-1)
  
  outer(x1, x2, Vectorize(f))
}