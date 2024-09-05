### ACT-3000
## numéro additionelle
## Jérémie Barde

### Fonctions
dcomo <- function(x1, x2, Fx1, Fx2){
  fx <- matrix(numeric(length(x1) * length(x2)), ncol = length(x2))
  pcomo <- function(x1, x2) {
    ifelse(x1 < 0 | x2 < 0, 0, pmin(Fx1[x1 + 1], Fx2[x2 + 1]))
  }
  
  for (i in x1){
    for (j in x2){
      fx[i + 1, j + 1] <-  pcomo(i, j) - pcomo(i-1, j) - pcomo(i, j-1) + pcomo(i-1, j-1)
    }
  }
  fx
}
danti <- function(x1, x2, Fx1, Fx2){
  fx <- matrix(numeric(length(x1) * length(x2)), ncol = length(x2))
  panti <- function(x1, x2) {
    ifelse(x1 < 0 | x2 < 0, 0, pmax(Fx1[x1 + 1] + Fx2[x2 + 1] - 1, 0))
  }
  
  for (i in x1){
    for (j in x2){
      
      fx[i + 1, j + 1] <-  panti(i, j) - panti(i-1, j) - panti(i, j-1) + panti(i-1, j-1)
    }
  }
  fx
}
dn <- function(k, f){
  
  fs <- function(s) sum(sapply(0:s, function(i) f[i + 1, s - i + 1]))
  sapply(k, fs)
  
}

### Variables
r <- c(3, 2)
q <- c(0.6, 0.8)
k <- 0:30

### Fonctions de répartition marginale
Fx1 <- pnbinom(k, r[1], q[1])
Fx2 <- pnbinom(k, r[2], q[2])

### (a)
## i)
pcomo <- function(x1, x2) {
  ifelse(x1 < 0 | x2 < 0, 0, pmin(Fx1[x1 + 1], Fx2[x2 + 1]))
}
pcomo(2, 3)

## ii)
fx12C <- dcomo(k, k, Fx1, Fx2)
fx12C[4, 2]

## iii)
EX12 <- sum(outer(k, k, '*')*fx12C)

## iv)
EX <- r*(1-q)/q
varX <- r*(1-q)/q^2
covX12 <- EX12 - prod(EX)
rho <-  covX12/sqrt(prod(varX))

### (b)
## i)
panti <- function(x1, x2) {
  ifelse(x1 < 0 | x2 < 0, 0, pmin(Fx1[x1 + 1], Fx2[x2 + 1]))
}

panti(2, 3)

## ii)
fx12A <- danti(k, k, Fx1, Fx2)
fx12A[2, 2]

## iii)
EX12 <- sum(outer(k, k, '*')*fx12A)

## iv)
EX <- r*(1-q)/q
varX <- r*(1-q)/q^2
covX12 <- EX12 - prod(EX)
rho <-  covX12/sqrt(prod(varX))










