### ACT-3000
## Exemple pour les cas comonotonne et animonotonne
## Jérémie Barde

### Fonction utile ###
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

#### Comonotone et antimonotone ####
lam <- c(2, 3)
m1 <- 0:50
m2 <- m1

Fm1 <- ppois(m1, lam[1])
Fm2 <- ppois(m2, lam[2])

fm12c <- dcomo(m1, m2, Fm1, Fm2)
Em12c <- sum(outer(m1, m2, "*") * fm12c)
covc <- Em12c - prod(lam)
ppc <- covc/sqrt(prod(lam))
cbind(Em12c, covc, ppc)

fm12a <- danti(m1, m2, Fm1, Fm2)
Em12a <- sum(outer(m1, m2, "*") * fm12a)
cova <- Em12a - prod(lam)
ppa <- cova/sqrt(prod(lam))
cbind(Em12a, cova, ppa)

fs <- dn(m1, fm12c)
sum(fs)
