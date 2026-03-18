### ACT-3000
## Exemple de base pour une v.a discrète bivariée
## Jérémie Barde

#### Fonction utile ####
repart <- function(x, y, fxy){
  f <- function(i, j) sum(fxy[0:i + 1, 0:j + 1])
  outer(x, y, Vectorize(f))
}

ds <- function(k, f){
  fs <- function(s) sum(sapply(0:s, function(i) f[i + 1, s - i + 1]))
  sapply(k, fs)
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
Fm12 <- repart(k, k, fm12)

## V.a. S = M1 + M2
kmax <- max(k)*2 + 1
fm12_0 <- matrix(0, ncol=kmax, nrow=kmax)
fm12_0[k + 1, k + 1] <- fm12

fs <- ds(0:4, fm12_0)
c("ES_test"=sum(0:4*fs), "ES"=Em1 + Em2)





