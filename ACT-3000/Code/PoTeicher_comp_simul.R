### ACT-3000
## Loi Poisson Teicher composée : simulation
## Jérémie Barde

##### Packages #####
library(actuar)

##### Fonction utiles #####
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

### Poisson Teicher composée simulation ###
n <- 1e5
lam <- c(2, 3)
a <- 1
m1 <- 0:20
m2 <- m1
k <- 0.95

### Densité de la loi Teicher
fm12 <- dPoTei(m1, m2, lam[1], lam[2], a)

### Simulation de la fréquence

M1 <- rpois(1, lam[1])
fm2c1 <- fm12[M1 + 1, ]/dpois(M1, lam[1])
Fm2c1 <- cumsum(fm2c1)
M2 <- m2[min(which(Fm2c1 >= runif(1)))]
c(M1, M2)

## Façon 1
simul_freq <- function(){
  M1 <- rpois(1, lam[1])
  fm2c1 <- fm12[M1 + 1, ]/dpois(M1, lam[1])
  Fm2c1 <- cumsum(fm2c1)
  M2 <- m2[min(which(Fm2c1 >= runif(1)))]
  c(M1, M2)
}
M12 <- t(replicate(n, simul_freq()))

# ## Façon 2
# U1 <- rpois(n, lam[1])
# fm2c1 <- sapply(U1, function(i) fm12[i + 1, ]/dpois(i, lam[1]))
# Fm2c1 <- apply(fm2c1, 2, cumsum)
# U2 <- apply(Fm2c1, 2, function(i) m2[min(which(i >= runif(1)))])
# U12 <- cbind(U1, U2)

### Simulation de la sévérité
X1 <- sapply(M12[, 1], function(i) sum(rgamma(i, 2, 0.1)))
X2 <- sapply(M12[, 2], function(i) sum(rpareto(i, 2, 10)))
X12 <- cbind(X1, X2)

### Simualtion de la somme X1 + X2
S <- rowSums(X12)

## Covariance
covX12 <- mean(X1 * X2) - mean(X1)*mean(X2)

## Correlation de pearson
rho <- covX12/(sd(X1) * sd(X2))

## Espérance et variance de S
espS <- mean(S)
varS <- var(S)  

## VaR de S
VaRS <- quantile(S, k)

## TVaR de S
TVaRS1 <- mean(S[S > VaRS])
TVaRS2 <- VaRS + 1/(1 - k)*mean(pmax(S - VaRS, 0))
cbind(TVaRS1, TVaRS2)  
  








