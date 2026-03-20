### ACT-3000
## Loi Poisson Teicher composée : simulation
## Jérémie Barde

#### Packages ####
library(actuar)

#### Fonction utiles ####
dPoTei <- function(k1, k2, lam1, lam2, a0){
  f <- function(i, j){
    k <- 0:min(i, j)
    sum(dpois(k, a0) * dpois(i - k, lam1 - a0) * dpois(j - k, lam2 - a0))
  }
  outer(k1, k2, Vectorize(f))
}

#### Poisson Teicher composée simulation ####
m <- 1e6
k <- 0:20
lam <- c(2, 3)
a0 <- 1.5
a <- 2
b <- c(0.1, 10)
u <- 0.99

## Densité de la loi Teicher
fm12 <- dPoTei(k, k, lam[1], lam[2], a0)

## Simulation de la fréquence
## Façon 1 -- Intuitive
simul_freq <- function(){
  M1 <- rpois(1, lam[1])
  fm2c1 <- fm12[M1 + 1, ]/dpois(M1, lam[1])
  Fm2c1 <- cumsum(fm2c1)
  M2 <- k[min(which(Fm2c1 >= runif(1)))]
  c(M1, M2)
}
U12 <- t(replicate(m, simul_freq()))

## Façon 1.5 -- Sans fonction
U1 <- rpois(m, lam[1])
fm2c1 <- sapply(U1, function(i) fm12[i + 1, ]/dpois(i, lam[1]))
Fm2c1 <- apply(fm2c1, 2, cumsum)
U2 <- apply(Fm2c1, 2, function(i) k[min(which(i >= runif(1)))])
U12 <- cbind(U1, U2)

## Façon 2 -- Simuler avec la fmp conjoite
comb <- expand.grid(k, k) |> as.matrix()
fm <- as.vector(fm12)
index <- sample(seq_along(fm), m, replace = TRUE, prob = fm)
M12 <- comb[index, ]

## Verif
CovM12 <- cov(M12)
cbind("a_test"=CovM12[1, 2], a0)
cbind("Vm1_test"=CovM12[1, 1], "Vm1"=lam[1],
      "Vm2_test"=CovM12[2, 2], "Vm2"=lam[2])

## Simulation de la sévérité
X1 <- sapply(M12[, 1], function(i) sum(rgamma(i, a, b[1])))
X2 <- sapply(M12[, 2], function(i) sum(rpareto(i, a, b[2])))
X12 <- cbind(X1, X2)

## Simualtion de la somme X1 + X2
S <- rowSums(X12)

# Covariance
CovX12 <- mean(X1 * X2) - mean(X1)*mean(X2)
cbind("CovX12_test"=CovX12, "CovX12"=a0*a/b[1]*b[2]/(a - 1))

# Correlation de pearson
rho <- CovX12/(sd(X1) * sd(X2))

# Espérance et variance
ES <- mean(S)
VS <- var(S)  

# VaR
VaRS <- quantile(S, u)

# TVaR
TVaRS1 <- mean(S[S > VaRS]) # Si continue
TVaRS2 <- VaRS + 1/(1 - u)*mean(pmax(S - VaRS, 0))
cbind(TVaRS1, TVaRS2)  
  
#### Loi discrete en dimension d simulation ####
fm <- array(c(0.083, 0.031, 0.054, 0.054, 0.050, 0.005, 0.035,
                 0.047, 0.051, 0.117, 0.044, 0.011, 0.114, 0.032,
                 0.025, 0.088, 0.028, 0.013, 0.005, 0.003, 0.009,
                 0.014, 0.001, 0.019, 0.024, 0.012, 0.031), dim=rep(3, 3))

# fm <- array(c(0.172125, 0.172125, 0.057375, 0.057375,
#                0.030375, 0.030375, 0.010125, 0.010125,
#                0.057375, 0.057375, 0.019125, 0.019125,
#                0.010125, 0.010125, 0.003375, 0.003375,
#                0.04303125, 0.04303125, 0.01434375, 0.01434375,
#                0.00759375, 0.00759375, 0.00253125, 0.00253125,
#                0.01434375, 0.01434375, 0.00478125, 0.00478125,
#                0.00253125, 0.00253125, 0.00084375, 0.10084375), dim = rep(2, 5))

m <- 1e6
d <- length(dim(fm))
k <- 0:(dim(fm)[1] - 1)

## Marginales
fmi <- sapply(1:d, function(i) apply(fm, i, sum))
Emi <- apply(fmi, 2, function(f) sum(k*f))
Es <- sum(Emi)

## Simulation du Vecteur M
comb <- arrayInd(seq_along(fm), dim(fm)) - 1
fm_vec <- as.vector(fm)
index <- sample(seq_along(fm_vec), m, replace = TRUE, prob = fm)
M <- comb[index, ]

# Repartion
# value <- c(0, 1, 1, 0, 1)
# mean(apply(M, 1, function(i) all(i <= value)))

## Simualtion de la v.a. S
S <- rowSums(M)
ES <- mean(S)
cbind(Es, ES)










