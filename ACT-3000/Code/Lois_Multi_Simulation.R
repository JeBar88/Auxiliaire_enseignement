### ACT-3000
## Simulation de v.a. multivariée et v.a. multivariée composée
## # Loi discrète quelconque
## # Loi Poisson de Teicher
## # Loi Bernoulli CGMR
## # Loi Binomial bivariée de Marshall-Olking
## # Loi gamma bivariée CRMM
## # Loi normle multivariée
## # Cas Comonotonne et antimonotonne
## Jérémie Barde

#### Packages ####
library(actuar)

#### Loi discrète quelconque ####
# fm <- matrix(c(0.5, 0.05, 0.05,
#                  0.05, 0.3, 0.01,
#                  0.01, 0.02, 0.01), ncol = 3, byrow = TRUE)

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

m <- 5e6
d <- length(dim(fm))
k <- 0:(dim(fm)[1] - 1)
a <- c(2, 2, 3)
b <- c(0.1, 0.5, 0.25)
u <- 0.90

### Marginales
fmi <- sapply(1:d, function(i) apply(fm, i, sum))
Fmi <- apply(fmi, 2, cumsum)
Emi <- apply(fmi, 2, function(f) sum(k*f))
Es <- sum(Emi)

### Simulation du Vecteur M (fréquence)
## Façon 1 -- Intuitive
## Code non flexible, on doit adapté pour la dimension
# simul_freq <- function(){
#   M1 <- k[min(which(Fmi[, 1] >= runif(1)))]
#   fm2c1 <- fm[M1 + 1, ]/fmi[M1 + 1, 1]
#   Fm2c1 <- cumsum(fm2c1)
#   M2 <- k[min(which(Fm2c1 >= runif(1)))]
#   c(M1, M2)
# }
# M <- t(replicate(m, simul_freq()))

## Façon 1.5 -- Intuitive Sans fonction
## Code non flexible, on doit adapté pour la dimension
# U <- runif(m)
# M1 <- sapply(U, function(u) k[min(which(Fmi[, 1] >= u))])
# fm2c1 <- sapply(M1, function(i) fm[i + 1, ]/fmi[i + 1, 1])
# Fm2c1 <- apply(fm2c1, 2, cumsum)
# M2 <- apply(Fm2c1, 2, function(i) k[min(which(i >= runif(1)))])
# M <- cbind(M1, M2)

## Façon 2
## Flexible, fonctionne en dimension d
comb <- arrayInd(seq_along(fm), dim(fm)) - 1
index <- sample(seq_along(fm), m, replace = TRUE, prob = fm)
M <- comb[index, ]

# Repartion
# value <- c(0, 1, 1, 0, 1)
# mean(apply(M, 1, function(i) all(i <= value)))

## Calcule covariance et correlation de Pearson
cov(M)
cor(M, method = "pearson")

## Simualtion de la v.a. S
S <- rowSums(M)

# Espérance et variance
ES <- mean(S)
VS <- var(S)  
cbind(Es, ES)

# VaR
VaRS <- quantile(S, u)

# TVaR
TVaRS1 <- mean(S[S > VaRS]) # Si continue, ici fonctionne pas
TVaRS2 <- VaRS + 1/(1 - u)*mean(pmax(S - VaRS, 0))
cbind(TVaRS1, TVaRS2)  

### Simulation du vecteur X (Sévérité)
## Si composée avec B~Bi~Ga(ai, bi). Ex. pour le cas 3d
Exi <- apply(fmi, 2, function(f) sum(k*f)) * a/b
Es <- sum(Exi)

X1 <- sapply(M[, 1], function(m) sum(rgamma(m, a[1], b[1])))
X2 <- sapply(M[, 2], function(m) sum(rgamma(m, a[2], b[2])))
X3 <- sapply(M[, 3], function(m) sum(rgamma(m, a[3], b[3])))
X <- cbind(X1, X2, X3)
# ou
# X <- sapply(1:d, function(i) sapply(M[, i], function(m) sum(rgamma(m, a[i], b[i]))))
# ou (Vraiment plus rapide, mais doit être Gamma)
# X <- matrix(0, m, d)
# for (i in 1:d) {
#   pos <- M[, i] > 0
#   X[pos, i] <- rgamma(sum(pos), M[pos, i]*a[i], b[i])
# }

## Calcule covariance et correlation de Pearson
cov(X)
cor(X, method = "pearson")

## Simualtion de la v.a. S
S <- rowSums(X)

# Espérance et variance
ES <- mean(S)
VS <- var(S)  
cbind(Es, ES)

# VaR
VaRS <- quantile(S, u)

# TVaR
TVaRS1 <- mean(S[S > VaRS]) # Si continue, ici fonctionne
TVaRS2 <- VaRS + 1/(1 - u)*mean(pmax(S - VaRS, 0))
cbind(TVaRS1, TVaRS2)  

#### Loi Poisson de Teicher ####
m <- 1e6
a0 <- 1.5
lam <- c(2, 3)
a <- 3
b <- c(0.1, 10)

Emi <- lam
Exi <- Emi*c(a/b[1], b[2]/(a - 1))
Es <- sum(Exi)


### Simulation du Vecteur M (Fréquence)
K1 <- rpois(m, lam[1] - a0)
K2 <- rpois(m, lam[2] - a0)
K0 <- rpois(m, a0)

M1 <- K1 + K0
M2 <- K2 + K0
M <- cbind(M1, M2)

### Verif
cbind("Emi_test"=apply(M, 2, mean), "Emi"=lam)
cbind("a0_test"=cov(M1, M2), a0)

### Simulation du Vecteur X (Sévérité)
X1 <- sapply(M[, 1], function(i) sum(rgamma(i, a, b[1])))
X2 <- sapply(M[, 2], function(i) sum(rpareto(i, a, b[2])))
X <- cbind(X1, X2)

### Simualtion de la somme X1 + X2
S <- rowSums(X)

### Covariance et correlation de Pearson
CovX <- cov(X)[1, 2]
rho <- cor(X, method = "pearson")[1, 2]

### Vérif
cbind("Exi_test"=apply(X, 2, mean), Exi)
cbind("CovX12_test"=CovX, "CovX12"=a0*a/b[1]*b[2]/(a - 1))
cbind("Es_test"=mean(S), Es)

#### Loi Bernouilli CGMR ####
m <- 1e6
d <- 3
q0 <- 0.1
q <- c(0.2, 0.4, 0.7)

Eii <- q

### Simulation du vecteur I
qt <- 1 - (1 - q)/(1 - q0)

J <- sapply(c(q0, qt), function(i) rbinom(m, 1, i))
I <- apply(J[, -1], 2, function(j) pmin(j + J[, 1], 1))

### Covariance
comb <- combn(d, 2)
Covi <- mapply(function(i, j) (1 - q0)*prod(qt[c(i, j)]) + q0 - prod(q[c(i, j)]), comb[1, ], comb[2, ])
names(Covi) <- paste0("Cov", comb[1, ], comb[2, ])

CovI <- cov(I)[upper.tri(cov(I))] 

### Simualtion de la v.a. S
S <- rowSums(I)

### Verif
cbind("Emi_test"=apply(I, 2, mean), "Emi"=Eii)
cbind("Covi_test"=CovI, Covi)
cbind("Es_test"=mean(S), "Es"=sum(q))

#### Loi Binomial bivariée de Marshall-Olking ####
m <- 1e6
n <- 5
q <- c(0.2, 0.35)
p11 <- 0.05

p01 <- q[2] - p11
p10 <- q[1] - p11
p00 <- 1 - p01 - p10 - p11
fi <- matrix(c(p00, p01,
               p10, p11), ncol = 2, byrow = TRUE)

### Simulation du vecteur I et In
In <- replicate(n, {
  comb <- arrayInd(seq_along(fi), dim(fi)) - 1
  index <- sample(seq_along(fi), m, replace = TRUE, prob = fi)
  I <- comb[index,]
}, simplify = FALSE)

### Simulation du vecteur M
M <- Reduce('+', In)

### Simualtion de la v.a. S
S <- rowSums(M)

### Verif
cbind("Eii_test"=apply(M, 2, mean), "Emi"=n*q)
cbind("Covi_test"=cov(M)[1, 2], "Covi"=n*(p11 - prod(q)))
cbind("Es_test"=mean(S), "Es"=n*sum(q))

#### Loi gamma CRMM ####
m <- 1e6
d <- 5
g0 <- 1
lam <- c(2, 3)
a <- c(1.5, 3, 2, 2.5, 5)
b <- c(0.1, 0.5, 1, 0.25, 0.8)

Exi <- a/b
Es <- sum(Exi)

### Simulation du Vecteur X (Fréquence)
Y0 <- rgamma(m, g0, 1) 
Y <- sapply(a, function(i) rgamma(m, i - g0, 1))
X <- sapply(1:d, function(i) 1/b[i]*(Y[, i] + Y0))

### Covariance
comb <- combn(d, 2)
Covx <- mapply(function(i, j) g0/prod(b[c(i, j)]), comb[1, ], comb[2, ])
names(Covx) <- paste0("Cov", comb[1, ], comb[2, ])

CovX <- mapply(function(i, j) cov(X)[i, j], comb[1, ], comb[2, ])
names(CovX) <- paste0("Cov", comb[1, ], comb[2, ])

### Simulation de la v.a. S
S <- rowSums(X)

### Verif
cbind("Ex_test"=apply(X, 2, mean), "Ex"=a/b)
cbind("Covx_test"=CovX, "Covx"=Covx)
cbind("Es_test"=mean(S), "Es"=sum(a/b))

#### Loi normle multivariée ####
m <- 1e6
d <- 3
mu <- c(15, 10, 16)
sig <- matrix(c(4, 2.4, -2,
              2.4, 9, -4.5,
              -2, -4.5, 25), ncol = 3, byrow = TRUE)

### Matrice de décompositoin de Choleski
B <- chol(sig)

### Simualtion du vecteur X
Z <- replicate(d, rnorm(m))
X <- sweep(Z %*% B, 2, mu, "+")

### Simulation de la v.a. S
S <- rowSums(X)

### Verif
cbind("Ex_test"=apply(X, 2, mean), "Ex"=mu)
list("Covx_test"=cov(X), "Covx"=sig)
cbind("Es_test"=mean(S), "Es"=sum(mu))



















#### Cas Comonotonne et antimonotonne ####
#### M1~Po(5) et M2~BinNeg(5, 0.2)
m <- 1e6
lam <- 5
r <- 5
q <- 0.2

### Espérance théorique
Emi <- c(lam, r*(1 - q)/q) 
Es <- sum(Emi)

### Cas comonotonne
## Simuler 1 million de réalisation du vecteur M et v.a. S
U <- runif(1e6)
M1 <- qpois(U, lam)
M2 <- qnbinom(U, r, q)
S <- M1 + M2

## Covariance et correlation de Pearson
cov(M1, M2)
cor(M1, M2, method = "pearson")/(sd(M1)*sd(M2))

### Cas antinomotonne
## Simuler 1 million de réalisation du vecteur M et v.a. S
U <- runif(1e6)
M1 <- qpois(U, lam)
M2 <- qnbinom(1 - U, r, q)
S <- M1 + M2

## Covariance et correlation de Pearson
cov(M1, M2)
cor(M1, M2, method = "pearson")/(sd(M1)*sd(M2))














