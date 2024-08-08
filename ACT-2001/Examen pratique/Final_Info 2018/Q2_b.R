# Examen informatique H2018

# Question 2

set.seed(2018)
U <- matrix(runif(300000), ncol =3, byrow =TRUE)
head(U)
X1 <- qgamma(U[, 1], 1/2, 0.05)
X2 <- qgamma(U[, 2], 1.5, 0.15)
X3 <- qgamma(U[, 3], 2.5, 0.25)

# cbind(X1, X2, X3)
u <- X1 + X2 + X3


# b) TVaR
M <- 100000
TVar <- function(x, k)
{
  X1[k*M] + (1/M)*(1/(1 - k))*(sum(X1[90001:M] - X1[k*M]))
}
TVar(X1, 0.9)

### Ajustement ###
TVaR2 <- function(x, k)
{
  sort(x)[k*M] + 1/(1-k) * mean(pmax(x - sort(x)[M*k], 0))
}

TVaR2(X1, 0.9)

### Démarche alternative dans le cas continue ###

VaRX1 <- quantile(X1, 0.9)
VaRX2 <- quantile(X2, 0.9)
VaRX3 <- quantile(X3, 0.9)

TVaR3 <- function(x, k)
{
  VaR <- quantile(x, k)
  mean(x[x > VaR])
}

TVaR3(X1, 0.9)
















