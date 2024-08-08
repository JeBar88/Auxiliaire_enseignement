### ACT-2001: Introduction à l'actuariat II
## Examen partiel informatique 2017 No 3

# b)
library(actuar)
m <- 1E5

lam <- c(1/50, 20, 100)
tau <- c(0.5, 2.5, 2)
a <- 2.5

set.seed(20160419)
matU <- matrix(runif(m * 3), ncol = 3, byrow = TRUE)
matU[c(1, 2, 3, 100000), ]

X <- matrix(numeric(100000 * 3), ncol = 3)
X[, 1] <- qweibull(matU[, 1], tau[1], 1/lam[1])
X[, 2] <- lam[2] * ((1/matU[, 2]) - 1)^(-1/tau[2])
X[, 3] <- lam[3]^(1/tau[3]) * (((1 - matU[, 3])^-(1/a)) - 1)^(1/tau[3])
X[c(1, 2, 3, 4, 100000), ]

# c)
VaRX1 <- quantile(X[,1], 0.9)
VaRX2 <- quantile(X[,2], 0.9)
VaRX3 <- quantile(X[,3], 0.9)

mean(X[,1][X[,1] > VaRX1])
mean(X[,2][X[,2] > VaRX2])
mean(X[,3][X[,3] > VaRX3])

# d)
S <- X[,1] + X[,2] + X[,3]

# e)
VaRS <- quantile(S, 0.9)
mean(S[S > VaRS])

# Vérif TVaR log-logistique
lam[2]/0.1 * gamma(1 + 1/tau[2]) * gamma(1 - 1/tau[2]) * pbeta(0.9, 1 + 1/tau[2], 1 - 1/tau[2], low = F)









