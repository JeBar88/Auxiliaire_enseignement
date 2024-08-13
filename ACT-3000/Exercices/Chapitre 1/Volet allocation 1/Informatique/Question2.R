##############################################
### Question 1 informatique volet allocation 2 ---------------------------------
##############################################
parP <- c(2.25, 125)
parG <- c(1/9, 1/900)
parLN <- c(3.453878, 1.517427)
m <- 1e5

set.seed(2019)
U <- matrix(runif(3e5), ncol = 3, byrow = T)

X1 <- qpareto(U[, 1], parP[1], parP[2])
X2 <- qgamma(U[, 2], parG[1], parG[2])
X3 <- qlnorm(U[, 3], parLN[1], parLN[2])

X1[c(1, m)]
X2[c(1, m)]
X3[c(1, m)]

S <- X1 + X2 + X3

Fs <- ecdf(S)

## c)
VaRS <- function(k) quantile(S, k)
sapply(c(0.9, 0.99, 0.999), VaRS)

CVaRX1 <- function(k) sum(X1 * I(S == sort(S)[k*m]))
sapply(c(0.9, 0.99, 0.999), CVaRX1)

CVaRX2 <- function(k) sum(X2 * I(S == sort(S)[k*m]))
sapply(c(0.9, 0.99, 0.999), CVaRX2)

CVaRX3 <- function(k) sum(X3 * I(S == sort(S)[k*m]))
sapply(c(0.9, 0.99, 0.999), CVaRX3)

CVaRX1(0.9) + CVaRX2(0.9) + CVaRX3(0.9)


## d)
TVaRS <- function(k) mean(S[S > VaRS(k)])
sapply(c(0.9, 0.99, 0.999), TVaRS)

CTVaRX1 <- function(k) mean(X1 * I(S > sort(S)[k*m]))/(1-k)
sapply(c(0.9, 0.99, 0.999), CTVaRX1)

CTVaRX2 <- function(k) mean(X2 * I(S > sort(S)[k*m]))/(1-k)
sapply(c(0.9, 0.99, 0.999), CTVaRX2)

CTVaRX3 <- function(k) mean(X3 * I(S > sort(S)[k*m]))/(1-k)
sapply(c(0.9, 0.99, 0.999), CTVaRX3)

CTVaRX1(0.999) + CTVaRX2(0.999) + CTVaRX3(0.999)













