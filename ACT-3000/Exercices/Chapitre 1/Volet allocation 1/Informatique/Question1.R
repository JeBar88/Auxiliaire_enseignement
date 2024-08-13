##############################################
### Question 1 informatique volet allocation 1 ---------------------------------
##############################################
### Question 1
m <- 1e5
a <- c(10, 4, 1)
b <- c(1/300, 1/500, 1/1000)

set.seed(2019)
U <- matrix(runif(3e5), ncol = 3, byrow = T)

X1 <- qgamma(U[, 1], a[1], b[1])
X2 <- qgamma(U[, 2], a[2], b[2])
X3 <- qgamma(U[, 3], a[3], b[3])

X1[c(1, m)]
X2[c(1, m)]
X3[c(1, m)]

S <- X1 + X2 + X3

## c)
Fs <- ecdf(S)

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







