##################################################
### ACT-2001 
### Calcule de contribution selon la régle Euleur
### Jérémie Barde
##################################################
m <- 1e6
set.seed(20221128)

## b)
X1 <- rpareto(m, 2.5, 150)
X2 <- rpareto(m, 3, 200)

c(X1[2], X2[2])

## c)
S <- X1 + X2
Fn <- ecdf(S)
1 - Fn(1000) 

u <- 0.95
VaRS <- sort(S)[u*m]
VaRS

CVaRX1 <- sum(X1 * I(S == VaRS))
CVaRX2 <- sum(X2 * I(S == VaRS))

CVaRX1 + CVaRX2

## d)
TVaRS <- function(k) mean(S[S > sort(S)[k * m]])
TVaRS(0.99)

## e)
CTVaRX1 <- function(k) mean(X1 * I(S > sort(S)[k*m]))/(1-k)
CTVaRX2 <- function(k) mean(X2 * I(S > sort(S)[k*m]))/(1-k)

CTVaRX2(0.99) + CTVaRX1(0.99)





