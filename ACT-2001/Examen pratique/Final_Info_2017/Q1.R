### ACT-2001 ###
### Examen final info 2017 : Question 1
m <- 1E5
set.seed(20160419)
U <- matrix(runif(3*m), ncol = 3, byrow = T)

X1 <- 1000 * (-log(U[,1]))^(-1/2)
X2 <- 1000 * (1/U[,2] -1)^(-1/2)
X3 <- -2000 * log(2 * (1 - 1/2^(1-U[,3])))

S <- X1 + X2 + X3

head(cbind(X1, X2, X3))

## Borne
mean(S > 20000) + c(-1, 1) * 1/sqrt(m) * sd(S > 20000) * qnorm(0.975)

## VaR
VaR <- quantile(S, 0.9999)

## TVaR
mean(S[S > VaR])
# où
VaR + 1/(1-0.9999) * mean(pmax(S - VaR,0))

# f)
mean(X2 * I(X1 + X3 <= 1500 | X1 + X3 > 150000))
