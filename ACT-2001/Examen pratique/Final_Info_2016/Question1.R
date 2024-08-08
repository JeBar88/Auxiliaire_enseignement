### ACT-2001 ###
### Examen final info 2016 : Question 1
library(actuar)
m <- 1e6
set.seed(20160419)
U <- matrix(runif(3*m), ncol = 3, byrow = T)
head(U)

X1 <- qpareto(U[, 1], 4.784422, 623.945768)
X2 <- qgamma(U[, 2], 0.581976707, 0.003529867)
X3 <- qlnorm(U[, 3], 4.60517, 1)

S <- X1 + X2 + X3

## Approximation 
mean(S > 1000) + c(-1, 1) * 1/sqrt(m) * sd(S > 1500) * qnorm(0.975)
