### ACT-3000
## Simualtion de variable discrètre multivarié
## Jérémie Barde

### Variables
m <- 1e6
a <- 3
lam <- 50
g <- 2.5
b <- 0.1
u <- 0.995

### Matrice de massse de probabilité
fm12 <- matrix(c(0.30, 0.20, 0.10,
                 0.10, 0.02, 0.03,
                 0.07, 0.10, 0.08), 3, 3)

m1 <- 0:2
m2 <- m1

### (a) marginal
fm1 <- apply(fm12, 1, sum)
fm2 <-  apply(fm12, 2, sum)

Fm1 <- cumsum(fm1)
Fm2 <- cumsum(fm2)

### (b) Esperance et variance
EM1 <- sum(m1*fm1)
EM2 <- sum(m2*fm2)
EB1 <- lam/(a - 1)
EB2 <- g/b
EX1 <- EM1*EB2
EX2 <- EM2*EB2
ES <- EX1 + EX2

VarM1 <- sum(m1^2*fm1) - EM1^2
VarM2 <- sum(m2^2*fm2) - EM2^2
VarB1 <- mpareto(2, a, lam) - EB1^2
VarB2 <- g/b^2
VarX1 <- EM1*VarB1 + VarM1*EB1^2
VarX2 <- EM2*VarB2 + VarM2*EB2^2 

CovM12 <- sum(outer(m1, m2,'*')*fm12) - EM1*EM2
CovX12 <- EB1*EB1*CovM12
VarS <- VarX1 + VarX2 + 2*CovX12

cbind(ES, VarS)

### (c) 
set.seed(21082024)
simul_freq <- function(){
  M1 <- m1[min(which(Fm1 >= runif(1)))]
  fm2c1 <- fm12[M1 + 1, ]/fm1[M1 + 1]
  Fm2c1 <- cumsum(fm2c1)
  M2 <- m2[min(which(Fm2c1 >= runif(1)))]
  c(M1, M2)
}

M12 <- t(replicate(m, simul_freq()))
M12[c(1, 2, m), ]

X1 <- sapply(M12[, 1], function(i) sum(rpareto(i, a, lam)))
X2 <- sapply(M12[, 2], function(i) sum(rgamma(i, g, b)))
X12 <- cbind(X1, X2)

X12[c(1, 2, m), ]

### (d) 
S <- X1 + X2
S[c(1, 2, m)]

### (e)
mean(S)
var(S)

### (f)
VaRS <- quantile(S, u)
TVaRS <- mean(S[S > VaRS])
cbind(VaRS, TVaRS)




