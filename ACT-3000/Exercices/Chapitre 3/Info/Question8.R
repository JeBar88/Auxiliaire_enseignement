### ACT-3000
## Chapitre 3 info : Question 8
## Jérémie Barde

### Variables
m <- 1e6
k <- 0:2
a <- 3
lam <- 50
gam <- 2.5
b <- 0.1

### a)
fm12 <- matrix(c(0.30, 0.20, 0.10,
                 0.10, 0.02, 0.03,
                 0.07, 0.10, 0.08), ncol=3, byrow = TRUE)
fm1 <- apply(fm12, 1, sum)
fm2 <- apply(fm12, 2, sum)

### b)
EM <- c(sum(k*fm1), sum(k*fm2))
VarM <- c(sum(k^2*fm1), sum(k^2*fm2)) - EM^2
CovM12 <- sum(outer(k, k)*fm12) - prod(EM)

EB <- c(mpareto(1, a, lam), gam/b)
VarB <- c(mpareto(2, a, lam) - EB[1]^2, gam/b^2)

EX <- EM*EB
VarX <- EM*VarB + VarM*EB^2
CovX12 <- prod(EB)*CovM12

ES <- sum(EX)
VarS <- sum(VarX) + 2*CovX12
  
### d)
set.seed(21082024)

simul <- function(){
  M1 <- k[min(which(cumsum(fm1) > runif(1)))]
  fm2c1 <- fm12[M1 + 1, ]/fm1[M1 + 1]
  M2 <-  k[min(which(cumsum(fm2c1) > runif(1)))]
  c(M1, M2)
}

M <- t(replicate(m, simul()))
M[c(1, 2, m), ]

X1 <- sapply(M[, 1], function(i) sum(rpareto(i, a, lam)))
X2 <- sapply(M[, 2], function(i) sum(rgamma(i, gam, b)))
X <- cbind(X1, X2)
X[c(1, 2, m), ]


# # ou
# perm <- expand.grid(0:2, 0:2)
# fm <- as.vector(fm12)
# simul <- sample(1:9, m, replace = TRUE, prob=fm)
# M1 <- perm[simul, 1]
# M2 <- perm[simul, 2]
# M <- cbind(M1, M2)









  