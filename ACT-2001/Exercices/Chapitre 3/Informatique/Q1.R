### ACT-2001
## Chapitre 3 informatique : Question 1
## Jérémie Barde

b <- c(0.04, 0.1)
simul <- function(n){
  U <- numeric(n)
  x0 <- 20150309
  a <- 41358
  m <- 2147483647
  U[1] <- ((a * x0) %% m)/m
  for (i in 1:(n-1)) {
      U[i+1] <- ((a * U[i]*m) %% m)/m
    }
  matrix(U, byrow = T, ncol = 2)
}
U <- simul(2000)

## a)
X1 <- qexp(U[, 1], b[2])
X2 <- qexp(U[, 2], b[1])
head(X1, 2)
head(X2, 2)

## b)
S <- X1 + X2
head(S, 2)
max(S)
min(S)

## c)
vk <- c(10, 50, 200)
# Fonction de répartition
Fm <- ecdf(S)
sapply(vk, Fm)
sapply(vk, function(k) mean(S <= k)) # Autre façon

# Stop-loss
sapply(vk, function(d) mean(pmax(S - d, 0)))

## d) Démarche à venir :)
# Répartition : façon 1
Fs <- function(k) pexp(k, b[1]) - 2/3 * (pexp(k, b[1], low = F) - pexp(k, b[2], low = F))
sapply(vk, Fs)
# Répartition : façon 2
k <- 0:1000
sum(dgeom(k, b[1]/b[2]) * pgamma(10, k+2, b[2]))
# stop-loss
SL <- function(k) 5/3 * mexp(1, b[1]) * pexp(k, b[1], low = F) - 2/3 * mexp(1, b[2]) * pexp(k, b[2], low = F)
sapply(vk, SL)

table2 <- cbind(sapply(vk, SL),sapply(vk, function(d) mean(pmax(S - d, 0))))
colnames(table2) <- c("Théorique", "Empirique")
table2











