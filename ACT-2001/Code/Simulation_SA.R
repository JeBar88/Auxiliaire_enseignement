###################################
### ACT-2001
### Simulation de somme aléatoire
### Jérémie Barde
###################################

### Param
m <- 1e5
a <- 1.5
b <- 0.1
lam <- 2

### Façon 1a
M <- numeric(m)
X <- numeric(m)
set.seed(2024)
for (i in 1:m) {
  M[i] <- qpois(runif(1), lam)
  
  if(M[i] > 0){
    X[i] <- sum(qgamma(runif(M[i]), a, b))
  }
  
}

### Façon 1b
M <- numeric(m)
Y <- numeric(m)
set.seed(2024)
for (i in 1:m) {
  M[i] <- rpois(1, lam)
  
  if(M[i] > 0){
  Y[i] <- sum(rgamma(M[i], a, b))
  }
  
}

#### Façon 2a
set.seed(2024)
M <- rpois(m, lam)
Z <- numeric(m)
for (i in 1:m) {
    Z[i] <- sum(rgamma(M[i], a, b))
}

#### Façon 2b
set.seed(2024)
M <- rpois(m, lam)
A <- sapply(1:m, function(i) sum(rgamma(M[i], a, b)))
mean(A)

#### Result
data <- cbind(X, Y, Z, A)
apply(data, 2, mean)
apply(data, 2, quantile, 0.99)








