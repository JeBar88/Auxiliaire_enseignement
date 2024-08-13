### ACT-3000
## Exemple de base pour une v.a discrète bivariée
## Jérémie Barde

##### Fonction utile
dPoTei <- function(m1, m2, lam1, lam2, a){
  f <- matrix(numeric(length(m1) * length(m2)), ncol = length(m2))
  
  for (i in m1){
    for (j in m2){
      
      k <- 0:min(i, j)
      
      f[i + 1, j + 1] <-  sum(dpois(k, a) * dpois(i - k, lam1 - a) * dpois(j - k, lam2 - a))
    }
  }
  f
}
repart <- function(x, y, fxy){
  f <- matrix(numeric(length(x)*length(y)), ncol = length(y))
  for (i in x) {
    for (j in y) {
      f[i + 1, j + 1] <- sum(fxy[0:i + 1, 0:j + 1])
    }
  }
  f
}
dn <- function(k, f){
  
  fs <- function(s) sum(sapply(0:s, function(i) f[i + 1, s - i + 1]))
  sapply(k, fs)
  
}
ds <- function(k, f){
  fx <- numeric(length(k))
  
  fs <- function(s){
    tot <- 0
    for (i in 0:s){
      tot <- tot + f[i + 1, s - i + 1]
    }
    tot
  }
  
  for(i in k){
    fx[i + 1] <- fs(i)
  }
  fx
} # même chose que dn

### Loi discrete ###
fxy <- matrix(c(0.5, 0.05, 0.05,
                0.05, 0.3, 0.01,
                0.01, 0.02, 0.01), ncol = 3, byrow = T)
x <- 0:2
y <- x

# 1)
fx <- apply(fxy, 1, sum)
fy <- apply(fxy, 2, sum)

espX <- sum(fx * x)
espY <- sum(fy * y)
EspXY <- sum(outer(0:2, 0:2, "*")*fxy)

# 2)
varX <- sum(fx * x^2) - espX^2
varY <- sum(fy * y^2) - espY^2
cov <- EspXY - espX * espY
pp <- cov/sqrt(varX*varY)

# 3)
Fxy <- repart(x, y, fxy)

# 4) S = X + Y
fxy0 <- matrix(0, ncol=5, nrow=5)
fxy0[1:nrow(fxy), 1:ncol(fxy)] <- fxy

fs <- dn(0:4, fxy0)
fs


