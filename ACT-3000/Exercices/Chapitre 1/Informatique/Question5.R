### ACT-3000
## Chapitre 1 informatique : Question 5
## Jérémie Barde

### Fonction
DirecteConvo <- function(f1, f2){
  smax <- length(f1)
  fs <- f1[1] * f2[1]
  f1 <- c(f1, numeric(smax))
  f2 <- c(f2, numeric(smax))
  for(i in 1:smax){
    j <- i+1
    fs <- c(fs, sum(f1[1:j] * f2[j:1]))
  }
  fs
}

### Variable
mu <- 8
s <- 1
a <- 4
lam <- 0.001
h <- 1000
k <- (0:500)*h

fx1 <- discretize(pgamma(x, a, lam), 0, 500000, step = 1000, method = "lower")
fx2 <- discretize(plnorm(x, mu, s), 0, 500000, step = 1000, method = "lower")
fs <- DirecteConvo(fx1, fx2)
fs[c(10000, 20000, 30000)/h + 1]
# ou
fs[which(k %in% c(10000, 20000, 30000))]


Fs <- cumsum(fs)

k[min(which(Fs >= 0.95))]
k[min(which(Fs >= 0.995))]
