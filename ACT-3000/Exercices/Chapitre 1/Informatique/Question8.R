### ACT-3000
## Chapitre 1 informartique : Question 8
## Jérémie Barde

### Fonction 
Panjer_Poisson <- function(lam, fb, smax){
  fs <- exp(lam * (fb[1] - 1))
  fb <- c(fb, numeric(smax))
  for (i in 1:smax) {
    j <- i+1
    fs <- c(fs, sum(lam * 1:i/i * fb[2:j] * fs[i:1]))
  }
  fs
}


### Variables
a <- 1:5
b <- 0.001

### (b) 
lam <- 0.6 - 0.1 * 1:5
lamN <- sum(lam)
pi <- lam/lamN

### (c)
Fc <- sapply(0:50000, function(x) sum(pi * pgamma(x, a, b)))
Fc[c(2000, 8000) + 1]

### (e) 
fk <- Panjer_Poisson(lamN, c(0, pi), 1000)
sum(fk)
fk[1:4]
# (e) --------------------------------------------------------------------------

Fs <- function(x){fk[1] + sum(fk[-1] * pgamma(x, 1:1000, b))}
sapply(c(0, 2000, 8000, 20000), Fs)



