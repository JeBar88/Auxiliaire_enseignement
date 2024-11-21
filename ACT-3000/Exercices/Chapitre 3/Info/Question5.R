### ACT-3000
## Chapitre 3 info : Question 5
## Jérémie Barde

### Package
library(Runuran)

### Variables
nfft <- 2^10
k <- 0:(nfft - 1)
b <- c(0.1, 0.2, 1/(1/0.1 + 1/0.2))
a <- c(2.5, 5)
gam0 <- 1
at <- sum(a) - gam0
q <- b[-which.max(b)]/max(b)

### b)
ES <- sum(a/b[1:2])
VarS <- sum(a/b[1:2]^2) + 2*gam0/prod(b[1:2])

### c)
j1 <- dnbinom(k, a[1] - gam0, q[1])
j2 <- dnbinom(k, gam0, q[2])
fkt <- fft(j1) * fft(j2)
fk <- Re(fft(fkt, T))/nfft

# Vérif
sum(fk*(at + k)/max(b))
sum(fk*(at + k)*(at + k + 1)/max(b)^2) - sum(fk*(at + k)/max(b))^2

## Répartition
Fs <- function(x) sum(fk * pgamma(x, at + k, max(b)))

## VaR
u <- 0.99
VaRS <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 400))$min
VaRS(0.99)

## TVaR
TVaRS <- function(u){
  VaR <- VaRS(u)
  1/(1 - u) * sum(fk*(at + k)/max(b)*pgamma(VaR, at + k + 1, max(b), low=FALSE))
}

TVaRS(u)
                                               
                                                                                    
                                                                                    
                                                                                    
                                                                                    
                                                                                    
                                                                                                                         0:1000 + 1, max(b), low = F))