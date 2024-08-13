### ACT-3000
## Loi Poisson Teicher bivariée
## Jérémie Barde

##### Fonction utile #####
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
      I <- matrix(numeric(length(x)*length(y)), ncol = length(y))
      I[0:i + 1, 0:j + 1] <- 1
      f[i + 1, j + 1] <- sum(I * fxy)
    }
  }
  f
}
ds <- function(k, f){
  
  fs <- function(s) sum(sapply(0:s, function(i) f[i + 1, s - i + 1]))
  sapply(k, fs)
  
}

##### Poisson Teicher #####
lam <- c(2, 3)
a <- 1
m1 <- 0:50
m2 <- m1
fm12 <- dPoTei(m1, m2, lam[1], lam[2], a)
Fm12 <- repart(m1, m2, fm12)

### Covariance et pearson
cov <- a
pp <- a/sqrt(prod(lam))

### Espérance conditionnelle
lam[1] - a + a/lam[2]*0:10
lam[2] - a + a/lam[1]*0:10

### Espérance de M1M2
Espm12 <- sum(outer(m1, m2, '*') * fm12)
cov <- Espm12 - prod(lam)
## E[min(M1, 7) * M2 * 1[M2 < = 6]]
sum(outer(pmin(m1, 7), m2 * I(m2 <= 6), '*')*fm12)

##### Poisson Teicher : loi de la somme #####
## Façon 1 : N = K1 + K2 + 2K0
nfft <- 2^6
n <- 0:(nfft - 1)
fk12 <- dpois(n, sum(lam - a))
f2k0 <- dpois(n/2, a)
fnt <- fft(fk12) * fft(f2k0)
fn <- Re(fft(fnt, T))/nfft
fn[c(1, 6, 11)]

## Façon 2 : utiliser directement fgp
fI <- c(0, 1, numeric(nfft - 2))
fnt <- exp((lam[1] - a)*(fft(fI) - 1)) * exp((lam[2] - a)*(fft(fI) - 1)) * exp(a*(fft(fI)^2 - 1))
fn2 <- Re(fft(fnt, T))/nfft
fn2[c(1, 6, 11)]

## Façon 3 : Connaitre la loi de N
lamN <- sum(lam) - a
fb <- c(0, (sum(lam) - 2*a)/lamN, a/lamN, numeric(nfft - 3))
fnt <- exp(lamN*(fft(fb) - 1))
fn3 <- Re(fft(fnt, T))/nfft
fn3[c(1, 6, 11)]

## Façon 4 : 'a la main'
fn4 <- ds(m1, fm12)
fn4[c(1, 6, 11)]

Fn <- cumsum(fn)
k <- 0.99
VaRN <- n[min(which(Fn > k))]
TVaRN <- VaRN + 1/(1-k)*sum(pmax(n - VaRN, 0)*fn)
cbind(VaRN, TVaRN)

