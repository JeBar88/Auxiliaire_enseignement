### 5 2016 info partiel
nfft <- 2^16
k <- 0:(nfft - 1)
n <- 10
q <- c(0.2, 0.3)
r <- 2
gam <- 0.8
lam <- c(1.5, 2)

## Espérances
Exi <- lam*r/r*n*q
Es <- sum(Exi)

### d)
## ii)
fb <- sapply(q, function(i) dbinom(k, n, i))
fbt <- mvfft(fb)

fgp_theta <- function(t1, t2){
  ((1 - gam)/((1 - (1 - gam)/r*t1)*(1 - (1 - gam)/r*t2)-gam))^r
}

fst <- fgp_theta(lam[1]*(fbt[, 1] - 1), lam[2]*(fbt[, 2] - 1))
fs <- Re(fft(fst, TRUE))/nfft

Fs <- cumsum(fs)
Fs[c(30, 35, 40) + 1]




