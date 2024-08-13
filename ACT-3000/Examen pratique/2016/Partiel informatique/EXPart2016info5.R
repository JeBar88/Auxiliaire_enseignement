## 5 2016 info partiel
nfft <- 2^16
n <- 10
q <- c(0.2, 0.3)
r <- 2
gam <- 0.8
lam <- c(1.5, 2)

f_I <- c(0, 1, numeric(nfft - 2))

fgp_theta <- function(t1, t2){
  ((1 - gam)/((1 - (1 - gam)/r*t1)*(1 - (1 - gam)/r*t2)-gam))^r
}

fgp_binom <- function(t, n, q){
  ((1 - q) + q*t)^n
}

fst <- fgp_theta(lam[1]*(fgp_binom(fft(f_I), n, q[1]) - 1), lam[2]*(fgp_binom(fft(f_I), n, q[2]) - 1))
fs <- Re(fft(fst, inverse = TRUE))/nfft
sum(fs)

Fs <- cumsum(fs)
Fs[c(31, 41)]




