####################################
### Examen partiel informatique 2022
####################################

### Question 3
nfft <- 2^20
mu <- log(3e3) - 0.32
s <- 0.8
r <- 2.5
q <- 1/3
x <- seq(0, 1e4, 1e3)

fx <- c(discretize(plnorm(x, mu, s), 0, 9e3, 1e3, "lower"), plnorm(9e3, mu, s, low = F), 
        numeric(nfft - length(x)))

fst <- (q/(1 - (1 - q) * fft(fx)))^r
fs <- Re(fft(fst, T))/nfft
sum(fs)
fs[16]

### Question 7
nfft <- 2^15
r <- c(2, 2.5)
p <- c(1/3, 1/5)
b <- c(0.1, 0.25)
q <- b[1]/b[2]

k <- 0:(nfft - 1)
fj <- c(0, q * (1 - q)^(1:(nfft - 1) - 1))
fi <- c(0, 1, numeric(nfft - 2))

vkt <- (p[1]/(1 - (1 - p[1]) * fft(fj)))^r[1] * (p[2]/(1 - (1 - p[2]) * fft(fi)))^r[2]
vk <- Re(fft(vkt, T))/nfft

Fs <- function(x) vk[1] + sum(vk[2:1002] * pgamma(x, 1:1001, b[2]))
Fs(300)









