############
# Question 10 informatique chapitre  1 -----------------------------------------
############
nFFT <- 2^10
b <- c(0.1, 0.5)
q <- b[1]/b[2]
a <- 5/8

gK <- c(0, a * q + (1 - a), a * q * (1 - q)^(2:(nFFT - 1) - 1))
sum(gK)
gK[1:5]

nKt <- fft(gK)^20
nK <- Re(fft(nKt, TRUE))/nFFT
sum(nK)
nK[c(51, 61, 71)]

Fs <- function(x) nK[1] + sum(nK[-1] * pgamma(x, 1:(nFFT - 1), b[2]))
sapply(c(100, 140, 200, 300), Fs)

nfft <- 2^10
x <- 0:(nfft - 1)
lam <- 1:2

fxi <- sapply(lam, function(i) dpois(x, i))
fxit <- apply(fxi, 2, fft)
fst <- apply(fxit, 1, prod)
fs <- Re(fft(fst, inverse = TRUE))/nfft


sum(fs * x)

à