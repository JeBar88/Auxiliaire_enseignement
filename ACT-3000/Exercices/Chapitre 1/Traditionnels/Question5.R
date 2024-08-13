######################################
### Question 5 traditionnel chapitre 1 -----------------------------------------
######################################
## a) b) c)
nFFT <- 2^15
n <- c(120, 80)
lam <- c(0.036, 0.054)
x <- seq(0, 1e3*1e4, 10000)

fb1 <- c(0, 0.4 * 0.6^(1:(nFFT-1) - 1))
fb2 <- c(0, 0.5 * 0.5^(1:(nFFT-1) - 1))

Eb1 <- sum(fb1[1:1001] * x)
Eb2 <- sum(fb2[1:1001] * x)

fst <- exp(lam[1] * (fft(fb1) - 1))^n[1] * exp(lam[2] * (fft(fb2) - 1))^n[2]
fs <- Re(fft(fst, T))/nFFT

fs[c(1, 2, 3, 4)]
fs[c(1, 11, 21, 31)]

ES <- sum(fs[1:1001] * x)

## d)
Fs <- cumsum(fs)
VaR <- function(k) x[min(which(Fs >= k))]
sapply(c(0.95, 0.99), VaR)

TVaR <- function(k) VaR(k) + sum(pmax(x - VaR(k), 0) * fs[1:1001])/(1-k)
sapply(c(0.95, 0.99), TVaR)

