##########################
### Section 1.2 question 3
##########################
nfft <- 2^20
a <- 1.5
lam <- 5
x <- seq(0, 5e4, 0.1)

fx <- c(discretise(ppareto(x, a, lam), 0, 5e4, 0.1, "lower"), numeric(nfft - length(x)))

fst <- fft(fx)^2
fs <- Re(fft(fst, T))/nfft
Fs <- cumsum(fs)

VaR <- function(k) (min(which(Fs > k)) - 1) * 0.1
sapply(c(0.9, 0.99, 0.999, 0.9999), VaR)
