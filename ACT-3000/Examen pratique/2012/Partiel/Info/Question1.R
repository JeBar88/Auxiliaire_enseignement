#######################################
### Examen partiel info 2012 question 1
#######################################
nFFT <- 2^20
mu <- 2.08
s <- 0.83
lam <- 0.64
x <- 0:1e5
fb <- c(discretise(plnorm(x, mu, s), 0, 1e5, 1, "lower"), numeric(nFFT - length(x)))
sum(fb)
length(fb) == nFFT

### a)
fst <- exp(lam * (fft(fb) - 1))
fs <- Re(fft(fst, T))/nFFT
sum(fs)

fs[c(1, 51)]

### b)
Fs <- cumsum(fs)
1 - Fs[41]
1 - Fs[51]

### c)
min(which(Fs > 0.95)) - 1
min(which(Fs > 0.99)) - 1
