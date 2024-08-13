######################################
### Question 7 informatique chapitre 1 -----------------------------------------
######################################

nFFT <- 2^16
r <- 0.01
q <- 0.25
a <- 2.5
lam <- 15000
k <- seq(0, 4e7, 1e3)

fb <- c(discretize(actuar::ppareto(x, a, lam), 0, 4e7, 1e3, "lower"), numeric(nFFT - length(k)))
sum(fb)

fst <- ((q/(1 - (1 - q) * fft(fb)))^r)^1000
fs <- Re(fft(fst, T))/nFFT

sum(k * fs[1:40001])
sum(fs[1:40001] * pmax(k - 2e6, 0))
