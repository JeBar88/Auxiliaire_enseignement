############
# Question 2 section 1.4 chapitre  1 -------------------------------------------
############
nFFT <- 2^20
lam <- 1.5
r <- 0.5
q <- 0.25
s <- 0.8
mu <- log(100) - (s^2)/2

fb <- c(discretize(plnorm(x, mu, s), 0, 1e5, 1, "lower"), numeric(nFFT - length(seq(0, 1e5, 1))))
sum(fb)

fxt <- (q^r * exp(lam * (fft(fb) - 1)))/(1 - (1 - q) * fft(fb))^r
fx <- Re(fft(fxt, T))/nFFT
sum(fx)

round(fx[c(1, 101, 201, 301, 401, 501)],6)

Fx <- cumsum(fx)

VaR <- function(k) min(which(Fx >= k)) - 1
round(sapply(c(0.5, 0.9, 0.99, 0.999, 0.9999), VaR), 6)

############
# Question 3 section 1.4 chapitre  1 -------------------------------------------
############
nFFT <- 2^20
lam <- 3
b <- 0.5
s <- 0.8
mu <- log(100) - (s^2)/2

fb <- c(discretize(plnorm(x, mu, s), 0, 1e5, 1, "lower"), numeric(nFFT - length(seq(0, 1e5, 1))))
sum(fb)

fxt <- exp(1/b * (1 - sqrt(1 - 2*b*lam*(fft(fb) - 1))))
fx <- Re(fft(fxt, T))/nFFT
sum(fx)

round(fx[c(1, 101, 201, 301, 401, 501)], 6)

Fx <- cumsum(fx)

VaR <- function(k) min(which(Fx >= k)) - 1
round(sapply(c(0.5, 0.9, 0.99, 0.999, 0.9999), VaR), 6)















