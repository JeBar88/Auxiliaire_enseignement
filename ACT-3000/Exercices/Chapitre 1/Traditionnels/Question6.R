######################################
### Question 6 traditionnel chapitre 1 -----------------------------------------
######################################
library(actuar)
nFFT <- 2^10
t <- 0.8
b <- 1/1000
r <- 0.4
q <- 2/3
k <- seq(1e3, 4e3, 1e3)
x <- 0:(nFFT - 1) * 1e3

fb <- c(discretize(pweibull(x, t, 1/b), 0, 4e3, 1e3, "lower"),
        pweibull(4000, t, 1/b, low = F),
        numeric(nFFT - 2 - length(k)))
sum(fb)

head(fb)

fxt <- (q/(1 - (1-q) * fft(fb)))^r
fx <- Re(fft(fxt, T))/nFFT

round(fx[1:8], 6)

sum(fx * x)

Fx <- cumsum(fx)
#x <- seq(0, 5e3, 1e3)
VaR <- x[min(which(Fx >= 0.99))]




