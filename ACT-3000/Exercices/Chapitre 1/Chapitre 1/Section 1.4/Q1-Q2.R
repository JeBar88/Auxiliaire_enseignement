##########################
### Section 1.4 question 1
##########################
nfft <- 2^10
lam <- 1.5
r <- 0.5
q <- 0.25

fb <- c(0, 1, numeric(nfft - 2))

fmt <- (q^r * exp(lam * (fft(fb) - 1)))/(1 - (1 - q) * fft(fb))^r
fm <- Re(fft(fmt, T))/nfft

fm[1:6]

##########################
### Section 1.4 question 2
##########################
nfft <- 2^15
lam <- 1.5
r <- 0.5
q <- 0.25
s <- 0.8
mu <- log(100) - s^2/2
x <- 0:(nfft - 1)

fb <- discretize(plnorm(x, mu, s), 0, nfft - 1, 1, "lower")

fxt <- (q^r * exp(lam * (fft(fb) - 1)))/(1 - (1 - q) * fft(fb))^r
fx <- Re(fft(fxt, T))/nfft
Fx <- cumsum(fx)


EM <- r*(1-q)/q + lam
EB <- sum(x*fb)
EX <- EM*EB

VarM <- r*(1-q)/q^2 + lam
VarB <- sum(x^2*fb) - EB^2
VarX <- EM*VarB + VarM*EB^2
# où
VarX <- sum(x^2*fx) - EX^2

VaR <- function(k) x[min(which(Fx > k))]
sapply(c(0.5, 0.9, 0.99, 0.999, 0.9999), VaR)














