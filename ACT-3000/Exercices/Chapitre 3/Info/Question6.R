### ACT-3000
## Chapitre 3 info : Question 6
## Jérémie Barde

### Variables
nfft <- 2^12
k <- (0:(nfft - 1))*500
lam <- c(2, 1)
a0 <- 0.5
n <- 4
q <- 0.6*1:2 - 0.4

fb1 <- c(0, dbinom(head(0:5, -1), n, q[1]))
fb1 <- c(fb1, numeric(nfft - length(fb1)))

sum(fb1)

fb2 <- c(0, dbinom(head(0:5, -1), n, q[2]))
fb2 <- c(fb2, numeric(nfft - length(fb2)))

fb <- cbind(fb1, fb2)

### a)
EB <- apply(fb, 2, function(i) sum(k*i))
EX <- lam*EB
ES <- sum(EX)

### b)
prod(EB)*a0

### c)
VarB <- apply(fb, 2, function(i) sum(k^2*i)) - EB^2
VarS <- sum(lam*VarB + lam*EB^2) + 2*prod(EB)*a0

### d)
fgpPoTei <- function(t1, t2){
  exp((lam[1] - a0)*(t1 - 1))*exp((lam[2] - a0)*(t2 - 1))*exp(a0*(t1*t2 - 1))
}

fst <- fgpPoTei(fft(fb[, 1]), fft(fb[, 2]))
fs <- Re(fft(fst, TRUE))/nfft

fs[0:5 + 1]

### e)
fx1t <- exp(lam[1]*(fft(fb[, 1]) - 1)) 
fx2t <- exp(lam[2]*(fft(fb[, 2]) - 1)) 

fx1 <- Re(fft(fx1t, TRUE))/nfft
fx2 <- Re(fft(fx2t, TRUE))/nfft

u <- 0.99
VaR <- function(u, f) k[min(which(cumsum(f) > u))]
VarX1 <- VaR(u, fx1)
VarX2 <- VaR(u, fx2)
VarS <- VaR(u, fs)
cbind(VarX1, VarX2, VarS)

TVaR <- function(u, f) VaR(u, f) + sum(pmax(k - VaR(u, f), 0)*f)/(1 - u)
TVarX1 <- TVaR(u, fx1)
TVarX2 <- TVaR(u, fx2)
TVarS <- TVaR(u, fs)
cbind(TVarX1, TVarX2, TVarS)

















