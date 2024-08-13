################################
### Chapitre 3 exercices trad 25
################################
nfft <- 2^15
lam <- c(0.003, 0.004)
g <- 0.001
lamS <- 500*(lam[1] - g) + 500*(lam[2] - g) + g
a <- c(2, 1)
b <- 1/1000
k <- 0:(nfft - 1)

### a)
fk <- c(0, 500*(lam[2] - g)/lamS, 500*(lam[1] - g)/lamS, numeric(1497) , g/lamS, numeric(nfft - 1501))
sum(fk)
length(fk)

fjt <- exp(lamS * (fft(fk) - 1))
fj <- Re(fft(fjt, T))/nfft

Fs <- function(x) fj[1] + sum(fj[-1] * pgamma(x, k[-1], b))
Fs(4000)

### b) Pour X1
fm <- dpois(0:10, lam[1])
fm[1]
# Pr(M = 0) = 0.9970045 donc VaR à 0.995 est 0

### c)
VaR995 <- optimize(function(x) abs(Fs(x) - 0.995), c(0, 4e4))$min
TVaR995 <- sum(fj[1:10001] * 0:10000/b * pgamma(VaR995, 0:10000 + 1, b, low = F))/(1-0.995)
cbind(VaR995, TVaR995)

