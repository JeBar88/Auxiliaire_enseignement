### ACT-3000
## Chapitre 1 informatique : Question 2
## Jérémie Barde

### Variable
nfft <- 2^13
lam <- 0.6 - 0.1 * 1:5
lamN <- sum(lam)
pi <- lam/lamN
a <- c(2.9, 2.7, 2.5, 2.3, 2.1)
eta <- c(1900, 1700, 1500, 1300, 1100)
x <- (0:(nfft - 1)) * 100

fbi <- sapply(1:5, function(i) c(discretize(ppareto(x, a[i], eta[i]), 0, 2e4, 100, method = "lower"), 
                                 numeric(nfft - length(seq(0, 2e4, 100)))))

fc <- colSums(pi * t(fbi))
fc[c(21, 81)]
sum(fc * x)

fst <- exp(lamN * 1.6*(fft(fc) - 1))
fs <- Re(fft(fst, T))/nfft
fs[c(1, 21, 81)]

Fs <- cumsum(fs)
Fs[c(1, 21, 81, 101)]

SL <- function(d) sum(pmax(x - d, 0) * fs)
sapply(c(0, 20, 80, 200)*100, SL)

k <- Fs[c(1, 21, 81, 101)]
VaR <- sapply(k, function(i) x[min(which(Fs >= i))])

TVaR <- VaR + 1/(1-k) * sapply(VaR, SL)


