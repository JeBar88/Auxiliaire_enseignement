### ACT-3000
## Chapitre 1 informatique : Question 7
## Jérémie Barde

### Variables
nfft <- 2^22
a <- 2.5
lam <- 15e3
r <- 0.01
q <- 0.25
h <- 1000
k <- (0:(nfft - 1))*h

fb <- discretize(ppareto(x, a, lam), 0, 4e6, step = h, method = "lower") 
fb <- c(fb, numeric(nfft - length(fb)))

fxt <- (q/(1 - (1 -  q)*fft(fb)))^r
fst <- fxt^1000
fs <- Re(fft(fst, TRUE))/nfft
Fs <- cumsum(fs)

ES <- 1000*r*(1-q)/q * lam/(a - 1)
ES_test <- sum(fs * k)
cbind(ES, ES_test)

SL <- sum(pmax(k - 2e6, 0)*fs)

Fs[c(5e5, 1e6, 2e6)/h + 1]




