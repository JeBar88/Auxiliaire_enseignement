### ACT-3000
## Chapitre 1 informatique : Question 6
## Jérémie Barde

### Variables
nfft <- 2^5
n <- c(10, 20)
q <- c(0.2, 0.3)
k <- 0:(nfft - 1)

fx1 <- dbinom(k, n[1], q[1])
fx2 <- dbinom(k, n[2], q[2])

fst <- fft(fx1)*fft(fx2)
fs <- Re(fft(fst, TRUE))/nfft

EX <- sum(n*q)
EX_test <- sum(fs*k)
cbind(EX, EX_test)
