### ACT-3000
## Aggrégattion avec FFT
## Jérémie Barde

##### Exemple 1 #####
### Xi - X - Po(5), i = {1, 10}, S = X1 + ... + X10
nfft <- 2^10
lam <- 5
x <- 0:(nfft - 1)

fx <- dpois(x, lam)
sum(fx)

fst <- fft(fx)^10
fs <- Re(fft(fst, TRUE))/nfft
sum(fs)
sum(fs*x)

##### Exemple 2 #####
### Si l'on connais seulement la fgp
nfft <- 2^10
lam <- 5
x <- 0:(nfft - 1)

fi <- c(0, 1, numeric(nfft - 2))

fxt <- exp(lam*(fft(fi) - 1))
fx <- Re(fft(fxt, TRUE))/nfft
sum(fx)
cbind(dpois(2, lam), fx[3])

##### Exemple 3 #####
### X1 - Po(5), X2 - Bin(10, 0.7), X3 - geo(0.8), S = X1 + X2 + X3
nfft <- 2^10
lam <- 5
n <- 10
q <- 0.7
p <- 0.8
x <- 0:(nfft - 1)

fx1 <- dpois(x, lam)
fx2 <- dbinom(x, n, q)
fx3 <- dgeom(x, p)
fx <- cbind(fx1, fx2, fx3)
apply(fx, 2, sum)

fxt <- mvfft(fx)
fst <- apply(fxt, 1, prod)
fs <- Re(fft(fst, TRUE))/nfft
sum(fs)
cbind(sum(fs*x), (1-p)/p + lam + n*q)

##### Exemple 4 #####

##### Exemple 1 #####
### Xi - Po(i/25), i = {1, 100}, S = X1 + ... + X100
nfft <- 2^10
i <- 1:100
lam <- i/25
x <- 0:(nfft - 1)

fx <- matrix(0, nfft, length(i))
for (i in 1:100) {
  fx[, i] <- dpois(x, lam[i])
}
# ou
fx <- sapply(lam, function(i) dpois(x, i))
apply(fx, 2, sum)

fxt <- mvfft(fx)
fst <- apply(fxt, 1, prod)
fs <- Re(fft(fst, TRUE))/nfft
sum(fs)
cbind(sum(fs*x), sum(lam))
cbind(fs[250 + 1], dpois(250, sum(lam)))





