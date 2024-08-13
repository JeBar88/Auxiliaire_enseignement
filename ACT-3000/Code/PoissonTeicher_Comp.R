### ACT-3000
## Loi poisson Teicher composée Expo/Gamma, bêta diff. 
## Jérémie Barde

### Variables utiles
nfft <- 2^10
k <- 0:(nfft - 1)
a <- c(2, 3)
b <- c(0.2, 0.3)
lam <- c(3, 2)
g0 <- 1
q <- b[1]/b[2]

### Variable Ji
fJ1 <- c(numeric(a[1]), dnbinom(head(k, - a[1]), a[1], q))
fJ2 <- c(numeric(a[2]), 1, numeric(nfft - 1 - a[2]))
fJ0t <- fft(fJ1) * fft(fJ2)

### Variable Wi
fW1t <- exp((lam[1] - g0)*(fft(fJ1) - 1))
fW2t <- exp((lam[2] - g0)*(fft(fJ2) - 1))
fW0t <- exp(g0*(fJ0t - 1))

### Variable L
fLt <- fW1t*fW2t*fW0t
fL <- Re(fft(fLt, T))/nfft

### Fonction de répartition de S
FS <- function(x) fL[1] + sum(fL[-1]*pgamma(x, k[-1], b[2]))









