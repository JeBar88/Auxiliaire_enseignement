### ACT-3000
## Exemple de code de base et FFT
## Jérémie Barde

### Ex 1 : X - Po(2)
lam <- 2
x <- 0:100
fx <- dpois(x, lam)
sum(fx)

## Repartition
Fx <- cumsum(fx)

## Espérance et variance
Esp <- sum(x * fx)
Var <- sum(x^2 * fx) - Esp^2
cbind(Esp, Var)

## VaR
k <- 0.9
VaR <- x[min(which(Fx >= k))]

## TVaR
sl <- sum(pmax(x - VaR, 0) * fx) # Stop-Loss
TVaR <- VaR + 1/(1 - k) * sl
cbind(VaR, sl, TVaR)

## Espérance tronquée
EspT <- sum(x * I(x > VaR) * fx)
TVaR2 <- 1/(1 - k) * (EspT + VaR * (Fx[VaR + 1] - k)) # Ne pas faire je vous en suplie !!!
cbind(TVaR, TVaR2)

### Ex 2 : produit de Convo : X - Po(10), Y - Bin(n = 10, q = 0.3)
lam <- 10
n <- 10
q <- 0.3
x <- 0:100
fx <- dpois(x, lam)  
fy <- dbinom(x, n, q)
apply(cbind(fx, fy), 2, sum) # sum(fx) et sum(fy)

s <- 0:100
fs <- sapply(s, function(s) sum(fx[1:(s+1)]*fy[(s+1):1]))
# Vérif
sum(fs)
Esp_theo <- lam + n*q
Esp <- sum(s * fs)  
cbind(Esp_theo, Esp)  

# Avec fft
nfft <- 2^6
x <- 0:(nfft - 1)
fx <- dpois(x, lam)  
fy <- dbinom(x, n, q)
apply(cbind(fx, fy), 2, sum) # sum(fx) et sum(fy)

fst <- fft(fx) * fft(fy)
fs <- Re(fft(fst, T))/nfft
sum(fs * x)

## Après même chose que Ex 1
# Repartition pour s = 10, 20, 30
Fs <- cumsum(fs)
Fs[c(11, 21, 31)]

# VaR k = 0.9
VaR <- s[min(which(Fs >= k))]

# TVaR k = 0.9
TVaR <- VaR + 1/(1 - k)*sum(pmax(s - VaR, 0) * fs) 
cbind(VaR, TVaR)






  
  


