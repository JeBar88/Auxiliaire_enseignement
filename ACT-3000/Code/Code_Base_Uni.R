### ACT-3000
## Exemple de code de base et FFT
## Jérémie Barde

#### Exemple 1 -- X - Po(2) ####
lam <- 10
k <- 0:1000
u <- 0.9

### Fontion de masse de probabilité
fx <- dpois(k, lam)
sum(fx)

### Repartition
Fx <- cumsum(fx)

### Espérance et variance
Esp <- sum(k*fx)
Var <- sum(k^2*fx) - Esp^2
cbind(Esp, Var)

### VaR
VaR <- k[min(which(Fx >= u))]

### TVaR
sl <- sum(pmax(k - VaR, 0) * fx)
TVaR <- VaR + 1/(1 - u) * sl # Formule avec Stop-Loss
cbind(VaR, sl, TVaR)

### Espérance tronquée
EspT <- sum(k * I(k > VaR) * fx)
TVaR2 <- 1/(1 - u) * (EspT + VaR * (Fx[VaR + 1] - u)) # Moins recommandé !
cbind(TVaR, TVaR2)

### Franchise et limite
# Y = min(max(X - 2, 0), 8)
k <- pmin(pmax(k - 2, 0), 8)
km <- unique(k)
fy <- tapply(fx, k, sum)
EY <- sum(km*fy)
EY2 <- sum(k*fx)
cbind(EY, EY2)

#### Exemple 2 -- produit de Convo : X - Po(10), Y - Bin(10, 0.3) ####
lam <- 10
n <- 10
q <- 0.3
k <- 0:100 # Choisir pour S en même temps
fm1 <- dpois(k, lam)  
fm2 <- dbinom(k, n, q)
sum(fm1); sum(fm2)

### Méthode brute
fs <- sapply(k, function(k) sum(fm1[0:k + 1]*fm2[k:0 + 1]))

# Vérif
sum(fs)
Es_test <- lam + n*q
Es <- sum(k * fs)  
cbind(Es_test, Es)  

### Avec fft
nfft <- 2^6
k <- 0:(nfft - 1)
fm1 <- dpois(k, lam)  
fm2 <- dbinom(k, n, q)
sum(fm1) ; sum(fm2)

# Vérif
fst <- fft(fm1) * fft(fm2)
fs <- Re(fft(fst, T))/nfft
sum(k*fs)

### Après même chose que Exemple 1
# Repartition pour s = 10, 20, 30
Fs <- cumsum(fs)
Fs[c(10, 20, 30) + 1]

# VaR k = 0.9
u <- 0.9
VaR <- k[min(which(Fs >= u))]

# TVaR k = 0.9
TVaR <- VaR + 1/(1 - u)*sum(pmax(k - VaR, 0) * fs) 
cbind(VaR, TVaR)




