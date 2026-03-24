### ACT-3000
## Exemple de code de base et FFT
## # Rappel univariée
## # Somme de v.a.s iid
## # Somme de v.a.s indépendantes
## # Fmp extraite de la fgp
## Jérémie Barde

#### Rappel univariée -- X - Po(2) ####
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

### Mesure entropique
rho <- 0.05
1/rho*log(sum(exp(k*rho)*fx))

### Franchise et limite
# Y = min(max(X - 2, 0), 8)
k <- pmin(pmax(k - 2, 0), 8)
km <- unique(k)
fy <- tapply(fx, k, sum)
EY <- sum(km*fy)
EY2 <- sum(k*fx)
cbind(EY, EY2)







#### Somme de v.a.s iid -- Xi~X~Po(5), i = {1, 10}, S = X1 + ... + X10 ####
nfft <- 2^10
n <- 10
lam <- 5
k <- 0:(nfft - 1)

fx <- dpois(x, lam)
sum(fx)

fst <- fft(fx)^n
fs <- Re(fft(fst, TRUE))/nfft
cbind("Es_test"=sum(k*fs), "Es"=n*lam)

#### Somme de v.a.s indépendantes ####
#### Exemple 1 : M1~Po(5), M2~Bin(10, 0.7), M3~geo(0.8)
nfft <- 2^10
n <- 10
lam <- 5
q <- c(0.7, 0.8)
k <- 0:(nfft - 1)

### Espéraces
Emi <- c(lam, n*q[1], (1 - q[2])/q[2])
Es <- sum(Emi)

### Fmp des v.a.s Mi
fm1 <- dpois(k, lam)
fm2 <- dbinom(k, n, q[1])
fm3 <- dgeom(k, q[2])

### V.a. S = M1 + ... + M10
fst <- fft(fm1)*fft(fm2)*fft(fm3)
fs <- Re(fft(fst, TRUE))/nfft
sum(fs)
cbind("Es_test"=sum(k*fs), Es)

#### Exemple 2 : Mi~Bern(qi), i=1,2,...,50
nfft <- 2^6
n <- 50
i <- 1:50
q <- (i + n/3)/(i + n)
k <- 0:(nfft - 1)

### Espérance
Emi <- q
Es <- sum(Emi)

### Fmp des v.a.s Mi
fmi <- sapply(q, function(i) dbinom(k, 1, i))
# Vérif
colSums(fmi)
apply(fmi, 2, function(f) sum(k*f))

### V.a S
fmit <- mvfft(fmi)
fst <- apply(fmit, 1, prod)
fs <- Re(fft(fst, TRUE))/nfft
cbind("Es_test"=sum(k*fs), Es)

#### Fmp extraite de la fgp -- M~BinNeg(2.5, 0.4) ####
nfft <- 2^10
r <- 2.5
q <- 0.4
k <- 0:(nfft - 1)

fi <- replace(numeric(nfft), 2, 1)

fmt <- (q/(1 - (1 - q)*fft(fi)))^r
fm <- Re(fft(fmt, TRUE))/nfft
cbind("Em_test"=sum(k*fm), "Em"=r*(1 - q)/q)









