### ACT-3000
## Exemple de code de base et FFT
## # Discrétisation avec h=1, Upper et Lower
## # Discrétisation avec h=0.1, Lower
## # V.a.s composées et discrétisation
## # Aggrégation et discretisation
## Jérémie Barde

#### Packages utiles ####
library(actuar)

#### Discrétisation avec h=1 -- X~Ga(2.5, 0.5) ####
k <- 0:1000
a <- 2.5
b <- 0.5
u <- 0.9

### Espérance, VaR et TVaR théorique
ExTh <- a/b 
VaRxTh <- qgamma(u, a, b)
TVaRxTh <- 1/(1 - u)*a/b*pgamma(VaRxTh, a + 1, b, low=FALSE)

#### Méthode Lower
### Fmp et cdf par discrétisation
fxl <- pgamma(k, a, b) - pgamma(k-1, a, b)
Fxl <- cumsum(fxl)
sum(fxl)

### Espérance approximative
Ex_lower <- sum(k*fxl)

### VaR approximative
VaRx_lower <- k[min(which(Fxl >= u))]

### TVaR approximative
TVaRx_lower <- VaRx_lower + sum(pmax(k - VaRx_lower, 0)*fxl)/(1 - u)

#### Méthode Upper
### Fmp et cdf par discrétisation
fxu <- pgamma(k + 1, a, b) - pgamma(k, a, b)
Fxu <- cumsum(fxu)
sum(fxu)

### Espérance approximative
Ex_upper <- sum(k*fxu)

### VaR approximative
VaRx_upper <- k[min(which(Fxu >= u))]

### TVaR approximative
TVaRx_upper <- VaRx_upper + sum(pmax(k - VaRx_upper, 0)*fxu)/(1 - u)

### Vérif
cbind(Ex_upper, ExTh, ExTh, Ex_lower)
cbind(VaRx_upper, VaRxTh, VaRx_lower)
cbind(TVaRx_upper, TVaRxTh, TVaRx_lower)

### Graphique
curve(pgamma(x, a, b), xlim = c(0, 15), lwd=2)
matplot(k, Fxl, add = TRUE, type = 's', col = "blue", lty=2)
matplot(k, Fxu, add = TRUE, type = 's', col = "darkred", lty=2)

#### Discrétisation avec h=1 -- X~Ga(2.5, 0.5) ####
h <- 0.1
k <- (0:1000)*h
a <- 2.5
b <- 0.5
u <- 0.9

### Espérance, VaR et TVaR théorique
ExTh <- a/b 
VaRxTh <- qgamma(u, a, b)
TVaRxTh <- 1/(1 - u)*a/b*pgamma(VaRxTh, a + 1, b, low=FALSE)

#### Méthode Lower
### Fmp et cdf par discrétisation
fx <- pgamma(k, a, b) - pgamma(k-h, a, b)
Fx <- cumsum(fx)

## Pr(X = 3) et Pr(X <= 3) approximation
fx[3/h + 1]
Fx[3/h + 1]

### Espérance approximative
Ex <- sum(k*fx)

### VaR approximative
VaRx <- k[min(which(Fx >= u))]

### TVaR approximative
TVaRx <- VaRx + sum(pmax(k - VaRx, 0)*fx)/(1 - u)

### Vérif
cbind(ExTh, ExTh, Ex)
cbind(VaRxTh, VaRx)
cbind(TVaRxTh, TVaRx)

### Graphique
curve(pgamma(x, a, b), xlim = c(0, 15), lwd=2)
matplot(k, Fx, add = TRUE, type = 's', col = "blue", lty=2)



#### V.a.s composées et discrétisation -- X~PoComp(5, Fb), B~Pa(3, 100) ####
nfft <- 2^20
h <- 0.1
k <- (0:(nfft - 1))*h
lam <- 5
a <- 3
b <- 100

### Espérance et variance théorique
EmTh <- lam
EbTh <- b/(a - 1)
ExTh <- EmTh*EbTh

VmTh <- lam
VbTh <- mpareto(2, a, b) - EbTh^2
VxTh <- EmTh*VbTh + VmTh*EbTh^2

### Fmp et cdf par discrétisation pour la v.a. B
fb <- ppareto(k, a, b) - ppareto(k - h, a, b)

### V.a. X
fxt <- exp(lam*(fft(fb) - 1))
fx <- Re(fft(fxt, TRUE))/nfft
fx[300/h + 1]

### Espérnace et variance par approximation
Eb <- sum(fb*k)
Vb <- sum(fb*k^2) - Eb^2

Ex <- sum(fx*k)
Vx <- sum(fx*k^2) - Ex^2

### Verif
cbind(EbTh, Eb, VbTh, Vb)
cbind(ExTh, Ex, VxTh, Vx)

#### Aggrégation et discretisation -- X~Ga() et Y~LN() ####
nfft <- 2^13
h <- 0.1
k <- (0:(nfft - 1))*h
lam <- 5
a <- 5
b <- 0.5
sig <- 0.6
mu <- log(10) - sig^2/2
mlnorm(1, mu, sig)

### Espérance théorique
ExTh <- a/b
EyTh <- mlnorm(1, mu, sig)
EsTh <- ExTh + EyTh

### Fmp et cdf par discrétisation
fx <- pgamma(k, a, b) - pgamma(k - h, a, b)
fy <- plnorm(k, mu, sig) - plnorm(k - h, mu, sig)

### V.a. S
fst <- fft(fx)*fft(fy)
fs <- Re(fft(fst, TRUE))/nfft
fs[20/h + 1]

### Espérnace et variance par approximation
Eb <- sum(fx*k)
Ex <- sum(fy*k)
Es <- sum(fs*k)

### Verif
cbind(ExTh, Eb)
cbind(EyTh, Ex)
cbind(EsTh, Es)






