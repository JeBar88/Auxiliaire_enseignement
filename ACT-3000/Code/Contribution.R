### ACT-3000
## Calcule de contriubiton au mesure VaR et TVaR
## Jérémie Barde

nfft <- 2^16
x <- 0:(nfft - 1)
n <- c(200, 300)
q <- c(0.15, 1/15)
k <- 0.9

fx <- dbinom(x, n[1], q[1])
fy <- dbinom(x, n[2], q[2])

espi <- n*q

fst <- fft(fx) * fft(fy)
fs <- Re(fft(fst, T))/nfft
Fs <- cumsum(fs)
Fs[201]

sum(fs*x)

VaRS <- x[min(which(Fs >= k))]
TVaRS <- VaRS + 1/(1-k)*sum(pmax(x - VaRS, 0) * fs)

gx <- x * fx/espi[1]
gy <- x * fy/espi[2]

faxt <- fft(gx) * fft(fy) 
fax <- Re(fft(faxt, T))/nfft

fayt <- fft(gy) * fft(fx)
fay <- Re(fft(fayt, T))/nfft

esp_aX <- espi[1] * fax/fs
esp_aY <- espi[2] * fay/fs

esp_aX[502]
esp_aY[502]

CVaRX <- esp_aX[VaRS + 1]
CVaRY <- esp_aY[VaRS + 1]
cbind(CVaRX, CVaRY, CVaRX+CVaRY, VaRS)

CTVaRX <- (sum(esp_aX[x > VaRS & x < 500] * fs[x > VaRS & x < 500]) + CVaRX * (Fs[VaRS + 1] - k))/(1-k) 
CTVaRY <- (sum(esp_aY[x > VaRS & x <= 500] * fs[x > VaRS & x <= 500]) + CVaRY * (Fs[VaRS + 1] - k))/(1-k)
cbind(CTVaRX, CTVaRY, CTVaRX+CTVaRY, TVaRS)



