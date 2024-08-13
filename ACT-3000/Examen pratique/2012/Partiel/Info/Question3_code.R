###############################################
### Examen informatique partiel 2012 question 3
###############################################
nfft <- 2^15
r <- 0.75*1:2
q <- 0.2 * 1:2
b <- 0.01
a <- 2
k <- 0:(nfft - 1)

fm1 <- dnbinom(k, r[1], q[1])
fi <- c(numeric(a), 1, numeric(nfft - a - 1))

vkt <- fft(fm1) * (q[2]/(1 - (1 - q[2]) * fft(fi)))^r[2]
vk <- Re(fft(vkt, T))/nfft
sum(vk)

vk[1:6]


Fs <- function(x) vk[1] + sum(vk[-1] * pgamma(x, 1:(nfft - 1), b))
Fs(1000)
