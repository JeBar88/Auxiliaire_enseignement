############
# Question 2 traditionnel chapitre  1 ------------------------------------------
############
nFFT <- 2^10
n <- c(40, 50, 30)
lam <- 0.04
r <- c(2, 3)
q <- c(0.97, 0.99)

fm1 <- dpois(0:(nFFT - 1), lam)
fm2 <- dnbinom(0:(nFFT - 1), r[1], q[1])
fm3 <- dnbinom(0:(nFFT - 1), r[2], q[2])

fnt <- fft(fm1)^n[1] * fft(fm2)^n[2] * fft(fm3)^n[3]
fn <- Re(fft(fnt, T))/nFFT

fn[1:4]
