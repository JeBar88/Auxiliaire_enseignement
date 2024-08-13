############
# Question 2 informatique chapitre  1 ------------------------------------------
############

x <- 0:100
fm <- 0.6 * I(x == 0) + 0.3 * dpois(x, 0.1) + 0.1 * dpois(x, 0.2)
sum(fm)


Em <- sum(fm * 0:100)
Em100 <- 100 * E
cbind(Em, Em100)

DePril <-function(ff,nn=5,smax=100){
  ll<-length(ff)
  ffs<-ff[1]^nn
  ff<-c(ff,rep(0,smax-ll+1))
  for (i in 1 :smax)
  {
    j<-i+1
    ffs<-c(ffs,(1/ff[1])*sum(ff[2 :j]*ffs[i :1]*((nn+1)*(1 :i)/i-1)))
  }
  return(ffs)
}

fn <- DePril(fm, 100, 100)
sum(fn)
En <- sum(fn * 0:100)

fn[c(1, 6, 11, 16)]

## Avec FFT
nFFT <- 2^10
lam <- c(0.1, 0.2)

fB <- c(0, 1, rep(0, nFFT-2))
fBt <- fft(fB)
fnt <- (0.6 +  0.3 * exp(lam[1] * (fBt - 1)) + 0.1 * exp(lam[2] * (fBt - 1)))^100
fn <- Re(fft(fnt, inverse = TRUE))/nFFT

sum(fn)
En <- sum(fn * 0:(nFFT-1))
fn[c(1, 6, 11, 16)]

radians(1)

