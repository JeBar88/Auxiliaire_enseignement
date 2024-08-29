### ACT-3000
## Chapitre 1 informatique : Question 1
## Jérémie Barde

DePril <- function(ff, nn=5, smax=100){
  ll <- length(ff)
  ffs <- ff[1]^nn
  ff <- c(ff, rep(0, smax - ll + 1))
  for(i in 1:smax){
    j <- i+1
    ffs <- c(ffs, 1/ff[1] * sum(ff[2:j] * ffs[i:1] * ((nn+1) * (1:i)/i-1)))
    
  }
  ffs
}

fm <- c(0.5, 0.1, 0.2, 0.15, 0.05)

fs <- DePril(fm, 10, 100)
sum(fs)

fs[c(1, 2, 11, 21)]

# ou
nFFT <- 2^5
fmt <- fft(c(fm, numeric(nFFT - 5)))
fst <- fmt^10
fs <- Re(fft(fst, TRUE))/nFFT
