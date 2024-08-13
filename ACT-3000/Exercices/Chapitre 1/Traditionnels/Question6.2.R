############
# Question 6 traditionelle chapitre  1 -----------------------------------------
############
nfft <- 2^10
r <- 0.4
q <- 2/3

fb <- c(discretise(pweibull(x, 0.8, 1000), 0, 4000, 1000, method = 'lower'),
        pweibull(4000, 0.8, 1000, low=F),
        numeric(nfft - 6))

fxt <- (q/(1 - (1 - q)*fft(fb)))^r
fx <- Re(fft(fxt, T))/nfft 
fx[1:8]

