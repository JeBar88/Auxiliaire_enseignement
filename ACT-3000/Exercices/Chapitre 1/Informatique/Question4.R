### ACT-3000
## Chapitre 1 informatique : Question 4
## Jérémie Barde

### Focntion
Panjer_BinNeg <- function(r, q, fb, smax){
  a <- 1 - q
  b <- a *(r-1)
  l <- length(fb)
  fx <- (q/(1 - a * fb[1]))^r
  fb <- c(fb, numeric(smax))
  for(i in 1:smax){
    j <- i + 1
    fx <- c(fx, sum((a + b * (1:i)/i) * fb[2:j] * fx[i:1]) * 1/(1 - a * fb[1]))
  }
  fx
}

### Variable
mu <- log(10) - 0.18
s <- 0.6
h <- 1
r <- 1.5
q <- 1/3
kmax <- 1000
k <- 0:kmax

### (a) 
fb_disc <- discretize(plnorm(x, mu, s), from = 0, to = 500, step = h, method = "lower")
sum(fb_disc)
fb_disc[16]

### (b) 
fx <- Panjer_BinNeg(r, q, fb_disc, kmax)
fx[16]

### (c)
Fx <- cumsum(fx)
Fx[61]
sum(pmax(k - 50, 0) * fx)
sum(fx * k)

### (e)
u <- Fx[51]
VaR <- k[min(which(Fx >= u))]
TVaR <- sum(pmax(k - VaR, 0) * fx)/(1 - u) + VaR
cbind(VaR, TVaR)





