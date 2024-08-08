### ACT-2001 ###
### Examen final info 2018 : Question 1
b <- 0.001
c <- 0.1

Fx <- function(x, n){
  pgamma(x, n*c, n*b)
}

Fx(300, 4)

VaRX <- function(k, n){
  qgamma(k, n*c, n*b)
}

VaRX(0.9999, 4)
VaRX(0.9999, 100)

TVaRX <- function(k, n){
  c/b * pgamma(VaRX(k, n), n*c + 1, b*n, low=F)/(1-k)
}

TVaRX(0.9999, 4)
TVaRX(0.9999, 100)
