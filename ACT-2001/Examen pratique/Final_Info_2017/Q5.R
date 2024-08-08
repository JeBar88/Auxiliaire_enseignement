### ACT-2001 ###
### Examen final info 2017 : Question 5
a <- c(1.5, 4.5)
b <- 1/100
k <- 0.9

Fs <- function(x) 2/4 * pgamma(x, a[1], b) + 2/4 * pgamma(x, a[2], b)

VaR <- optimise(function(x) abs(Fs(x) - k), c(0, 10000))$min

TVaR <- 2/4 * 1/(1-0.9) * a[1]/b * pgamma(VaR, a[1] + 1, b, low=F) + 
  2/4 * 1/(1-0.9) * a[2]/b * pgamma(VaR, a[2] + 1, b, low=F)
cbind(VaR, TVaR)