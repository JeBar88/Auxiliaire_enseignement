############################
### Examen partiel info 2015
############################

### question 2)
## a)
a <- 1.5 + 0.1 * 1:10
lam <- 100*(a - 1)

VaR <- function(k) sum(sapply(1:10, function(i) qpareto(k, a[i], lam[i])))

Fs <- function(y) optimize(function(x) abs(VaR(x) - y), c(0, 1))$min
Fs(2000)

### Question 4)
nfft <- 2^10

fk <- c(0.6, 0, 0.1, 0.1, 0, 0.2, numeric(nfft - 6))
length(fk)
sum(fk)

fnt <- fft(fk)^20
fn <- Re(fft(fnt, T))/nfft
sum(fn)

fn[c(1, 6, 21)]

10^6
