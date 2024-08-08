### ACT-2001 ###
### Examen info 2018 : Question 2
m <- 1E5
a <- c(0.5, 1.5, 2.5)
b <- c(0.05, 0.15, 0.25)

set.seed(2018)
U <- matrix(runif(3*m), ncol = 3, byrow = T)
head(U, 2);tail(U, 1)


x <- qgamma(U[,1:3], a, b)
head(x)
tail(x, 1)

x[4, 2]

x1 <- qgamma(U[,1:3], a[1], b[1])
head(x1, 2)
x1 <- qgamma(U[,1:3], a[1], b[1])
head(x1, 2)

x <- matrix(sapply(1:3, function(i) qgamma(U[,i], a[i], b[i])), ncol = 3)
head(x, 4)

qgamma(U[3,1], a[2], b[2])
