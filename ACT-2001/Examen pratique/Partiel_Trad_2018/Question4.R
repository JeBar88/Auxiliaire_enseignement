#################################
### Examen partiel trad 2018 : Q4
#################################
fX1X2 <- matrix(c(9/180, 0, 0, 71/180, 20/180, 71/180, 0, 0, 9/180), ncol = 3)
fY1Y2 <- matrix(c(4/180, 1/180, 4/180, 72/180, 18/180, 72/180, 4/180, 1/180, 4/180), ncol = 3)

fX1 <- apply(fX1X2, 1, sum)
fX2 <- apply(fX1X2, 2, sum)
fY1 <- apply(fY1Y2, 1, sum)
fY2 <- apply(fY1Y2, 2, sum)

FX1 <- cumsum(fX1)
FX2 <- cumsum(fX2)
FY1 <- cumsum(fY1)
FY2 <- cumsum(fY2)

X1 <- c(0, 10, 20)
X2 <- c(0, 20, 40)
Y1 <- c(0, 10, 20)
Y2 <- c(0, 20, 40)

EX1 <- sum(fX1 * X1)
EX2 <- sum(fX2 * X2)

VarX1 <- sum(fX1 * X1^2) - EX1^2
VarX2 <- sum(fX2 * X2^2) - EX2^2
VarY1 <- sum(fY1 * Y1^2) - sum(fY1 * Y1)^2
VarY2 <- sum(fY2 * Y2^2) - sum(fY2 * Y2)^2
cbind(VarX1, VarX2, VarY1, VarY2)

phik <- dnorm(qnorm(0.99))/(1-0.99)

rhoX1 <- EX1 + sqrt(VarX1) * phik
rhoX2 <- EX2 + sqrt(VarX2) * phik                      
cbind(rhoX1, rhoX2)
