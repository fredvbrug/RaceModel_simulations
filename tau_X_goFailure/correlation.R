n     <- 20                    # length of vector
rho   <- 0.3                   # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle


x1    <- rnorm(n, 1, 1)        # fixed given data
x2    <- rnorm(n, 2, 0.5)      # new random data
X12     <- cbind(x1, x2)         # matrix
X12ctr  <- scale(X12, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(X12ctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% X12ctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(X12ctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x2b <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
cor(x1, x2b)    




x3    <- rnorm(n, 3, 1)      # new random data
X23     <- cbind(x2b, x3)         # matrix
X23ctr  <- scale(X23, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(X23ctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% X23ctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(X23ctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x3b <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
cor(x2b, x3b)   


cor(x1, x3b)   