# square.root.matrix() is a function that calculates the square root of a matrix of a PD matrix 
# note: due to rounding, identical(B%*%B, og.mat) will most likely not return TRUE (where B is the sq root matrix)
#, but all.equal(B%*%B, og.mat) will return TRUE 
# the p.root.matrix() calculates the pth root of the matrix


#p = the root you're interest in 
# A = the PD matrix you want to find the pth root of
# if A is not symmetric and you are SURE that is PD, switch to symmetric is true to bypass the error catching  
p.root.matrix <- function(p, A, symmetric = TRUE) {
# error catching for a symmetric matrix 
if(symmetric == TRUE){
if(all(eigen(A)$values > 0) == FALSE){
stop("symmetric matricies must have positive eigenvalues in order to by positive definite")
}}

E <- eigen(A) 
V <- E$vectors
U <- solve(V) 
D <- diag(E$values) 
B <- V %*% D^(1/p) %*% U

return(B)
}

square.root.matrix <- function(A, symmetric = TRUE){

B <- p.root.matrix(p=2, A=A, symmetric = symmetric)

return(B)
}
