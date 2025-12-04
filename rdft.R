# rdft.R

# This source file's function can be used to make an orthogonal matrix with real discrete fourier transform structure (From "Fourier-Structured Tensor-Variate Distributions and PCA" by Carlos Llosa-Vite and Ranjan Maitra)

# function name: make.rdft.matrix()
# input:
## m = number of variables 
# output:
## rdft.mat - the real discrete fourier transformation matrix for in an m-variate situation

# the formula is: (1/sqrt(m))[1_(mx1 vector), sqrt(2) Mc, +-1_(mx1 vector), sqrt(2)Ms)

# Mc matrix, dimension of m x floor((m-1)/2), (i,j)th element = cos(2pi j(i-1)/m)

# Ms matrix, dimension of m x floor((m-1)/2), (i,j)th = sin(2 pi j(i-1)/m)

# note: if m is odd, rdft.mat has no +-1_(mx1 vector)

# makes the Ms and Mc matricies 
make.MsMc <- function(m) {
# if m <= 2, Ms and Mc do not exist
if(m <= 2) {Ms <- c(); Mc <- c()
} else{ 
        indexes <- cbind(rep(1:m, each=floor((m-1)/2)), rep(1:floor((m-1)/2), times=m))
        Mc <- matrix(apply(indexes, MAR = 1, function(y){cos(2 * pi * y[2] * (y[1] - 1) / m)}), 
	nrow=m, ncol=floor((m-1)/2), byrow=TRUE)
        Ms <- matrix(apply(indexes, MAR = 1, function(y){sin(2 * pi * y[2] * (y[1] - 1) / m)}),
        nrow=m, ncol=floor((m-1)/2), byrow=TRUE)
}
return(list("Mc" = Mc, "Ms" = Ms))
}

make.rdft.matrix <- function(m)
{
# error statement
if(m %% 1 != 0 || m < 1) {
stop("m must be an natural number")
}
one.vec <- rep(1, m)

MsMc.list <- make.MsMc(m)
Ms <- MsMc.list$Ms
Mc <- MsMc.list$Mc

# even numbers of variables have a vector of (1, -1, 1, -1, ...) in the middle
if(m %% 2 == 0){
pm.one.vec <- rep(c(1,-1), m/2)
rdft.mat <- (1/sqrt(m)) * cbind(one.vec, sqrt(2) * Mc, pm.one.vec, sqrt(2) * Ms) 
} else{
rdft.mat <- (1/sqrt(m)) * cbind(one.vec, sqrt(2) * Mc, sqrt(2) * Ms)
}
return(unname(rdft.mat))
}

# old version with for loops
#make.Mc <- function(m){
#if(m <= 2) {Mc <- c() # Mc doesn't exist if m <= 2
#} else{
#Mc <- matrix(0, nrow=m, ncol = floor((m-1)/2))
#for(i in 1:m){
#       for(j in 1:floor((m-1)/2)){
#       Mc[i,j] <- cos(2*pi*j*(i - 1)/m)
#}}
#}
#return(Mc) 
#}

# old version with for loops 
#make.Ms <- function(m){
#if(m <= 2){Ms <- c() # Mc doesn't exist if m <= 2
#} else{
#Ms <- matrix(0, nrow=m, ncol = floor((m-1)/2))
#for(i in 1:m){
#        for(j in 1:floor((m-1)/2)){
 #       Ms[i,j] <- sin(2*pi*j*(i - 1)/m)
#}}
#}
#return(Ms)
#}

