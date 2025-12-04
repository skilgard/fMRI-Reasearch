# simulation1.R
# created on 08/29/2025, last edited 9/2/2025
source("/Users/stk/Downloads/grad_school/fMRI research/source.code/rdft.R") 
# this source file uses make.rdft.matrix(m=number of variable) to make a real discrete fourier tranform 
source("/Users/stk/Downloads/grad_school/fMRI research/source.code/square.root.mat.R")
# this source file uses square.root.matrix(A=matrix, symmetric==TRUE,by default) to find a square root matrix for A 


# seed = 319 (generated from floor(abs(rnorm(1, mean=100, sd=100)))
seed <- 319
set.seed(319)

# 1. generate a fourier vector variate time series, AND save the beta1 and beta0 used to generate this data 

m <- 100

# if Sigma_m has a fourier vector variate covariance structure, Sigma_m = Gamma_m Delta t(Gamma_m), where Gamma_m = rdftmatrix, and Delta = eigenvalues of Sigma_m in a diagonal matrix 

Delta <- diag(runif(n=m, min=0.1, max=100))
Gamma_m <- make.rdft.matrix(m)
Sigma_m <- Gamma_m %*% Delta %*% t(Gamma_m)

# beta_0 and beta_1 are in this vector
#beta.vec <- runif(n=2, min=-10, max=10)
 beta.vec <- c(2, -4)

X.vec <- runif(n=m, min=0, max=200)

Z.vec <- rnorm(n=m, mean=0, sd = 1) # to be used for generating data as in (16) on page 10 of fourier covariance paper 

sq.root.Sigma_m <- square.root.matrix(A = Sigma_m)

# this is the simulated data 
sim.y <-  rep(x=beta.vec[1], times=m) + beta.vec[2] * X.vec + sq.root.Sigma_m %*% Z.vec

# 2) Use paper to find beta[2] = beta_1 again (and calculate SE)

rdft.lm <- lm(formula = Sigma_m %*% sim.y ~ Sigma_m %*% X.vec)

summary(rdft.lm)

#Call:
#lm(formula = Sigma_m %*% sim.y ~ Sigma_m %*% X.vec)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-825.73 -369.37   -2.11  284.10 1275.50 
#
#Coefficients:
#                   Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)        90.00475  128.31349    0.701    0.485    
#Sigma_m %*% X.vec  -3.99147    0.01455 -274.395   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 470.2 on 98 degrees of freedom
#Multiple R-squared:  0.9987,	Adjusted R-squared:  0.9987 
#F-statistic: 7.529e+04 on 1 and 98 DF,  p-value: < 2.2e-16

# or using equation from page 11, (21)

rdft.beta2 <- c(0,0)
# the first column if for b_0 and the second one for b_1
X.mat <- cbind(rep(0, m), X.vec)
for(j in 1:2){
	rdft.beta2[j] <- solve(sum(apply(X.mat, 1, function(x){x %*% t(x)}))) %*% sum(apply(cbind(X.mat, Z.vec), 2, function(x){x[1:2] * x[3]}))
}

# 3) Use previous methods to find beta[2] = beta_1 again (and calculate SE) 

traditional.lm <- lm(formula = sim.y ~ X.vec)

summary(traditional.lm)
#Call:
#lm(formula = sim.y ~ X.vec)

#Residuals:
#    Min       1Q   Median       3Q      Max 
# -17.4793  -4.6671  -0.2164   4.0732  18.8762 

# Coefficients:
           # Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)  1.50346    1.29113    1.164    0.247    
#X.vec       -3.99535    0.01086 -368.061   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 6.775 on 98 degrees of freedom
#Multiple R-squared:  0.9993,	Adjusted R-squared:  0.9993 
#F-statistic: 1.355e+05 on 1 and 98 DF,  p-value: < 2.2e-16


