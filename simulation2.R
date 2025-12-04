# created on 9/6/2025
# last edited on 9/25/2025
# simulation2.R
# this source file uses make.rdft.matrix(p=number of time points) to make a real discrete fourier tranform
source("/Users/stk/Downloads/grad_school/fMRI.research/source.code/rdft.R") 
source("/Users/stk/Downloads/grad_school/fMRI.research/source.code/Refft.R")


# last time i ran it was: simulation.fn(reps = 100, beta.vec = c(0.9, -1.8), seed = 319)
# seed = 319 (generated from floor(abs(rnorm(1, mean=100, sd=100)))
# beta.vec is fixed 2 entry vector of y-int and slope, in that order 
# n = sample size
simulation.fn <- function(reps, beta.vec, seed = 319, n){
if(length(beta.vec) != 2){
stop("beta.vec must have a length of 2, beta.vec[1] = yint and beta.vec[2] = slope")
}

set.seed(seed)

# time series data 
X <- scan("/Users/stk/Downloads/grad_school/fMRI.research/data/f.true")
p <- length(X) # p = number of time points in our time series 

results <- matrix(nrow=reps, ncol=8)

colnames(results) <- c("rdft b0", "trad b0", "rdft b0 SE", "trad b0 SE", "rdft b1", "trad b1", "rdft b1 SE", "trad b1 SE")

for(i in 1:reps){
### PART 1: simulate data 

# error: add some N(0,1) noise to the simulated data 
epsilon <- rnorm(mean=0, sd=1, n=p)

sim.y <- rep(beta.vec[1], p) + beta.vec[2] * X + ReFFTf(epsilon)
### PART 2: linear regression using traditional lm()
# trad.lm <- lm(sim.y ~ cbind(rep(1, p),X))
# the cbind() does not work 
trad.lm <- lm(sim.y ~ X)

trad.b0 <- unname(trad.lm$coefficients[1])
trad.b1 <- unname(trad.lm$coefficients[2])
trad.b1.se <- unname(summary(trad.lm)$coefficients[2,2])
trad.b0.se <- unname(summary(trad.lm)$coefficients[1,2])

### PART 3: transform into fourier space using C and then use lm() 
# the cbind() does not work 
# rdft.lm <- lm(ReFFTf(sim.y) ~ ReFFTf2(cbind(rep(1, p),X)))
# rdft.lm <- lm(ReFFTf(sim.y) ~ ReFFTf(X))
# rdft.lm <- lm(ReFFTf(sim.y) ~ 0 + ReFFTf(rep(1, p)) + ReFFTf(X) )# is this what was meant?
 rdft.lm <- lm(ReFFTi(sim.y) ~ 0 + ReFFTi(rep(1, p)) + ReFFTi(X) )# is this what was meant?

rdft.b0 <- unname(rdft.lm$coefficients[1])
rdft.b1 <- unname(rdft.lm$coefficients[2])
rdft.b1.se <- unname(summary(rdft.lm)$coefficients[2,2])
rdft.b0.se <- unname(summary(rdft.lm)$coefficients[1,2])

#results[i, ] <- c(rdft.b0, rdft.b1, rdft.b1.se, rdft.b0, rdft.b1, trad.b1.se)
results[i, ] <- c(rdft.b0, trad.b0, rdft.b0.se, trad.b0.se, rdft.b1, trad.b1, rdft.b1.se, trad.b1.se)
}

summary.means <- apply(results, 2, mean)
summary.output <- data.frame(rbind(c(beta.vec[1], NA, beta.vec[2], NA), summary.means[c(1, 3, 5, 7)], 
summary.means[c(2, 4, 6, 8)]), row.names=c("true", "rdft mean", "trad lm mean"))
colnames(summary.output) <- c("b0", "b0 SE", "b1", "b1 SE")
return(list(results = results, beta.vec = beta.vec, summary = summary.output))
}
