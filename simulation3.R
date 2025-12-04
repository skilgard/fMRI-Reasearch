source("/Users/stk/Downloads/grad_school/fMRI.research/source.code/Refft.R")

simulation.fn <- function(reps, beta.vec, n, seed){
  if(length(beta.vec) != 2){
    stop("beta.vec must have a length of 2, beta.vec[1] = yint and beta.vec[2] = slope")
  }
  
  if(n < 1 || n %% 1 != 0){
    stop("n must be a natural number")
  }
  
set.seed(seed)



# time series data
X <- scan("/Users/stk/Downloads/grad_school/fMRI.research/data/f.true", quiet = TRUE)
X <- rep(X, 3)
p <- length(X) # p = number of time points in our time series

# initialize vectors to store data
trad.b0 <- trad.b1 <- trad.b0.se <- trad.b1.se <- matrix(NA, nrow = reps, ncol = n)
rdft.b0 <- rdft.b1 <- rdft.b0.se <- rdft.b1.se <- matrix(NA, nrow = reps, ncol = n)

for(i in 1:reps) {
not.trans.error <- replicate(n = n, expr = sapply(X, function(x){rnorm(mean = 0, sd = (x + 0.005)*10, 1)}))
Error <- apply(not.trans.error, 2, function(x)ReFFTf(x))

y.sim.mat <- t(matrix(data = rep(beta.vec[1], p) + beta.vec[2] * X, byrow=TRUE, nrow = n, ncol = p)) + Error

# PART 2: linear regression using traditional lm()
#trad.lm <- lm(y.sim.mat ~ X, weights = (1/(X + 0.005))^2)
trad.lm <- lm(y.sim.mat ~ X)
trad.b0[i, ] <- unname(trad.lm$coefficients[1,])
trad.b1[i, ] <- unname(trad.lm$coefficients[2,])
for(j in 1:n){
  trad.b0.se[i, j] <- summary(trad.lm)[[j]]$coefficients[1,2]
  trad.b1.se[i, j] <- summary(trad.lm)[[j]]$coefficients[2,2]
}


# PART 3: linear regression into fourier space and then use lm()
#rdft.lm <- lm(apply(y.sim.mat, 2, function(x){ReFFTi(x)}) ~ 0 + ReFFTi(rep(1, p)) + ReFFTi(X), weights = (1/(X + 0.005))^2) # this essentially performs linear regression on each "subject" 

rdft.lm <- lm(apply(y.sim.mat, 2, function(x){ReFFTi(x)}) ~ 0 + ReFFTi(rep(1, p)) + ReFFTi(X)) # this essentially performs linear regression on each "subject" 
rdft.b0[i, ] <- unname(rdft.lm$coefficients[1,])
rdft.b1[i, ] <- unname(rdft.lm$coefficients[2,])
for(j in 1:n){
  rdft.b0.se[i, j] <- summary(rdft.lm)[[j]]$coefficients[1,2]
  rdft.b1.se[i, j] <- summary(rdft.lm)[[j]]$coefficients[2,2]
}

}

results <- list(trad.b0 = trad.b0, trad.b1 = trad.b1, trad.b0.se = trad.b0.se, trad.b1.se = trad.b1.se,
                rdft.b0 = rdft.b0, rdft.b1 = rdft.b1, rdft.b0.se = rdft.b0.se, rdft.b1.se = rdft.b1.se)
summary.output <- data.frame(rbind(c(beta.vec[1], NA, beta.vec[2], NA), 
                             c(mean(trad.b0), mean(trad.b0.se), mean(trad.b1), mean(trad.b1.se)), 
                             c(mean(rdft.b0), mean(rdft.b0.se), mean(rdft.b1), mean(rdft.b1.se))), 
                             row.names = c("True Values", "Using lm()", "Using lm() + RDFT"))
colnames(summary.output) <- c("y-int", "SE of y-int", "slope", "SE of slope")

# maybe add Error = Error 
return(list(y.sim.mat = y.sim.mat, results = results, 
            summary = summary.output, Error = Error))
}

# df1 <- simulation.fn(reps = 1, beta.vec = c(15, 47), n= 3, seed = 794398)
# df2 <- simulation.fn(reps = 1, beta.vec = c(29, 47), n= 3, seed = 794398)
# df3 <- simulation.fn(reps = 1, beta.vec = c(1.5, 4.7), n= 3, seed = 794398)
# df2 <- simulation.fn(reps = 3, beta.vec = c(29, 47), n= 3, seed =90722457)
# print(list("take one: traditional", df1[[1]], "take two: traditional", df2[[1]], 
#            "take one: rdft", df1[[2]], "take two: rdft", df2[[2]]))
# 
# df3 <- simulation.fn(reps = 1, beta.vec = c(29, 47), n= 13, seed = 794398)
# df4 <- simulation.fn(reps = 1, beta.vec = c(29, 47), n= 13, seed = 9835789)

