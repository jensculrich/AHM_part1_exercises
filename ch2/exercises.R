
##--------------------------------------------------------------
# 1. bayes rule: x | y = 
# [y|x](conditional of y)*[x](marginal of x) / [y](marginal of y)

# marginal of x = sum of probabilities of y for x. We will need the joint distr first

# joint distr:
y <- 0:5
x <- 0:5
# p is the inv logit transformation of the realized values of x
p <- plogis(-1.2 + 2*x)
p
# joint distribution [y, x]
lambda <- 0.4
joint <- matrix(NA, length(y), length(x))
rownames(joint) <- paste ("y=", y, sep="")
colnames(joint) <- paste ("x=", x, sep="")

for(i in 1:length(y)) {
  joint[,i] <- dbinom(y, 5, p[i]) * dpois(x[i], lambda)
}

round(joint, 3)

# marginal of x
margX <- colSums(joint)
# marginal of y 
margY <- rowSums(joint)


# conditional likelihood of x|y =
# [Y, X](joint distribution) / [Y] (marginal distr)
YgivenX <- joint / matrix(margX, nrow(joint), ncol(joint), byrow = TRUE)
round(YgivenX, 2)

# multiply each row by margX
temp <- YgivenX*rep(margX,each=nrow(YgivenX))
conditionalX <- temp / margY
(round(conditionalX, 2))

##--------------------------------------------------------------
# 3. 

n <- 1000

mean <- 0 
sd <- 1
sim_data_logit <- runif(n, mean, sd)
hist(sim_data_logit)
hist(plogis(sim_data_logit), xlim = c(0,1))

mean <- 0 
sd <- 1
sim_data_logit <- rnorm(n, mean, sd)
hist(sim_data_logit)
hist(plogis(sim_data_logit), xlim = c(0,1))

mean <- 0 
sd <- 1.5
sim_data_logit <- rnorm(n, mean, sd)
hist(sim_data_logit)
hist(plogis(sim_data_logit), xlim = c(0,1))

mean <- 0 
sd <- 2
sim_data_logit <- rnorm(n, mean, sd)
hist(sim_data_logit)
hist(plogis(sim_data_logit), xlim = c(0,1))

mean <- 0 
sd <- 5
sim_data_logit <- rnorm(n, mean, sd)
hist(sim_data_logit)
hist(plogis(sim_data_logit), xlim = c(0,1))

##--------------------------------------------------------------
# 2.
# original
# ------------------------------------------------------------------------

sim.data <- function(beta0 = -3, beta1 = 2, p = 0.6, x=NULL){
  # Function allows input of covariate "x", or simulates new
  
  M <- 100
  if(is.null(x)) # x - covariate - a vector of veg heights
    vegHt <- runif(M, 1, 3) # uniform from 1 to 3
  
  # Suppose that occupancy probability increases with vegHt
  # The relationship is described (default) by an intercept of -3 and
  #    a slope parameter of 2 on the logit scale
  # plogis is the inverse-logit (constrains us back to the [0-1] scale)
  psi <- plogis(beta0 + beta1*vegHt)
  
  # Now we simulated true presence/absence for 100 sites
  z <- rbinom(M, 1, psi)
  
  # Now generate observations
  J <- 3 # sample each site 3 times
  y <- rbinom(M,J,p*z)
  
  list(y=y, J=J, vegHt=vegHt)
}

# This is the negative log-likelihood based on the marginal distribution
# of y. It is the pmf of a zero-inflated binomial random variable.
#
negLogLikeocc <- function(beta, y, x, J) {
  beta0 <- beta[1]
  # beta1 <- beta[2]
  p <- plogis(beta[2])
  # p <- plogis(logitp)
  # psi <- plogis(beta0 + beta1*x)
  # let's try a "bad" model below
  psi <- plogis(beta0)
  # marginal likelihood is likelihood of detecting on y occassions out of J given p.
  # TIMES averaged over all possible values of x (which we may or may not have info about)
  marg.likelihood <- dbinom(y, J, p) * psi + ifelse(y==0, 1, 0) * (1-psi)
  return(-sum(log(marg.likelihood)))
}

data <- sim.data()        # Generate a data set

# Let's minimize it
# starting.values <- c(beta0=0, beta1=0, logitp=0)
starting.values <- c(beta0=0, logitp=0) # bad model
opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x=data$vegHt,
                 J=data$J, hessian=TRUE)
(mles <- opt.out$par)

# Make a table with estimates, SEs, and 95% CI
mle.table <- data.frame(Est=mles,
                        SE = sqrt(diag(solve(opt.out$hessian))))
mle.table$lower <- mle.table$Est - 1.96*mle.table$SE
mle.table$upper <- mle.table$Est + 1.96*mle.table$SE
mle.table


# Define a fit statistic
fitstat <- function(y, Ey){
  sum((sqrt(y) - sqrt(Ey)))
} # distance from expected value
# Compute it for the observed data
# ~~~ 3 lines of code added to ensure we are using output from sim.data(), see Errata 2021-10-09
y <- data$y
J <- data$J
vegHt <- data$vegHt
# T.obs <- fitstat(y, J*plogis(mles[1] + mles[2]*vegHt)*plogis(mles[3]))
T.obs <- fitstat(y, J*plogis(mles[1])*plogis(mles[2]))

# Get bootstrap distribution of fit statistic
T.boot <- rep(NA, 100)
for(i in 1:100){
  # Simulate a new data set and extract the elements. Note we use
  # the previously simulated "vegHt" covariate
  data <- sim.data(beta0=mles[1],p=plogis(mles[2]),x=vegHt)
  # Next we fit the model
  starting.values <- c(0,0)
  opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x= data$vegHt, J=data$J, hessian=TRUE)
  (parms <- opt.out$par)
  # Obtain the fit statistic
  T.boot[i]<- fitstat(y, J*plogis(parms[1])*plogis(parms[2]) )
}

(T.obs)

summary(T.boot)
