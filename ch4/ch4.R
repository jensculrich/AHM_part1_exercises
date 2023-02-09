# Simulate some point count data

M <- 267 # number of sites
J <- 3 # Number of temporal replicates (counts)

# initialize random number generator
set.seed(24)
# repeatedly running the code without the above seed will give you a better feel
# for the stochastic process that generates data sets and for the sampling error.
# Use the seed to recover values presented in the textbook.

# Covariate data
elev <- runif(n = M, -1, 1)
forest <- runif(n = M, -1, 1)
wind <- array(runif(n = M*J, -1, 1), dim = c(M, J))

# mean abundance when all other covariates are at 0, (in terms of actual abundance)
mean.lambda <- 2
# and log abundnace
beta0 <- log(mean.lambda)

beta1 <- -2 # effect (slope) of elevation
beta2 <- 2 # effect of forest cover
beta3 <- 1 # interaction effect of elev and forest

log.lambda <- beta0 + beta1*elev + beta2*forest + beta3*elev*forest
lambda <- exp(log.lambda)

# plot everything
par(mfrow = c(2,2), mar = c(5,4,2,2), cex.main=1)
curve(exp(beta0+beta1*x), -1, 1, col="red", frame.plot = FALSE, ylim = c(0, 18),
      xlab = "elevation", ylab = "lambda", lwd = 2)
text(-0.9, 17, "A", cex = 1.5)
# B has variation because forest cover (another covariate) is also changing at the same time
# but we can see a general downward trend across elevation. Similar pattern with plot D
plot(elev, lambda, frame.plot = FALSE, ylim = c(0, 38), xlab = "Elevation", ylab = "")
text(-0.9, 36, "B", cex = 1.5)
curve(exp(beta0 + beta2*x), -1, 1, col = "red", frame.plot = FALSE,
      ylim = c(0, 18), xlab = "Forest cover", ylab = "lambda", lwd = 2)
text(-0.9, 17, "C", cex = 1.5)
plot(forest, lambda, frame.plot = FALSE, ylim = c(0, 38), xlab = "Forest cover", ylab = "")
text(-0.9, 36, "D", cex = 1.5)
par(op)

# Compute expected abundance for a grid of elevation and forest cover
cov1 <- seq(-1, 1,,100)                       # Values for elevation
cov2 <- seq(-1,1,,100)                        # Values for forest cover
lambda.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix,
# for every combination of values of elevation and forest cover
for(i in 1:100){
  for(j in 1:100){
    lambda.matrix[i, j] <- exp(beta0 + beta1 * cov1[i] + beta2 * cov2[j] +
                                 beta3 * cov1[i] * cov2[j])
  }
}

op <- par(mfrow = c(1, 2), mar = c(5,4,3,2), cex.main = 1.6)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x = cov1, y = cov2, z = lambda.matrix, col = mapPalette(100),
      xlab = "Elevation", ylab = "Forest cover", cex.lab = 1.2)
contour(x = cov1, y = cov2, z = lambda.matrix, add = TRUE, lwd = 1)
title(main = "A")
matpoints(elev, forest, pch="+", cex=0.8)
par(op)

# build stocahsticity into the relationships by describimng the random variation around the expected value lambda
N <- rpois(n = M, lambda = lambda)   # Realised abundance
sum(N)                               # Total population size at M sites
table(N)                             # Frequency distribution of tit abundance

# now move towards detection (observation process and its outcome)
mean.detection <- 0.3            # Mean expected detection
alpha0 <- qlogis(mean.detection) # same on logit scale (intercept)
alpha1 <- 1                      # Effect (slope) of elevation
alpha2 <- -3                     # Effect (slope) of wind speed
alpha3 <- 0                      # Interaction effect (slope) of elev and wind

logit.p <- alpha0 + alpha1 * elev + alpha2 * wind + alpha3 * elev * wind
p <- plogis(logit.p)             # Inverse link transform
mean(p)                          # average per-individual p is 0.38

op <- par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", frame.plot = FALSE,
      ylim = c(0, 1.1), xlab = "Elevation", ylab = "p", lwd = 2)
text(-0.9, 1.05, "A", cex = 1.5)
# matplot plots the matrix by colour?
matplot(elev, p, pch = "*", frame.plot = FALSE, ylim = c(0, 1.1),
        xlab = "Elevation", ylab = "")
text(-0.9, 1.05, "B", cex = 1.5)
curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", frame.plot = FALSE,
      ylim = c(0, 1.1), xlab = "Wind speed", ylab = "p", lwd = 2)
text(-0.9, 1.05, "C", cex = 1.5)
matplot(wind, p, pch = "*", frame.plot = FALSE, ylim = c(0, 1.1),
        xlab = "Wind speed", ylab = "p")
text(-0.9, 1.05, "D", cex = 1.5)
par(op)

# Compute expected detection probability for a grid of elevation and wind speed
cov1 <- seq(-1, 1,,100)                  # Values of elevation
cov2 <- seq(-1,1,,100)                   # Values of wind speed
p.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix which combines every value in cov 1 with every other in cov2
for(i in 1:100){
  for(j in 1:100){
    p.matrix[i, j] <- plogis(alpha0 + alpha1 * cov1[i] + alpha2 * cov2[j] + alpha3 * cov1[i] * cov2[j])
  }
}
image(x = cov1, y = cov2, z = p.matrix, col = mapPalette(100),
      xlab = "Elevation", ylab = "Wind speed", cex.lab = 1.2)
contour(x = cov1, y = cov2, z = p.matrix, add = TRUE, lwd = 1)
title(main = "B")
matpoints(elev, wind, pch="+", cex=0.7, col = "black")

# now actually simulate some counts given our stochatsically generated abundances
C <- matrix(NA, nrow = M, ncol = J)      # Prepare array for counts
for (i in 1:J){                          # Generate counts
  C[,i] <- rbinom(n = M, size = N, prob = p[,i])
}

head(cbind("True N" = N, "1st count" = C[,1], "2nd count" = C[,2],
           "3rd count" = C[,3]), 10)                    # First 10 rows (= sites)

table(C)

op <- par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
matplot(elev, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38),
        xlab = "Elevation", ylab = "Count (C)")
text(-0.9, 36, "A", cex = 1.5)
matplot(forest, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38),
        xlab = "Forest cover", ylab = "Count (C)")
text(-0.9, 36, "B", cex = 1.5)
matplot(wind, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38),
        xlab = "Wind speed", ylab = "Count (C)")
text(-0.9, 36, "C", cex = 1.5)
hist(C, breaks = 50, col = "grey", ylim = c(0, 460), main = "", xlab = "Count (C)")
text(3, 450, "D", cex = 1.5)
par(op)

real <- sum(N)                   # True total abundance (all sites)
obs <- sum(apply(C, 1, max))    # 'Observed' total abundance (all sites)
# combined estimation error is an underestimate of pop size by: _%
(estimation_error <- (real - obs) / real)

real <- sum(N>0)                 # True number of occupied sites
obs <- sum(apply(C, 1, max)>0)  # 'Observed' number of occupied sites
# combined estimation error is an underestimate of occurence by: _%
(estimation_error <- (real - obs) / real)
