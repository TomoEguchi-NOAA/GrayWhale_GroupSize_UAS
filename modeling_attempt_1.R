# data exploration
# Look at Trevor's data files:
# GROUPS_BEST.csv
# GROUPS_MIN.csv
# SIGHTINGS_GROUPS.csv
# 
# SIGHTING_ID- Concatenation of Date and Sighting Number within each day
# GROUP_SIZE_LAST - Last (and presumably best?) visual group size estimate
# Group_Size_Min - minimum UAS group size estimate (every group flown over will have one of these)
# Group_Size - best UAS group size estimate (only groups where I had a good idea of "true" group size were assigned abest group size estimate)
# Proximity - Ordinal score (0 = near perfect alignment in space and time, 1 = not perfect but an acceptable distance forinclusion in the analysis, 2 = UAS group too far from visual group to be included in analysis); I've only merged groups withscores 0-1 in these tables
# Other.Confounding.Groups - Ordinal score (0 = no other groups around that could confound linkage, 1 = other possiblegroups around but high enough confidence in linkage to be included in analysis, 2 = too many potential alternative groupsto be included in analysis); I've only merged groups with scores 0-1 in these tables
# DIST_MIN- minimum radial distance from observer trailer to the nearest visual observation of the SIGHTING_ID(calculated from bearings and reticles - so a fair amount of error)
# BEAUFORT_LAST - last sea state observation associated with a group
# VISIBILITY_LAST - last visibility code associated with a group
# 
# 

library(tidyverse)
library(ggplot2)
library(jagsUI)
library(bayesplot)
library(loo)

compute.LOOIC <- function(loglik.array, data.array, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  # Convert the log-likelihood and data arrays into vectors
  loglik.vec <- c(loglik.array)
  data.vec <- c(data.array)
  
  # Convert the log-likelihood vector into a 2D array
  n.data <- length(data.vec)
  loglik.mat <- array(loglik.vec, 
                      c((n.per.chain * MCMC.params$n.chains), n.data))
  
  # remove the log-likelihood values that correspond to NA data points
  loglik.mat <- loglik.mat[, !is.na(data.vec)]
  
  # Also, some of data points (0s) were unobserved and no log likelihood values
  # exist
  colsums.loglik <- colSums(loglik.mat)
  loglik.mat <- loglik.mat[, !is.na(colsums.loglik)]
  
  # remove NAs in the data vector
  #data.vec <- data.vec[!is.na(data.vec)]
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  loo.out <- rstanarm::loo(loglik.mat, 
                           r_eff = Reff, 
                           cores = 4, k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}


# Best and Min are subsets of Sightings. Sightings df has all the information
Best <- read.csv(file = "data/GROUPS_BEST.csv") %>%
   rename(Group_Size_Best = Group_Size)
Min <- read.csv(file = "data/GROUPS_MIN.csv")
Sightings <- read.csv(file = "data/SIGHTINGS_GROUPS.csv")

ggplot(Sightings) + 
  #geom_point(aes(x = Group_Size, y = GROUP_SIZE_LAST)) +
  geom_jitter(aes(x = Group_Size, y = GROUP_SIZE_LAST))

ggplot(Sightings) +
  geom_freqpoly(aes(x = DIST_MIN), bins = 20)

ggplot(Sightings) +
  geom_jitter(aes(x = DIST_MIN, y = GROUP_SIZE_LAST)) +
  geom_jitter(aes(x = DIST_MIN, y = Group_Size),
              color = "red", size = 2, alpha = 0.6)

MCMC.params <- list(n.samples = 10000,
                    n.thin = 50,
                    n.burnin = 5000,
                    n.chains = 5)

# MCMC.params <- list(n.samples = 2500,
#                     n.thin = 10,
#                     n.burnin = 200,
#                     n.chains = 5)

jags.params.v1 <- c("GS.UAS", "mu", "B0", "B1", "B2",
                    "B3", "B4", "log.lkhd")

jags.data.v1 <- list(n.obs = nrow(Sightings),
                     GS.Vis = Sightings$GROUP_SIZE_LAST,
                     GS.UAS = Sightings$Group_Size,
                     Dist = Sightings$DIST_MIN,
                     Bft = Sightings$BEAUFORT_LAST,
                     Vis = Sightings$VISIBILITY_LAST)

jags.params.v2 <- c("mu", "B0", "B1", "B2",
                    "B3", "B4", "GS.mu", "log.lkhd")

pred.data <- Sightings %>% filter(is.na(Group_Size))
jags.data.v2 <- list(n.obs = nrow(Best),
                     GS.Vis = Best$GROUP_SIZE_LAST,
                     GS.UAS = Best$Group_Size,
                     Dist = Best$DIST_MIN,
                     Bft = Best$BEAUFORT,
                     Vis = Best$VISIBILITY,
                     GS.pred = pred.data$GROUP_SIZE_LAST,
                     n.pred = nrow(pred.data),
                     Dist.pred = pred.data$DIST_MIN,
                     Bft.pred = pred.data$BEAUFORT_LAST,
                     Vis.pred = pred.data$VISIBILITY_LAST)

jm <- jags(data = jags.data.v2,
           parameters.to.save = jags.params.v2,
           model.file = "models/model_Pois_logMu_v2.txt",
           n.chains = MCMC.params$n.chains,
           n.burnin = MCMC.params$n.burnin,
           n.iter = MCMC.params$n.samples)

LOOIC <- compute.LOOIC(loglik.array = jm$sims.list$log.lkhd,
                         data.array = jags.data$GS.UAS,
                         MCMC.params = MCMC.params)

mcmc_trace(jm$samples,
           c("B0", "B1", "B2", "B3", "B4"))

GS.pred <- data.frame(Pred = jm$mean$mu.pred,
                      Vis = pred.data$GROUP_SIZE_LAST)

ggplot() + geom_point(data = Sightings,
                      aes(x = GROUP_SIZE_LAST,
                          y = Group_Size),
                      color = "green", size = 2) +
  geom_jitter(data = GS, aes(x = Vis, y = UAS), 
             color = "red", alpha = 0.3)
