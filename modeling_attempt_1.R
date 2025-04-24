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

rm(list = ls())
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

# ggplot(Sightings) + 
#   #geom_point(aes(x = Group_Size, y = GROUP_SIZE_LAST)) +
#   geom_jitter(aes(x = Group_Size, y = GROUP_SIZE_LAST))
# 
# ggplot(Sightings) +
#   geom_freqpoly(aes(x = DIST_MIN), bins = 20)
# 
# ggplot(Sightings) +
#   geom_jitter(aes(x = DIST_MIN, y = GROUP_SIZE_LAST)) +
#   geom_jitter(aes(x = DIST_MIN, y = Group_Size),
#               color = "red", size = 2, alpha = 0.6)

MCMC.params <- list(n.samples = 10000,
                    n.thin = 50,
                    n.burnin = 5000,
                    n.chains = 5)

# MCMC.params <- list(n.samples = 2500,
#                     n.thin = 10,
#                     n.burnin = 200,
#                     n.chains = 5)

#Sightings <- filter(Sightings, GROUP_SIZE_LAST < 7)
# For testing - keep only sightings with both visual and UAS data
# Sightings %>%
#   filter(!is.na(Group_Size)) -> Sightings

# Index for whether or not there is a UAS-based group size
min.group.size <- Sightings$Group_Size_Min
GS.I <- vector(mode = "numeric", length = length(min.group.size))
GS.I[!is.na(Sightings$Group_Size)] <- 1
min.group.size[is.na(min.group.size)] <- 1

jags.data <- list(n.obs = nrow(Sightings),
                  GS.Vis = Sightings$GROUP_SIZE_LAST,
                  GS.UAS = Sightings$Group_Size,
                  Dist = Sightings$DIST_MIN,
                  Bft = Sightings$BEAUFORT_LAST,
                  Vis = Sightings$VISIBILITY_LAST,
                  GS.min = min.group.size,
                  GS.I = GS.I)

model.ver <- "v6"
model.file <- switch(model.ver,
                     "v1" = "models/model_Pois_logMu_v1.txt",
                     "v1-1" = "models/model_Pois_logMu_v1-1.txt",
                     "v2" = "models/model_Pois_Pois_logMu_v2.txt",
                     "v2-1" = "models/model_Pois_Pois_logMu_v2-1.txt",
                     "v3" = "models/model_Pois_Gam_logMu_v3.txt",
                     "v3-1"= "models/model_Pois_Gam_logMu_v3-1.txt",
                     "v4" = "models/model_Pois_Pois_logitP_v4.txt",
                     "v5" = "models/model_Pois_Gam_logitP_v5.txt",
                     "v6" = "models/model_Pois_Gam_logitP_v6.txt")

jags.params <- switch(model.ver,
       "v1" = c("GS.UAS", "mu", "B0", "B1", "B2",
                "B3", "B4", "log.lkhd"),
       "v1-1" = c("GS.UAS", "mu", "B0", "B1", "B2",
                "B3", "B4", "B5", "log.lkhd"),
       "v2" = c("B0", "B1", "B2", "B3", "B4", 
                "GS.UAS", "mu.UAS", "log.lkhd"),
       "v2-1" = c("B0", "B1", "B2", "B3", "B4", "B5",
                "GS.UAS", "mu.UAS", "log.lkhd"),
       "v3" = c("B0", "B1", "B2", "B3", "B4", 
                "GS.UAS", "alpha.UAS", "beta.UAS",
                "log.lkhd"),
       "v3-1" = c("B0", "B1", "B2", "B3", "B4", "B5",
                "GS.UAS", "alpha.UAS", "beta.UAS",
                "log.lkhd"),
       
       "v4" = c("B0", "B1", "B2", "B3", "p.Vis",
                "GS.UAS", "mu.GS", "log.lkhd"),
       "v5" = c("B0", "B1", "B2", "B3", "p.Vis",
                "GS.UAS", "alpha.UAS", "beta.UAS",
                "log.lkhd"),
       "v6" = c("B0", "B1", "B2", "B3", "p.Vis", "p.UAS",
                "B0.uas", "B1.uas", "B2.uas", "B3.uas",
                "GS.", "GS.UAS", "GS.Vis", "alpha.", "beta.",
                "log.lkhd"))

jm <- jags(data = jags.data,
           parameters.to.save = jags.params,
           model.file = model.file,
           n.chains = MCMC.params$n.chains,
           n.burnin = MCMC.params$n.burnin,
           n.iter = MCMC.params$n.samples,
           parallel = T,
           DIC = T)

LOOIC <- compute.LOOIC(loglik.array = jm$sims.list$log.lkhd,
                       data.array = jags.data$GS.UAS,
                       MCMC.params = MCMC.params)

Rmax <- lapply(jm$Rhat, FUN = max, na.rm = T)

params.to.plot <- switch(model.ver,
                         "v1" = c("B0", "B1", "B2", "B3", "B4"),
                         "v1-1" = c("B0", "B1", "B2", "B3", "B4", "B5"),
                         "v2" = c("B0", "B1", "B2", "B3", "B4", "mu.UAS"),
                         "v2-1" = c("B0", "B1", "B2", "B3", "B4", "B5", "mu.UAS"),
                         "v3" = c("B0", "B1", "B2", "B3", "B4", 
                                  "alpha.UAS", "beta.UAS"),
                         "v3-1"= c("B0", "B1", "B2", "B3", "B4", "B5", 
                                   "alpha.UAS", "beta.UAS"),
                         "v4" = c("B0", "B1", "B2", "B3", "mu.GS"),
                         "v5" = c("B0", "B1", "B2", "B3", "alpha.UAS", "beta.UAS"),
                         "v6" = c("B0", "B1", "B2", "B3", "B0.uas", 
                                  "B1.uas", "B2.uas", "B3.uas", "alpha.", "beta."))

mcmc_trace(jm$samples,
           params.to.plot)

mcmc_dens(jm$samples,
           params.to.plot)

if (model.ver == "v5" | model.ver == "v3" | model.ver == "v3-1"){
  GS.UAS.df <- data.frame(x = seq(0, 15, by = 0.01)) %>%
  mutate(y = dgamma(x, jm$mean$alpha.UAS, jm$mean$beta.UAS))

  ggplot(GS.UAS.df) + 
    geom_line(aes(x = x, y = y)) +
    labs(x = "Group size", y = "Density")
} elseif (mode.ver == "v6"){
  GS.df <- data.frame(x = seq(0, 15, by = 0.01)) %>%
    mutate(y = dgamma(x, jm$mean$alpha., jm$mean$beta.))
  
  ggplot(GS.df) + 
    geom_line(aes(x = x, y = y)) +
    labs(x = "Group size", y = "Density")
}


if (model.ver == "v6"){
  GS.UAS.pred <- data.frame(GS.mean = jm$mean$GS.,
                            GS.low = jm$q2.5$GS.,
                            GS.high = jm$q97.5$GS.,
                            p.Vis.mean = jm$mean$p.Vis)
  
} else{
  GS.UAS.pred <- data.frame(GS.mean = jm$mean$GS.UAS,
                            GS.low = jm$q2.5$GS.UAS,
                            GS.high = jm$q97.5$GS.UAS,
                            p.Vis.mean = jm$mean$p.Vis)
  
}
# for v6 testing

Sightings.pred <- cbind(Sightings, GS.UAS.pred)

ggplot(Sightings.pred) + 
  geom_point(aes(x = GROUP_SIZE_LAST,
                 y = GS.mean),
             color = "blue", size = 2) +
  geom_errorbar(aes(x = GROUP_SIZE_LAST,
                    ymin = GS.low,
                    ymax = GS.high),
                color = "blue", alpha = 0.5) +
  geom_jitter(aes(x = GROUP_SIZE_LAST,
                 y = GS.mean),
             color = "blue", size = 2) +
  geom_jitter(aes(x = GROUP_SIZE_LAST, y = Group_Size),
             color = "red", alpha = 0.5,
             height = 0) +
  geom_rug(aes(x = GROUP_SIZE_LAST)) +
  geom_abline(slope = 1.0)


ggsave(filename = paste0("figures/GroupSizePredictions_", model.ver, ".png"),
         device = "png", dpi = 600)
