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
library(ggridges)

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
Sightings.1 <- read.csv(file = "data/SIGHTINGS_GROUPS_with_ObsIDs.csv")

# Remove influential data points:
Sightings.1 %>%
  filter(GROUP_SIZE_LAST < 7) -> Sightings.1

uniq.obs.df <- data.frame(OBS_1 = unique(Sightings.1$OBS_1)) %>%
  mutate(Obs.ID = 1:length(OBS_1))

Sightings.1 %>% 
  left_join(uniq.obs.df, by = "OBS_1") -> Sightings

# Create a new column of minimum group sizes
Sightings %>%
  mutate(Group_Size_Min_Vis = GROUP_SIZE_LAST - 2) -> Sightings

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
#min.group.size <- Sightings$Group_Size_Min
min.group.size.Vis <- Sightings$Group_Size_Min_Vis

GS.I <- vector(mode = "numeric", length = length(min.group.size.Vis))
GS.I[!is.na(Sightings$Group_Size)] <- 1

# Adjust the minimum group size. Minimum group size should be less than
# the observed 
#min.group.size[is.na(min.group.size)] <- 1
min.group.size.Vis[min.group.size.Vis < 1] <- 1

min.group.size.UAS <- Sightings$Group_Size_Min
min.group.size.UAS[is.na(min.group.size.UAS)] <- min.group.size.Vis[is.na(min.group.size.UAS)]

min.group.size.df <- data.frame(UAS = min.group.size.UAS,
                                Vis = Sightings$GROUP_SIZE_LAST)

min.group.size <- apply(min.group.size.df, FUN = max, MARGIN = 1, na.rm = T)

# group size categories
GS.cat <- vector(mode = "numeric", length = length(min.group.size))
GS.cat[Sightings$GROUP_SIZE_LAST < 3] <- 1
GS.cat[Sightings$GROUP_SIZE_LAST > 2 & Sightings$GROUP_SIZE_LAST < 5] <- 2
GS.cat[Sightings$GROUP_SIZE_LAST > 4] <- 3

GS.min <- min.group.size.UAS
GS.min[which((min.group.size.UAS - min.group.size) < 0)] <- Sightings$GROUP_SIZE_LAST[which((min.group.size.UAS - min.group.size) < 0)]

jags.data <- list(n.grp = nrow(Sightings),
                  GS.Vis = Sightings$GROUP_SIZE_LAST,
                  GS.UAS = Sightings$Group_Size,
                  Dist = Sightings$DIST_MIN,
                  Bft = Sightings$BEAUFORT_LAST,
                  Vis = Sightings$VISIBILITY_LAST,
                  obs = Sightings$Obs.ID,
                  GS.min.Vis = min.group.size.Vis,
                  GS.min.UAS = min.group.size.UAS,
                  GS.min = min.group.size,
                  GS.I = GS.I,
                  GS.max = max(Sightings$Group_Size_Min, na.rm = T) + 5,
                  GS.cat = GS.cat,
                  n.obs = nrow(uniq.obs.df))

# -1 models don't seem to work so well... 2025-04-30
#models<-c("v1", "v1-1", "v2", "v2-1", "v3", "v3-1", "v4", "v5", "v6")
model.ver <- "v10"
model.file <- switch(model.ver,
                     "v1" = "models/model_Pois_logMu_v1.txt",
                     "v1-1" = "models/model_Pois_logMu_v1-1.txt",
                     "v1-2" = "models/model_Pois_logMu_v1-2.txt",
                     "v2" = "models/model_Pois_Pois_logMu_v2.txt",  
                     "v2-1" = "models/model_Pois_Pois_logMu_v2-1.txt",
                     "v3" = "models/model_Pois_Gam_logMu_v3.txt",
                     "v3-1"= "models/model_Pois_Gam_logMu_v3-1.txt",
                     "v4" = "models/model_Pois_Pois_logitP_v4.txt",
                     "v5" = "models/model_Pois_Gam_logitP_v5.txt",
                     "v6" = "models/model_Pois_Gam_logitP_v6.txt",
                     "v7" = "models/model_Binom_Gam_logitP_v7.txt",
                     "v8" =  "models/model_Pois_pois_logitP_v8.txt",
                     "v9" = "models/model_Gam_Gam_logMu_v9.txt",
                     "v10" = "models/model_Binom_Pois_logitP_v10.txt")

jags.params <- switch(model.ver,
                      "v1" = c("GS.UAS", "B0", "B1", "B2",
                               "B3", "B4", "mu.Vis", "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v1-1" = c("GS.UAS", "B0", "B1", "B2", "B3", "B4", "B5", 
                                 "mu.Vis", "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v1-2" = c("GS.UAS", "B0", "B1",
                                 "mu.Vis", "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v2" = c("B0", "B1", "B2", "B3", "B4", "mu.Vis",
                               "GS.UAS", "mu.UAS","Obs.RF", "sigma.Obs", "log.lkhd"),
                      "v2-1" = c("B0", "B1", "B2", "B3", "B4", "B5", "mu.Vis",
                                 "GS.UAS", "mu.UAS", "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v3" = c("B0", "B1", "B2", "B3", "B4", "mu.Vis",
                               "GS.UAS", "alpha.UAS", "beta.UAS",
                               "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v3-1" = c("B0", "B1", "B2", "B3", "B4", "B5", "mu.Vis",
                                 "GS.UAS", "alpha.UAS", "beta.UAS",
                                 "Obs.RF", "sigma.Obs","log.lkhd"),
                      
                      "v4" = c("B0", "B1", "B2", "B3", "p.Vis", "mu.Vis",
                               "GS.UAS", "mu.GS", "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v5" = c("B0", "B1", "B2", "B3", "p.Vis", "mu.Vis",
                               "GS.UAS", "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v6" = c("B0", "B1", "B2", "B3", "B4", "p.Vis", "mu.Vis",
                               "GS.", "GS.UAS", "GS.Vis", 
                               "Obs.RF", "sigma.Obs", "log.lkhd"),
                      "v7" = c("B0", "B1", "B2", "B3", "B4", "p.Vis", "B0.uas",
                               "GS.", "GS.UAS", "GS.Vis", "p.UAS",
                               "Obs.RF", "alpha.", "beta.", "sigma.Obs", "log.lkhd"),
                      "v8" = c("B0", "B1", "B2", "B3", "B4", "p.Vis", "p.UAS",
                               "mu.Vis",
                               "B0.uas", #"B1.uas", "B2.uas", "B3.uas",
                               "GS.", "GS.UAS", "GS.Vis", "GS.mean",
                               "Obs.RF", "sigma.Obs", "log.lkhd"),
                      "v9" = c("B0", "B1", "B2", "B3", "B4", "mu.Vis",
                               "GS.UAS", 
                               "Obs.RF", "sigma.Obs","log.lkhd"),
                      "v10" = c("B0", "B1", "B2", "B3", "B4", "p.Vis", "B0.uas",
                                "GS.", "GS.UAS", "GS.Vis", "p.UAS",
                                "Obs.RF", "GS.mean", "sigma.Obs", "log.lkhd"))

out.file.name <- paste0("RData/jm_out_", model.ver, "_max6.rds")
if (!file.exists(out.file.name)){
  jm <- jags(data = jags.data,
             parameters.to.save = jags.params,
             model.file = model.file,
             n.chains = MCMC.params$n.chains,
             n.burnin = MCMC.params$n.burnin,
             n.iter = MCMC.params$n.samples,
             parallel = T,
             DIC = T)  
  
  jm.out <- list(jm = jm,
                 jags.data = jags.data,
                 jags.params = jags.params,
                 MCMC.params = MCMC.params)
  saveRDS(jm.out, file = out.file.name)
  
  
} else {
  jm.out <- readRDS(out.file.name)  
  jm <- jm.out$jm
}

LOOIC <- compute.LOOIC(loglik.array = jm$sims.list$log.lkhd,
                       data.array = jags.data$GS.UAS,
                       MCMC.params = MCMC.params)

max.Rhat <- lapply(jm$Rhat, FUN = max, na.rm = T) 
max.max.Rhat <- max(unlist(max.Rhat))

params.to.plot <- switch(model.ver,
                         "v1" = c("B0", "B1", "B2", "B3", "B4", "sigma.Obs"),
                         "v1-1" = c("B0", "B1", "B2", "B3", "B4", "B5", "sigma.Obs"),
                         "v1-2" = c("B0", "B1",  "sigma.Obs"),
                         "v2" = c("B0", "B1", "B2", "B3", "B4", "mu.UAS", "sigma.Obs"),
                         "v2-1" = c("B0", "B1", "B2", "B3", "B4", "B5", "mu.UAS", "sigma.Obs"),
                         "v3" = c("B0", "B1", "B2", "B3", "B4", 
                                  "alpha.UAS", "beta.UAS", "sigma.Obs"),
                         "v3-1"= c("B0", "B1", "B2", "B3", "B4", "B5", 
                                   "alpha.UAS", "beta.UAS", "sigma.Obs"),
                         "v4" = c("B0", "B1", "B2", "B3", "mu.GS", "sigma.Obs"),
                         "v5" = c("B0", "B1", "B2", "B3", "sigma.Obs"),
                         "v6" = c("B0", "B1", "B2", "B3", "sigma.Obs"),
                         "v7" = c("B0", "B1", "B2", "B3", "B4", "B0.uas", "sigma.Obs",
                                  "alpha.", "beta."),
                         "v8" = c("B0", "B1", "B2", "B3", "GS.mean",
                                  "sigma.Obs"),
                         "v9" = c("B0", "B1", "B2", "B3", "B4", 
                                  "sigma.Obs"),
                         "v10" = c("B0", "B1", "B2", "B3", "B4", "B0.uas", "sigma.Obs",
                                   "GS.mean"))

# The following two lines don't run for "v1-1" for some reason... 
mcmc_trace(jm$samples,
           params.to.plot)

mcmc_dens(jm$samples,
          params.to.plot)



if (model.ver == "v3" | model.ver == "v3-1"){
  GS.UAS.df <- data.frame(x = seq(0, 15, by = 0.01)) %>%
    mutate(y = dgamma(x, jm$mean$alpha.UAS, jm$mean$beta.UAS))
  
  ggplot(GS.UAS.df) + 
    geom_line(aes(x = x, y = y)) +
    labs(x = "Group size", y = "Density")
} else if (model.ver == "v6"){
  # GAM(2.952, 0.760) comes from fitting the gamma to UAS group sizes
  alpha. <- 2.952 #1.5
  beta. <- 0.760 #0.5
  GS.df <- data.frame(x = seq(0, 15, by = 0.01)) %>%
    mutate(y = dgamma(x, alpha.,beta.))
  
  ggplot(GS.df) + 
    geom_line(aes(x = x, y = y)) +
    labs(x = "Group size", y = "Density")
} else if (model.ver == "v7"){
  beta. <- 1.0
  GS.df <- data.frame(x = seq(0, 15, by = 0.01)) %>%
    mutate(y = dgamma(x, jm$mean$alpha., jm$mean$beta.))
  
  ggplot(GS.df) + 
    geom_line(aes(x = x, y = y)) +
    labs(x = "Group size", y = "Density")
  
} else if (model.ver == "v10"){
  GS.df <- data.frame(x = seq(0, 15, by = 1)) %>%
    mutate(y = dpois(x, jm$mean$GS.mean))
  
  ggplot(GS.df) + 
    geom_line(aes(x = x, y = y)) +
    labs(x = "Group size", y = "Density")
  
}

if (model.ver == "v6" | model.ver == "v7" | model.ver == "v8" | model.ver == "v10"){
  GS.pred <- data.frame(GS.mean = jm$mean$GS.,
                        GS.low = jm$q2.5$GS.,
                        GS.high = jm$q97.5$GS.)
  
}

GS.UAS.pred <- data.frame(GS.UAS.mean = jm$mean$GS.UAS,
                          GS.UAS.low = jm$q2.5$GS.UAS,
                          GS.UAS.high = jm$q97.5$GS.UAS)



if (model.ver == "v6" | model.ver == "v7" | model.ver == "v8" | model.ver == "v10"){
  
  Sightings.pred <- cbind(Sightings, GS.UAS.pred, GS.pred)
  p.pred <- ggplot(Sightings.pred) + 
    # geom_point(aes(x = GROUP_SIZE_LAST,
    #                y = GS.mean),
    #            color = "blue", size = 2) +
    geom_errorbar(aes(x = GROUP_SIZE_LAST,
                      ymin = GS.low,
                      ymax = GS.high),
                  color = "blue", alpha = 0.5) +
    geom_jitter(aes(x = GROUP_SIZE_LAST,
                    y = GS.mean),
                color = "blue", size = 2,
                height = 0) +
    geom_jitter(aes(x = GROUP_SIZE_LAST, 
                    y = Group_Size),
                color = "red", alpha = 0.5,
                height = 0) +
    geom_rug(aes(x = GROUP_SIZE_LAST)) +
    geom_abline(slope = 1.0)
  ggsave(p.pred,
         filename = paste0("figures/GroupSizePredictions_", model.ver, "_max6.png"),
         device = "png", dpi = 600)
  
} else {
  GS.Vis.mean <- data.frame(GS.Vis.mean = jm$mean$mu.Vis,
                            GS.Vis.low = jm$q2.5$mu.Vis,
                            GS.Vis.high = jm$q97.5$mu.Vis)
  
  Sightings.pred <- cbind(Sightings, GS.UAS.pred, GS.Vis.mean)
  p.pred.Vis <- ggplot(Sightings.pred) + 
    # geom_point(aes(x = GROUP_SIZE_LAST,
    #                y = GS.UAS.mean),
    #            color = "blue", size = 2) +
    geom_errorbar(aes(x = GROUP_SIZE_LAST,
                      ymin = GS.Vis.low,
                      ymax = GS.Vis.high),
                  color = "blue", alpha = 0.5) +
    geom_jitter(aes(x = GROUP_SIZE_LAST,
                    y = GS.Vis.mean),
                color = "blue", size = 2,
                height = 0) +
    geom_jitter(aes(x = GROUP_SIZE_LAST, 
                    y = Group_Size),
                color = "red", alpha = 0.5,
                height = 0) +
    geom_rug(aes(x = GROUP_SIZE_LAST)) +
    geom_abline(slope = 1.0)
  
  ggsave(p.pred.Vis,
         filename = paste0("figures/Vis_GroupSizePredictions_", model.ver, "_max6.png"),
         device = "png", dpi = 600)
}

p.pred.UAS <- ggplot(Sightings.pred) + 
  # geom_point(aes(x = GROUP_SIZE_LAST,
  #                y = GS.UAS.mean),
  #            color = "blue", size = 2) +
  geom_errorbar(aes(x = GROUP_SIZE_LAST,
                    ymin = GS.UAS.low,
                    ymax = GS.UAS.high),
                color = "blue", alpha = 0.5) +
  geom_jitter(aes(x = GROUP_SIZE_LAST,
                  y = GS.UAS.mean),
              color = "blue", size = 2,
              height = 0) +
  geom_jitter(aes(x = GROUP_SIZE_LAST, 
                  y = Group_Size),
              color = "red", alpha = 0.5,
              height = 0) +
  geom_rug(aes(x = GROUP_SIZE_LAST)) +
  geom_abline(slope = 1.0)

ggsave(p.pred.UAS,
       filename = paste0("figures/UAS_GroupSizePredictions_", model.ver, "_max6.png"),
       device = "png", dpi = 600)


plot.ridges.2 <- function(x, levels = c("1", "2", "3", "4", 
                                        "5", "6", "7", "8", 
                                        "9", "12"),
                          bandwidth = 0.12){
  jm.out <- readRDS(x)
  
  n.UAS <- jm.out$jm$sims.list$GS.
  n.Vis <- jm.out$jags.data$GS.Vis
  
  Sightings. <- data.frame(UAS = jm.out$jags.data$GS.UAS,
                           Vis = jm.out$jags.data$GS.Vis) %>%
    na.omit() %>%
    mutate(Vis.f = factor(Vis, levels = levels))
  
  n. <- data.frame(UAS = as.vector(n.UAS),
                   Vis = rep(n.Vis, each = nrow(n.UAS))) %>%
    mutate(Vis.f = factor(Vis, levels = levels))
  
  p.ridges <- ggplot() +
    # geom_density_ridges(data = Sightings.2,
    #                     aes(x = UAS, y = Vis.f)) +
    geom_density_ridges(data = n.,
                        aes(x = UAS, y = Vis.f, fill = Vis.f),
                        bandwidth = bandwidth) +
    geom_jitter(data = Sightings.,
                aes(x = UAS, y = Vis.f, color = Vis.f),
                width = 0) +
    geom_abline(slope = 1.0) +
    theme(legend.position = "none") +
    ylab("Visual group size") +
    xlab("Group size")
  
  out <- list(plot = p.ridges,
              sightings = Sightings.,
              n = n.)  
}

out. <- plot.ridges.2(out.file.name,
                      levels = c("1", "2", "3", "4", 
                                          "5", "6"),
                      bandwidth = 0.3)

# # Extract posterior samples from jagsUI output
# # zm is the object that comes back from jags in jagsUI
# # Do not include the index, e.g., [1], [2], etc.
# extract.samples.jagsUI <- function(varname, jm){
#   par.names <- unlist(dimnames(jm$samples[[1]])[2])
#   col.idx <- grep(varname, par.names)    
#   samples.list <- list()
#   
#   samples <- lapply(jm$samples, FUN = function(x) x[, col.idx])
#   # for (k in 1:length(col.idx)){
#   #   samples.list[[k]] <- unlist(lapply(samples, FUN = function(x) x[,k]))
#   # }
#   
#   return(samples)
# }
# 
# if (model.ver == "v6"){
#   GS.samples.list <- extract.samples.jagsUI("GS.\\[", jm)  
# } else { #if (model.ver == "v1" | model.ver == "v1-1"){
#   GS.samples.list <- extract.samples.jagsUI("GS.UAS\\[", jm)  
# }
# 
# k <- 7
# GS.Vis.samples <- list()
# GS.plots <- list()
# GS.median <- vector(mode = "numeric")
# for (k in 1:9){
#   GS.Vis.samples[[k]] <- lapply(GS.samples.list, 
#                                 FUN = function(x){
#                                   x[, which(jags.data$GS.Vis == k)]}) %>% unlist()
#   
#   #tmp <- GS.Vis.samples[[k]] %>% unlist()  
#   
#   GS.samples <- data.frame(GS = GS.Vis.samples[[k]])
#   GS.median[k] <- median(GS.Vis.samples[[k]])
#   GS.plots[[k]] <- ggplot(GS.samples) + 
#     geom_histogram(aes(x = GS),
#                    binwidth = 1) +
#     labs(title = paste0("Visual group size = ", k))
# }
# 
# GS.median.df <- data.frame(Vis.GS = c(1:9),
#                            GS.median = GS.median)
