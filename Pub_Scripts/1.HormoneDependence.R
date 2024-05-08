#!/usr/bin/env Rscript
######################################################################
# Puppy Hormone Analyses: OT-F Dependence
# GEG
# Thu Mar  9 21:30:41 2023 ------------------------------ GEG
# Wed May  8 09:31:40 2024 ------------------------------ GEG
# To run: 
# Rscript ../Pub_Scripts/1.HormoneDependence.R 
# Model formula is specified below, not by command line arguments.
# Note: cannot run without any covariates, since prior is then not used and throws an error
#####################################################################
library(here)
library(brms)
library(rstan) # Needed to specify rstan_options
library(tidyverse)
#####################################################################
# Prepared file with data and a matrix
load(here("Pub_Data", "PuppyHormone_HeritabilityPrep_Anon.RData")) 
# Note: this version has two more individuals than "DCDBModelPrep", which requires cognitive data.

# BRMS parameters to tweak, if necessary.
tot <- 200000
thn <- 10
burn <- 2000
adelta <- 0.999
rstan_options(auto_write = TRUE)
seed_n <- 619

######################################################################
# Specify model formula, type, and matching model name for files.
mod_formula <- "Corr.F.log.z ~ 1 + Time.Collect.F.z + Corr.OT.log.z + (1|Plate.F) + (1|Plate.OT) + (1|p|gr(Name, cov = A))"
mod_name <- "F-OT_Dependence"
mod_type <- "gaussian"

# Create corresponding output directory
mod_dir <- here("Output_Heritability", paste(Sys.Date(), mod_name, sep = "_"))
dir.create(file.path(mod_dir), recursive = TRUE, showWarnings = FALSE)

######################################################################  
# Get everything ready
######################################################################
# Drop individuals without data on both analytes.
analytes <- c("Corr.F.log.z", "Corr.OT.log.z")
mod_data <- all_data %>% drop_na(all_of(analytes))
 
# Filter A-matrix to match. This also sorts in the same order, so Z can be an identity matrix.
a_mat_filt <- a_mat[mod_data$Name, mod_data$Name] 
# Note: despite inconsistencies in documentation, we do want A matrix (relatedness), not kinship.
rm(a_mat) # Remove unfiltered a_matrix to save space.

# Weakly informative prior, with mean 0 and stdev 1.
b_prior <- prior(normal(0,1), class = "b")
# Note: there are versioning incompatibilities with prior, so I've switched to the prior
# being specified in this script, instead of in the workspace prep.

######################################################################
# Run model
######################################################################
mod <- brm(mod_formula,
           data = mod_data, 
           data2 = list(A = a_mat_filt), 
           prior = b_prior,
           family = mod_type,
           chains = 4, cores = 4, thin = thn, iter = tot, warmup = burn, 
           silent = F, sample_prior = TRUE, seed = seed_n,
           control = list(adapt_delta = adelta))

summary(mod)
######################################################################
# Pull out data
######################################################################
# Function to get full distribution of heritabilities
get_heritability <- function(brms_model){
  tmp <- as_draws_df(mod, variable = c("sd_Name__Intercept","sigma"))
  tmp <- tmp[,1:2]
  colnames(tmp) <- c('sd_g', 'sd_e')
  tmp <- mutate(tmp, h2 = (sd_g^2) / ((sd_g^2) + (sd_e^2)))
  h2 <- tmp
  return(h2)
}

h2 <- get_heritability(mod)
print(paste0("Mean h2 = ", mean(h2$h2)))

png(here(mod_dir, paste0(mod_name, "_heritability_hist.png")), 
    width = 6, height = 4, units = "in", res = 300 )
hist(h2$h2, main = mod_formula)
dev.off()

# Save  h2 information as part of the model
mod$h2 <- h2
save(mod,	h2, 
     file = here(mod_dir, paste0(mod_name, "_brms_heritability.RData")))

# Function to read all other brms data
brms_sum <- function(brms_model){
  # Bayesian model
  brms_smry <- summary(brms_model)
  # Fixed effects
  brms_fixed <- brms_smry$fixed
  colnames(brms_fixed) <- make.names(colnames(brms_fixed))
  brms_fixed$par <- rownames(brms_fixed)
  # Sample size
  brms_fixed$n <- brms_smry$ngrps[[1]]
  # Number of iterations and divergences
  brms_fixed$iter <- brms_smry$iter
  n_divergences <- get_num_divergent(brms_model$fit)
  brms_fixed$divergences <- n_divergences
  brms_fixed$diverge.p <- n_divergences/brms_smry$iter
  return(brms_fixed)
}

mod_sum <- brms_sum(mod)
print(mod$formula)
print(paste0("% Divergences = ", mod_sum$diverge.p[1]))

######################################################################
# PP plot as diagnostic
######################################################################
check_pp_plot <- function(brms_model){
  name <- paste0("PP_", mod_name, ".png")
  pp_plot <- pp_check(brms_model, type = "dens_overlay", ndraws = 20) +
    ggtitle(brms_model$formula) +
    theme(axis.text = element_text(size = 16))
  ggsave(here(mod_dir, name), pp_plot, width = 12, height = 6, units = "in")
  return(pp_plot)
}
pp <- check_pp_plot(mod)

######################################################################
# Write both heritability and model summary stats
######################################################################
write_csv(h2, here(mod_dir, paste0(Sys.Date(), "_Puppy_h2_", mod_name, '.csv')))
write_csv(mod_sum, here(mod_dir, paste0(Sys.Date(), "_Puppy_modsum_", mod_name,'.csv')))
