#!/usr/bin/env Rscript
######################################################################
# Puppy Hormone-Cognition Analyses (Additive)
# GEG
# Mon May  1 20:27:21 2023 ------------------------------ GEG
# Wed May  8 10:00:58 2024 ------------------------------ GEG
# To run:
# Rscript ../Pub_Scripts/4.PuppyDCDBHormones_QAdditive.R 
# The model itself is specified below.
#####################################################################
library(tidyverse)
library(here)
library(brms)
library(rstan)
#####################################################################
# Load prepared file with data but also prior, a matrix, etc.
load(here("Pub_Data", "PuppyHormone_DCDBModelPrep_Anon.RData")) 

# Specify analyte columns, model name, and create output directory.
analytes <- c("Corr.OT.log.z", "Corr.F.log.z")
mod_name <- "Additive"
out_path <- here("Output_Cognitive", paste(Sys.Date(), mod_name, sep = "_"))
dir.create(out_path)

######################################################################
# BRMS parameters to tweak, if necessary.
tot <- 100000
thn <- 10
burn <- 1000
adelta <- 0.999
rstan_options(auto_write = TRUE)
seed_n <- 719

######################################################################
# Set trait and prepare data and model
######################################################################
# Set n of trait based on the array number 
trait_n <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
trait <- phenotypes_focal[trait_n]
print(trait)

# Determine gaussian or bernoulli based on variable naming
if(grepl(".rn", trait)|grepl(".z", trait)){
  mod_type <- "gaussian"} else{
    if(grepl(".bin", trait)){
      mod_type <- "bernoulli"
    } else {
      print("ERROR in model specification")
      quit(save = "no")
    }
  }
print(mod_type)

# Drop individuals without data on this trait.
mod_data <- all_data %>% drop_na(all_of(trait), all_of(analytes))
 
# Filter A-matrix to match. This also sorts in the same order, so Z can be an identity matrix.
a_mat_filt <- a_mat[mod_data$Name, mod_data$Name] 
# Note: despite inconsistencies in documentation, we do want A matrix (relatedness), not kinship.
rm(a_mat) # Remove unfiltered a_matrix to save space.

# Create the full model formula
Plate <- " + (1|Plate.OT) + (1|Plate.F)"
mod_formula <- paste0(trait, " ~ 1 + sex + plab.fac + experimenter_code + puptestage.z + whelp.location + Time.Collect.F.z + ", 
                      analytes[1], "+ I(", analytes[1], "^2)", " + ", 
                      analytes[2], "+ I(", analytes[2], "^2)",
                      Plate, " + (1|p|gr(Name, cov = A))")
print(mod_formula)

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
    
  
# Save trait and locus information as part of the model
mod$trait <- trait
save(mod,	file = here(out_path, paste(mod_name, trait, "brms.RData", sep = "_")))
