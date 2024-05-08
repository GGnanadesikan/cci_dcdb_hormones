#!/usr/bin/env Rscript
######################################################################
# Puppy Hormone-Cognition Analyses
# GEG
# Thu Mar  9 21:30:41 2023 ------------------------------ GEG
# Wed May  8 09:54:14 2024 ------------------------------ GEG
# To run (oxytocin):
# Rscript ../Pub_Scripts/3.PuppyDCDBHormones_Quad.R ""Corr.OT.log.z" "OT_Quadratic"
# To run (cortisol): 
# Rscript ../Pub_Scripts/3.PuppyDCDBHormones_Quad.R "Corr.F.log.z" "F_Quadratic"
# This is designed to be run on a high performance computing cluster using an array job
# and four nodes per job. Comments below should help adjust to run locally and interactively.
#####################################################################
library(tidyverse)
library(here)
library(brms)
library(rstan)
#####################################################################
# Read in prepared file with data, prior & a matrix
load(here("Pub_Data", "PuppyHormone_DCDBModelPrep_Anon.RData")) 

# Get some model info from the call (analyte and a name for files)
args <- commandArgs(trailingOnly = TRUE) 
# Or if running interactively, rather than from a command line, specify either:
# args <- c("Corr.F.log.z", "F_Quadratic")
# args <- c("Corr.OT.log.z", "F_Quadratic")
print(args)
analyte <- args[1] # OT or F, and which transformation, i.e. "Corr.OT.log.z"
mod_name <- args[2] # Use this to make file names.

# Create output directory
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
# Set n of trait based on the array number, or set as any number 1:20
trait_n <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Or can set trait manually here by replacing with the column name.
trait <- phenotypes_focal[trait_n]
print(trait) # Print to check trait.

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

# Drop individuals without data on this trait or analyte.
mod_data <- all_data %>% drop_na(all_of(trait), all_of(analyte))
 
# Filter A-matrix to match. This also sorts in the same order, so Z can be an identity matrix.
a_mat_filt <- a_mat[mod_data$Name, mod_data$Name] 
# Note: despite inconsistencies in documentation, we do want A matrix (relatedness), not kinship.
rm(a_mat) # Remove unfiltered a_matrix to save space.

# Models for OT and F are a little different, since F depends on collect time, but OT doesn't.
if (grepl("OT", analyte)){
  Plate <- " + (1|Plate.OT)"
  mod_formula <- paste0(trait, " ~ 1 + sex + plab.fac + experimenter_code + puptestage.z + whelp.location +", 
                        analyte, "+ I(", analyte, "^2)", Plate, " + (1|p|gr(Name, cov = A))")
  } 
if (grepl("F", analyte)){
  Plate <- "+ (1|Plate.F)"
  mod_formula <- paste0(trait, " ~ 1 + sex + plab.fac + experimenter_code + puptestage.z + whelp.location + Time.Collect.F.z + ", 
                        analyte, "+ I(", analyte, "^2)", Plate, " + (1|p|gr(Name, cov = A))")
  }

print(mod_formula) # Print the model formula to check.

# Weakly informative prior, with mean 0 and stdev 1.
b_prior <- prior(normal(0,1), class = "b")
# Note: there are versioning incompatibilities with prior, so I've switched to the prior
# being specified in this script, instead of in the workspace prep.

######################################################################
# Run model
######################################################################
trait_mod_bayes <- brm(mod_formula,
                       data = mod_data, 
                       data2 = list(A = a_mat_filt), 
                       prior = b_prior,
                       family = mod_type,
                       chains = 4, cores = 4, thin = thn, iter = tot, warmup = burn, 
                       silent = F, sample_prior = TRUE, seed = seed_n, 
                       control = list(adapt_delta = adelta))
    
summary(trait_mod_bayes)  
# Save trait and locus information as part of the model
trait_mod_bayes$trait <- trait
trait_mod_bayes$analyte <- analyte
save(trait_mod_bayes,	file = here(out_path, paste(analyte, trait, "brms.RData", sep = "-")))
