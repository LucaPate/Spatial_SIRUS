# Guided example for running S-SIRUS
# Authors: L. Patelli, N. Golini, R. Ignaccolo, M. Cameletti
# Last update: 09/08/2024
# For details see the companion paper: Arxiv????


rm(list = ls())

# The following packages are required
library(RandomForestsGLS)
library(tidyverse) 
library(glmnet) 
library(BRISC) 
library(sirus)

# Set the working directory to your folder 
# from where all the following R files will be sourced

# Load the necessary functions ---------------------------------
source("s.sirus.R") # This script contains the main S-SIRUS functions
source("s.sirus_rfgls.R") 
source("s.sirus_formula.R")
source("s.sirus_RcppExports.R") 
source("s.sirus_utility.R")
#source("utility.R")


# Load the simulated data ---------------------------------
# See Scenario A in the companion paper
load("SimulatedData.RData")
str(data)

# Save the response variable in a separate vector
y = data$y
data$y = NULL

# In the following code the first 400 observations 
# are used as training data and the last 100 observations
# as test data. 
# Please, note this setting is quite LONG to run
# For a quick overview of the results,
# use a smaller number of observations 
# (e.g., 100 training data and 25 test data)

n_tr = 100 # number of training data
n_te = 25  # number of test data
seed = 1234

# S-SIRUS cross-validation --------------------------------------------------------------------------------
# This is required for tuning the p0 parameter
s.sirus.cv_output <- s.sirus.cv(data = data[1:n_tr,],
                                y = y[1:n_tr],
                                incl.coords = FALSE, 
                                type = "reg",
                                nfold = 2, # number of K CV folds
                                ncv = 2,   # number of CV repetitions
                                num.rule.max = 25,
                                q = 10, # for continuous predictors discretization
                                num.trees.step = 1000, # b: number of trees per step
                                alpha = 0.05, # 
                                max.depth = 2,
                                seed = seed)

# This function returns the 2 CV plots: 
# error and stability versus the number of rules when p0 varies.
s.sirus.plot.cv(s.sirus.cv_output) 

# S-SIRUS fit ----------------------------------------------------------------------------
# This is required for the fit of the S-SIRUS rules
s.sirus.fit.p0_output = s.sirus.fit(data = data[1:n_tr,],
                                    y = y[1:n_tr],
                                    incl.coords = FALSE,
                                    type = "reg",
                                    max.depth = 2,
                                    p0 = s.sirus.cv_output$p0.stab,
                                    q = 10,
                                    num.trees.step = 1000, # b: number of trees per step
                                    seed = seed)

# This function returns the S-SIRUS list of rules
# with the corresponding weights
s.sirus.print(s.sirus.fit.p0_output, digits = 3)


# S-SIRUS predictions --------------------------------------------------------------------
# This is required for predicting the response variable 
# in any new location using S-SIRUS-RK

# Large Scale (LS) prediction for the training data
# (coordinates are not provided as predictors) 
LS_hat_s.sirus_tr = s.sirus.predict(s.sirus.fit.p0_output,
                                    data[1:n_tr,-c(1:2)]) #omit coordinates
# Compute residuals
residuals_s.sirus_tr = y[1:n_tr] - LS_hat_s.sirus_tr

# Prepare the coordinates matrix separately for training and test data
coords_tr = as.matrix(data[1:n_tr, 1:2])
coords_te = as.matrix(data[(n_tr+1):(n_tr+n_te), 1:2])

# Estimation of the spatial parameters
# using the training data residuals and BRISC functions
set.seed(seed)
est_theta_tr_res = BRISC_estimation(coords_tr,
                                        x = matrix(1, nrow(coords_tr), 1), 
                                        y = residuals_s.sirus_tr,
                                        verbose = FALSE)
est_theta_tr_res$Theta
  
# Small scale prediction for the test data
# using BRISC functions
set.seed(seed)
pred_smallscale_te <- BRISC_prediction(est_theta_tr_res,
                                       coords_te,
                                       X.0 = matrix(1, nrow(coords_te), 1),
                                       verbose = FALSE)
length(pred_smallscale_te$prediction)

# Large Scale (LS) prediction for the test data
# (coordinates are not provided as predictors) 
LS_hat_s.sirus_te = s.sirus.predict(s.sirus.fit.p0_output,
                                    data[(n_tr+1):(n_tr+n_te),-c(1:2)])

# Response variable prediction for the test data
# as sum of the large scale and the small scale
y_hat_s.sirus_te = LS_hat_s.sirus_te + pred_smallscale_te$prediction


