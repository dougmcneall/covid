# GenerateSensDesign.R
# Generate a design for a simple sensitivity analysis of MetaWards.

# Choose some limits for incubation period, rzero, infectious period etc.

# From a MetaWards error message:
# Put the result through Ellen's code to generate a design matrix
# Available variables are ['beta', 'progress', 'too_ill_to_move', 'contrib_foi', 'user', 'length_day', 
# 'plength_day', 'UV', 'initial_inf', 'static_play_at_home', 'dyn_play_at_home', 'data_dist_cutoff', 
# ''dyn_dist_cutoff', 'play_to_work', 'work_to_play', 'local_vaccination_thesh', 'global_detection_thresh', 
# ''daily_ward_vaccination_capacity', 'neighbour_weight_threshold', 'daily_imports'],
# ' or to set a user parameter 'user.parameter' or '.parameter'. 
# ' To set an index use 'parameter[index]', e.g. 'beta[2]'"

library(lhs)
library(MASS)
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")


# helper functions for calculating the outputs we want
get_beta <- function(rzero,IP){
  beta <- rzero/IP
  beta
}

get_dt <- function(IP, incu){
  b <- -(1/IP + 1/incu)
  a <- (beta/incu - (1/incu)*(1/IP))
  c <- -1
  dt <- log(2)*(-b+sqrt((b^2-4*a*c)))/(2*a)
  dt
}

get_IP2 <- function(IP, IP1){
  IP2 <- IP - IP1
  IP2
}

exptdir <- "experiments/2020-05-12-sensitivity-analysis/"

# Number of ensemble members in the design
n = 50
varnames <- c('incu', 'IP', 'rzero')

lhs = maximinLHS(n = n, k = length(varnames), method = "build", dup = 1, eps = 0.05,
           maxIter = 100, optimize.on = "grid", debug = FALSE)

# minima and maxima for the parameters
param_mins  <- c(4, 2, 2.5)
param_maxes <- c(6, 4, 4)

# rescale the lhs to the parameter maxes and mins
design <- unnormalize(
  lhs,
  un.mins = param_mins,
  un.maxes = param_maxes
)

colnames(design) <- varnames
write.matrix(design, file = paste0(exptdir,'lhs.csv'))

# Calculate the model inputs (and extras)
beta <- get_beta(rzero = design[,'rzero'], IP = design[,'IP'])
dt <- get_dt(IP = design[,'IP'], incu = design[,'incu'])
rzero <- design[,'rzero']
IP <- design[,'IP']

IP1 <- rep(1,length(beta))
IP2 <- get_IP2(IP = IP, IP1 = IP1)

Progress2 = 1/design[,'incu']
Progress3 = 1/IP1
Progress4 = 1/IP2

parameters3 = cbind(beta, beta, Progress2, Progress3, Progress4, rzero, dt)
colnames(parameters3) <- c("beta[3]" , "beta[4]", "progress[2]", 'progress[3]', 'progress[4]', 'r0', 'dt')


write.matrix(as.matrix(parameters3), paste0(exptdir,"design.csv"))
write.matrix(as.matrix(parameters3[,1:5]), paste0(exptdir,"ncovparams.csv"))





