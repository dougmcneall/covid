# GenerateSensDesign.R

# Choose some limits for incubation period, rzero, infectious period etc.

# Put the result through Ellen's code to generate a design matrix

library(lhs)
library(MASS)
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")

# Suggested ranges

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

exptdir <- "experiments/2020-05-12-sensitivity-analysis/"

n = 3
varnames <- c('incu', 'IP', 'rzero')

lhs = maximinLHS(n = n, k = length(varnames), method = "build", dup = 1, eps = 0.05,
           maxIter = 100, optimize.on = "grid", debug = FALSE)


param_mins  <- c(4, 2, 2)
param_maxes <- c(6, 4, 4)

design <- unnormalize(
  lhs,
  un.mins = param_mins,
  un.maxes = param_maxes
)


colnames(design) <- varnames
#write.matrix(design, file = 'SA_design.csv')

#head(design)

beta <- get_beta(rzero = design[,'rzero'], IP = design[,'IP'])
dt <- get_dt(IP = design[,'IP'], incu = design[,'incu'])
rzero <- design[,'rzero']
IP <- design[,'IP']

IP1 <- rep(1,length(beta))
IP2 <- IP - IP1

Progress2 = 1/design[,'incu']
Progress3 = 1/IP1
Progress4 = 1/IP2
## "It is not possible to adjust the variable beta3 to equal 0.598. 
#Available variables are ['beta', 'progress', 'too_ill_to_move', 'contrib_foi', 'user', 'length_day', 
#'plength_day', 'UV', 'initial_inf', 'static_play_at_home', 'dyn_play_at_home', 'data_dist_cutoff', 
#''dyn_dist_cutoff', 'play_to_work', 'work_to_play', 'local_vaccination_thesh', 'global_detection_thresh', 
#''daily_ward_vaccination_capacity', 'neighbour_weight_threshold', 'daily_imports'],
#' or to set a user parameter 'user.parameter' or '.parameter'. 
#' To set an index use 'parameter[index]', e.g. 'beta[2]'"

parameters3 = cbind(beta, beta, Progress2, Progress3, Progress4, rzero, dt)

colnames(parameters3) <- c("beta[3]" , "beta[4]", "progress[2]", 'progress[3]', 'progress[4]', 'r0', 'dt')

#parameters3  <-  data.frame("beta[3]" = round(beta,3),
#                           "beta[4]" = round(beta,3),
#                           "progress[2]" = round(Progress2,3),
#                           "progress[3]"=round(Progress3,3),
#                           "progress[4]"=round(Progress4,3),
#                       #                       HospitalRate = hosp,
                       #                       CriticalRate = crit,
#                           "r0"=round(rzero,3),
#                           "dt"=round(dt,3))

#write.table(as.matrix(parameters3[,1:5]),paste0(exptdir,"ncovparams.csv"),row.names=F,col.names =F,sep=",")

#write.csv(parameters3,paste0(exptdir, "ncovparams_r0dt.csv"))

write.matrix(as.matrix(parameters3), paste0(exptdir,"design.csv"))
write.matrix(as.matrix(parameters3[,1:5]), paste0(exptdir,"ncovparams.csv"))





