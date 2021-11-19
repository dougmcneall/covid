## FFBS for early stage data
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(coda)

## source FF function
sourceCpp("iFFBS.cpp")

## read parameters
disease <- read_delim("disease.dat", delim = " ") %>%
    select(starts_with("beta"), ends_with("_8"))
colnames(disease) <- gsub("_8", "", colnames(disease))
colnames(disease) <- gsub("\\.", "", colnames(disease))

## read data
seeds <- read_csv("time_seeds.csv", col_names = FALSE) %>%
    filter(X2 == X2[1]) %>%
    select(-X2) %>%
    complete(X1 = seq(as.Date("2020-01-01"), as.Date("2020-03-13"), by = 1),
        fill = list(X3 = 0))
DH_prime <- seeds$X3
DI_prime <- rep(0, length(DH_prime))

ntest <- 30
## set parameters
test_max <- matrix(nrow = ntest, ncol = 11)

#test_i <- runif(30, min = 1, max = 250)

for(i in 1:length(test_i)){
    
parind <- test_i[i]
pH <- disease$pH[parind]
pHD <- disease$pHD[parind]
pI1 <- disease$pI1[parind]
pI1D <- disease$pI1D[parind]
pI1H <- disease$pI1H[parind]
pI2 <- disease$pI2[parind]
pP <- disease$pP[parind]
pE <- disease$pE[parind]
pEP <- disease$pEP[parind]
pA <- disease$pA[parind]
beta <- disease$`beta[3]`[parind]
betaA <- disease$`beta[6]`[parind]
#Npop <- 56082077
Npop <- 100000
Nposs <- 100
pini <- 1 / (Nposs * length(DI_prime))
#pini <- 1000 / (Nposs * length(DI_prime))

#parcheck <- c(pH, pHD, pI1, pI1D, pI1H, pI2, pP, pE, pEP, pA, beta , betaA, Npop/1e8, Nposs/100, pini )
#print(parcheck)

### set data
#DI_prime <- c(rep(0, 30), 0, 3)
#DH_prime <- c(rep(0, 30), 1, 1)
#Npop <- 100000

### set parameters
#pH <- 0.1
#pHD <- 0.1
#pI1 <- 0.2
#pI1D <- 0.2    
#pI1H <- 0.5
#pI2 <- 0.3
#pP <- 0.3
#pE <- 0.3
#pEP <- 0.8
#pA <- 0.8
#beta <- 1
#betaA <- 0.5


## run iFFBS algorithm
test <- iFFBS(DI_prime, DH_prime, 1:length(DI_prime), Npop, Nposs, 1000,
           pH, pHD, pI1, pI1D, pI1H, pI2, pP, pE, pEP, pA, beta, betaA, pini)

## trace plots
colnames(test) <- apply(expand.grid(c("S", "E", "A", "RA", "P", "I1", "I2", "RI", "DI", "H", "RH", "DH"), 1:length(DI_prime)), 1, function(x) paste0(x[1], "_", as.numeric(x[2])))
test <- as.mcmc(test)
#plot(test)
           
## plot marginal densities over time
test <- as.matrix(test) %>%
    as_tibble() %>%
    mutate(iter = 1:n()) %>%
    pivot_longer(!iter, names_to = "var", values_to = "n") %>%
    separate(var, c("var", "t"), sep = "_") %>%
    filter(var != "S") %>%
    mutate(t = as.numeric(t))

#p <- ggplot(test) +
#        geom_violin(aes(x = t, y = n, group = t)) +
#        facet_wrap(~var) +
#        xlab("Days from 1st January 2020") + 
#        ylab("Cumulative counts")
#ggsave("iFFBS2.pdf", p, width = 20, height = 10)

test73 <- filter(test, t==73, iter > 700)

test_out <- test73 %>%
    group_by(var) %>%
    summarise(mean = mean(n), min = min(n), max = max(n))

test_max[i, ] <- test_out$max

}

test_max_pini_1 <- test_max
colnames(test_max_pini_1) <- test_out$var

save(test_max_pini_1, test_max_pini_10, test_max_pini_10,test_max_pini_100, test_max_pini_1000, test_i, file = 'test_pini.rdata' )

## check counts
pivot_wider(test, names_from = var, values_from = n) %>%
    group_by(iter, t) %>%
    transmute(
        EAP = (E >= A + P),
        ARA = (A >= RA),
        PI1 = (P >= I1),
        I1I2HDI = (I1 >= I2 + H + DI),
        I2RI = (I2 >= RI),
        HRHDH = (H >= RH + DH)
    ) %>%
    mutate(all = all(EAP, ARA, PI1, I1I2HDI, I2RI, HRHDH)) %>%
    summary()
pivot_wider(test, names_from = var, values_from = n) %>%
    group_by(iter) %>%
    mutate(across(E:DH, ~{. - lag(.)})) %>%
    summary()
