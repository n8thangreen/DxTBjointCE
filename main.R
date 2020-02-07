#
# Nathan Green
# IDEA data
# Dx pathway joint costs, times


library(dplyr)
library(tidyr)


## readin data

data_dir <- "C:/Users/ngreen1/Google Drive/TB/IDEA/R/packages/IDEAdectree"
load(paste0(data_dir, "/data/TBdata_clinical_cleaned.RData"))


#############
# data prep #
#############

# in case of division by 0 cost, time
##TODO: test sensitivity
epsilon <- 0.0001

dat <-
  data %>%   
  transmute(Dosanjh,
            EPTBorPTB,
            Smear,
            start.to.diag,
            testDiagCon_diff,
            testDrug_diff,
            totalcost) %>% 
  drop_na(start.to.diag) %>% 
  mutate(log_cost = log(totalcost + epsilon),      #log-transform for normal distns
         log_time = log(start.to.diag + epsilon))

plot(dat$start.to.diag, dat$totalcost, xlim = c(0, 100))


# stratified
# jags_dat <- 
#   list(cost =
#          split(dat$log_cost, list(dat$EPTBorPTB, dat$Smear)),
#        time = 
#          split(dat$log_time, list(dat$EPTBorPTB, dat$Smear)))
# 
# ##TODO: to test...
# jags_dat <-
#   list(n = length(jags_dat$cost$PTB.POSITIVE),
#        c = jags_dat$cost$PTB.POSITIVE,
#        e = jags_dat$time$PTB.POSITIVE)

# covariate adjustment
jags_dat <- 
  list(n = nrow(dat),
       c = dat$log_cost,
       e = dat$log_time,
       # EPTBorPTB = dat$EPTBorPTB,
       status = as.numeric(dat$EPTBorPTB))

# 1: blank
# 2: EPTB;PTB
# 3: PTB


############
# run jags #
############

library(R2jags)
library(coda)


# jagsfile <- "jags.txt"
jagsfile <- "jags_covariates.txt"
params <- c("mu.c", "mu.e",
            "sigma.c", "sigma.e",
            "gamma.c", "gamma.e",
            "c_pred", "e_pred",
            "phi.c_pred", "phi.e_pred",
            "beta")
inits <- c(mu.c = 6,
           mu.e = 1,
           beta = 0,
           sigma.c = 1,
           sigma.e = 1,
           gamma.c = c(0,0,0),
           gamma.e = c(0,0,0))
n.iter <- 2000
n.burnin <- 100
n.thin <- 10

out <- 
  jags(data = jags_dat,
       parameters.to.save = params,
       # inits = inits,
       model.file = jagsfile,
       n.chains = 2,
       n.iter = n.iter,
       n.burnin = n.burnin,
       n.thin = n.thin,
       DIC = TRUE)


out.ls <- as.mcmc(out$BUGSoutput)

## mean costs on natural scale
## lognormal transformation

# m.c <- exp(out.ls$sims.list$mu.c + 0.5*out.ls$sims.list$sigma.c^2)
# m.e <- exp(out.ls$sims.list$mu.e + 0.5*out.ls$sims.list$sigma.e^2)

m.c <- exp(out.ls[[1]][, 'mu.c'] + 0.5*out.ls[[1]][, 'sigma.c']^2)
m.e <- exp(out.ls[[1]][, 'mu.e'] + 0.5*out.ls[[1]][, 'sigma.e']^2)

m.c_pred <- exp(out.ls[[1]][, 'phi.c_pred[2]'] + 0.5*out.ls[[1]][, 'sigma.c']^2)
m.e_pred <- exp(out.ls[[1]][, 'phi.e_pred[2]'] + 0.5*out.ls[[1]][, 'sigma.e']^2)



## save output


#########
# plots #
#########

## posteriors
library(MCMCvis)
library(mcmcplots)


MCMCvis::MCMCplot(out, params = c("beta", "mu.c", "mu.e", "sigma.c", "sigma.e", "gamma.c", "gamma.e"))
mcmcplots::mcmcplot(out)

plot(density(m.e))#, xlim = c(0, 5000))
abline(v = median(exp(jags_dat$e)), col= "red")
plot(density(m.c), xlim = c(0, 1500))
abline(v = mean(exp(jags_dat$c)), col= "red")

# TB status
plot(density(out.ls[[1]][, 'gamma.c[1]']), xlim = c(-1,1), main = "cost") #neither
lines(density(out.ls[[1]][, 'gamma.c[2]']), col = "red")   #EPTB;PTB
lines(density(out.ls[[1]][, 'gamma.c[3]']), col = "blue")  #PTB

plot(density(out.ls[[1]][, 'gamma.e[1]']), xlim = c(-2,2), main = "time") #neither
lines(density(out.ls[[1]][, 'gamma.e[2]']), col = "red")   #EPTB;PTB
lines(density(out.ls[[1]][, 'gamma.e[3]']), col = "blue")  #PTB

plot(as.numeric(m.e),
     as.numeric(m.c))#, xlim = c(0,2000))
plot(as.numeric(m.e_pred),
     as.numeric(m.c_pred))


## results table

library(reshape2)
library(tidyr)

##TODO:
# pivot_wider(data = tab)#, id_cols = "status")

tab <-
  MCMCvis::MCMCsummary(out) %>% 
  data.frame(param = rownames(.)) %>% 
  filter(grepl("phi", param)) %>% 
  select("param", "mean", "sd") %>% 
  separate(param, c("val", "status"), sep = "_") %>% 
  split(.$val)

merge(x = tab[[1]][, -1],
      y =  tab[[2]][, -1],
      by = "status",
      suffixes = c(".c",".e")) %>% 
  mutate_if(is.numeric, round, 2)


## model checking?




