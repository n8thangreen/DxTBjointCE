#
# Nathan Green
# IDEA data
# Diagnostic pathway (joint) costs, times


library(dplyr)
library(tidyr)
library(purrr)


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
            SiteID,
            EPTBorPTB,
            Smear,
            # start.to.diag,
            # testDiagCon_diff,
            # SmearSputumRes, #or Smear? whats the difference? Smear is more complete
            Smear,
            testDrug_diff,
            totalcost) %>% 
  # drop_na(start.to.diag) %>% 
  drop_na(testDrug_diff) %>% 
  mutate(
    testDrug_diff = as.numeric(testDrug_diff),
    testDrug_diff = if_else(testDrug_diff == 0, 1, testDrug_diff), # at least a single day before starting treatment
    log_cost = log(totalcost + epsilon),       #log-transform for normal distns
    # log_time = log(start.to.diag + epsilon)
    # log_time = log(testDrug_diff + epsilon)
    log_time = log(testDrug_diff)
  ) %>% 
  filter(Dosanjh %in% c(1,2))

# plot(x = dat$start.to.diag,
#      y = dat$totalcost, xlim = c(0, 100))

hist(dat$totalcost, breaks = 30)
hist(dat$testDrug_diff, breaks = 80)

hist(dat$log_cost, breaks = 20)
hist(dat$log_time, breaks = 20)


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
       ssm = as.numeric(dat$Smear),
       tbstatus = as.numeric(dat$EPTBorPTB))

# 1: blank
# 2: EPTB;PTB
# 3: PTB

# 1: Not taken
# 2: NEGATIVE
# 3: POSITIVE


############
# run jags #
############

library(R2jags)
library(coda)


# jagsfile <- "jags.txt"
# jagsfile <- "jags_covariates.txt"
jagsfile <- "jags_covariates_iid.txt"

params <-
  c("mu.c", "mu.e",
    "gamma.c", "gamma.e",
    "sigma.c", "sigma.e",
    "c_pred", "e_pred",
    "phi.c_pred", "phi.e_pred"#,
    # "beta"
  )

inits <- c(mu.c = 6,
           mu.e = 1,
           beta = 0,
           sigma.c = 1,
           sigma.e = 1,
           gamma.c = c(0,0,0),
           gamma.e = c(0,0,0))

n.iter <- 20000
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

save(out.ls, file = "data/out.ls.RData")

## mean costs on natural scale
## lognormal transformation

# m.c <- exp(out.ls$sims.list$mu.c + 0.5*out.ls$sims.list$sigma.c^2)
# m.e <- exp(out.ls$sims.list$mu.e + 0.5*out.ls$sims.list$sigma.e^2)
nchain <- 1

m.c <- exp(out.ls[[nchain]][, 'mu.c'] + 0.5*out.ls[[nchain]][, 'sigma.c']^2)
m.e <- exp(out.ls[[nchain]][, 'mu.e'] + 0.5*out.ls[[nchain]][, 'sigma.e']^2)

# m.c_pred <- exp(out.ls$sims.list$phi.c_pred[,2] + 0.5*out.ls$sims.list$sigma.c^2)
# m.e_pred <- exp(out.ls$sims.list$phi.c_pred[,2] + 0.5*out.ls$sims.list$sigma.e^2)

m.c_pred <-
  rep( list(list()), 3) %>%
  setNames(c("", "E;PTB", "PTB"))

m.e_pred <-
  rep( list(list()), 3) %>%
  setNames(c("", "E;PTB", "PTB"))

for (i in 1:3) {
  for (j in 1:3) {
    
    phi.c_pred <- paste0("phi.c_pred[", i, ",", j, "]")
    phi.e_pred <- paste0("phi.e_pred[", i, ",", j, "]")
    
    m.c_pred[[i]][[j]] <- exp(out.ls[[nchain]][, phi.c_pred] + 0.5*out.ls[[nchain]][, 'sigma.c']^2)
    m.e_pred[[i]][[j]] <- exp(out.ls[[nchain]][, phi.e_pred] + 0.5*out.ls[[nchain]][, 'sigma.e']^2)
  }
}


## save output


#########
# plots #
#########

library(MCMCvis)
library(mcmcplots)


MCMCvis::MCMCplot(out,
                  params =
                    c(
                      # "beta",
                      "mu.c", "mu.e",
                      "sigma.c", "sigma.e",
                      "gamma.c", "gamma.e"))

mcmcplots::mcmcplot(out)

plot(density(m.e), xlim = c(0, 5000))
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


## for tbstatus only model
par(mfrow = c(2, 2))
hist(as.numeric(m.e_pred$`E;PTB`), breaks = 30, xlim = c(10, 50))
hist(as.numeric(m.c_pred$`E;PTB`), breaks = 40, xlim = c(500, 1400))

hist(as.numeric(m.e_pred$PTB), breaks = 30, xlim = c(10, 50))
hist(as.numeric(m.c_pred$PTB), breaks = 30, xlim = c(500, 1400))

## for ssm & tbstatus  model
par(mfrow = c(2, 2))

for (i in 1:3) {
  
  hist(as.numeric(m.e_pred$`E;PTB`[[i]]), breaks = 30, xlim = c(10, 50), main = "")
  hist(as.numeric(m.c_pred$`E;PTB`[[i]]), breaks = 40, xlim = c(500, 1400), main = "")
  
  hist(as.numeric(m.e_pred$PTB[[i]]), breaks = 30, xlim = c(10, 50), main = "")
  hist(as.numeric(m.c_pred$PTB[[i]]), breaks = 30, xlim = c(500, 1400), main = "")
  
  mtext(levels(data$Smear)[i], side = 3, line = -1, outer = TRUE)
}


## results tables

library(reshape2)
library(tidyr)

##TODO:
# pivot_wider(data = tab)#, id_cols = "status")

# tab <-
#   MCMCvis::MCMCsummary(out) %>% 
#   data.frame(param = rownames(.)) %>% 
#   filter(grepl("phi", param)) %>% 
#   select("param", "mean", "sd") %>% 
#   separate(param, c("val", "status"), sep = "_") %>% 
#   split(.$val)
# 
# merge(x = tab[[1]][, -1],
#       y =  tab[[2]][, -1],
#       by = "status",
#       suffixes = c(".c",".e")) %>% 
#   mutate_if(is.numeric, round, 2)


# mean cost, time by grouping

## tbstatus only model
tribble(~name,   ~'mean cost',           ~'sd cost',           ~'mean time',           ~'sd time',
        "E;PTB", mean(m.c_pred$`E;PTB`), sd(m.c_pred$`E;PTB`), mean(m.e_pred$`E;PTB`), sd(m.e_pred$`E;PTB`),
        "PTB",   mean(m.c_pred$`PTB`),   sd(m.c_pred$`PTB`),   mean(m.e_pred$`PTB`),   sd(m.e_pred$`PTB`)) 

## smm & tbstatus
tribble(~ ssm,       ~'tb status', ~'mean cost',                ~'sd cost',                ~'mean time',                ~'sd time',
        "Not taken", "E;PTB",      mean(m.c_pred$`E;PTB`[[1]]), sd(m.c_pred$`E;PTB`[[1]]), mean(m.e_pred$`E;PTB`[[1]]), sd(m.e_pred$`E;PTB`[[1]]),
        "Not taken", "PTB",        mean(m.c_pred$`PTB`[[1]]),   sd(m.c_pred$`PTB`[[1]]),   mean(m.e_pred$`PTB`[[1]]),   sd(m.e_pred$`PTB`[[1]]), 
        "Negative",  "E;PTB",      mean(m.c_pred$`E;PTB`[[2]]), sd(m.c_pred$`E;PTB`[[2]]), mean(m.e_pred$`E;PTB`[[2]]), sd(m.e_pred$`E;PTB`[[2]]),
        "Negative",  "PTB",        mean(m.c_pred$`PTB`[[2]]),   sd(m.c_pred$`PTB`[[2]]),   mean(m.e_pred$`PTB`[[2]]),   sd(m.e_pred$`PTB`[[2]]),
        "Positive",  "E;PTB",      mean(m.c_pred$`E;PTB`[[3]]), sd(m.c_pred$`E;PTB`[[3]]), mean(m.e_pred$`E;PTB`[[3]]), sd(m.e_pred$`E;PTB`[[3]]),
        "Positive",  "PTB",        mean(m.c_pred$`PTB`[[3]]),   sd(m.c_pred$`PTB`[[3]]),   mean(m.e_pred$`PTB`[[3]]),   sd(m.e_pred$`PTB`[[3]])) 



#################
# pathway probs #
#################

# raw pathway probabilities
p_tbl <-
  dat %>%
  count(Smear, EPTBorPTB) %>%
  mutate(p = n/sum(n))

knitr::kable(p_tbl)


############
# run jags #
############

jagsfile <- "jags_pathprobs.txt"

jags_dat <- 
  list(N = nrow(dat),
       X = p_tbl$n)

# 1: blank
# 2: EPTB;PTB
# 3: PTB

# 1: Not taken
# 2: NEGATIVE
# 3: POSITIVE

params <- c("prob")

# inits <- c(mu.c = 6,
#            sigma.e = 1)

n.iter <- 20000
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


out_pathprobs <- as.mcmc(out$BUGSoutput)

save(out_pathprobs, file = "data/out.ls_pathprobs.RData")
