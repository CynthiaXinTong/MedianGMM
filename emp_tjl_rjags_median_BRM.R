library(loo)
library(rjags)
library(mvtnorm)

options(scipen=999)
rm(list=ls())

Niter = 60000
burnIn = 30000

#-------Load Data-----------------------------------------------------------------
load(file = paste0('data/avg_tri_1819_n1000.Rdata'))

y = resp1
N = nrow(y)
Time = ncol(y)
Lambda = matrix(c(rep(1,Time), (c(1:Time)-1)), nrow = Time)

tau<-0.5

#------------------------------------------
# GCM
dat <- list("N" = N, "y" = y, "tau" = tau, "Time" = ncol(y))

# Set initial values
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 331)

jags.gmm.md1 <- jags.model( file = "JAGS/median/model_gcm_md.txt", data=dat, inits= initial, n.chains=1, n.adapt=1000 )
params <- c("par")
samps.gmm.md1 <- coda.samples(jags.gmm.md1, params, n.iter = Niter)

smp.md1 = window(samps.gmm.md1[[1]], start = burnIn)
geweke.md1 = apply(smp.md1, 2, function(x) geweke.diag(x)$z)
geweke.md1
summary(smp.md1)

#------------------------------------------
# GMM (G=2)
Niter = 60000
burnIn = 30000

load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
y = resp1
N = nrow(y)
Time = ncol(y)
Lambda = matrix(c(rep(1,Time), (c(1:Time)-1)), nrow = Time)
tau<-0.5

# Set initial values
source(file = paste0('sources/initial_BUGS_gmm_G2.R'))
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 88788,
             muLS = muLS.init)

dat.gmm.md2 <- list("N" = N, "y" = y, "Time" = ncol(y), "tau" = tau,  "alpha" = rep(10,2))
jags.gmm.md2 <- jags.model( file = paste0("JAGS/median/model_gmm_md_G2_diffvar.txt"), data=dat.gmm.md2,  n.chains=1, inits= initial, n.adapt=1000 )

params <- c("par")
samps.gmm.md2 <- coda.samples(jags.gmm.md2, params, n.iter = Niter)

smp.md2 = window(samps.gmm.md2[[1]], start = burnIn)
geweke.md2 = apply(smp.md2, 2, function(x) geweke.diag(x)$z)
geweke.md2
summary(smp.md2)
plot(smp.md2)


#------------------------------------------
# GMM (G=3)
Niter = 60000
burnIn = 30000

load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
y = resp1
N = nrow(y)
Time = ncol(y)
Lambda = matrix(c(rep(1,Time), (c(1:Time)-1)), nrow = Time)
tau<-0.5

# Set initial values
source(file = paste0('sources/initial_BUGS_gmm_G3.R'))
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 110,
             muLS = muLS.init)

dat.gmm.md3 <- list("N" = N, "y" = y, "Time" = ncol(y), "tau" = tau, "alpha" = rep(10,3))

jags.gmm.md3 <- jags.model( file = "JAGS/median/model_gmm_md_G3_diffvar.txt", data=dat.gmm.md3,  n.chains=1, inits= initial, n.adapt=1000 )
params <- c("par")
samps.gmm.md3 <- coda.samples(jags.gmm.md3, params, n.iter = Niter)

smp.md3 = window(samps.gmm.md3[[1]], start = burnIn)
geweke.md3 = apply(smp.md3, 2, function(x) geweke.diag(x)$z)
geweke.md3
summary(smp.md3)

plot(smp.md3)



#------------------------------------------
# GMM (G=4)
Niter = 60000
burnIn = 30000

load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
y = resp1
N = nrow(y)
Time = ncol(y)
Lambda = matrix(c(rep(1,Time), (c(1:Time)-1)), nrow = Time)
tau<-0.5


source(file = paste0('sources/initial_BUGS_gmm_G4.R'))
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 1551)
dat.gmm.md4 <- list("N" = N, "y" = y, "Time" = ncol(y), "tau" = tau, "alpha" = rep(10,4))

jags.gmm.md4 <- jags.model( file = "JAGS/median/model_gmm_md_G4_diffvar.txt", data=dat.gmm.md4, inits=initial, n.chains=1, n.adapt=1000 )
params <- c("par")
samps.gmm.md4 <- coda.samples(jags.gmm.md4, params, n.iter = Niter)

smp.md4 = window(samps.gmm.md4[[1]], start = burnIn)
geweke.md4 = apply(smp.md4, 2, function(x) geweke.diag(x)$z)
geweke.md4
summary(smp.md4)

plot(smp.md4)

