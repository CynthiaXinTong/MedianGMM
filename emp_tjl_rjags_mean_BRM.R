library(mvtnorm)
library(rjags)
library(loo)

rm(list=ls())

Niter = 60000
burnIn = 30000

#-------Load Data-----------------------------------------------------------------
load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
y = resp1
N = nrow(y)
Time = ncol(y)


#------------------------------------------
# GCM
dat <- list("N" = N, "y" = y, "Time" = ncol(y))

# Set initial values
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 40)

jags.gmm.nm1 <- jags.model( file = "JAGS/traditional/model_gcm_nm.txt", data=dat, inits= initial, n.chains=1, n.adapt=1000 )
params <- c("par")
samps.gmm.nm1 <- coda.samples(jags.gmm.nm1, params, n.iter = Niter)

smp.nm1 = window(samps.gmm.nm1[[1]], start = burnIn)
geweke.nm1 = apply(smp.nm1, 2, function(x) geweke.diag(x)$z)
geweke.nm1
summary(smp.nm1)
plot(smp.nm1)


#------------------------------------------
# GMM (G=2)
Niter = 60000
burnIn = 30000

load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
y = resp1
N = nrow(y)
Time = ncol(y)

# Set initial values
source(file = paste0('sources/initial_BUGS_gmm_G2.R'))
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 7747,
             muLS = muLS.init)
dat.gmm.nm2 <- list("N" = N, "y" = y, "Time" = ncol(y), "alpha" = rep(10,2))

jags.gmm.nm2 <- jags.model( file = paste0("JAGS/traditional/model_gmm_nm_G2_diffvar.txt"), data=dat.gmm.nm2, n.chains=1, inits= initial, n.adapt=1000 )
params <- c("par")
samps.gmm.nm2 <- coda.samples(jags.gmm.nm2, params, n.iter = Niter)

smp.nm2 = window(samps.gmm.nm2[[1]], start = burnIn)
geweke.nm2 = apply(smp.nm2, 2, function(x) geweke.diag(x)$z)
geweke.nm2
summary(smp.nm2)
plot(smp.nm2)


#------------------------------------------
# GMM (G=3)
Niter = 60000
burnIn = 30000

load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
y = resp1
N = nrow(y)
Time = ncol(y)

# Set initial values
source(file = paste0('sources/initial_BUGS_gmm_G3.R'))
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 90012,
             muLS = muLS.init)

dat.gmm.nm3 <- list("N" = N, "y" = y, "Time" = ncol(y), "alpha" = rep(10,3))
jags.gmm.nm3 <- jags.model( file = paste0("JAGS/traditional/model_gmm_nm_G3_diffvar.txt"), data=dat.gmm.nm3, n.chains=1, inits= initial, n.adapt=1000 )
params <- c("par")
samps.gmm.nm3 <- coda.samples(jags.gmm.nm3, params, n.iter = Niter)

smp.nm3 = window(samps.gmm.nm3[[1]], start = burnIn)
geweke.nm3 = apply(smp.nm3, 2, function(x) geweke.diag(x)$z)
geweke.nm3
summary(smp.nm3)

plot(smp.nm3)




#------------------------------------------
# GMM (G=4)
Niter = 60000
burnIn = 30000

load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
y = resp1
N = nrow(y)
Time = ncol(y)

# Set initial values
source(file = paste0('sources/initial_BUGS_gmm_G4.R'))
initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 3888,
             muLS = muLS.init)

dat.gmm.nm4 <- list("N" = N, "y" = y, "Time" = ncol(y), "alpha" = rep(10,4))

jags.gmm.nm4 <- jags.model( file = paste0("JAGS/traditional/model_gmm_nm_G4_diffvar.txt"), data=dat.gmm.nm4, n.chains=1, inits= initial, n.adapt=1000 )
params <- c("par")
samps.gmm.nm4 <- coda.samples(jags.gmm.nm4, params, n.iter = Niter)

smp.nm4 = window(samps.gmm.nm4[[1]], start = burnIn)
geweke.nm4 = apply(smp.nm4, 2, function(x) geweke.diag(x)$z)
geweke.nm4
summary(smp.nm4)

plot(smp.nm4)
