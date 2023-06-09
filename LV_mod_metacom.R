#setwd("~/Dropbox/Projects/095_OldField2022/src/LV_metacom/")

rm(list=ls())

require(deSolve)
require(viridis)

system("R CMD SHLIB LVmod_metacom_log.c")


n0mu = 0.01
n0sd = n0mu/2
Kmu = 1
Ksd = 0
rmu = 1
rsd = 0
Amu = -0.5
Asd = 0.2
cmu = 0.5
csd = 0.2

minval = -999
M = 50  # number of sites
N0 = 20 # number of species (in global pool)
D = N0 # dimensionality of interaction matrix

Amat = matrix(nrow = N0, ncol = N0)
Amat[row(Amat)!=col(Amat)] = rnorm(N0^2-N0, Amu, Asd)
diag(Amat) = -1

Pars  <- c(N   = N0,                         # number of species
           M   = M,                          # number of patches
           minval = minval,                  # minimum population size to track
           K   = pmax(0, rnorm(N0, Kmu, Ksd)),   # vector of carrying capacities
           r  = rnorm(N0, rmu, rsd),         # vector of initial growth rates
           A  = c(Amat),                     # interaction coefficients
           d_c = 0.1,                        # proc. noise c
           d_z = 1.2,                        # proc. noise z
           d_n = 0,                          # proc. noise nugget
           d_m = 0,                          # proc. noise mean
           d_w = 2,                          # disturbance waiting time
           cv  = pmax(0, rnorm(N0, cmu, csd)),   # dispersal rate
           n0 = abs(rep(n0mu, N0)))          # initial abundance post-colonization

# initial abundances
nini  <- abs(rnorm(N0*M, n0mu, n0sd))
#nini[4:6] = minval
#nini[N0+1:5-1] = 0

times <- seq(0, 1e2, by = 1)

# log transform
lnini = log(nini)
lnini[!is.finite(lnini)] = minval

# add value for disturbance
lnini = c(lnini, 1, runif(1))

# add values for dispersal
lnini = c(lnini, rep(1, N0), runif(N0, 0.5, 1))

# check t_dist and t_dist
tdist0 = -log(1-lnini[c(N0*M+2)])*Pars["d_w"]
tdisp0 = -log(1-lnini[c(N0*M+2)])/Pars["d_w"]

#dyn.load("LVmod_metacom_log.so")
#DLLfunc(y = lnini, times = times, func = "derivs",
#        parms = Pars, dllname = "LVmod_metacom_log", initfunc = "initmod",
#        nout = N0)
#dyn.unload("LVmod_metacom_log.so")

dyn.load("LVmod_metacom_log.so")
system.time(out_C <- ode(y = lnini, times = times, func = "derivs",
                         parms = Pars, dllname = "LVmod_metacom_log", initfunc = "initmod",
                         events = list(func = "event", root = TRUE), rootfun = "myroot", nout = N0, nroot = N0+1))
dyn.unload("LVmod_metacom_log.so")

#times = out_C[,1]
Mout = array(dim =c(nrow(out_C), N0, M))
Mout[] = exp(out_C[,1:(N0*M)+1])

# plot average global dynamics
par(mfrow=c(1,1), mar=c(4,4,2,2))
matplot(times, apply(Mout, 1:2, mean), lty = 1, type = "l", lwd = 2,
        xlab = "time", ylab = "Species Abundance", col = viridis(N0))
abline(h=0, lty=3)

# plot local dynamics
par(mfrow=c(2,1), mar=c(4,4,2,2))
for(i in 1:M) {
  matplot(times, Mout[,,i], lty = 1, type = "l", lwd = 2,
          xlab = "time", ylab = "n", col = viridis(N0))
  abline(h=0, lty=3)
}

Pout = out_C[,ncol(out_C)+1-(N0:1)]
par(mfrow=c(1,1), mar=c(4,4,2,2))
matplot(times, Pout, lty = 1, type = "l", lwd = 2,
        xlab = "time", ylab = "p", col = viridis(N0))
abline(h=0, lty=3)

# check p vs. occupency
tmp = apply(Mout, 1:2, function(x) mean(x>0))
plot(c(Pout), c(tmp)); abline(a=0, b=1, lty=2)

# check events
# show events
par(mfrow=c(1,1), mar=c(4,4,2,2))
matplot(out_C[,1], out_C[,c((N0*M)+2, (N0*M)+2+1+1:N0)], type = "l", col = c(1, viridis(N0)), lty = 1)


# TODO:
# add spatial covariance in disturbances
# add species-level covariance in disturbances
# add spatially explicit dispersal?
# add back in Jacobian
# add niche dimensionality to alpha


