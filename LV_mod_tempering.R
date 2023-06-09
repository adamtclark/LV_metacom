############################################
# try tempering experiment
############################################
M = 1  # number of sites
N0 = 100 # number of species (in global pool)
D = N0 # dimensionality of interaction matrix
minval = -99

Amu = -0.6
Asd = 0.3

n0mu = 0.05
n0sd = 0

Amat = matrix(nrow = N0, ncol = N0)
Amat[row(Amat)!=col(Amat)] = rnorm(N0^2-N0, Amu, Asd)
diag(Amat) = -1

Pars  <- c(N   = N0,                         # number of species
           M   = M,                          # number of patches
           minval = minval,                  # minimum population size to track
           K   = pmax(0, rnorm(N0, Kmu, Ksd)),   # vector of carrying capacities
           r  = rnorm(N0, rmu, rsd),         # vector of initial growth rates
           A  = c(Amat),                     # interaction coefficients
           d_c = 0.05,                        # proc. noise c
           d_z = 2.2,                        # proc. noise z
           d_n = 0,                          # proc. noise nugget
           d_m = 0,                          # proc. noise mean
           d_w = 2,                          # disturbance waiting time
           cv  = pmax(0, rnorm(N0, cmu, csd)),   # dispersal rate
           n0 = abs(rep(n0mu, N0)))          # initial abundance post-colonization

# plot disturbance function
n = seq(1e-3, 1, length=100)
dsd = sqrt(Pars["d_c"]*(n^Pars["d_z"]))
plot(n, dsd, log = "xy", type = "b"); abline(a=0, b=1, lty=2)
abline(v=exp(minval), lty =3)

plot(n, pnorm(0, n, dsd), type = "b")

# initial abundances
nini  <- abs(rnorm(N0*M, n0mu, n0sd))
#nini[as.logical(rbinom(length(nini), 1, 0.9))] = 0

# log transform
lnini = log(nini)
lnini[!is.finite(lnini)] = minval
lnini[as.logical(rbinom(N0,1,0.8))] = minval

# add value for disturbance
lnini = c(lnini, 1, runif(1))

# add values for dispersal
lnini = c(lnini, rep(1, N0), runif(N0, 0.5, 1))


niter = 1000
simout = matrix(nrow=niter, ncol = 5)
colnames(simout) = c("CV", "prod", "abund_mean", "rich", "aij_mean")

i = 1
while(i <= niter) {
  if(sum(which(lnini==minval))>0) {
    newps = sample(which(lnini==minval),1)
    lnini[newps] = log(abs(rnorm(1, n0mu, n0sd)))
  }
  
  dyn.load("LVmod_metacom_log.so")
  out_C <- ode(y = lnini, times = times, func = "derivs",
               parms = Pars, dllname = "LVmod_metacom_log", initfunc = "initmod",
               events = list(func = "event", root = TRUE), rootfun = "myroot", nout = N0, nroot = N0+1)
  dyn.unload("LVmod_metacom_log.so")
  
  Mout = array(dim =c(nrow(out_C), N0, M))
  Mout[] = exp(out_C[,1:(N0*M)+1])
  Mout[Mout<=exp(minval)] = 0
  
  if(FALSE) {
    par(mfrow=c(1,1), mar=c(4,4,2,2))
    matplot(times, apply(Mout, 1:2, mean), lty = 1, type = "l", lwd = 2,
            xlab = "time", ylab = "Species Abundance", col = viridis(N0))
    abline(h=0, lty=3)
  }
  
  ts_tmp = apply(Mout, 1, sum)[(max(times)/2):max(times)]
  simout[i,"CV"] = sd(ts_tmp)/mean(ts_tmp)
  simout[i,"prod"] = sum(Mout[length(times),,])
  alivesp = Mout[length(times),,]>0
  simout[i,"rich"] = sum(alivesp)
  
  simout[i,"abund_mean"] = mean(Mout[length(times),alivesp,])
  
  tmp = Amat[which(alivesp),,drop=FALSE][,which(alivesp),drop=FALSE]
  simout[i,"aij_mean"] = mean(tmp[row(tmp)!=col(tmp)])
  
  # get new starting conditions
  lnini = out_C[nrow(out_C),-1][1:length(lnini)]
  lnini[N0+1:2] = c(1, runif(1))
  lnini[N0+2+1:(N0*2)] = c(rep(1, N0), runif(N0, 0.5, 1))
  
  if(i/10 == floor(i/10)) {
    print(i/niter)
  }
  
  #simout[i,]
  i = i+1
}

plot(1:niter, simout[,"rich"], type = "l", ylab = "rich", xlab = "iteration")
plot(1:niter, simout[,"CV"], type = "l", log="y", ylab = "CV", xlab = "iteration")
plot(1:niter, simout[,"prod"], type = "l", ylab = "prod", xlab = "iteration")
plot(1:niter, simout[,"abund_mean"], type = "l", ylab = "abund_mean", xlab = "iteration")
plot(1:niter, simout[,"aij_mean"], type = "l", ylab = "aij_mean", xlab = "iteration")


## maybe add Eigenvalues?
## see how this changes with z?