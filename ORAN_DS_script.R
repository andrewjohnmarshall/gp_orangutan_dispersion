#------------------------------------------#
#----Script to run orangutan BUGS model----#
#----by Matthew Farr-----------------------#
#------------------------------------------#

# AJM: All code confirmed functional in R 4.1.0 on 2021-06-16 

#---------------#
#-Load Packages-#
#---------------#

library(jagsUI)

#-----------#
#-Load Data-#
#-----------#

DSdata <- dget("DSdata.R")

#-------------#
#-Attach Data-#
#-------------#

search()
attach(DSdata)

#------------#
#-BUGS MODEL-#
#------------#

sink("ORAN_DS.txt")
cat("
    model{

##Priors

#Transect effect
tau.j ~ dgamma(0.01, 0.01)
sig.j <- 1/sqrt(tau.j)

#Overdispersion priors
sig <- 1/sqrt(tau)
tau ~ dgamma(0.01, 0.01)

#Covariate parameters
for(k in 1:nhab){
alpha[k] ~ dnorm(mu.a, tau.a)
betaf[k] ~ dnorm(mu.f, tau.f)
betar[k] ~ dnorm(mu.r, tau.r)
betaH[k] ~ dnorm(mu.H, tau.H)
betaL[k] ~ dnorm(mu.L, tau.L)
}

mu.a ~ dnorm(0, 0.1)
mu.f ~ dnorm(0, 0.1)
mu.r ~ dnorm(0, 0.1)
mu.H ~ dnorm(0, 0.1)
mu.L ~ dnorm(0, 0.1)

tau.a ~ dgamma(0.1, 0.1)
tau.f ~ dgamma(0.1, 0.1)
tau.r ~ dgamma(0.1, 0.1)
tau.H ~ dgamma(0.1, 0.1)
tau.L ~ dgamma(0.1, 0.1)

#Detection parameter
asig ~ dnorm(0, 0.01)

#Zero inflation parameter
omega ~ dunif(0.0001, 1)

for(t in 1:nreps){

#Random effect of time
time[t] ~ dnorm(0, tau)

}

##Likelihood

for(j in 1:nsites){

psi[j] ~ dnorm(0, tau.j)

for(t in 1:nreps){

#Detection function scale parameter
sigma[t,j] <- exp(asig +log(offset.S[t,j]))

#Construct cell probabilities for nD cells using numerical integration
#Sum of the area (rectangles) under the detection function
for(k in 1:nD){

#Half-normal detection function at midpoint (length of rectangle)
p[k,t,j] <- exp(-midpt[k]*midpt[k]/(2*sigma[t,j]*sigma[t,j]))

#Probability of x in each interval (width of rectangle) for both sides of the transect
#B is half width of the transect and v is the width of each interval
pi[k,t,j] <- v/B 

#Detection probability for each interval (area of each rectangle)
f[k,t,j] <- p[k,t,j] * pi[k,t,j]

#Conditional detection probability (scale to 1)
fc[k,t,j] <- f[k,t,j]/pcap[t,j]

}#end k loop

#Detection probability at each site (sum of rectangles)
pcap[t,j] <- sum(f[1:nD,t,j])

#Observed population @ each t,j
y[t,j] ~ dbin(pcap[t,j], N[t,j])  ##Part 2 of HM

#Abundance @ each t,j
N[t,j] ~ dpois(lambda.eff[t,j])  ##Part 3 of HM

#Zero inflation component
lambda.eff[t,j] <- z[t,j]*lambda[t,j]
z[t,j] ~ dbern(omega)

#Linear predictor for mean abundance with offset for transect length and overdispersion parameter
lambda[t,j] <- exp(alpha[habID[j]] + 
                   betaf[habID[j]] * fruit[t,j] + betar[habID[j]] * rain[t,j] +
                   betaH[habID[j]] * tempH[t,j] + betaL[habID[j]] * tempL[t,j] +
                   log(offset.A[j]) + psi[j] + time[t])    

#Density; Area = site length * transect width * m to km converstion (1e+06)
Den[t,j] <- N[t,j] / (site.length[j] * 100 / 1e+06)

}#end t loop

}#end j loop

#Multinomial component to determine detection probability
for(i in 1:nind){

#Conditional detection probability at each nD interval
dclass[i] ~ dcat(fc[1:nD, rep[i], site[i]]) ## Part 1 of HM

}#end i loop

##Derived quantities

#Mean Sigma
mean.sig <- mean(sigma[,])

for(j in 1:nsites){

for(t in 1:nreps){

#Estimate fruit covariate NAs
fruit[t,j] ~ dnorm(0, 1)

} #end t loop

} #end j loop

for(tt in nmon){

#Derived latent density
D[tt,1] <- mean(Den[tt:(tt+1),1:4])
D[tt,2] <- mean(Den[tt:(tt+1),5:7])
D[tt,3] <- mean(Den[tt:(tt+1),8:12])
D[tt,4] <- mean(Den[tt:(tt+1),13:16])
D[tt,5] <- mean(Den[tt:(tt+1),17:21])
D[tt,6] <- mean(Den[tt:(tt+1),22:23])
D[tt,7] <- mean(Den[tt:(tt+1),24:27])

#Expected abundance given linear predictor
ED[tt,1] <- mean(lambda[tt:(tt+1),1:4])
ED[tt,2] <- mean(lambda[tt:(tt+1),5:7])
ED[tt,3] <- mean(lambda[tt:(tt+1),8:12])
ED[tt,4] <- mean(lambda[tt:(tt+1),13:16])
ED[tt,5] <- mean(lambda[tt:(tt+1),17:21])
ED[tt,6] <- mean(lambda[tt:(tt+1),22:23])
ED[tt,7] <- mean(lambda[tt:(tt+1),24:27])

}#end tt

}#end model
    ",fill=TRUE)
sink()

#-----------------------------#
#-Compile Data for BUGS Model-#
#-----------------------------#

str(data <- list(nD = nD, v = v, site = site, rep = rep, 
                 site.length = site.length,
                 y = y, B = B, midpt = midpt, 
                 nind = nind, nmon = nmon, 
                 dclass = dclass, nsites = nsites, 
                 nreps = nreps, nhab = nhab, 
                 offset.S = offset.S, offset.A = offset.A, 
                 fruit = fruit, habID = habID, 
                 rain = rain, tempL = tempL, tempH = tempH))


#-----------------------#
#-Create Initial Values-#
#-----------------------#

Nst <- y + 1

inits <- function(){list(N = Nst, alpha = runif(nhab, -1, 1),
                         asig = runif(1, 3, 4))}

#-----------------------#
#-Parameters to Monitor-#
#-----------------------#

params <- c("tau.j", "sig.j", "tau", "alpha",
            "betaf", "betar", "betaH", "betaL", 
            "mu.a", "mu.f", "mu.r", "mu.H", "mu.L", 
            "tau.a", "tau.f", "tau.r", "tau.H", "tau.L", 
            "asig", "omega", "D")

#---------------#
#-MCMC Settings-#
#---------------#

nc <- 3
ni <- 120000
nb <- 100000
nt <- 10

#----------------#
#-Run BUGS Model-#
#----------------#

ORAN_DS <- jags(data = data, inits = inits, parameters.to.save = params, 
                model.file = "ORAN_DS.txt", 
                n.chains = nc, n.iter = ni, n.burnin = nb, 
                n.thin = nt, store.data = TRUE, parallel = TRUE)

#-------------#
#-Detach Data-#
#-------------#

detach(DSdata)
search()

#-------------#
#-Save Output-#
#-------------#

save(ORAN_DS, file = "ORAN_OUTPUT.Rdata")
