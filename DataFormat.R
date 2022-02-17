#----------------------------------------------------------------------------------#
#----Data formatting for orangutan BUGS model. Harvesting from loaded matrices.----#
#----Last updated on 3/28/2017 by Matthew Farr-------------------------------------#
#----------------------------------------------------------------------------------#

## AJM: All code confirmed functional in R 4.1.0 on 2021-06-16 

#-------------#
#-Import Data-#
#-------------#

## (AJM:changed to import strings as factors, so compatible with R 4.x)

M1 <- read.csv("Matrix1_OHcounts.csv", 
               stringsAsFactors = TRUE,
               header = TRUE)

M2 <- read.csv("Matrix2_StemsMRfruits.csv",  
               stringsAsFactors = TRUE,
               header = TRUE)

M3 <- read.csv("Matrix3_ObservationCovariates.csv",  
               stringsAsFactors = TRUE,
               header = TRUE)

M4 <- read.csv("Matrix4_SiteCovariates.csv",  
               stringsAsFactors = TRUE,
               header = TRUE)

M5 <- read.csv("actual_m_walked.csv",  
               stringsAsFactors = TRUE,
               header = TRUE)

M6 <- read.csv("Ntimes_sitelength_walked.csv",  
               stringsAsFactors = TRUE,
               header = TRUE)

M7 <- read.csv("Matrix_rain.csv",  
               stringsAsFactors = TRUE,
               header = TRUE)

M8 <- read.csv("Matrix_tmax.csv",
               stringsAsFactors = TRUE,
               header = TRUE)

M9 <- read.csv("Matrix_tmin.csv",  
               stringsAsFactors = TRUE,
               header = TRUE)


#--------------#
#-Harvest Data-#
#--------------#

#Number of times transects were walked
nx <- t(data.matrix(M6[-1]))

#Offset for detection
offset.S <- nx
for(i in 1:length(nx[,1])){
  for(j in 1:length(nx[1,])){
    if(nx[i,j] == 0)
      offset.S[i,j] <- 1
  }
}

#Observed Count of orangutans
y <- t(data.matrix(M1[-1]))
for(i in 1:length(nx[,1])){
  for(j in 1:length(nx[1,])){
    if(nx[i,j] == 0)
      y[i,j] <- NA
  }
}

#Habitat ID of each site
habitat <- as.numeric(M4$Habitat)

#Observed fruit availibility
fruit <- t(data.matrix(M2[-1]))
mu <- NULL
for(i in 1:7){mu[i] <- mean(fruit[,habitat==i], na.rm = TRUE)}
for(j in 1:27){for(t in 1:198){if(is.na(fruit[t,j])){fruit[t,j] <- mu[habitat[j]]}}}
fruit <- (fruit - mean(fruit, na.rm = TRUE))/sd(fruit, na.rm = TRUE)

colnames(fruit) <- NULL
rownames(fruit) <- NULL

#Observed precipitation
rain <- t(data.matrix(M7[-1]))
rain <- (rain - mean(rain, na.rm = TRUE))/sd(rain, na.rm = TRUE)
mu <- NULL
for(i in 1:7){mu[i] <- mean(rain[,habitat==i], na.rm = TRUE)}
for(j in 1:27){for(t in 1:198){if(is.na(rain[t,j])){rain[t,j] <- mu[habitat[j]]}}}

colnames(rain) <- NULL
rownames(rain) <- NULL

#Temperature highs
tempH <- t(data.matrix(M8[-1]))
tempH[is.infinite(tempH)] <- NA 
mu <- NULL
for(i in 1:7){mu[i] <- mean(tempH[,habitat==i], na.rm = TRUE)}
for(j in 1:27){for(t in 1:198){if(is.na(tempH[t,j])){tempH[t,j] <- mu[habitat[j]]}}}
tempH <- (tempH - mean(tempH, na.rm = TRUE))/sd(tempH, na.rm = TRUE)

colnames(tempH) <- NULL
rownames(tempH) <- NULL

#Temperature lows
tempL <- t(data.matrix(M9[-1]))
tempL[is.infinite(tempL)] <- NA
mu <- NULL
for(i in 1:7){mu[i] <- mean(tempL[,habitat==i], na.rm = TRUE)}
for(j in 1:27){for(t in 1:198){if(is.na(tempL[t,j])){tempL[t,j] <- mu[habitat[j]]}}}
tempL <- (tempL - mean(tempL, na.rm = TRUE))/sd(tempL, na.rm = TRUE)

colnames(tempL) <- NULL
rownames(tempL) <- NULL

#Distance of each observation
dst <- as.numeric(M3$Obs.dist)

#Site ID of each observation
site <- M3$Obs.Site
levels(site) <- c(levels(site), levels(M1[,1]))
site <- factor(site, levels = levels(M1[,1]))
site <- as.numeric(site)

#Replicate ID of each observation
rep <- M3$Rep.ID
levels(rep) <- c(levels(rep), sub("X", "", colnames(M1[-1])))
rep <- factor(rep, levels = sub("X", "", colnames(M1[-1])))
rep <- as.numeric(rep)

#Length of each site
site.length <- as.numeric(M4$Site.Length)

#Offset for abundance
offset.A <- site.length/mean(site.length)

#Elevation of each site
elev <- as.numeric(M4$Elevation)

#----------------#
#-Create Indices-#
#----------------#

#Width of distance categories
v <- 10

#Transect half width
B <- 50            

#Mid point of each distance category
midpt <- seq(5, 45, 10)       

#Number of individuals/observations
nind <- length(dst)

#Number of sites
nsites <- length(offset.A)

#Number of replicates
nreps <- length(y[,1])

#Number of months
nmon <- seq(1, nreps, 2)

#Number of distance categories
db <- seq(0, 50, 10)		  
nD <- length(db)-1

#Number of habitats
nhab <- length(levels(M4$Habitat))

#---------------#
#-Reformat Data-#
#---------------#

#Distance class (bins)
dclass <- array(NA, dim = nind)

dst[is.na(dst)] <- 999

for(i in 1:nind){
  if(dst[i] <= 10)
    dclass[i] <- 1
  if((10 < dst[i]) & (dst[i] <= 20))
    dclass[i] <- 2
  if((20 < dst[i]) & (dst[i] <= 30))
    dclass[i] <- 3
  if((30 < dst[i]) & (dst[i] <= 40))
    dclass[i] <- 4
  if((40 < dst[i]) & (dst[i] <= 50))
    dclass[i] <- 5
  if(dst[i] == 999)
    dclass[i] <- NA
}

#--------------#
#-Compile Data-#
#--------------#

DSdata <- list(y, fruit, rain, tempH, tempL, dclass, site, rep, elev, habitat, nhab, offset.S,
               offset.A, nind, nsites, nreps, nmon, nD, site.length, v, B, midpt)
names <- c("y", "fruit", "rain", "tempH", "tempL", "dclass", "site", "rep", "elev", "habID", "nhab", "offset.S",
           "offset.A", "nind", "nsites", "nreps", "nmon", "nD", "site.length", "v", "B", "midpt")
DSdata <- setNames(DSdata, nm = names)


dput(DSdata, file ="DSdata.R")
