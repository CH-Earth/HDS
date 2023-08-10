# define the time series
ndays <- 90
nhour <- 1
dt <- 1.0 / nhour
t <- seq(0, (nhour*ndays-1)*dt, by=dt) + 0.5

# define the seasonal statistics
pop  =   0.4
qTot =  75.0
pTot = 250.0
eMax =   7.0



# Define the synthetic forcing
nTot <- length(t)

# Define precipitation
uRand <- runif(nTot)
pDays <- ceiling(pop - uRand)
ixPOP <- which(pDays == 1)
lambda <- dt * length(ixPOP) / pTot
pRate <- rep(0, nTot)
if (length(ixPOP) > 0) {
  pRate[ixPOP] <- -log(runif(length(ixPOP))) / lambda
}

# Define the seasonal runoff
mu <- 1.5
sig <- 0.5
fSeas <- (dt / (t * sig * sqrt(2 * pi))) * exp(-((log(t) - mu)^2) / (2 * (sig^2)))
qSeas <- qTot * fSeas / dt

# Define ET
dStart <- 31 + 28 + 31 + 30 + 10  # time since the winter solstice
dSolar <- t + dStart
xSolar <- 0.5 - 0.5 * cos(2 * pi * dSolar / 365)
etPond <- xSolar * eMax

#create output forcing df
forcing <- data.frame(t, qSeas, pRate, etPond)
#write output to csv
write.csv(forcing, 'syntheticForcing.csv', row.names = F)
