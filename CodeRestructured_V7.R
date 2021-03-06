################  CONSTANT DECLARATION, LOAD PACKAGES, LOAD DATA  ################  

OUTPUT_PLOTS <- FALSE

library(scatterplot3d)
library("DEoptim")
library("plotly")
library("reshape2")

install.packages('quantmod')
library('quantmod')

load(file = "fm408_exam_data.RData")
epsilon<-0.0001
mean<-0
var<-0
varHelp<-0

## parameters for optimisation function DEoptim
## l: vector of lower bound of parameters
## u: vector of upper bound of parameters
l <- c(-20,0,-0.99,-20 ,1e-5)
u <- c(20,50,0.99,5,20)


# model assumptions
#should we replace with something?
kappa=1 # rate of reversion to "current target variance theta(t)"
thetaConst=0.1 # long-term variance
sigma1=0.1 # volatility of volatility 
sigma2=0.1 # volatility of theta(t)
rho=0.2 # weight factor
eta=0.2 # rate of reversion to long-term target variance thetaConst


# Option Assumptions
#S=100 # stock price today
#S=exp(-raw$InterestRate[1]+raw$DividendYield[1])*raw$Forward[1]
V=0.5 # volatility of Future
#Rf=0.5 # risk-free rate; Assumed constant
#months=12 # months until maturity
#Div=0 #dividend yield
#F=S*(1+Rf)^(months/12)

# Initialize variables
#dt is monthly

dt=1/12 

################  Utility Functions  ################  

#Taken from PS7 solution
BSprice<-function(pc, S, k, vol, d, r, t)
{
  #pc  put/call indicator call=1, put=-1
  #S   Stock price at 0
  #K   strike
  #vol volatility
  #d   dividend yield
  #r   riskless rate
  #t   time to maturity
  
  
  d1 = (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  d2 = d1 - vol * sqrt(t)
  
  BSprice = pc * exp(-d * t) * S * 
    pnorm(pc * d1) - pc * k * exp(-r * t) * pnorm(pc * d2)
  return(BSprice)
}

#Derive a function f two times at point x (using approximation formula, x is the delta used)
secondDerivative = function(f, x, delta) {
  res = (f(x - delta) - 2*f(x) + f(x + delta)) / (delta^2)
  return (res)
}


#Like vfunc, but gives the implied vol according to the SVI curve instead of the squared error
## p: vector of svi parameters in the following order
##    (a, b, rho, m, sigma)
## k: log moneyness (log(K/F))
## maturity: options' time to maturity in years
implVol = function(p, k, maturity) {
  a=p[1]
  b=p[2]
  rho=p[3]
  m=p[4]
  sigma2=p[5]^2
  
  kk=(k-m)^2+sigma2
  totalIV=a+b*(rho * (k-m)+sqrt(kk))
  return(sqrt(totalIV/maturity))
}


## p: vector of svi parameters in the following order
##    (a, b, rho, m, sigma)
## k: log moneyness (log(K/F))
## maturity: options' time to maturity in years
## v: observed implied volatility
vfunc=function(p,k, maturity, v){
  a=p[1]
  b=p[2]
  rho=p[3]
  m=p[4]
  sigma2=p[5]^2
  
  kk=(k-m)^2+sigma2
  totalIV=a+b*(rho * (k-m)+sqrt(kk))
  
  if(!is.numeric(totalIV)) return(1000)
  if(min(totalIV)<0) return(1000)
  
  ## students should add restrictions here
  # e.g. if(some conditions give TRUE/FALSE as output) return(1000)
  # restrictions from task description
  if(!is.numeric(a)) return(1000)
  if(b<0) return(1000)
  if(abs(rho)>=1) return(1000)
  if(!is.numeric(m)) return(1000) 
  if(p[5]<=0) return(1000) #p[5] = sigma
  if((a + b*p[5]*sqrt(1-rho^2)) < 0) return(1000)
  
  # sum of squares
  res=sum((totalIV-maturity*((v)^2))^2)
  
  return(res)
}

calcImplDensity <- function(date, term) {
  var<--1
  
  while(var < 0) {
    #Table of options that are used to fit the SVI curve (i.e. all options that match a given t and T (but have different Strikes))
    optionsToFit<-raw[raw$ValuationDate == date, ]
    optionsToFit<-optionsToFit[optionsToFit$Term == term, ]
    
    lIntegrBound <- 100
    uIntegrBound <- 10000
    
    #Vectors used for the implied vol surface plot
    moneynessPlotVec <- c()
    impliedVolPlotVec <- c()
    termPlotVec <- c()
    
    moneynessLogPlotVec <- c()
    impliedLogVolPlotVec <- c()
    
    #Note: The SVI interpolation is non-deterministic because the optimization routine is random - results therefore vary
    var<--1
    
    func = function(p) {
      sum = 0
      for (i in 1:nrow(optionsToFit)) {
        #Minimize the squared error (vfunc calculates squared error, sum gives the total squared error from all options used to calculate the SVI curve)
        sum = sum + vfunc(p, optionsToFit$Moneyness[i], optionsToFit$TTM[i], optionsToFit$ImpliedVol[i])
      }
      return(sum)
    }
    
    #Then write the optimization function applied to vfunc
    outDEoptim <- DEoptim(func, l, u, DEoptim.control(VTR = 0.5e-3, itermax =  1e4))
    summary(outDEoptim)
    
    #The SVI curve. Note:k is the log moneyness, not the strike
    bestFit = function(k) {implVol(outDEoptim$optim$bestmem, k, optionsToFit$TTM[1])}
    bestFitLog = function(k) {bestFit(log(k))} 
    
    #Code for Implied Vol plot
    for(i in 1:3000) {
      r <- (i/3000)*4 - 2
      moneynessPlotVec <- c(moneynessPlotVec, r)
      impliedVolPlotVec <- c(impliedVolPlotVec, bestFit(r))
      termPlotVec <- c(termPlotVec, term)
      
      r2 <- (i/3000)*(exp(2) - exp(-2)) + exp(-2)
      moneynessLogPlotVec <- c(moneynessLogPlotVec, r2)
      impliedLogVolPlotVec <- c(impliedLogVolPlotVec, bestFitLog(r2))
    }
    
    #Use this to plot the SVI curve
    if (isTRUE(OUTPUT_PLOTS)) {
      png(filename = paste("SVILOG", date, "_", term, ".png", sep = ""))
      curve(bestFit, from=-2, to=2, xlab="Moneyness [log(K/F)]", ylab="Implied Vol")
      title(main = paste("Implied Vol Curve - ", date, "Term:", term))
      dev.off()
      
      png(filename = paste("SVI", date, "_", term, ".png", sep = ""))
      curve(bestFitLog, from=exp(-2), to=exp(2), xlab="Moneyness [K/F]", ylab="Implied Vol")
      title(main = paste("Implied Vol Curve - ", date, "Term:", term))
      dev.off()
    }
    
    C = function(K) {BSprice(1, optionsToFit$Forward[1]*exp(-(optionsToFit$InterestRate[1] - optionsToFit$DividendYield[1])*optionsToFit$TTM[1]), K, bestFit(log(K/optionsToFit$Forward[1])), optionsToFit$DividendYield[1], optionsToFit$InterestRate[1], optionsToFit$TTM[1])}
    
    #The implied density as a function of the strike. Note: Not the implied return density
    implDensity = function(K) {
      Dt = exp(-optionsToFit$TTM[1]*(optionsToFit$InterestRate[1]))
      secDeriv = secondDerivative(C, K, 0.5)
      return(secDeriv / Dt)
    }

    integral = tryCatch(integrate(implDensity, lower = lIntegrBound, upper = 10000), error = function(e) e)  
    if(inherits(integral, "error")) next
    
    implDensityScaled = function(K) {return(implDensity(K) / integral$value)}
    
    #Use this to plot implied densities    
    if(isTRUE(OUTPUT_PLOTS)) {
      png(filename = paste("RNDST", date, "_", term, ".png", sep = ""))
      curve(implDensityScaled, from=200, to=2500, xlab="ST", ylab="Risk Neutral Density")
      title(main = paste("Risk Neutral Density - ", date, "Term:", term))
      dev.off()
    }
    
    #Function only used for plotting
    #Input: Return - ouput: RND
    returnRND = function(R) {
      ST = exp(R)*optionsToFit$Forward[1]
      return(implDensityScaled(ST))
    }
    
    integralReturn = tryCatch(integrate(returnRND, lower = log(lIntegrBound/optionsToFit$Forward[1]), upper = log(10000/optionsToFit$Forward[1])), error = function(e) e)  
    if(inherits(integralReturn, "error")) next
    
    returnRNDScaled = function(R) {return(returnRND(R) / integralReturn$value)}

    #Use this to plot implied densities
    if (isTRUE(OUTPUT_PLOTS)) {
      png(filename = paste("RND", date, "_", term, ".png", sep = ""))
      curve(returnRNDScaled, from=log(200/optionsToFit$Forward[1]), to=log(2500/optionsToFit$Forward[1]), xlab="Log Holding Period Excess Return R", ylab="Risk Neutral Density")    
      title(main = paste("Risk Neutral Density -", date, "Term:", term))
      dev.off()
    }
    
    #Output: "return times probability of the return"
    returnImplDensity = function(ST) {
      res<-log(ST/optionsToFit$Forward[1])*implDensityScaled(ST)
      return(res)
    } 
    
    #Output: "return^2 times probability of the return"
    #cannot just square returnImplDensity (otherwise you would square the probability as well)
    returnImplDensitySquared = function(ST) {
      res<-log(ST/optionsToFit$Forward[1])^2*implDensityScaled(ST)
      return(res)
    } 
    #returnImplDensity = function(r) {return(returnImplDensityUnscaled(r)/integrate(returnImplDensityUnscaled, lower=lIntegrBound,upper=uIntegrBound)$value)}
      
    #curve(implDensity, from=lIntegrBound, to=5000, xlab="R", ylab="Implied Density")
    #curve(implDensityScaled, from=lIntegrBound, to=5000, xlab="R", ylab="Implied Log Returns")
      
      
    #calculate the mean of the function by integration, and the var by using the usual formula
    mean = tryCatch(integrate(function(x) {return(returnImplDensity(x))}, lIntegrBound, uIntegrBound), error = function(e) e)
    if(inherits(mean, "error")) next
      
    varHelp = tryCatch(integrate(function(x) {return(returnImplDensitySquared(x))}, lIntegrBound, uIntegrBound), error = function(e) e)
    if(inherits(varHelp, "error")) next
    
    mean = mean$value
    var = (varHelp$value - mean^2)/(term/12)  
  
    res <- data.frame(mean, var, moneynessPlotVec, impliedVolPlotVec, termPlotVec, moneynessLogPlotVec, impliedLogVolPlotVec)
    
  }
  return(res)
}

average <- function(list) {
  sum <- 0
  for (i in 1:length(list)) {
    sum = sum + list[i]
  }
  return(sum/length(list))
}

monteCarlo <- function(date, term, iter) {  
  
  observedData <- raw[raw$ValuationDate == date, ]
  observedData <- observedData[observedData$Term == term, ]
  
  #TODO: Dividend yield??
  Rf <- observedData$InterestRate[1]
  Div <- observedData$DividendYield[1]
  S <- observedData$Forward[1]*exp(-(Rf - Div)*observedData$TTM[1])
  Forward <- observedData$Forward[1]
  
  simulatedReturn = rep(1,iter)
  simulatedVol = rep(1,iter)
  
  for (i in 1:iter) {
    lnS=log(S)
    St=S
    Vt<-V
    #not sure it is the right inizialization
    #I made a new variable TT to prevent 
    thetaT<-thetaConst
    TT<-thetaT  
    
    for (j in 1:term) {
      
      #brownian of stock
      BrownS=qnorm(runif(1,0,1),0,1)
      #Brownian of variance
      BrownV=rho*BrownS+sqrt(1-rho^2)*qnorm(runif(1,0,1),0,1)
      #brownian of theta
      BrownT=sigma2*qnorm(runif(1,0,1),0,1)
      
      lnS=lnS+(Rf-Div-0.5*Vt)*dt+sqrt(Vt)*BrownS
      St=exp(lnS)

      Vt=Vt+kappa*(thetaT-Vt)*dt+sigma1*sqrt(Vt)*BrownV
      if(Vt<0) Vt=-Vt
      
      TT=TT+eta*(thetaConst-TT)*dt+BrownT*sqrt(TT)
      if(TT<0) TT=-TT
      thetaT=TT
    }
    
    #TODO: Vol?! Just final value?!
    simulatedReturn[i] = log(St/Forward)
    simulatedVol[i] = Vt
  }
  res <- data.frame(simulatedReturn, simulatedVol)
  return(res)
}  

plotMatrix <- function(xvec, matrix, xlab, ylab, maintitle, subtitle, plotVIX = FALSE) { 
  #determine min and max for yaxis
  min <- matrix[1,1]
  max <- matrix[1,1]
  
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      if (matrix[i,j] > max) {max = matrix[i,j]}
      if (matrix[i,j] < min) {min = matrix[i,j]}
    }
  }
  
  #need that for loess interpolation (plot Q20)
  numbers <- 1:length(xvec)
  
  if (isTRUE(plotVIX)) {
    for (i in 1:length(VIX)) {
      if (VIX[i] > max) {max = VIX[i]}
      if (VIX[i] < min) {min = VIX[i]}
    }
  }
  
  if (max >= 0) {max = max * 1.05} else {max = max * 0.95}
  if (min >= 0) {min = min * 0.95} else {min = min * 1.05}
  
  if (missing(maintitle)) {
    plot(xvec, matrix[1, ], type = "n", xlab = xlab, ylab = ylab, ylim = c(min, max))
    
  } else if (missing(subtitle)) {
    plot(xvec, matrix[1, ], type = "n", xlab = xlab, ylab = ylab, ylim = c(min, max), main = maintitle)
    
  } else {
    plot(xvec, matrix[1, ], type = "n", xlab = xlab, ylab = ylab, ylim = c(min, max), main = maintitle, sub = subtitle)
  }

  for (i in 1:nrow(matrix)) {
    lw <- loess(matrix[i, ] ~ numbers)
    #Alpha = 1: dark
    lines(xvec, lw$fitted, col = rgb(0,0,0, alpha = (i / nrow(matrix))), lwd = 2)
  }
  
  if(isTRUE(plotVIX)) {
    lines(xvec, VIX, col = rgb(1,0,0), lwd = 2)
  }
  
}

################  Q3  ################  

#We add 3 columns to the given data:
# - implied vols
# - time to maturity (expressed in years)
# - moneyness (ln(Strike/Forward))
volvec<-vector()
ttmvec<-vector()
moneynessvec<-vector()
length(volvec) = nrow(raw)
length(ttmvec) = nrow(raw)
length(moneynessvec) = nrow(raw)

for (i in 1:nrow(raw)) {
  #Important: Discount Forward back (using interest rate and dividend yield) to get Spot (just using Forward does not work)
  #Define f as the difference between the BS price with implied vol and the observed price
  f<-function (x) BSprice(1, raw$Forward[i]*exp(-(raw$InterestRate[i]-raw$DividendYield[i])*(raw$Term[i]/12)), raw$Strike[i], x, raw$DividendYield[i], raw$InterestRate[i], raw$Term[i]/12) - raw$CallPrice[i]
  
  #Find the root of the function f (i.e. the implied vol that makes the BS price equal to the observed price)
  #vol cannot be lower than 0 (lower boundary) - and is most likely not higher than 3 (arbitrarily chosen)
  voli<-uniroot(f, lower=0, upper=3)[[1]]
  ttm<-raw$Term[i]/12
  moneyness<-log(raw$Strike[i]/raw$Forward[i])
  
  #assign results to vectors
  volvec[i]<-voli
  ttmvec[i]<-ttm
  moneynessvec[i]<-moneyness
}

#Enrich given data with new columns
raw["ImpliedVol"]<-volvec
raw["TTM"]<-ttmvec
raw["Moneyness"]<-moneynessvec

#Plot implied vol surface (w/o inter/extrapolation) for a given date (here: 2006-01-31)
#install.packages("plotly")
plotdata <- raw[raw[1] == "2009-01-30", ]
p<-plot_ly(plotdata, x=plotdata$TTM, y=plotdata$Strike, z=plotdata$ImpliedVol) %>%
  add_markers() %>%
  layout(title = "Implied Vol Surface - 2009-01-30",
         scene = list(xaxis = list(title = 'Time to Maturity'),
                      yaxis = list(title = 'Strike'),
                      zaxis = list(title = 'Implied Volatility')))
p
Sys.setenv("plotly_username"="tneuber")
Sys.setenv("plotly_api_key"="7CknaAatVziORktIj116")
#chart_link = plotly_POST(p, filename="ImpliedVol")
#chart_link

################  Q4+6+9  ################  

#Vectors for Q4 plot
moneynessPlotVec <- c()
impliedVolPlotVec <- c()
termPlotVec <- c()

moneynessLogPlotVec <- c()
impliedLogVolPlotVec <- c()

#Choose 4 dates for Q9: 2006-01-31, 2007-01-31, 2008-01-31, 2009-01-30
#These vectors are used for Q9:
# - meanvec is a vector of means of log excess returns that were calculated in Q6 
# - varvec is a vector of annualised var of log ex returns (needed for Q9)
# - termvec is a vector of time to maturities in months
mean06vec<-vector()
var06vec<-vector()
mean07vec<-vector()
var07vec<-vector()
mean08vec<-vector()
var08vec<-vector()
mean09vec<-vector()
var09vec<-vector()
termvec<-vector()

#TODO: still needed?!
#Indicates whether the calculated implied density curve is valid, if not, the SVI interpolation is redone
valid<-FALSE

#Looping over different times to maturities. Idea: calculate the SVI curve (Q4) + RN densities (Q6) for different maturities
#Choose t fix (here: 2006-01-31)
#Note: annualisation is done in calcImplDensity function
for(term in c(1, 3, 6, 12, 24, 36, 48, 60, 84, 120)) {
  res <- calcImplDensity("2006-01-31", term)
  mean06vec<-c(mean06vec, res$mean)
  var06vec<-c(var06vec, res$var)
  
  moneynessPlotVec <- c(moneynessPlotVec, res$moneynessPlotVec)
  impliedVolPlotVec <- c(impliedVolPlotVec, res$impliedVolPlotVec)
  termPlotVec <- c(termPlotVec, res$termPlotVec)
  
  moneynessLogPlotVec <- c(moneynessLogPlotVec, res$moneynessLogPlotVec)
  impliedLogVolPlotVec <- c(impliedLogVolPlotVec, res$impliedLogVolPlotVec)
  
  #res <- calcImplDensity("2007-01-31", term)
  #mean07vec<-c(mean07vec, res$mean)
  #var07vec<-c(var07vec, res$var)
  
  #res <- calcImplDensity("2008-01-31", term)
  #mean08vec<-c(mean08vec, res$mean)
  #var08vec<-c(var08vec, res$var)
  
  #res <- calcImplDensity("2009-01-30", term)
  #mean09vec<-c(mean09vec, res$mean)
  #var09vec<-c(var09vec, res$var)
  
  #termvec<-c(termvec, term)
}

#Careful: LogPlotVecs do NOT give log moneyness - other way around!
q4frame <- data.frame(termPlotVec, moneynessPlotVec, impliedVolPlotVec)
scatterplot3d(q4frame, cex.symbols = 0.3, xlab = "Time to Maturity (months)", ylab = "Moneyness [log(K/F)]", zlab = "Implied Vol", main = "Implied Vols - 2006-01-31", highlight.3d = TRUE, box = FALSE)

q4Logframe <- data.frame(termPlotVec, moneynessLogPlotVec, impliedLogVolPlotVec)
scatterplot3d(q4Logframe, cex.symbols = 0.3, xlab = "Time to Maturity (months)", ylab = "Moneyness [K/F]", zlab = "Implied Vol", main = "Implied Vols - 2006-01-31", highlight.3d = TRUE, box = FALSE)


q9frame<-data.frame(termvec, mean06vec, var06vec, mean07vec, var07vec, mean08vec, var08vec, mean09vec, var09vec)
ggplot(q9frame, aes(termvec, y = value, color = variable)) +
  xlab("Time to Maturity (months)") +
  ylab("Value") +
  geom_smooth(aes(y = mean06vec, col = "Returns 2006-01-31"), se = FALSE) +
  geom_smooth(aes(y = var06vec, col = "Vol 2006-01-31"), se = FALSE) +
  geom_smooth(aes(y = mean07vec, col = "Returns 2007-01-31"), se = FALSE) +
  geom_smooth(aes(y = var07vec, col = "Vol 2007-01-31"), se = FALSE) +
  geom_smooth(aes(y = mean08vec, col = "Returns 2008-01-31"), se = FALSE) +
  geom_smooth(aes(y = var08vec, col = "Returns 2008-01-31"), se = FALSE) +
  geom_smooth(aes(y = mean09vec, col = "Returns 2009-01-30"), se = FALSE) +
  geom_smooth(aes(y = var09vec, col = "Returns 2008-01-30"), se = FALSE) 
  #geom_point(aes(y = mean06vec, col = "Returns 2006-01-31")) +
  #geom_point(aes(y = var06vec, col = "Vol 2006-01-31")) +
  #geom_point(aes(y = mean07vec, col = "Returns 2007-01-31")) +
  #geom_point(aes(y = var07vec, col = "Vol 2007-01-31")) +
  #geom_point(aes(y = mean08vec, col = "Returns 2008-01-31")) +
  #geom_point(aes(y = var08vec, col = "Vol 2008-01-31")) +
  #geom_point(aes(y = mean09vec, col = "Returns 2008-01-30")) +
  #geom_point(aes(y = var09vec, col = "Vol 2009-01-30")) +
  #geom_smooth(aes(y = mean06vec, col = "Returns 2006-01-31"), se = FALSE) +
  #geom_smooth(aes(y = var06vec, col = "Vol 2006-01-31"), se = FALSE) +
  #geom_smooth(aes(y = mean07vec, col = "Returns 2007-01-31"), se = FALSE) +
  #geom_smooth(aes(y = var07vec, col = "Vol 2007-01-31"), se = FALSE) +
  #geom_smooth(aes(y = mean08vec, col = "Returns 2008-01-31"), se = FALSE) +
  #geom_smooth(aes(y = var08vec, col = "Returns 2008-01-31"), se = FALSE) +
  #geom_smooth(aes(y = mean09vec, col = "Returns 2009-01-30"), se = FALSE) +
  #geom_smooth(aes(y = var09vec, col = "Returns 2008-01-30"), se = FALSE) 


################  Q11  ################  

#downloading VIX

t<-c(unique(raw$ValuationDate))

get(getSymbols("^VIX",src="yahoo",from="2006-01-31",to="2009-12-31"))
VIX<-VIX[,4] #getting closing prices
# to get the values of VIX for the 48 dates from raw_data, use this
VIX<-as.numeric(VIX[which((index(VIX)) %in% t)]/100)
# where t is the vector of the 48 dates

#These vectors are used for Q11:
# - meanTTM12vec is a vector of means of log excess returns with a fixed TTM of 12m for all 48 dates
# - varvec is a vector of annualised vars of log ex returns with a fixed TTM of 12m for all 48 dates
meanTTM12vec<-vector()
varTTM12vec<-vector()
integral<-0
valid<-FALSE

#Looping over different dates 
#Choose term fix as 12
for(date in t) {
  print(date)
  sink(tempfile())
  res <- calcImplDensity(date, 12)
  sink()
  meanTTM12vec<-c(meanTTM12vec, res$mean)
  varTTM12vec<-c(varTTM12vec, res$var)
}

varTTM12vecSCALED<-vector()
meanTTM12vecSCALED<-vector()

#meanTTM12vec<-c(meanTTM12vec, 0)
#varTTM12vec<-c(varTTM12vec, 0)

for(i in 1:48) {
  varTTM12vecSCALED[i] = varTTM12vec[i]*exp(1)
  meanTTM12vecSCALED[i] = meanTTM12vec[i]*(-10)
}

TTM12frame<-data.frame(t, VIX, meanTTM12vec, varTTM12vec)
ggplot(TTM12frame, aes(t,y = value, color = variable)) + 
  geom_point(aes(y = VIX, col = "VIX")) + 
  geom_point(aes(y = meanTTM12vec, col = "mean")) + 
  geom_point(aes(y = varTTM12vec, col = "var")) +
  geom_smooth(aes(y = VIX, col = "VIX")) +
  geom_smooth(aes(y = meanTTM12vec, col = "mean")) + 
  geom_smooth(aes(y = varTTM12vec, col = "var")) 

TTM12frameSCALED<-data.frame(t, VIX, meanTTM12vecSCALED, varTTM12vecSCALED)
ggplot(TTM12frame, aes(t,y = value, color = variable)) + 
  geom_point(aes(y = VIX, col = "VIX")) + 
  geom_point(aes(y = meanTTM12vecSCALED, col = "Return (times (-10))")) + 
  geom_point(aes(y = varTTM12vecSCALED, col = "Vol")) +
  geom_smooth(aes(y = VIX, col = "VIX")) +
  geom_smooth(aes(y = meanTTM12vecSCALED, col = "Return (times (-10))")) + 
  geom_smooth(aes(y = varTTM12vecSCALED, col = "Vol")) 


################  Q13  ################  

varmatrix<-matrix(nrow=48, ncol=0)

for(term in c(1, 3, 6, 12, 24, 36, 48, 60, 84, 120)) {
  meanvec<-vector()
  varvec<-vector()
  for(date in t) {
    print(term)
    sink(tempfile())
    res <- calcImplDensity(date, term)
    sink()
    meanvec<-c(meanvec, res$mean)
    varvec<-c(varvec, res$var)
  }
  varmatrix<-cbind(varmatrix, varvec)
}

## here you will need to construct calculated variance from risk neutral density in the following format
## column : time to maturity; row: valuationdate
## if Q12 is the resulting matrix, to do the PCA, use
pz <- prcomp(varmatrix, tol = 0.1,na.rm=T)
summary(pz)

################  Q17 + Q18 ################  

t<-c(unique(raw$ValuationDate))
terms <- c(unique(raw$Term))

#min 100
nsims <- 100

#meanTermStructureSeries <- matrix(nrow = length(t), ncol = length(terms))
#volTermStructureSeries <- matrix(nrow = length(t), ncol = length(terms))
meanTermStructureSeries <- matrix(ncol = length((terms)))
volTermStructureSeries <- matrix(ncol = length((terms)))

mean1PercQuantTermStructureSeries <- matrix(ncol = length((terms)))
vol1PercQuantTermStructureSeries <- matrix(ncol = length((terms)))

mean5PercQuantTermStructureSeries <- matrix(ncol = length((terms)))
vol5PercQuantTermStructureSeries <- matrix(ncol = length((terms)))

for(dateIndex in 1:length(t)) {
  date <- t[dateIndex]
  
  meanTermStructure <- vector()
  volTermStructure <- vector()
  
  mean1PercTermStructure <- vector()
  mean5PercTermStructure <- vector()
  
  vol1PercTermStructure <- vector()
  vol5PercTermStructure <- vector()
  
  termvec <- vector()
  print(paste(dateIndex,"out of", length(t)))
  
  for(term in terms) {
    res <- monteCarlo(date, term, nsims)
    
    meanTermStructure <- c(meanTermStructure, average(res$simulatedReturn))
    volTermStructure <- c(volTermStructure, average(res$simulatedVol))
    
    sortedReturn <- sort(res$simulatedReturn)
    return1PercQuantile <- sortedReturn[round(nsims*0.01, digits = 0)]
    return5PercQuantile <- sortedReturn[round(nsims*0.05, digits = 0)]
    
    sortedVol <- sort(res$simulatedVol)
    vol1PercQuantile <- sortedVol[round(nsims*0.01, digits = 0)]
    vol5PercQuantile <- sortedVol[round(nsims*0.05, digits = 0)]
    
    mean1PercTermStructure <- c(mean1PercTermStructure, return1PercQuantile)
    mean5PercTermStructure <- c(mean5PercTermStructure, return5PercQuantile)
    
    vol1PercTermStructure <- c(vol1PercTermStructure, vol1PercQuantile)
    vol5PercTermStructure <- c(vol5PercTermStructure, vol5PercQuantile)
    
    termvec <- c(termvec, term)
  }
  
  meanTermStructureSeries <- rbind(meanTermStructureSeries, meanTermStructure)
  volTermStructureSeries <- rbind(volTermStructureSeries, volTermStructure)
 
  mean1PercQuantTermStructureSeries <- rbind(mean1PercQuantTermStructureSeries, mean1PercTermStructure)
  vol1PercQuantTermStructureSeries <- rbind(vol1PercQuantTermStructureSeries, vol1PercTermStructure)
  
  mean5PercQuantTermStructureSeries <- rbind(mean5PercQuantTermStructureSeries, mean5PercTermStructure)
  vol5PercQuantTermStructureSeries <- rbind(vol5PercQuantTermStructureSeries, vol5PercTermStructure) 
}

#Delete first NA row
meanTermStructureSeries <- meanTermStructureSeries[2:nrow(meanTermStructureSeries),]
volTermStructureSeries <- volTermStructureSeries[2:nrow(volTermStructureSeries),]
mean1PercQuantTermStructureSeries <- mean1PercQuantTermStructureSeries[2:nrow(mean1PercQuantTermStructureSeries),]
vol1PercQuantTermStructureSeries <- vol1PercQuantTermStructureSeries[2:nrow(vol1PercQuantTermStructureSeries),]
mean5PercQuantTermStructureSeries <- mean5PercQuantTermStructureSeries[2:nrow(mean5PercQuantTermStructureSeries),]
vol5PercQuantTermStructureSeries <- vol5PercQuantTermStructureSeries[2:nrow(vol5PercQuantTermStructureSeries),]

#Plot for Q17
#TODO: add empirical data
plotMatrix(termvec, meanTermStructureSeries, "Time to Maturity", "Log Holding Period Excess Return")

################  Q19 ################  

plotMatrix(termvec, mean1PercQuantTermStructureSeries, "Time to Maturity", "Log Holding Period Excess Return", maintitle = "Question 19 - Return 0.01 Quantile", subtitle = "(More recent term structures are depicted darker)")
plotMatrix(termvec, vol1PercQuantTermStructureSeries, "Time to Maturity", "Volatility", maintitle = "Question 19 - Volatility 0.01 Quantile", subtitle = "(More recent term structures are depicted darker)")

plotMatrix(termvec, mean5PercQuantTermStructureSeries, "Time to Maturity", "Log Holding Period Excess Return", maintitle = "Question 19 - Return 0.05 Quantile", subtitle = "(More recent term structures are depicted darker)")
plotMatrix(termvec, vol5PercQuantTermStructureSeries, "Time to Maturity", "Volatility", maintitle = "Question 19 - Volatility 0.05 Quantile", subtitle = "(More recent term structures are depicted darker)")

################  Q20 ################  

plotMatrix(t, t(meanTermStructureSeries), "Date", "Mean Log Holding Period Excess Return")
#plotMatrix(t, t(volTermStructureSeries), "Date", "Volatility")

plotMatrix(t, t(mean1PercQuantTermStructureSeries), "Date", "Log Holding Period Excess Return", maintitle = "Question 20 - Return 0.01 Quantile")
plotMatrix(t, t(vol1PercQuantTermStructureSeries), "Date", "Volatility", maintitle = "Question 20 - Volatility 0.01 Quantile", plotVIX = TRUE)

plotMatrix(t, t(mean5PercQuantTermStructureSeries), "Date", "Log Holding Period Excess Return", maintitle = "Question 20 - Return 0.05 Quantile")
plotMatrix(t, t(vol5PercQuantTermStructureSeries), "Date", "Volatility", maintitle = "Question 20 - Volatility 0.05 Quantile", plotVIX = TRUE)


################  Q15  ################  

#% updated 15/03/2017

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coefficients for Term 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T2 = (kappa - rho*sigma1 + exp(kappa*(t - T))*(rho*sigma1 + kappa*(-1 + rho*sigma1*(-t + T))))/kappa^2

T3 = -(((-1 + exp(kappa*(t - T)))*rho*sigma1 + kappa*(kappa - rho*sigma1)*(t - T))/kappa^2)

T4 = (exp(kappa*T)*(eta - kappa)^2*(kappa - rho*sigma1) + exp(eta*t - eta*T + kappa*T)*kappa^2*(eta - kappa + rho*sigma1) + exp(kappa*t)*eta*(kappa*(-eta + kappa) + rho*sigma1*(eta + kappa^2*(t - T) + kappa*(-2 - eta*t + eta*T))))/(exp(kappa*T)*eta*(eta - kappa)^2*kappa)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (non-zero) coefficients for Term 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ttilde2 = (sigma1^2*(-exp(2*kappa*t) + exp(2*kappa*T) + 2*exp(kappa*(t + T))*kappa*(t - T)))/(2*exp(2*kappa*T)*kappa^3)

Ttilde3 = -((3 - 4*exp(eta*(t - T)) + exp(2*eta*(t - T)))*kappa^6*sigma2^2 + eta*kappa^5*sigma2^2*(-1 + exp(2*eta*(t - T)) + 2*kappa*t - 2*kappa*T) + eta^6*sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^5*kappa*sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^4*kappa^2*(sigma1 - sigma2)*(sigma1 + sigma2)*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - 2*eta^2*kappa^4*sigma2^2*(2 - 2*(exp(eta*(t - T)) + exp(kappa*(t - T)) - exp((eta + kappa)*(t - T))) + kappa*t - kappa*T) + eta^3*kappa^3*(sigma1^2*(3 - 4*exp(kappa*(t - T)) + exp(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) + sigma2^2*(-1 + exp(2*kappa*(t - T)) - 2*kappa*t + 2*kappa*T)))/(4*eta^3*(eta - kappa)^2*kappa^3*(eta + kappa))

Ttilde4 = (exp(-2*eta*t - 3*kappa*t - 7*eta*T - 6*kappa*T)*(2*exp(3*eta*t + 4*kappa*t + 6*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa^2*sigma2^2 - exp(4*eta*t + 3*kappa*t + 5*eta*T + 6*kappa*T)*(eta - 2*kappa)*kappa^4*sigma2^2 - exp(2*eta*t + 5*kappa*t + 7*eta*T + 4*kappa*T)*eta^3*((eta - kappa)^2*sigma1^2 - kappa^2*sigma2^2) + exp(2*eta*t + 3*kappa*t + 7*eta*T + 6*kappa*T)*(eta - 2*kappa)*(eta - kappa)^2*(eta^2*sigma1^2 + kappa^2*sigma2^2) + 2*exp(2*eta*t + 4*kappa*t + 7*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa*(-(kappa*sigma2^2) + eta*sigma1^2*(1 + (eta - kappa)*(t - T))) - 2*exp(3*(eta + kappa)*(t + 2*T))*eta*kappa^2*(-(eta*kappa*sigma1^2) + sigma2^2*(eta^2 + 2*kappa^3*(t - T) + eta*kappa*(-2 + eta*t - eta*T) + kappa^2*(2 - 3*eta*t + 3*eta*T)))))/(2*eta^3*(eta - 2*kappa)*(eta - kappa)^2*kappa^2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TSigmaSquared = T2 + 0.5*Ttilde2
TThetaBar = T3 + 0.5*Ttilde3
TTheta = T4 + 0.5*Ttilde4

Et = (sigma^2 - thetaBar)*TSigmaSquared + thetaBar*TThetaBar + (theta - thetaBar)*TTheta

################  Q16  ################  

# use Milstein scheme to disretize SDE for theta
for (i in 2:(n+1)){
  theta[,i] = (sqrt(theta[,i-1]) + sigma2/2 * sqrt(Delta) * Z3[,i-1])^2 
  - eta * (theta[,i-1] - thetabar) * Delta - sigma2^2/4 * Delta
}

# use Milstein scheme to disretize SDE for v
for (i in 2:(n+1)){
  v[,i] = (sqrt(v[,i-1]) + sigma1/2 * sqrt(Delta) * Z2[,i-1])^2 
  - kappa * (v[,i-1] - theta[,i-1]) * Delta - sigma1^2/4 * Delta
}

# create return time series using simple discretization of SDE
for (i in 2:(n+1)){
  x[,i] = x[,i-1] - v[,i-1]/2 * Delta + sqrt(Delta * v[,i-1]) * Z1[,i-1]
}
