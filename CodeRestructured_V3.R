################  CONSTANT DECLARATION, LOAD PACKAGES, LOAD DATA  ################  

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
#l <- c(-2,0,-0.99,-2 ,1e-5)
#u <- c(2,2,.99,1,2)


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
    #date<-"2006-01-31"
    #term<-12
    #Table of options that are used to fit the SVI curve (i.e. all options that match a given t and T (but have different Strikes))
    optionsToFit<-raw[raw$ValuationDate == date, ]
    optionsToFit<-optionsToFit[optionsToFit$Term == term, ]
    
    
    #TODO: uIntegrBound = Inf does not work
    uIntegrBound <- 10000
    
    #TODO: Bug: the calculated var is sometimes negative. Repeat calcuation if that is the case until var is positive.
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
    outDEoptim <- DEoptim(func, l, u, DEoptim.control(VTR = 1e-4, itermax =  2e4))
    summary(outDEoptim)
    
    #The SVI curve. Note:k is the log moneyness, not the strike
    bestFit = function(k) {implVol(outDEoptim$optim$bestmem, k, optionsToFit$TTM[1])}
    
    #curve(bestFit, from=-1, to=3, xlab="Moneyness [log(K/F)]", ylab="Implied Vol")
    
    C = function(K) {BSprice(1, optionsToFit$Forward[1]*exp(-(optionsToFit$InterestRate[1] - optionsToFit$DividendYield[1])*optionsToFit$TTM[1]), K, bestFit(log(K/optionsToFit$Forward[1])), optionsToFit$DividendYield[1], optionsToFit$InterestRate[1], optionsToFit$TTM[1])}
    
    #The implied density as a function of the strike. Note: Not the implied return density
    implDensity = function(K) {
      Dt = exp(-optionsToFit$TTM[1]*(optionsToFit$InterestRate[1]))
      secDeriv = secondDerivative(C, K, 0.5)
      return(secDeriv / Dt)
    }
    #curve(implDensity, from=0.5, to=5000, xlab="K", ylab="Implied Density")
    
    integral = tryCatch(integrate(implDensity, lower = 0.5 + 1, upper = 10000), error = function(e) e)  
    if(inherits(integral, "error")) next
    
    implDensityScaled = function(K) {return(implDensity(K) / integral$value)}
    
    
    returnImplDensity = function(ST) {
      res<-log(ST/optionsToFit$Forward[1])*implDensity(ST)
      return(res)
    } 
    #returnImplDensity = function(r) {return(returnImplDensityUnscaled(r)/integrate(returnImplDensityUnscaled, lower=lIntegrBound,upper=uIntegrBound)$value)}
      
    #curve(returnImplDensity, from=lIntegrBound, to=3, xlab="R", ylab="Implied Density")
    #curve(returnImplDensityUnscaled, from=0.5, to=5000, xlab="R", ylab="Implied Log Returns")
      
      
    #calculate the mean of the function by integration, and the var by using the usual formula
    mean = tryCatch(integrate(function(x) {return(returnImplDensity(x))}, 0.5, 10000), error = function(e) e)
    if(inherits(mean, "error")) next
      
    varHelp = tryCatch(integrate(function(x) {return(returnImplDensity(x)^2)}, 0.5, 10000), error = function(e) e)
    if(inherits(varHelp, "error")) next
    
    mean = mean$value
    var = varHelp$value/(term/12) - mean^2  
  
  }
  return(c(mean, var))
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
plotdata <- raw[raw[1] == "2006-01-31", ]
p<-plot_ly(plotdata, x=plotdata$TTM, y=plotdata$Moneyness, z=plotdata$ImpliedVol) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Time to Maturity'),
                      yaxis = list(title = 'Moneyness [log(K/F)]'),
                      zaxis = list(title = 'Implied Volatility')))
p
Sys.setenv("plotly_username"="tneuber")
Sys.setenv("plotly_api_key"="7CknaAatVziORktIj116")
#chart_link = plotly_POST(p, filename="ImpliedVol")
#chart_link

################  Q4+6  ################  

#These vectors are used for Q9:
# - meanvec is a vector of means of log excess returns that were calculated in Q6 
# - varvec is a vector of annualised var of log ex returns (needed for Q9)
# - termvec is a vector of time to maturities in months
meanvec<-vector()
varvec<-vector()
termvec<-vector()

#Indicates whether the calculated implied density curve is valid, if not, the SVI interpolation is redone
valid<-FALSE

#Looping over different times to maturities. Idea: calculate the SVI curve (Q4) + RN densities (Q6) for different maturities
#In the current version: choose t fix (here: 2006-01-31)
#for(term in c(1, 3, 6, 12, 24, 36, 48, 60, 84, 120)) {
for(term in c(1)) {
  res <- calcImplDensity("2006-01-31", term)
  meanvec<-c(meanvec, res[1])
  varvec<-c(varvec, res[2])
  termvec<-c(termvec, term)
}

varframe<-data.frame(termvec, varvec)
ggplot(varframe, aes(termvec,varvec)) + geom_point() + geom_smooth()

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
  meanTTM12vec<-c(meanTTM12vec, res[1])
  varTTM12vec<-c(varTTM12vec, res[2])
}

meanTTM12vec<-c(meanTTM12vec, 0)
varTTM12vec<-c(varTTM12vec, 0)

for(i in 1:48) {
  varTTM12vec[i] = varTTM12vec[i]*100
}

TTM12frame<-data.frame(t, VIX, meanTTM12vec, varTTM12vec)
ggplot(TTM12frame, aes(t,y = value, color = variable)) + 
  geom_point(aes(y = VIX, col = "VIX")) +
  geom_point(aes(y = meanTTM12vec, col = "mean")) +
  geom_point(aes(y = varTTM12vec, col = "var")) 


################  Q13  ################  

varmatrix<-matrix(nrow=48, ncol=0)

for(term in c(1, 3, 6, 12, 24, 36, 48, 60, 84, 120)) {
  meanvec<-vector()
  varvec<-vector()
  for(date in t) {
    res <- calcImplDensity(date, term)
    meanvec<-c(meanvec, res[1])
    varvec<-c(varvec, res[2])
  }
  varmatrix<-cbind(varmatrix, varvec)
}

## here you will need to construct calculated variance from risk neutral density in the following format
## column : time to maturity; row: valuationdate
## if Q12 is the resulting matrix, to do the PCA, use
pz <- prcomp(varmatrix, tol = 0.1,na.rm=T)
summary(pz)


################  Q15  ################  

#% updated 15/03/2017

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coefficients for Term 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T2 = (kappa - rho*sigma1 + E^(kappa*(t - T))*(rho*sigma1 + kappa*(-1 + rho*sigma1*(-t + T))))/kappa^2

T3 = -(((-1 + E^(kappa*(t - T)))*rho*sigma1 + kappa*(kappa - rho*sigma1)*(t - T))/kappa^2)

T4 = (E^(kappa*T)*(eta - kappa)^2*(kappa - rho*sigma1) + E^(eta*t - eta*T + kappa*T)*kappa^2*(eta - kappa + rho*sigma1) + E^(kappa*t)*eta*(kappa*(-eta + kappa) + rho*sigma1*(eta + kappa^2*(t - T) + kappa*(-2 - eta*t + eta*T))))/(E^(kappa*T)*eta*(eta - kappa)^2*kappa)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (non-zero) coefficients for Term 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ttilde2 = (sigma1^2*(-E^(2*kappa*t) + E^(2*kappa*T) + 2*E^(kappa*(t + T))*kappa*(t - T)))/(2*E^(2*kappa*T)*kappa^3)

Ttilde3 = -((3 - 4*E^(eta*(t - T)) + E^(2*eta*(t - T)))*kappa^6*sigma2^2 + eta*kappa^5*sigma2^2*(-1 + E^(2*eta*(t - T)) + 2*kappa*t - 2*kappa*T) + eta^6*sigma1^2*(3 - 4*E^(kappa*(t - T)) + E^(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^5*kappa*sigma1^2*(3 - 4*E^(kappa*(t - T)) + E^(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - eta^4*kappa^2*(sigma1 - sigma2)*(sigma1 + sigma2)*(3 - 4*E^(kappa*(t - T)) + E^(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) - 2*eta^2*kappa^4*sigma2^2*(2 - 2*(E^(eta*(t - T)) + E^(kappa*(t - T)) - E^((eta + kappa)*(t - T))) + kappa*t - kappa*T) + eta^3*kappa^3*(sigma1^2*(3 - 4*E^(kappa*(t - T)) + E^(2*kappa*(t - T)) + 2*kappa*t - 2*kappa*T) + sigma2^2*(-1 + E^(2*kappa*(t - T)) - 2*kappa*t + 2*kappa*T)))/(4*eta^3*(eta - kappa)^2*kappa^3*(eta + kappa))

Ttilde4 = (E^(-2*eta*t - 3*kappa*t - 7*eta*T - 6*kappa*T)*(2*E^(3*eta*t + 4*kappa*t + 6*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa^2*sigma2^2 - E^(4*eta*t + 3*kappa*t + 5*eta*T + 6*kappa*T)*(eta - 2*kappa)*kappa^4*sigma2^2 - E^(2*eta*t + 5*kappa*t + 7*eta*T + 4*kappa*T)*eta^3*((eta - kappa)^2*sigma1^2 - kappa^2*sigma2^2) + E^(2*eta*t + 3*kappa*t + 7*eta*T + 6*kappa*T)*(eta - 2*kappa)*(eta - kappa)^2*(eta^2*sigma1^2 + kappa^2*sigma2^2) + 2*E^(2*eta*t + 4*kappa*t + 7*eta*T + 5*kappa*T)*eta^2*(eta - 2*kappa)*kappa*(-(kappa*sigma2^2) + eta*sigma1^2*(1 + (eta - kappa)*(t - T))) - 2*E^(3*(eta + kappa)*(t + 2*T))*eta*kappa^2*(-(eta*kappa*sigma1^2) + sigma2^2*(eta^2 + 2*kappa^3*(t - T) + eta*kappa*(-2 + eta*t - eta*T) + kappa^2*(2 - 3*eta*t + 3*eta*T)))))/(2*eta^3*(eta - 2*kappa)*(eta - kappa)^2*kappa^2)



