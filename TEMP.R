################  Q3  ################  

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

load(file = "fm408_exam_data.RData")

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
library("plotly")
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
epsilon<-0.05

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

## parameters for optimisation function DEoptim
## l: vector of lower bound of parameters
## u: vector of upper bound of parameters
## itermax: number of iterations
## VTR: minimum sum of squares to be reached

library("DEoptim")

l <- c(-20,0,-0.99,-20 ,1e-5)
u <- c(20,50,0.99,5,20)



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
for(term in c(120)) {
#for(term in c(12)) {
  
  valid<-FALSE
  var<--1
  while(!isTRUE(valid) || var < 0) {
    
    #TODO: make selection variable based
    #Table of options that are used to fit the SVI curve (i.e. all options that match a given t and T (but have different Strikes))
    optionsToFit<-raw[raw$ValuationDate == "2006-01-31", ]
    optionsToFit<-optionsToFit[optionsToFit$Term == term, ]
    
    #TODO: uIntegrBound = Inf does not work
    lIntegrBound <- optionsToFit$Forward[1]
    uIntegrBound <- 5000
    
    #TODO: Bug: the calculated var is sometimes negative. Repeat calcuation if that is the case until var is positive.
    #Note: The SVI interpolation is non-deterministic because the optimization routine is random - results therefore vary
    var<--1
    
    func = function(p) {
      sum = 0
      for (i in 1:nrow(optionsToFit)) {
        #Minimize the squared error (vfunc calculates squared error, sum gives the total squared error from all options used to calculate the SVI curve)
        sum =+ vfunc(p, optionsToFit$Moneyness[i], optionsToFit$TTM[i], optionsToFit$ImpliedVol[i])
      }
      return(sum)
    }
    
    #Then write the optimization function applied to vfunc
    outDEoptim <- DEoptim(func, l, u, DEoptim.control(VTR = 1e-4, itermax =  1e4))
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
    curve(implDensity, from=0.5, to=5000, xlab="K", ylab="Implied Density")
    #curve(bestFit2, from=-2, to=2, xlab="Moneyness [log(K/F)]", ylab="Implied Vol")
    
    
    #p <- c(-2.9773*(10^(-2)),0.1739,-0.84412,-0.3486967,3.1999622*(10^(-1)),7.74761*(10^(-6)))
    #bestFit2 = function(k) {implVol(p, k, optionsToFit$TTM[1])}
    
    integral = tryCatch(integrate(implDensity, lower = 0.5 + 1, upper = 10000), error = function(e) e)  
    if(inherits(integral, "error")) next
    if(integral$value - integral$abs.error - epsilon <= 1 && 1 <= integral$value + integral$abs.error + epsilon) {
      valid<-TRUE
      #TODO: Check this.
      #use the implied density to get the implied log excess return density. 
      #Problem: doing this leads to function that does not yield 1 when integrated form -Inf to Inf (implDensity does not do that either)
      #Dirty trick: scale the resulting function so that the probabilities sum up to 1
      #Problem 2: Integrating from -Inf to Inf leads to strange results because the function behaves strange from extrem values
      #Even dirtier trick: only integrate from lIntegrBound to uIntegrBound. Chose this boundaries as the min/max of log moneyness we have data for (does that make sense?)
      returnImplDensityUnscaled = function(ST) {
        res<-log(ST/optionsToFit$Forward[1])*implDensity(ST)
        return(res)
      } 
      #returnImplDensity = function(r) {return(returnImplDensityUnscaled(r)/integrate(returnImplDensityUnscaled, lower=lIntegrBound,upper=uIntegrBound)$value)}
    
      #curve(returnImplDensity, from=0.5, to=5000, xlab="R", ylab="Implied Density")
      curve(returnImplDensityUnscaled, from=0.5, to=5000, xlab="R", ylab="Implied Log Returns")
      
      
      #calculate the mean of the function by integration, and the var by using the usual formula
      mean = tryCatch(integrate(function(x) {return(returnImplDensityUnscaled(x))}, 0.5, 10000), error = function(e) e)
      if(inherits(mean, "error")) next
      
      varHelp = tryCatch(integrate(function(x) {return(returnImplDensityUnscaled(x)^2)}, 0.5, 10000), error = function(e) e)
      if(inherits(varHelp, "error")) next
      
      mean = mean$value
      var = varHelp$value/(term/12) - mean^2  
    } else {
      valid<-FALSE
    }
    
  }
  meanvec<-c(meanvec, mean)
  varvec<-c(varvec, var)
  termvec<-c(termvec, term)
}

varframe<-data.frame(termvec, varvec)
ggplot(varframe, aes(termvec,varvec)) + geom_point() + geom_smooth()
################  Q11  ################  

