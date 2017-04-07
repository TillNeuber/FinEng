################  CONSTANT DECLARATION, LOAD PACKAGES, LOAD DATA  ################  

# model assumptions
#should we replace with something?
kappa=1 # rate of reversion to "current target variance theta(t)"
thetaConst=0.1 # long-term variance
sigma1=0.1 # volatility of volatility 
sigma2=0.1 # volatility of theta(t)
rho=0.2 # weight factor
eta=0.2 # rate of reversion to long-term target variance thetaConst


# Option Assumptions
S=100 # stock price today
#S=exp(-raw$InterestRate[1]+raw$DividendYield[1])*raw$Forward[1]
V=0.5 # volatility of Future
K=90 # option strike price
#K=raw$Strike[1]
Rf=0.5 # risk-free rate; Assumed constant
months=12 # months until maturity
Div=0 #dividend yield
F=S*(1+Rf)^(months/12)

# Initialize variables
#dt is monthly

iter=5000
allS=rep(1,months)
dt=1/12 

################  Utility Functions  ################  

average <- function(list) {
  sum <- 0
  for (i in 1:length(list)) {
    sum = sum + list[i]
  }
  return(sum/length(list))
}

monteCarlo <- function() {     
  for (j in 2:months) {
    
    #brownian of stock
    BrownS=qnorm(runif(1,0,1),0,1)
    #Brownian of variance
    BrownV=rho*BrownS+sqrt(1-rho^2)*qnorm(runif(1,0,1),0,1)
    #brownian of theta
    BrownT=sigma2*qnorm(runif(1,0,1),0,1)
    
    lnS=lnS+(Rf-Div-0.5*Vt)*dt+sqrt(Vt)*BrownS
    St=exp(lnS)
    allS[j]=St
    
    V=V+kappa*(thetaT-Vt)*dt+sigma1*sqrt(Vt)*BrownV
    if(V<0) V=-V
    Vt=V
    
    TT=TT+eta*(thetaConst-TT)*dt+BrownT*sqrt(TT)
    if(TT<0) TT=-TT
    thetaT=TT
  }
  return(c(log(St/F), Vt))
}  

################  Q17  ################  

simulatedReturn = rep(1,iter)
simulatedVol = rep(1,iter)

### Calculation
for (i in 1:iter) {
  lnS=log(S)
  St=S
  Vt<-V
  #not sure it is the right inizialization
  #I made a new variable TT to prevent 
  thetaT<-thetaConst
  TT<-thetaT  
  res <- monteCarlo()
  simulatedReturn[i] = res[1]
  simulatedVol[i] = res[2]
}

mean <- average(simulatedReturn)
vol <- average(simulatedVol)

################  Q18  ################  

sortedReturn <- sort(simulatedReturn)
return1PercQuantile <- sortedReturn[round(iter*0.01, digits = 0)]
return5PercQuantile <- sortedReturn[round(iter*0.05, digits = 0)]

sortedVol <- sort(simulatedVol)
vol1PercQuantile <- sortedVol[round(iter*0.01, digits = 0)]
vol5PercQuantile <- sortedVol[round(iter*0.05, digits = 0)]



####### Testcode ######

sum = 0
for (i in 1:iter) {
  sum = sum + simulatedVol[i]
}
sum/iter

ModelM=simulpath[iter]/iter
ModelM


