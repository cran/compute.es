##=========================== computeES ==============================##

# d to es

des <- function(d, n.1, n.2) {
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  n <- n.1 + n.2
  r <- d/sqrt((d^2) + a) # to compute r from d
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- d*(pi/sqrt(3))
  var.lor <- var.d*(pi^2/3)
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}


# Formulas for computing effect sizes (d, r, log odds ratio)
# in designs with independent groups.

# Computing effect sizes d and g, independent groups
# (0) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# sd.1 (treatment standard deviation at post-test), sd.2 (comparison 
# standard deviation at post-test), n.1 (treatment), n.2 (comparison/control).


mes <- function(m.1,m.2,sd.1,sd.2,n.1, n.2) {
  s.within<-sqrt(((n.1-1)*sd.1^2+(n.2-1)*sd.2^2)/(n.1+n.2-2))
  d<-(m.1-m.2)/s.within
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}

# (1) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control).

mes2 <- function(m.1,m.2,s.pooled,n.1, n.2) {
  d<-(m.1-m.2)/s.pooled
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}

# (2) Study reported: 
# t (t-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).

tes <- function(t, n.1, n.2) { #If only total n reported just split .5
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  n <- n.1 + n.2
  r <- sqrt((t^2)/(t^2 + n-2))
  var.r <- ((1-r^2)^2)/(n-1)
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}

# (3) Study reported: 
# f (F-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).


fes <- function(f,n.1, n.2) {
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  t <- sqrt(f)
  n <- n.1 + n.2
  #r<- sqrt((t^2)/(t^2 + n-2)) # to compute r from f with 1 df
  r <- d/sqrt((d^2) + a) # to compute r from d
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}

# (4) Study reported: 
# p-value,  n.1 (treatment), n.2 (comparison/control), tail (one or two tailed?).

pes <- function(p, n.1, n.2, tail = "two") {
  n <- n.1 + n.2
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  a <- ((n.1 + n.2)^2)/(n.1*n.2) 
  if(tail == "one") {
    pxtwo<-p*2
    TINV<-qt((1-pxtwo/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
    var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
    
  }
  if(tail == "two") {
    TINV<-qt((1-p/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
    var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}



#  Pearson r to effect size.
res <- function(r, var.r = NULL, n ) { # If var.r not reported use n
  if(missing(var.r)){
    r <- r
    var.r <-(1-r^2)^2/(n-1) 
    #d<-2*r*sqrt((n-1)/(n*(1-r^2)))*abs(r)/r 
    d <- (2*r)/(sqrt(1-r^2))
    var.d <- 4*var.r/(1-r^2)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  r <- r
  var.r <- var.r
  d <- 2*r/sqrt(1-r^2)
  var.d <- 4*var.r/(1-r^2)^3 
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  out<-list(MeanDifference= c(d = d,var.d = var.d), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))  
  return(out)
}  


# Formulas for computing effect sizes in designs with independent groups
# using ANCOVA. 
# Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# s.adj (adjusted standard deviation), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.mes <- function(m.1.adj,m.2.adj,sd.adj,n.1, n.2, R, q) {
   s.within<-sd.adj/sqrt(1-R^2)
   d<-(m.1.adj-m.2.adj)/s.within
   var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
   df<- (n.1+n.2)-2 - q
   j<-1-(3/(4*df-1))
   g<-j*d
   var.g<-j^2*var.d
   a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
   r <- d/sqrt((d^2) + a)
   var.r <- (a^2*var.d)/(d^2 + a)^3
   lor <- pi*d/sqrt(3)
   var.lor <- pi^2*var.d/3
   n <- n.1 + n.2
   z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
   var.z <- 1/(n-3) 
   #z.score <- r/sqrt(var.r)
   out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
            Correlation= c(r= r, var.r = var.r),
            Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
            Fishers_z = c(z = z, var.z = var.z),
            TotalSample = c(n= n))
   return(out)
}

# Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.mes2 <- function(m.1.adj, m.2.adj, s.pooled, n.1, n.2, R, q) {
  d<-(m.1.adj-m.2.adj)/s.pooled
  var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}


# Study reported: 
# t (t-test value from ANCOVA), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.tes <- function(t, n.1, n.2, R, q) {
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)

}

# Study reported: 
# f (F-test value from ANCOVA) with independent groups, n.1 (treatment),
# n.2 (comparison/control),R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.fes<-function(f,n.1, n.2, R, q) {
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)

}

# Study reported: 
# p-value (for ONE-tailed test, from ANCOVA), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.pes <- function(p, n.1, n.2, R, q, tail = "two") {
  n <- n.1 + n.2
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  a <- ((n.1 + n.2)^2)/(n.1*n.2) 
  if(tail == "one") {
    pxtwo<-p*2
    TINV<-qt((1-pxtwo/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
    var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  if(tail == "two") {
    TINV<-qt((1-p/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
    var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}

# computing es from log odds ratio

lores <- function(lor, var.lor, n.1, n.2) { # n.1 =  tmt grp
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  df<- (n.1+n.2)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
 return(out)
}  

# compute or from proportions

propes <- function(p1, p2, n.ab, n.cd) {
  or <-(p1*(1-p2))/(p2*(1-p1))
  lor <- log(or)
  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  df<- (n.ab+n.cd)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.ab + n.cd)^2)/(n.ab*n.cd)  # not sure if this is appropriate*
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.ab+n.cd
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))
  return(out)
}


# Odds Ratio to es: if have info for 'failure' in both conditions 
# (B = # tmt failure; D = # non-tmt failure) and the sample size
# for each group (n.1 & n.0 respectively):

failes <- function(B, D, n.1, n.0) {
  A <- n.1 - B  # tmt success
  B <- B        # tmt failure
  C <- n.0 - D  # non-tmt success
  D <- D        # non-tmt failure
  p1 <- A/n.1   # proportion 1 
  p2 <- C/n.0   # proportion 2
  n.ab <-  A+B  # n of A+B
  n.cd <-  C+D  # n of C+D        
  or <- (p1 * (1 - p2))/(p2 * (1 - p1))  # odds ratio
  lor <- log(or)  # log odds ratio
  var.lor <-  1/A + 1/B + 1/C + 1/D  # variance of log odds ratio
  #var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor * sqrt(3)/pi  # conversion to d
  var.d <- 3 * var.lor/pi^2  # variance of d
  df<- (n.1 + n.0)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.0)^2)/(n.1*n.0)  # not sure if this is appropriate*
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.0
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  #z.score <- r/sqrt(var.r)
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))

  return(out)
}

# Converting Chi-squared statistic with 1 df to es

chies <- function(chi.sq,  n) {
  r <- sqrt(chi.sq/n)
  var.r <-(1-r^2)^2/(n-1) 
  d<-2*r*sqrt((n-1)/(n*(1-r^2)))*abs(r)/r 
  var.d <- 4*var.r/(1-r^2)^3
  df<- (n)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
           Correlation= c(r= r, var.r = var.r),
           Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
           Fishers_z = c(z = z, var.z = var.z),
           TotalSample = c(n= n))

  return(out)
}




