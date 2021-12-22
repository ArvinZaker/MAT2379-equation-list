#load this package in R and enjoy!
#does a t.test if datapoints are known
t.test2 <- function(m1,m2,s1,s2,n1,n2,mu=0,var.equal=FALSE)
{
  if( var.equal==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-mu)/se 
  dat <- c(m1-m2, se, t, df, 2*pt(-abs(t),df), pt(-abs(t),df) )    
  names(dat) <- c("Difference of means", "Std Error", "t", "df", "p-value two tail", "p-value one tail")
  return(dat) 
}
#sinds pooled sd if mean and t are known
spcalc.t.mean <- function(t,x1,x2,n1,n2)
{
  return((x1-x2)/(t*(sqrt(1/n1+1/n2))))
}
#calculates sp for you if s1 and s2 are known
spcalc <- function(s1,s2,n1,n2)
{
  if(n1<30 && n2<30)
  {
    HH <- ((n1-1)*(s1^2))+((n2-1)*(s2^2))
    return(sqrt(HH/(n1+n2-2)))
  } else
  {
    return("calculation of Sp not applicable! sample is too large!")
  }
}

# calculates v for you: welsh t.test degree of freedom
vcalc <- function(s1,s2,n1,n2)
{
  if(n1<31 && n2<31)
  {
    U=((s1^2/n1)+(s2^2/n2))^2
    D=((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
    V=U/D
    return(V)
  } else
  {
    return("calculation of v not applicable! sample is too large!")
  }
}

#calc CI for z.test, use CI value between 0 and 1
cicalc.ztest <- function(x,s,n,CI=1)
{
  o=(1-CI)/2
  a=-abs(qnorm(o))
  SE=s/sqrt(n)
  dat <- c(CI*100,x+a*SE,x-a*SE)
  names(dat) <- c("CI %", "the interval")
  if(n>30)
  {
    return(dat)
  } else
  {
    names(dat) <- c("CI %", "the interval", "N IS TOO SMALL!")
    return(dat)
  }
}

#calc CI for t.test, use CI value between 0 and 1
cicalc.ttest <- function(x,s,n,CI=1)
{
  o=(1-CI)/2
  b=-abs(qt(o,n-1))
  SE=s/sqrt(n)
  dat <- c(CI*100,x+b*SE,x-b*SE)
  names(dat) <- c("CI %", "the interval")
  if(n<30)
  {
    return(dat)
  } else
  {
    names(dat) <- c("CI %", "the interval", "N IS TOO LARGE!")
    return(dat)
  }
}

#calc CI for proportions, use CI value between 0 and 1
cicalc.prop <- function(p,n,CI=1)
{
  o=(1-CI)/2
  a=-abs(qnorm(o))
  SE=sqrt(p*(1-p)/n)
  dat <- c(CI*100,p+a*SE,p-a*SE)
  names(dat) <- c("CI %", "the interval")
  if(n*p>10 && n*(1-p)>10)
  {
    return(dat)
  } else
  {
    names(dat) <- c("CI %", "the interval", "N IS TOO SMALL!")
    return(dat)
  }
}
#calc t.test
t.test1 <- function(x,s,n,mu=0)
{
  t0=(x-mu)/(s/sqrt(n))
  pv=c(pt(t0,n-1),pt(-t0,n-1))
  dat <- c(t0, n-1, min(pv), 2*min(pv))
  names(dat) <- c("t0", "df", "p-value one sided","p-value two sided")
  return(dat)
}
#calc z.test
z.test1 <- function(x,s,n,mu=0)
{
  z0=(x-mu)/(s/sqrt(n))
  pv=c(pnorm(z0),pnorm(-z0))
  dat <- c(z0, min(pv), 2*min(pv))
  names(dat) <- c("z0","p-value one sided","p-value two sided")
  return(dat)
}
#calc prop
p.test <- function(p,n,p0=0)
{
  z0=(p-p0)/(sqrt(p0*(1-p0)/n))
  pv=c(pnorm(z0),pnorm(-z0))
  dat <- c(z0, min(pv), 2*min(pv))
  names(dat) <- c("z0","p-value one sided","p-value two sided")
  return(dat)
}
#calc 2 variable z test if values known
z.test2 <- function(x1,x2,s1,s2,n1,n2,mu=0)
{
  sdt=sqrt((s1^2)/n1+(s2^2)/n2)
  z0=(x1-x2-mu)/sdt
  pv=c(pnorm(z0),pnorm(-z0))
  dat <- c(z0, min(pv), 2*min(pv))
  names(dat) <- c("z0","p-value one sided","p-value two sided")
  return(dat)
}
#confidence interval 2 value t.test
cicalc.ttest2 <- function(m1,m2,s1,s2,n1,n2,CI=1,var.equal=FALSE)
{
  o=(1-CI)/2
  if( var.equal==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  b=abs(qt(o,df)) 
  dat <- c(df,CI*100,m1-m2-b*se,m1-m2+b*se)    
  names(dat) <- c("df", "CI", "the interval")
  return(dat) 
}
#calc ci of 2 variable z test
cicalc.ztest2 <- function(x1,x2,s1,s2,n1,n2,CI=1)
{
  o=(1-CI)/2
  a=abs(qnorm(o))
  SE=sqrt( (s1^2/n1) + (s2^2/n2) )
  dat <- c(CI*100,x1-x2-a*SE,x1-x2+a*SE)
  names(dat) <- c("CI %", "the interval")
  if(n1>30 && n2>30)
  {
    return(dat)
  } else
  {
    names(dat) <- c("CI %", "the interval", "N IS TOO SMALL!")
    return(dat)
  }
}