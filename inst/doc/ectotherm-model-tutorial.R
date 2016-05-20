## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ------------------------------------------------------------------------
library(NicheMapR)

## ------------------------------------------------------------------------
a=getLoadedDLLs()
if(is.loaded("microclimate", "MICROCLIMATE", type = "FORTRAN")==TRUE){
  dyn.unload(a$MICROCLIMATE[[2]])
dyn.unload(a$ECTOTHERM[[2]])
}
library(NicheMapR)
micro<-micro_global(loc = "Townsville, Queensland", runmoist=1)
ecto<-ectotherm()

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(ecto$environ[,1:15], 12), digits = 2)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(ecto$enbal, 12), digits = 2)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(ecto$masbal[,1:10], 12))
knitr::kable(head(ecto$masbal[,11:19], 12))

## ---- fig.width=7, fig.height=5, fig.show = "hold"-----------------------
# retrieve output
metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
environ<-as.data.frame(ecto$environ) # activity, Tb and environment
enbal<-as.data.frame(ecto$enbal) # energy balance values
masbal<-as.data.frame(ecto$masbal) # mass balance value (note most missing if DEB 
#model not running)

# append dates
days<-rep(seq(1,12),24)
days<-days[order(days)]
dates<-days+metout$TIME/60/24-1 # dates for hourly output
dates2<-seq(1,12,1) # dates for daily output
metout<-cbind(dates,metout)
environ<-cbind(dates,environ)
masbal<-cbind(dates,masbal)
enbal<-cbind(dates,enbal)

# Hourly Tb (black), activity (orange, 5=bask, 10=forage), depth (brown, m) and shade 
# (green, %/10)
with(environ, plot(TC~dates,ylab="Tb, depth, activity and shade", xlab="month of year",
ylim=c(-20,70),type = "l", main = 
    "Fig. 1, Hourly Tb, depth, activity level and shade, 90% max shade"))
with(environ, points(ACT*5~dates,type = "l",col="orange"))
with(environ, points(SHADE/10~dates,type = "l",col="green"))
with(environ, points(DEP/10~dates,type = "l",col="brown"))
abline(ecto$VTMAX,0,lty=2,col='red')
abline(ecto$VTMIN,0,lty=2,col='blue')

# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
forage<-subset(environ,ACT==2)
bask<-subset(environ,ACT==1)
night<-subset(metout,ZEN==90)
day<-subset(metout,ZEN!=90)
with(night,plot(TIME/60~JULDAY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col=
    'dark blue', main = "Fig. 2 Annual activity window, 90% max shade"))
# nighttime hours
with(forage,points((TIME-1)~JULDAY,pch=15,cex=2,col='orange')) # foraging Tbs
with(bask,points((TIME-1)~JULDAY,pch=15,cex=2,col='light blue')) # basking Tbs


## ------------------------------------------------------------------------
ecto<-ectotherm(maxshades = rep(50,12))

## ---- echo=FALSE, fig.width=7, fig.height=5, fig.show = "hold"-----------
# retrieve output
metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
environ<-as.data.frame(ecto$environ) # activity, Tb and environment
enbal<-as.data.frame(ecto$enbal) # energy balance values
masbal<-as.data.frame(ecto$masbal) # mass balance value (note most missing if DEB
# model not running)

# append dates
days<-rep(seq(1,12),24)
days<-days[order(days)]
dates<-days+metout$TIME/60/24-1 # dates for hourly output
dates2<-seq(1,12,1) # dates for daily output
metout<-cbind(dates,metout)
environ<-cbind(dates,environ)
masbal<-cbind(dates,masbal)
enbal<-cbind(dates,enbal)

# Hourly Tb (black), activity (orange, 5=bask, 10=forage), depth (brown, m) and shade (green,
# %/10)
with(environ, plot(TC~dates,ylab="Tb, depth, activity and shade", xlab="month of year",
ylim=c(-20,70),type = "l", main = "Fig. 3, Hourly Tb, depth, activity level and shade, 50% max shade"))
with(environ, points(ACT*5~dates,type = "l",col="orange"))
with(environ, points(SHADE/10~dates,type = "l",col="green"))
with(environ, points(DEP/10~dates,type = "l",col="brown"))
abline(ecto$VTMAX,0,lty=2,col='red')
abline(ecto$VTMIN,0,lty=2,col='blue')

# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
forage<-subset(environ,ACT==2)
bask<-subset(environ,ACT==1)
night<-subset(metout,ZEN==90)
day<-subset(metout,ZEN!=90)
with(night,plot(TIME/60~JULDAY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col=
    'dark blue', main = "Fig. 4 Annual activity window, 50% max shade"))
with(forage,points((TIME-1)~JULDAY,pch=15,cex=2,col='orange')) # foraging Tbs
with(bask,points((TIME-1)~JULDAY,pch=15,cex=2,col='light blue')) # basking Tbs


## ---- fig.width=7, fig.height=5, fig.show = "hold"-----------------------
# run the microclimate model daily for 5 years
timeinterval<-365
nyears<-1
micro<-micro_global(loc = "Townsville, Queensland", timeinterval = timeinterval, nyears = nyears,
  runmoist = 1)

# run the ectotherm model with the DEB model turned on and in viviparous mode, simulating the
# Eastern Water Skink, Eulamprus quoyii
ecto<-ectotherm(DEB = 1, viviparous = 1)

metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
environ<-as.data.frame(ecto$environ) # activity, Tb and environment
masbal<-as.data.frame(ecto$masbal) # activity, Tb and environment
debout<-as.data.frame(ecto$debout) # activity, Tb and environment

# append dates
days<-rep(seq(1,timeinterval*nyears),24)
days<-days[order(days)]
dates<-(days+metout$TIME/60/24-1)/365 # dates for hourly output
metout<-cbind(dates,metout)
environ<-cbind(dates,environ)
masbal<-cbind(dates,masbal)
debout<-cbind(dates,debout)

# Hourly Tb (black), activity (orange, 5=bask, 10=forage), depth (brown, m) and shade 
# (green, %/10)
with(environ, plot(TC~dates,ylab="Tb, depth, activity and shade", xlab="year",
ylim=c(-20,70),type = "l", main = "Hourly Tb, depth, activity level and shade"))
with(environ, points(ACT*5~dates,type = "l",col="orange"))
with(environ, points(SHADE/10~dates,type = "l",col="green"))
with(environ, points(DEP/10~dates,type = "l",col="brown"))
abline(ecto$VTMAX,0,lty=2,col='red')
abline(ecto$VTMIN,0,lty=2,col='blue')

# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
forage<-subset(environ,ACT==2)
bask<-subset(environ,ACT==1)
night<-subset(cbind(metout,debout$DAY),ZEN==90)
day<-subset(cbind(metout,debout$DAY),ZEN!=90)
colnames(night)[20]<-"DAY"
colnames(day)[20]<-"DAY"
with(night,plot(TIME/60~DAY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col='dark blue',
  main = "Annual activity windows"))
with(bask,points((TIME-1)~DAY,pch=15,cex=2,col='light blue')) # basking Tbs
with(forage,points((TIME-1)~DAY,pch=15,cex=2,col='orange')) # foraging Tbs

with(debout, plot(WETMASS~dates,ylab="wet mass, g", xlab="year",
type = "l", main = "wet mass through time"))


