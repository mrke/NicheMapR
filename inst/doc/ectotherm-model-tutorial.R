## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ---- message=FALSE, warnings=FALSE--------------------------------------
library(NicheMapR)
longlat <- c(146.77, -19.29) # Townsville, northern Australia
micro <- micro_global(loc = longlat)
ecto <- ectotherm() # uses default settings (the Eastern Water Skink)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(ecto$environ[,c(4:15,25:27)], 13), digits = 2)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(ecto$enbal, 13), digits = 2)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(ecto$masbal[,c(1:5, 15:17)], 12))

## ---- warning=FALSE, message=FALSE---------------------------------------
Ww_g <- 40        # wet weight of animal (g)
pct_wet <- 0.2    # % of surface area acting as a free-water exchanger
alpha_min <-0.85  # minimum solar absorbtivity (dec %)
alpha_max <- 0.85 # maximum solar absorbtivity (dec %)
shape <- 3        # lizard shape
T_RB_min <- 17.5  # min Tb at which they will attempt to leave retreat
T_B_min <- 17.5   # min Tb at which leaves retreat to bask
T_F_min <- 24     # minimum Tb at which activity occurs
T_F_max <- 34     # maximum Tb at which activity occurs
T_pref <- 30      # preferred Tb (will try and regulate to this)
CT_max <- 40      # critical thermal minimum (affects choice of retreat)
CT_min <- 6       # critical thermal maximum (affects choice of retreat)
mindepth <- 2     # min depth (node, 1-10) allowed
maxdepth <- 10    # max depth (node, 1-10) allowed
shade_seek <- 1   # shade seeking?
burrow <- 1       # can it burrow?
climb <- 0        # can it climb to thermoregulate?
nocturn <- 0      # nocturnal activity
crepus <- 0       # crepuscular activity
diurn <- 1        # diurnal activity
minshades <- rep(0, 12)   # min available shade?
maxshades <- micro$MAXSHADES # max available shade?

## ---- warning=FALSE, message=FALSE---------------------------------------
ecto <- ectotherm(Ww_g = Ww_g, alpha_max = alpha_max, alpha_min = alpha_min, pct_wet = pct_wet, T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min, CT_max = CT_max, CT_min = CT_min, T_pref = T_pref, mindepth = mindepth, maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb, minshades = minshades, nocturn = nocturn, diurn = diurn, crepus = crepus, maxshades = maxshades)

# retrieve output
environ <- as.data.frame(ecto$environ) # behaviour, Tb and environment
enbal <- as.data.frame(ecto$enbal) # heat balance outputs
masbal <- as.data.frame(ecto$masbal) # mass balance outputs
metout <- as.data.frame(micro$metout) # above ground microclimate
environ <- cbind(environ,metout$SOLR) # add solar radiation for activity window plots
colnames(environ)[ncol(environ)] <- "Solar"

# append dates
days <- rep(seq(1,12),24)
days <- days[order(days)]
dates <- days+metout$TIME/60/24-1 # dates for hourly output
dates2 <- seq(1,12,1) # dates for daily output
metout <- cbind(dates,metout)
environ <- cbind(dates,environ)
masbal <- cbind(dates,masbal)
enbal <- cbind(dates,enbal)

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
with(environ, plot(TC ~ dates, ylab = "T_b, depth, activity and shade", 
                   xlab = "month of year", ylim = c(0, 50), type = "l", 
                   main = "Fig. 1 Body temperature, depth, shade and activity, 
                   90% max shade"))
with(environ, points(ACT * 5 ~ dates, type = "l", col = "orange"))
with(environ, points(SHADE / 10 ~ dates, type = "h", col = "dark green"))
with(environ, points(DEP / 10 ~ dates, type = "l",col = "brown"))
abline(T_F_max, 0, lty = 2,col = 'red')
abline(T_F_min, 0, lty = 2,col = 'blue')
abline(T_pref, 0, lty = 2,col = 'orange')
text(0.1, T_F_max + 2, "T_F_max", col = 'red')
text(0.1, T_F_min - 2, "T_F_min", col = 'blue')
text(0.1, T_pref + 2, "T_pref", col = 'orange')
legend(x = 5, y = T_F_max + 17, legend = c("T_b (°C)", "depth (cm/10)", "activity, 0 5 or 10)", "shade (%/10)"), 
       col = c("black", "brown", "orange", "dark green"), lty = rep(1, 4), bty = "n")

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
forage<-subset(environ,ACT==2)
bask<-subset(environ,ACT==1)
night<-subset(environ,Solar==0)
day<-subset(environ,Solar==0)
with(night,plot(TIME ~ DOY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col=
    'dark blue', main = "Fig. 2 Annual activity window, 90% max shade"))
# nighttime hours
with(forage,points(TIME~DOY,pch=15,cex=2,col='orange')) # foraging Tbs
with(bask,points(TIME~DOY,pch=15,cex=2,col='light blue')) # basking Tbs

## ---- echo=FALSE, fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
micro <- micro_global(loc = longlat, maxshade = 10)
maxshades <- micro$MAXSHADES

ecto <- ectotherm(Ww_g = Ww_g, alpha_max = alpha_max, alpha_min = alpha_min, pct_wet = pct_wet, T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min, CT_max = CT_max, CT_min = CT_min, T_pref = T_pref, mindepth = mindepth, maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb, minshades = minshades, nocturn = nocturn, diurn = diurn, crepus = crepus, maxshades = maxshades)

# retrieve output
environ <- as.data.frame(ecto$environ) # behaviour, Tb and environment
enbal <- as.data.frame(ecto$enbal) # heat balance outputs
masbal <- as.data.frame(ecto$masbal) # mass balance outputs
metout <- as.data.frame(micro$metout) # above ground microclimate
environ <- cbind(environ,metout$SOLR) # add solar radiation for activity window plots
colnames(environ)[ncol(environ)] <- "Solar"

# append dates
days <- rep(seq(1,12),24)
days <- days[order(days)]
dates <- days+metout$TIME/60/24-1 # dates for hourly output
dates2 <- seq(1,12,1) # dates for daily output
metout <- cbind(dates,metout)
environ <- cbind(dates,environ)
masbal <- cbind(dates,masbal)
enbal <- cbind(dates,enbal)

with(environ, plot(TC ~ dates, ylab = "T_b, depth, activity and shade", xlab = "month of year", ylim = c(0, 50), type = "l", main = "Fig. 3 Body temperature, depth, shade and activity, 10% max shade"))
with(environ, points(ACT * 5 ~ dates, type = "l", col = "orange"))
with(environ, points(SHADE / 10 ~ dates, type = "h", col = "dark green"))
with(environ, points(DEP / 10 ~ dates, type = "l",col = "brown"))
abline(T_F_max, 0, lty = 2,col = 'red')
abline(T_F_min, 0, lty = 2,col = 'blue')
abline(T_pref, 0, lty = 2,col = 'orange')
text(0.1, T_F_max + 2, "T_F_max", col = 'red')
text(0.1, T_F_min - 2, "T_F_min", col = 'blue')
text(0.1, T_pref + 2, "T_pref", col = 'orange')
legend(x = 5, y = T_F_max + 17, legend = c("T_b (°C)", "depth (cm/10)", "activity, 0 5 or 10)", "shade (%/10)"), col = c("black", "brown", "orange", "dark green"), lty = rep(1, 4), bty = "n")

# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
forage<-subset(environ,ACT==2)
bask<-subset(environ,ACT==1)
night<-subset(environ,Solar==0)
day<-subset(environ,Solar==0)
with(night,plot(TIME ~ DOY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col=
    'dark blue', main = "Fig. 4 Annual activity window, 10% max shade"))
# nighttime hours
with(forage,points(TIME~DOY,pch=15,cex=2,col='orange')) # foraging Tbs
with(bask,points(TIME~DOY,pch=15,cex=2,col='light blue')) # basking Tbs

