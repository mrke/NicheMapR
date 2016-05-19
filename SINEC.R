microdaily=0
iday=1
TMIN=5 # daily minimum air temperature
TMAX=40 # daily maxmum air temperature
TMIN2=2 # daily minimum air temperature of next day (for 365 days of year calcs)
TMAX2=41 # daily maxmum air temperature of next day (for 365 days of year calcs)
TIMTMX=1300 # time of maximum (from SOLRAD)
TIMSR=500 # time of sunrise (from SOLRAD)
TIMSS=2000 # time of sunset (from SOLRAD)
ITAIR=0 # initial air temperature

TMIN=TMIN+273
TMAX=TMAX+273
TMIN2=TMIN2+273
TMAX2=TMAX2+273
A=(TMAX-TMIN)/2 # average temp
TSR=TMIN # temperature at sunrise
TREF=(TIMTMX-TIMSR)/2+TIMSR # reference time, half way between sunrise and time of the maximum value
SS=360*(TIMSS-TREF)/(2*(TIMTMX-TIMSR)) # hour angle for
SY=SS/57.29577
ZS=sin(SY)
TSS=A*ZS+TMIN+A
TAU=3/((2400-TIMSS)+TIMSR)

TIMARY<-rep(0,25)
TAIRRY<-TIMARY
# THIS LOOP COMPUTES TEMPERATURE EACH HOUR OF THE DAY
for(I in 1:24){
  J=I+1
  TIME=I*100
  if(TIME-TIMSR==0){
    #sunrise
    T=TMIN
    T=T-273
  }else{
    if(TIME-TIMSR<0){
      TI=(2400-TIMSS)+TIME
      if((microdaily==1)&(iday>1)){
        T=((TMIN-273-ITAIR)/TIMSR)*TIME+ITAIR
      }else{
        E=TI*TAU
        T=(TSS-TSR)/exp(E)+TSR
        T=T-273
      }
    }else{
      if(TIME-TIMSS>0){
        #  after sunset
        TI=(2400-TIMSS)-(2400-TIME)
        E=TI*TAU
        if(microdaily==1){
          T=(TSS-TMIN2)/exp(E)+TMIN2
        }else{
          T=(TSS-TSR)/exp(E)+TSR
        }
        T=T-273
      }else{
        # before or at sunset
        X=360*(TIME-TREF)/(2.*(TIMTMX-TIMSR))
        Y=X/57.29577
        Z=sin(Y)
        T=A*Z+TMIN+A
        T=T-273
      }
    }}
  #  CONVERTING FROM MILITARY TIME TO BIOME TIME
  ITIME=ceiling(TIME)/100
  FRMIN=(TIME/100)-ITIME
  ITIME=ITIME*60
  FRMIN=FRMIN*100
  TIMEC=ITIME+FRMIN
  TIMARY[J] = TIMEC
  TAIRRY[J] = T
}

TIMARY[1] = 0.0
if((microdaily==1)&(iday>1)){
  TAIRRY[1] = ITAIR
}else{
  TAIRRY[1] = T
}

ITAIR=TAIRRY[25]

plot(TAIRRY~TIMARY)
