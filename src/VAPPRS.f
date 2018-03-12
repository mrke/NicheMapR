         function vapprs(db)

C    COPYRIGHT 1997  WARREN P. PORTER,  ALL RIGHTS RESERVED
      implicit none
         real*4 loge
         real t,db,estar,vapprs
          t=db+273.16 
         if (t .le. 273.16)go to 20 
         loge=-7.90298*(373.16/t-1.)+5.02808*   
     x   alog10(373.16/t)-1.3816E-07 
     x   *(10.**(11.344*(1.-t/373.16))-1.)+8.1328E-03 
     x   *(10.**(-3.49149* 
     x   (373.16/t-1.))-1.)+alog10(1013.246)
         estar=(10.**loge)*100  
         go to 30   
20       loge=-9.09718*(273.16/t-1.)-3.56654*   
     x   alog10(273.16/t)+.876793*  
     x   (1.-t/273.16)+alog10(6.1071)   
         estar=(10.**loge)*100. 
30       continue   
         vapprs=estar   
         return 
         end
