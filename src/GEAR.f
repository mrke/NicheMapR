C	  THIS IS THE GEAR ADAMS PREDICTOR CORRECTOR ALGORITHM - LSODE
      SUBROUTINE LSODE (F, NEQ, Y, TT, TOUT, ITOL, RTOL, ATOL, ITASK,    
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, TRANCT,
     &            NDAYY)  
      use aacommondat
      EXTERNAL F, JAC                                                   
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF, NN   

      Character*1 TranIn       
      DOUBLE PRECISION Y, TT, TOUT, RTOL, ATOL, RWORK, Timout                    
      DIMENSION NEQ(1), Y(1), RTOL(1), ATOL(1), RWORK(LRW), IWORK(LIW)

C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      EXTERNAL PREPJ, SOLSY  
      double precision Enberr,Tprint 
      INTEGER I, I1, I2, IER, IFLAG, ILLIN, IMXER, INIT, IOWNS,         
     1   JSTART, KFLAG, KGO, L, LACOR, LENIW, LENRW, LENWM, LEWT, LF0,  
     2   LIWM, LSAVF, LWM, LYH, MAXORD, METH, MITER, ML, MORD, MU,      
     3   MXHNIL, MXHNL0, MXSTEP, MXSTP0, N, NFE, NHNIL, NJE, NQ, NQU,   
     4   NSLAST, NST, NTREP, NYH 
      INTEGER TRANCT,NDAYY,IT                                       
      DOUBLE PRECISION TRET, ROWNS, EL0, H, HMIN, HMXI, HU, TN, UROUND, 
     1   ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI, TCRIT, TDIST, 
     2   TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0, D1MACH, VNORM            
      DOUBLE PRECISION TPROD 
      DIMENSION MORD(2)                                                 
      LOGICAL IHIT                                                      
C-----------------------------------------------------------------------
C THE FOLLOWING INTERNAL COMMON BLOCK CONTAINS                          
C (A) VARIABLES WHICH ARE LOCAL TO ANY SUBROUTINE BUT WHOSE VALUES MUST 
C     BE PRESERVED BETWEEN CALLS TO THE ROUTINE (OWN VARIABLES), AND    
C (B) VARIABLES WHICH ARE COMMUNICATED BETWEEN SUBROUTINES.             
C THE STRUCTURE OF THE BLOCK IS AS FOLLOWS..  ALL REAL VAIABLES ARE     
C LISTED FIRST, FOLLOWED BY ALL INTEGERS.  WITHIN EACH TYPE, THE        
C VARIABLES ARE GROUPED WITH THOSE LOCAL TO SUBROUTINE LSODE FIRST,     
C THEN THOSE LOCAL TO SUBROUTINE STODE, AND FINALLY THOSE USED          
C FOR COMMUNICATION.  THE BLOCK IS DECLARED IN SUBROUTINES              
C LSODE, INTDY, STODE, PREPJ, AND SOLSY.  GROUPS OF VARIABLES ARE       
C REPLACED BY DUMMY ARRAYS IN THE COMMON DECLARATIONS IN ROUTINES       
C WHERE THOSE VARIABLES ARE NOT USED.                                   
C-----------------------------------------------------------------------
      COMMON /LS0001/ TRET, ROWNS(210),                                 
     1   EL0, H, HMIN, HMXI, HU, TN, UROUND,                            
     2   ILLIN, INIT, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,               
     3   MXSTEP, MXHNIL, NHNIL, NTREP, NSLAST, NYH, IOWNS(6),           
     4   IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE,   
     5   NJE, NQU 
c	Output print interval in from Sub. Itaday
      Common/Usrop2/Enberr,tprint   
      common/outsub/IT  
      COMMON/TRANINIT/tranin
      common/daystorun/nn
C     ALLOCATE (XP(nn*25),YP(nn*25),ZP1(nn*25),ZP2(nn*25),ZP3(nn*25))
C     ALLOCATE (ZP4(nn*25),ZP5(nn*25),ZP6(nn*25),ZP7(nn*25))      
C                                                                       
      DATA  MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
c	Initializing first output time
      Timout = Tprint
      TN = TT
      if((TRANIN .eq. 'n') .or. (TRANIN .eq. 'N'))then
        ISTATE = 1    
      ENDIF
C-----------------------------------------------------------------------
C BLOCK A.                                                              
C THIS CODE BLOCK IS EXECUTED ON EVERY CALL.                            
C IT TESTS ISTATE AND ITASK FOR LEGALITY AND BRANCHES APPROPIATELY.     
C IF ISTATE .GT. 1 BUT THE FLAG INIT SHOWS THAT INITIALIZATION HAS      
C NOT YET BEEN DONE, AN ERROR RETURN OCCURS.                            
C IF ISTATE = 1 AND TOUT = TT, JUMP TO BLOCK G AND RETURN IMMEDIATELY.   
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601                   
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602                     
      IF (ISTATE .EQ. 1) GO TO 10                                       
      IF (INIT .EQ. 0) GO TO 603                                        
      IF (ISTATE .EQ. 2) GO TO 200                                      
      GO TO 20                                                          
 10   INIT = 0                                                          
      IF (TOUT .EQ. TT) GO TO 430                                        
 20   NTREP = 0                                                         
C-----------------------------------------------------------------------
C BLOCK B.                                                              
C THE NEXT CODE BLOCK IS EXECUTED FOR THE INITIAL CALL (ISTATE = 1),    
C OR FOR A CONTINUATION CALL WITH PARAMETER CHANGES (ISTATE = 3).       
C IT CONTAINS CHECKING OF ALL INPUTS AND VARIOUS INITIALIZATIONS.       
C                                                                       
C FIRST CHECK LEGALITY OF THE NON-OPTIONAL INPUTS NEQ, ITOL, IOPT,      
C MF, ML, AND MU.                                                       
C-----------------------------------------------------------------------
      IF (NEQ(1) .LE. 0) GO TO 604                                      
      IF (ISTATE .EQ. 1) GO TO 25                                       
      IF (NEQ(1) .GT. N) GO TO 605                                      
 25   N = NEQ(1)                                                        
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606                       
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607                       
      METH = MF/10                                                      
      MITER = MF - 10*METH                                              
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608                       
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608                     
      IF (MITER .LE. 3) GO TO 30                                        
      ML = IWORK(1)                                                     
      MU = IWORK(2)                                                     
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609                           
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610                           
 30   CONTINUE                                                          
C NEXT PROCESS AND CHECK THE OPTIONAL INPUTS. --------------------------
      IF (IOPT .EQ. 1) GO TO 40                                         
      MAXORD = MORD(METH)                                               
      MXSTEP = MXSTP0                                                   
      MXHNIL = MXHNL0                                                   
      IF (ISTATE .EQ. 1) H0 = 0.0D0                                     
      HMXI = 0.0D0                                                      
      HMIN = 0.0D0                                                      
      GO TO 60                                                          
 40   MAXORD = IWORK(5)                                                 
      IF (MAXORD .LT. 0) GO TO 611                                      
      IF (MAXORD .EQ. 0) MAXORD = 100                                   
      MAXORD = MIN0(MAXORD,MORD(METH))                                  
      MXSTEP = IWORK(6)                                                 
      IF (MXSTEP .LT. 0) GO TO 612                                      
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0                                
      MXHNIL = IWORK(7)                                                 
      IF (MXHNIL .LT. 0) GO TO 613                                      
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0                                
      IF (ISTATE .NE. 1) GO TO 50                                       
      H0 = RWORK(5)                                                     
      IF ((TOUT - TT)*H0 .LT. 0.0D0) GO TO 614                           
 50   HMAX = RWORK(6)                                                   
      IF (HMAX .LT. 0.0D0) GO TO 615                                    
      HMXI = 0.0D0                                                      
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX                            
      HMIN = RWORK(7)                                                   
      IF (HMIN .LT. 0.0D0) GO TO 616                                    
C-----------------------------------------------------------------------
C SET WORK ARRAY POINTERS AND CHECK LENGTHS LRW AND LIW.                
C POINTERS TO SEGMENTS OF RWORK AND IWORK ARE NAMED BY PREFIXING L TO   
C THE NAME OF THE SEGMENT.  E.G., THE SEGMENT YH STARTS AT RWORK(LYH).  
C-----------------------------------------------------------------------
 60   LYH = 21                                                          
      IF (ISTATE .EQ. 1) NYH = N                                        
      LWM = LYH + (MAXORD + 1)*NYH                                      
      IF (MITER .EQ. 0) LENWM = 0                                       
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENWM = N*N + 2               
      IF (MITER .EQ. 3) LENWM = N + 2                                   
      IF (MITER .GE. 4) LENWM = (2*ML + MU + 1)*N + 2                   
      LEWT = LWM + LENWM                                                
      LSAVF = LEWT + N                                                  
      LACOR = LSAVF + N                                                 
      LENRW = LACOR + N - 1                                             
      IWORK(17) = LENRW                                                 
      LIWM = 1                                                          
      LENIW = 20 + N                                                    
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 20                    
      IWORK(18) = LENIW                                                 
      IF (LENRW .GT. LRW) GO TO 617                                     
      IF (LENIW .GT. LIW) GO TO 618                                     
C CHECK RTOL AND ATOL FOR LEGALITY. ------------------------------------
      RTOLI = RTOL(1)                                                   
      ATOLI = ATOL(1)                                                   
      DO 80 I = 1,N                                                     
        IF (ITOL .GE. 3) RTOLI = RTOL(I)                                
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)               
        IF (RTOLI .LT. 0.0D0) GO TO 619                                 
        IF (ATOLI .LT. 0.0D0) GO TO 620                                 
 80     CONTINUE                                                        
      IF (ISTATE .EQ. 1) GO TO 100                                      
C IF ISTATE = 3, SET FLAG TO SIGNAL PARAMETER CHANGES TO STODE. --------
      JSTART = -1                                                       
      IF (NQ .LE. MAXORD) GO TO 85                                      
      DO 83 I = 1,N                                                     
  83    RWORK(I+LSAVF-1) = RWORK(I+LWM-1)                               
  85  RWORK(LWM) = DSQRT(UROUND)                                        
      IF (N .EQ. NYH) GO TO 200                                         
C NEQ WAS REDUCED.  ZERO PART OF YH TO AVOID UNDEFINED REFERENCES. -----
      I1 = LYH + L*NYH                                                  
      I2 = LYH + (MAXORD + 1)*NYH - 1                                   
      IF (I1 .GT. I2) GO TO 200                                         
      DO 90 I = I1,I2                                                   
 90     RWORK(I) = 0.0D0                                                
      GO TO 200                                                         
C-----------------------------------------------------------------------
C BLOCK C.                                                              
C THE NEXT BLOCK IS FOR THE INITIAL CALL ONLY (ISTATE = 1).             
C IT CONTAINS ALL REMAINING INITIALIZATIONS, THE INITIAL CALL TO F,     
C AND THE CALCULATION OF THE INITIAL STEP SIZE.                         
C-----------------------------------------------------------------------
 100  UROUND = D1MACH(4)                                                
      TN = TT                                                            
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110                    
      TCRIT = RWORK(1)                                                  
      IF ((TCRIT - TOUT)*(TOUT - TT) .LT. 0.0D0) GO TO 625               
      IF (H0 .NE. 0.0D0 .AND. (TT + H0 - TCRIT)*H0 .GT. 0.0D0)           
     1   H0 = TCRIT - TT                                                 
 110  JSTART = 0                                                        
      RWORK(LWM) = DSQRT(UROUND)                                        
      NHNIL = 0                                                         
      NST = 0                                                           
      NJE = 0                                                           
      NSLAST = 0                                                        
      HU = 0.0D0                                                        
      NQU = 0                                                           
C INITIAL CALL TO F.  (LF0 POINTS TO YH(*,2). --------------------------
      LF0 = LYH + NYH                                                   
      CALL F (NEQ, TT, Y, RWORK(LF0))                                    
      NFE = 1                                                           
C LOAD THE INITIAL VALUE VECTOR IN YH. ---------------------------------
      DO 115 I = 1,N                                                    
 115    RWORK(I+LYH-1) = Y(I)                                           
C LOAD THE ERROR WEIGHT VECTOR EWT.  (H IS TEMPORARILY SET TO 1.0.) ----
      NQ = 1                                                            
      H = 1.0D0                                                         
      CALL EWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))         
      DO 120 I = 1,N                                                    
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621                       
 120  CONTINUE                                                          
C-----------------------------------------------------------------------
C THE CODING BELOW COMPUTES THE STEP SIZE, H0, TO BE ATTEMPTED ON THE   
C FIRST STEP, UNLESS THE USER HAS SUPPLIED A VALUE FOR THIS.            
C FIRST CHECK THAT TOUT - TT DIFFERS SIGNIFICANTLY FROM ZERO.            
C A SCALAR TOLERANCE QUANTITY TOL IS COMPUTED, AS MAX(RTOL(I))          
C IF THIS IS POSITIVE, OR MAX(ATOL(I)/ABS(Y(I))) OTHERWISE, ADJUSTED    
C SO AS TO BE BETWEEN 100*UROUND AND 1.0E-3.                            
C THEN THE COMPUTED VALUE H0 IS GIVEN BY..                              
C                                      NEQ                              
C   H0**2 = TOL / ( W0**-2 + (1/NEQ) * SUM ( F(I)/YWT(I) )**2  )        
C                                       1                               
C WHERE   W0     = MAX ( ABS(TT), ABS(TOUT) ),                           
C         F(I)   = I-TH COMPONENT OF INITIAL VALUE OF F,                
C         YWT(I) = EWT(I)/TOL  (A WEIGHT FOR Y(I)).                     
C THE SIGN OF H0 IS INFERRED FROM THE INITIAL VALUES OF TOUT AND TT.     
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180                                      
      TDIST = DABS(TOUT - TT)                                            
      W0 = DMAX1(DABS(TT),DABS(TOUT))                                    
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622                         
      TOL = RTOL(1)                                                     
      IF (ITOL .LE. 2) GO TO 140                                        
      DO 130 I = 1,N                                                    
 130    TOL = DMAX1(TOL,RTOL(I))                                        
 140  IF (TOL .GT. 0.0D0) GO TO 160                                     
      ATOLI = ATOL(1)                                                   
      DO 150 I = 1,N                                                    
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)               
        AYI = DABS(Y(I))                                                
        IF (AYI .NE. 0.0D0) TOL = DMAX1(TOL,ATOLI/AYI)                  
 150    CONTINUE                                                        
 160  TOL = DMAX1(TOL,100.0D0*UROUND)                                   
      TOL = DMIN1(TOL,0.001D0)                                          
      SUM = VNORM (N, RWORK(LF0), RWORK(LEWT))                          
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2                              
      H0 = 1.0D0/DSQRT(SUM)                                             
      H0 = DMIN1(H0,TDIST)                                              
      H0 = DSIGN(H0,TOUT-TT)                                             
C ADJUST H0 IF NECESSARY TO MEET HMAX BOUND. ---------------------------
 180  RH = DABS(H0)*HMXI                                                
      IF (RH .GT. 1.0D0) H0 = H0/RH                                     
C LOAD H WITH H0 AND SCALE YH(*,2) BY H0. ------------------------------
      H = H0                                                            
      DO 190 I = 1,N                                                    
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)                              
      GO TO 270                                                         
C-----------------------------------------------------------------------
C BLOCK D.                                                              
C THE NEXT CODE BLOCK IS FOR CONTINUATION CALLS ONLY (ISTATE = 2 OR 3)  
C AND IS TO CHECK STOP CONDITIONS BEFORE TAKING A STEP.                 
C-----------------------------------------------------------------------
 200  NSLAST = NST                                                      
      GO TO (210, 250, 220, 230, 240), ITASK                            
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250                           
      CALL INTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)                   
      IF (IFLAG .NE. 0) GO TO 627                                       
      TT = TOUT                                                          
      GO TO 420                                                         
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)                             
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623                           
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250                           
      GO TO 400                                                         
 230  TCRIT = RWORK(1)                                                  
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624                          
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625                        
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245                           
      CALL INTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)                   
      IF (IFLAG .NE. 0) GO TO 627                                       
      TT = TOUT                                                          
      GO TO 420                                                         
 240  TCRIT = RWORK(1)                                                  
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624                          
 245  HMX = DABS(TN) + DABS(H)                                          
      IHIT = DABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX                   
      IF (IHIT) GO TO 400                                               
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)                             
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250                       
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)                           
      IF (ISTATE .EQ. 2) JSTART = -2                                    
C-----------------------------------------------------------------------
C BLOCK E.                                                              
C THE NEXT BLOCK IS NORMALLY EXECUTED FOR ALL CALLS AND CONTAINS        
C THE CALL TO THE ONE-STEP CORE INTEGRATOR STODE.                       
C                                                                       
C THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.                    
C                                                                       
C FIRST CHECK FOR TOO MANY STEPS BEING TAKEN, UPDATE EWT (IF NOT AT     
C START OF PROBLEM), CHECK FOR TOO MUCH ACCURACY BEING REQUESTED, AND   
C CHECK FOR H BELOW THE ROUNDOFF LEVEL IN TT.                            
C-----------------------------------------------------------------------
 250  CONTINUE                                                          
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500                           
      CALL EWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))         
      DO 260 I = 1,N                                                    
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510                       
 260    CONTINUE                                                        
 270  TOLSF = UROUND*VNORM (N, RWORK(LYH), RWORK(LEWT))                 
      IF (TOLSF .LE. 1.0D0) GO TO 280                                   
      TOLSF = TOLSF*2.0D0                                               
      IF (NST .EQ. 0) GO TO 626                                         
      GO TO 520                                                         
 280  IF ((TN + H) .NE. TN) GO TO 290                                   
      NHNIL = NHNIL + 1  

      IF (NHNIL .GT. MXHNIL) then 
      Go to 290
      else                                 
        CALL XERRWV(101, 1, 0, 0, 0, 2, TN, H)
      Endif

      IF (NHNIL .LT. MXHNIL)then
      GO TO 290
      else                  
        CALL XERRWV(102, 1, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      endif                     
 290  CONTINUE                                                          
C-----------------------------------------------------------------------
C     CALL STODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,PREPJ,SOLSY)
C-----------------------------------------------------------------------
      CALL STODE (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),     
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),           
     2   F, JAC, PREPJ, SOLSY)                                          
      KGO = 1 - KFLAG                                                   
      GO TO (300, 530, 540), KGO                                        
C-----------------------------------------------------------------------
C BLOCK F.                                                              
C THE FOLLOWING BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN FROM THE  
C CORE INTEGRATOR (KFLAG = 0).  TEST FOR STOP CONDITIONS.               
C-----------------------------------------------------------------------
 300  INIT = 1                                                          
      GO TO (310, 400, 330, 340, 350), ITASK                            
C ITASK = 1.  IF TOUT HAS BEEN REACHED, INTERPOLATE. -------------------
 310  Continue
        If ((TN - Timout)*H .LT. 0.0D0)then
c         Try again
          GO TO 250
         else 
c         Integrator time greater than the print interval: interpolate to get
c         a value at the print interval                         
        CALL INTDY (Timout, 0, RWORK(LYH), NYH, Y, IFLAG)
c         Intdy sets time to zero because it thinks things are over. Correcting this.
          TT = TOUT
C         IF THE PRINT TIME INTERVAL IS SMALL, e.g. 10 min AND THE TIME CONSTANT IS LONG (LARGE 
C         ANIMAL), e.g. 90 min, THE INTEGRATOR TIME INTERVAL, H, GETS LARGE AND FREQUENT CALLS 
C         TO NERR = 52 RESULT STATING THAT THE TIME IS NOT IN THE INTERVAL.  Some
C         of this can be ameliorated by inserting the code just below that continually
c         downsizes H. It slows execution some, but eliminates a lot of the error 
c         messages.
          IF(H.GT.TPRINT)THEN
                H = TPRINT/2.
          ENDIF
c         Call output subroutine; successful output for this interval
          Call Osub2(Timout,Y,TRANCT,NDAYY)
c       Checking to see if the day is done
          if(TN.lt. Tout)then
c           Incrementing output by another hour if not the end of the day
            Timout = Timout + Tprint
c           Reinitialize counter for problems
            Nst = 0
           else
            go to 420
          endif
          Go to 250
        Endif
                                                       
      GO TO 420                                                         
C ITASK = 3.  JUMP TO EXIT IF TOUT WAS REACHED. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400                           
      GO TO 250                                                         
C ITASK = 4.  SEE IF TOUT OR TCRIT WAS REACHED.  ADJUST H IF NECESSARY. 
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345                           
      CALL INTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)                   
      TT = TOUT                                                          
      GO TO 420                                                         
 345  HMX = DABS(TN) + DABS(H)                                          
      IHIT = DABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX                   
      IF (IHIT) GO TO 400                                               
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)                             
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250                       
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)                           
      JSTART = -2                                                       
      GO TO 250                                                         
C ITASK = 5.  SEE IF TCRIT WAS REACHED AND JUMP TO EXIT. ---------------
 350  HMX = DABS(TN) + DABS(H)                                          
      IHIT = DABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX                   
C-----------------------------------------------------------------------
C BLOCK G.                                                              
C THE FOLLOWING BLOCK HANDLES ALL SUCCESSFUL RETURNS FROM LSODE.        
C IF ITASK .NE. 1, Y IS LOADED FROM YH AND TT IS SET ACCORDINGLY.        
C ISTATE IS SET TO 2, THE ILLEGAL INPUT COUNTER IS ZEROED, AND THE      
C OPTIONAL OUTPUTS ARE LOADED INTO THE WORK ARRAYS BEFORE RETURNING.    
C IF ISTATE = 1 AND TOUT = TT, THERE IS A RETURN WITH NO ACTION TAKEN,   
C EXCEPT THAT IF THIS HAS HAPPENED REPEATEDLY, THE RUN IS TERMINATED.   
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N                                                    
 410    Y(I) = RWORK(I+LYH-1)                                           
      TT = TN                                                            
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420                    
      IF (IHIT) TT = TCRIT                                               
 420  ISTATE = 2                                                        
      ILLIN = 0                                                         
      RWORK(11) = HU                                                    
      RWORK(12) = H                                                     
      RWORK(13) = TN                                                    
      IWORK(11) = NST                                                   
      IWORK(12) = NFE                                                   
      IWORK(13) = NJE                                                   
      IWORK(14) = NQU                                                   
      IWORK(15) = NQ                                                    
      RETURN                                                            
C                                                                       
 430  NTREP = NTREP + 1                                                 
      IF (NTREP .LT. 5)then
      RETURN 
      else                                         
        CALL XERRWV(301, 1, 0, 0, 0, 1, TT, 0.0D0)
      endif                              
      GO TO 800                                                         
C-----------------------------------------------------------------------
C BLOCK H.                                                              
C THE FOLLOWING BLOCK HANDLES ALL UNSUCCESSFUL RETURNS OTHER THAN       
C THOSE FOR ILLEGAL INPUT.  FIRST THE ERROR MESSAGE ROUTINE IS CALLED.  
C IF THERE WAS AN ERROR TEST OR CONVERGENCE TEST FAILURE, IMXER IS SET. 
C THEN Y IS LOADED FROM YHIS, TT IS SET TO TN, AND THE ILLEGAL INPUT     
C COUNTER ILLIN IS SET TO 0.  THE OPTIONAL OUTPUTS ARE LOADED INTO      
C THE WORK ARRAYS BEFORE RETURNING.                                     
C-----------------------------------------------------------------------
C THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE REACHING TOUT. -----
 500  CALL XERRWV(201, 1, 1, MXSTEP, 0, 1, TN, 0.0D0)                        
      ISTATE = -1                                                       
      GO TO 580                                                         
C EWT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM). -----------
 510  EWTI = RWORK(LEWT+I-1)                                            
      CALL XERRWV(202, 1, 1, I, 0, 2, TN, EWTI)                              
      ISTATE = -6                                                       
      GO TO 580                                                         
C TOO MUCH ACCURACY REQUESTED FOR MACHINE PRECISION. --------------
 520  CALL XERRWV(203, 1, 0, 0, 0, 2, TN, TOLSF)                             
      RWORK(14) = TOLSF                                                 
      ISTATE = -2                                                       
      GO TO 580                                                         
C KFLAG = -1.  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H) = HMIN. 
 530  CALL XERRWV(204, 1, 0, 0, 0, 2, TN, H)                                 
      ISTATE = -4                                                       
      GO TO 560                                                         
C KFLAG = -2.  CONVERGENCE FAILED REPEATEDLY OR WITH ABS(H) = HMIN. 
 540  CALL XERRWV(205, 1, 0, 0, 0, 2, TN, H)                                 
      ISTATE = -5                                                       
C COMPUTE IMXER IF RELEVANT. -------------------------------------------
 560  BIG = 0.0D0                                                       
      IMXER = 1                                                         
      DO 570 I = 1,N                                                    
        SIZE = DABS(RWORK(I+LACOR-1)/RWORK(I+LEWT-1))                   
        IF (BIG .GE. SIZE) GO TO 570                                    
        BIG = SIZE                                                      
        IMXER = I                                                       
 570    CONTINUE                                                        
      IWORK(16) = IMXER                                                 
C SET Y VECTOR, TT, ILLIN, AND OPTIONAL OUTPUTS. ------------------------
 580  DO 590 I = 1,N                                                    
        Y(I) = RWORK(I+LYH-1)
 590  CONTINUE                                      
      TT = TN                                                            
      ILLIN = 0                                                         
      RWORK(11) = HU                                                    
      RWORK(12) = H                                                     
      RWORK(13) = TN                                                    
      IWORK(11) = NST                                                   
      IWORK(12) = NFE                                                   
      IWORK(13) = NJE                                                   
      IWORK(14) = NQU                                                   
      IWORK(15) = NQ                                                    
      RETURN                                                            
C-----------------------------------------------------------------------
C BLOCK I.                                                              
C THE FOLLOWING BLOCK HANDLES ALL ERROR RETURNS DUE TO ILLEGAL INPUT    
C (ISTATE = -3), AS DETECTED BEFORE CALLING THE CORE INTEGRATOR.        
C FIRST THE ERROR MESSAGE ROUTINE IS CALLED.  THEN IF THERE HAVE BEEN   
C 5 CONSECUTIVE SUCH RETURNS JUST BEFORE THIS CALL TO THE SOLVER,       
C THE RUN IS HALTED.                                                    
C-----------------------------------------------------------------------
C ISTATE = -3, ILLEGAL INPUT
 601  CALL XERRWV(1, 1, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)                       
      GO TO 700                                                         
C602  CALL XERRWV(30HLSODE--  ITASK (I1) ILLEGAL   ,                    
 602  CALL XERRWV(2, 1, 1, ITASK, 0, 0, 0.0D0, 0.0D0)                        
      GO TO 700                                                         
C603  CALL XERRWV(50HLSODE--  ISTATE .GT. 1 BUT LSODE NOT INITIALIZED  ,
C    1   50, 3, 1, 0, 0, 0, 0, 0.0D0, 0.0D0)                            
 603  CALL XERRWV(3, 1, 0, 0, 0, 0, 0.0D0, 0.0D0)                            
      GO TO 700                                                         
C604  CALL XERRWV(30HLSODE--  NEQ (=I1) .LT. 0     ,                    
 604  CALL XERRWV(4, 1, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)                       
      GO TO 700                                                         
C605  CALL XERRWV(50HLSODE--  ISTATE = 3 AND NEQ INCREASED (I1 TO I2)  ,
C    1   50, 5, 1, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)                       
 605  CALL XERRWV(5, 1, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)                       
      GO TO 700                                                         
C606  CALL XERRWV(30HLSODE--  ITOL (I1) ILLEGAL    ,                    
 606  CALL XERRWV(6, 1, 1, ITOL, 0, 0, 0.0D0, 0.0D0)                         
      GO TO 700                                                         
C607  CALL XERRWV(30HLSODE--  IOPT (I1) ILLEGAL    ,                    
 607  CALL XERRWV(7, 1, 1, IOPT, 0, 0, 0.0D0, 0.0D0)                         
      GO TO 700                                                         
C608  CALL XERRWV(30HLSODE--  MF (I1) ILLEGAL      ,                    
 608  CALL XERRWV(8, 1, 1, MF, 0, 0, 0.0D0, 0.0D0)                           
      GO TO 700                                                         
C609  CALL XERRWV(9, 1, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)
C 52HLSODE--  ML (I1) ILLEGAL.. .LT. 0 OR .GE. NEQ (I2) 
 609  CALL XERRWV(9, 1, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)                      
      GO TO 700                                                         
C610  CALL XERRWV(50HLSODE--  MU (I1) ILLEGAL.. .LT. 0 OR .GE. NEQ (I2),
C    1   50, 10, 1, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)                     
 610  CALL XERRWV(10, 1, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)                     
      GO TO 700                                                         
C611  CALL XERRWV(30HLSODE--  MAXORD (I1) .LT. 0   ,                    
 611  CALL XERRWV(11, 1, 1, MAXORD, 0, 0, 0.0D0, 0.0D0)                      
      GO TO 700                                                         
C612  CALL XERRWV(30HLSODE--  MXSTEP (I1) .LT. 0   ,                    
 612  CALL XERRWV(12, 1, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)                      
      GO TO 700                                                         
C613  CALL XERRWV(30HLSODE--  MXHNIL (I1) .LT. 0   ,                    
 613  CALL XERRWV(13, 1, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)                      
      GO TO 700                                                         
C     CALL XERRWV(50H      INTEGRATION DIRECTION IS GIVEN BY H0 (R1)   ,
C    1   50, 14, 1, 0, 0, 0, 1, H0, 0.0D0)
 614  CALL XERRWV(14, 1, 0, 0, 0, 2, TOUT, TT)                                
      GO TO 700                                                         
C615  CALL XERRWV(30HLSODE--  HMAX (R1) .LT. 0.0   ,                    
 615  CALL XERRWV(15, 1, 0, 0, 0, 1, HMAX, 0.0D0)                            
      GO TO 700                                                         
C616  CALL XERRWV(30HLSODE--  HMIN (R1) .LT. 0.0   ,                    
 616  CALL XERRWV(16, 1, 0, 0, 0, 1, HMIN, 0.0D0)                            
      GO TO 700                                                         
 617  CALL XERRWV(17, 1, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)                     
      GO TO 700                                                         
 618  CALL XERRWV(18, 1, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)                     
      GO TO 700                                                         
 619  CALL XERRWV(19, 1, 1, I, 0, 1, RTOLI, 0.0D0)                           
      GO TO 700                                                         
 620  CALL XERRWV(20, 1, 1, I, 0, 1, ATOLI, 0.0D0)                           
      GO TO 700                                                         
 621  EWTI = RWORK(LEWT+I-1)                                            
      CALL XERRWV(21, 1, 1, I, 0, 1, EWTI, 0.0D0)                            
      GO TO 700                                                         
 622  CALL XERRWV(22, 1, 0, 0, 0, 2, TOUT, TT)                                
      GO TO 700                                                         
 623  CALL XERRWV(23, 1, 1, ITASK, 0, 2, TOUT, TP)                           
      GO TO 700                                                         
 624  CALL XERRWV(24, 1, 0, 0, 0, 2, TCRIT, TN)                              
      GO TO 700                                                         
 625  CALL XERRWV(25, 1, 0, 0, 0, 2, TCRIT, TOUT)                            
      GO TO 700                                                         
C626  CALL XERRWV(50HLSODE--  AT START OF PROBLEM, TOO MUCH ACCURACY   ,
C    1   50, 26, 1, 0, 0, 0, 0, 0.0D0, 0.0D0)                           
 626  CALL XERRWV(26, 1, 0, 0, 0, 0, 0.0D0, 0.0D0)                           
      RWORK(14) = TOLSF                                                 
      GO TO 700                                                         
C627  CALL XERRWV(50HLSODE--  TROUBLE FROM INTDY. ITASK = I1, TOUT = R1,
C    1   50, 27, 1, 1, ITASK, 0, 1, TOUT, 0.0D0)                        
 627  CALL XERRWV(27, 1, 1, ITASK, 0, 1, TOUT, 0.0D0)                        
C                                                                       
 700  IF (ILLIN .EQ. 5) GO TO 710                                       
      ILLIN = ILLIN + 1                                                 
      ISTATE = -3                                                       
      RETURN                                                            
C710  CALL XERRWV(50HLSODE--  REPEATED OCCURRENCES OF ILLEGAL INPUT    ,
C    1   50, 302, 1, 0, 0, 0, 0, 0.0D0, 0.0D0)                          
 710  CALL XERRWV(302, 1, 0, 0, 0, 0, 0.0D0, 0.0D0)                          
C                                                                       
C800  CALL XERRWV(50HLSODE--  RUN ABORTED.. APPARENT INFINITE LOOP     ,
C    1   50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)                          
 800  CALL XERRWV(303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)                          
      RETURN                                                            
C----------------------- END OF SUBROUTINE LSODE -----------------------
      END                                                               
      SUBROUTINE INTDY (TT, K, YH, NYH, DKY, IFLAG)                      
      INTEGER K, NYH, IFLAG, I, IC, IER, IOWND, IOWNS, J, JB, JB2,      
     1   JJ, JJ1, JP1, JSTART, KFLAG, L, MAXORD, METH, MITER, N, NFE,   
     2   NJE, NQ, NQU, NST                                              
      DOUBLE PRECISION TT, YH, DKY, TPRES, TNEX, TPROD,                                    
     1   ROWND, ROWNS, EL0, H, HMIN, HMXI, HU, TN, UROUND,              
     2   C, R, S, TP                                                    
C     DIMENSION YH(NYH,1), DKY(1)     
c     kearney changed dimensions below to fix array bound error
      DIMENSION YH(NYH,2), DKY(1)      
      COMMON /LS0001/ ROWND, ROWNS(210),                                
     1   EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(14), IOWNS(6),       
     4   IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE,   
     5   NJE, NQU                                                       
C-----------------------------------------------------------------------
C INTDY COMPUTES INTERPOLATED VALUES OF THE K-TH DERIVATIVE OF THE      
C DEPENDENT VARIABLE VECTOR Y, AND STORES IT IN DKY.                    
C THIS ROUTINE IS CALLED BY LSODE WITH K = 0 AND TT = TOUT, BUT MAY      
C ALSO BE CALLED BY THE USER FOR ANY K UP TO THE CURRENT ORDER.         
C (SEE DETAILED INSTRUCTIONS IN LSODE USAGE DOCUMENTATION.)             
C-----------------------------------------------------------------------
C THE COMPUTED VALUES IN DKY ARE GOTTEN BY INTERPOLATION USING THE      
C NORDSIECK HISTORY ARRAY YH.  THIS ARRAY CORRESPONDS UNIQUELY TO A     
C VECTOR-VALUED POLYNOMIAL OF DEGREE NQCUR OR LESS, AND DKY IS SET      
C TO THE K-TH DERIVATIVE OF THIS POLYNOMIAL AT TT.                       
C THE FORMULA FOR DKY IS..                                              
C              Q                                                        
C  DKY(I)  =  SUM  C(J,K) * (TT - TN)**(J-K) * H**(-J) * YH(I,J+1)       
C             J=K                                                       
C WHERE  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.  
C THE QUANTITIES  NQ = NQCUR, L = NQ+1, N = NEQ, TN, AND H ARE          
C COMMUNICATED BY COMMON.  THE ABOVE SUM IS DONE IN REVERSE ORDER.      
C IFLAG IS RETURNED NEGATIVE IF EITHER K OR TT IS OUT OF BOUNDS.         
C-----------------------------------------------------------------------
      IFLAG = 0                                                         
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80                             
      TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
C       TP = PRIOR CURRENT TIME, TN = CURRENT INDEPENDENT VARIABLE (TIME, USUALLY)
C       FROM STODE BELOW: H = HU = THE STEP SIZE TO BE ATTEMPTED ON THE NEXT STEP.              
C          H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING THE       
C          PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE, BUT ITS      
C          SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM. 
      TPRES = TT-TP
      TNEX = TT-TN
        TPROD = TPRES*TNEX                       
      IF (TPROD .GT. 0.0D0) GO TO 90                              
C                                                                       
      S = (TT - TN)/H                                                    
      IC = 1                                                            
      IF (K .EQ. 0) GO TO 15                                            
      JJ1 = L - K                                                       
      DO 10 JJ = JJ1,NQ                                                 
 10     IC = IC*JJ                                                      
 15   C = FLOAT(IC)                                                    
      DO 20 I = 1,N                                                     
 20     DKY(I) = C*YH(I,L)                                              
      IF (K .EQ. NQ) GO TO 55                                           
      JB2 = NQ - K                                                      
      DO 50 JB = 1,JB2                                                  
        J = NQ - JB                                                     
        JP1 = J + 1                                                     
        IC = 1                                                          
        IF (K .EQ. 0) GO TO 35                                          
        JJ1 = JP1 - K                                                   
        DO 30 JJ = JJ1,J                                                
 30       IC = IC*JJ                                                    
 35     C = FLOAT(IC)                                                  
        DO 40 I = 1,N                                                   
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)                               
 50     CONTINUE                                                        
      IF (K .EQ. 0) RETURN                                              
 55   R = H**(-K)                                                       
      DO 60 I = 1,N                                                     
 60     DKY(I) = R*DKY(I)                                               
      RETURN                                                            
C                                                                       
C80   CALL XERRWV(30HINTDY--  K (I1) ILLEGAL       ,                    
 80   CALL XERRWV(51, 1, 1, K, 0, 0, 0.0D0, 0.0D0)                           
      IFLAG = -1                                                        
      RETURN                                                            
C90   CALL XERRWV(30HINTDY--  TT (R1) ILLEGAL       ,                    
 90   CALL XERRWV(52, 1, 0, 0, 0, 1, TT, 0.0D0)                               
c      CALL XERRWV(                                                      
c     1  60H      TT NOT IN INTERVAL TCUR - HU (= R1) TO TCUR (R2)       ,
c     1   60, 52, 1, 0, 0, 0, 2, TP, TN)                                 
      IFLAG = -2                                                        
      RETURN                                                            
C----------------------- END OF SUBROUTINE INTDY -----------------------
      END                                                               
      SUBROUTINE STODE (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,          
     1   WM, IWM, F, JAC, PJAC, SLVS)                                   
      EXTERNAL F, JAC, PJAC, SLVS                                       
      INTEGER NEQ, NYH, IWM, I, I1, IALTH, IER, IOWND, IREDO, IRET,     
     1   IPUP, J, JB, JSTART, KFLAG, L, LMAX, M, MAXORD, MEO, METH,     
     2   MITER, N, NCF, NEWQ, NFE, NJE, NQ, NQNYH, NQU, NST, NSTEPJ     
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM,                 
     1   ROWND, CONIT, CRATE, EL, ELCO, HOLD, RC, RMAX, TESCO,          
     2   EL0, H, HMIN, HMXI, HU, TN, UROUND,                            
     3   DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,              
     4   R, RH, RHDN, RHSM, RHUP, TOLD, VNORM                           
c      DIMENSION NEQ(1), Y(1), YH(NYH,1), YH1(1), EWT(1), SAVF(1),       
c     1   ACOR(1), WM(1), IWM(1)
c     kearney changed dimensions below to solve array bound error
      DIMENSION NEQ(1), Y(1), YH(NYH,2), YH1(1), EWT(1), SAVF(1),       
     1   ACOR(1), WM(3), IWM(21)     
      COMMON /LS0001/ ROWND, CONIT, CRATE, EL(13), ELCO(13,12),         
     1   HOLD, RC, RMAX, TESCO(3,12),                                   
     2   EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(14),                 
     3   IALTH, IPUP, LMAX, MEO, NQNYH, NSTEPJ,                         
     4   IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE,   
     5   NJE, NQU                                                       
C-----------------------------------------------------------------------
C STODE PERFORMS ONE STEP OF THE INTEGRATION OF AN INITIAL VALUE        
C PROBLEM FOR A SYSTEM OF ORDINARY DIFFERENTIAL EQUATIONS.              
C NOTE.. STODE IS INDEPENDENT OF THE VALUE OF THE ITERATION METHOD      
C INDICATOR MITER, WHEN THIS IS .NE. 0, AND HENCE IS INDEPENDENT        
C OF THE TYPE OF CHORD METHOD USED, OR THE JACOBIAN STRUCTURE.          
C COMMUNICATION WITH STODE IS DONE WITH THE FOLLOWING VARIABLES..       
C                                                                       
C Y      = AN ARRAY OF LENGTH .GE. N USED AS THE Y ARGUMENT IN          
C          ALL CALLS TO F AND JAC.                                      
C NEQ    = INTEGER ARRAY CONTAINING PROBLEM SIZE IN NEQ(1), AND         
C          PASSED AS THE NEQ ARGUMENT IN ALL CALLS TO F AND JAC.        
C YH     = AN NYH BY LMAX ARRAY CONTAINING THE DEPENDENT VARIABLES      
C          AND THEIR APPROXIMATE SCALED DERIVATIVES, WHERE              
C          LMAX = MAXORD + 1.  YH(I,J+1) CONTAINS THE APPROXIMATE       
C          J-TH DERIVATIVE OF Y(I), SCALED BY H**J/FACTORIAL(J)         
C          (J = 0,1,...,NQ).  ON ENTRY FOR THE FIRST STEP, THE FIRST    
C          TWO COLUMNS OF YH MUST BE SET FROM THE INITIAL VALUES.       
C NYH    = A CONSTANT INTEGER .GE. N, THE FIRST DIMENSION OF YH.        
C YH1    = A ONE-DIMENSIONAL ARRAY OCCUPYING THE SAME SPACE AS YH.      
C EWT    = AN ARRAY OF N ELEMENTS WITH WHICH THE ESTIMATED LOCAL        
C          ERRORS IN YH ARE COMPARED.                                   
C SAVF   = AN ARRAY OF WORKING STORAGE, OF LENGTH N.                    
C ACOR   = A WORK ARRAY OF LENGTH N, USED FOR THE ACCUMULATED           
C          CORRECTIONS.  ON A SUCCESSFUL RETURN, ACOR(I) CONTAINS       
C          THE ESTIMATED ONE-STEP LOCAL ERROR IN Y(I).                  
C WM,IWM = REAL AND INTEGER WORK ARRAYS ASSOCIATED WITH MATRIX          
C          OPERATIONS IN CHORD ITERATION (MITER .NE. 0).                
C PJAC   = NAME OF ROUTINE TO EVALUATE AND PREPROCESS JACOBIAN MATRIX   
C          IF A CHORD METHOD IS BEING USED.                             
C SLVS   = NAME OF ROUTINE TO SOLVE LINEAR SYSTEM IN CHORD ITERATION.   
C H      = THE STEP SIZE TO BE ATTEMPTED ON THE NEXT STEP.              
C          H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING THE       
C          PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE, BUT ITS      
C          SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM.            
C HMIN   = THE MINIMUM ABSOLUTE VALUE OF THE STEP SIZE H TO BE USED.    
C HMXI   = INVERSE OF THE MAXIMUM ABSOLUTE VALUE OF H TO BE USED.       
C          HMXI = 0.0 IS ALLOWED AND CORRESPONDS TO AN INFINITE HMAX.   
C          HMIN AND HMXI MAY BE CHANGED AT ANY TIME, BUT WILL NOT       
C          TAKE EFFECT UNTIL THE NEXT CHANGE OF H IS CONSIDERED.        
C TN     = THE INDEPENDENT VARIABLE. TN IS UPDATED ON EACH STEP TAKEN.  
C JSTART = AN INTEGER USED FOR INPUT ONLY, WITH THE FOLLOWING           
C          VALUES AND MEANINGS..                                        
C               0  PERFORM THE FIRST STEP.                              
C           .GT.0  TAKE A NEW STEP CONTINUING FROM THE LAST.            
C              -1  TAKE THE NEXT STEP WITH A NEW VALUE OF H, MAXORD,    
C                    N, METH, MITER, AND/OR MATRIX PARAMETERS.          
C              -2  TAKE THE NEXT STEP WITH A NEW VALUE OF H,            
C                    BUT WITH OTHER INPUTS UNCHANGED.                   
C          ON RETURN, JSTART IS SET TO 1 TO FACILITATE CONTINUATION.    
C KFLAG  = A COMPLETION CODE WITH THE FOLLOWING MEANINGS..              
C               0  THE STEP WAS SUCCESFUL.                              
C              -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED.           
C              -2  CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED.         
C          A RETURN WITH KFLAG = -1 OR -2 MEANS EITHER                  
C          ABS(H) = HMIN OR 10 CONSECUTIVE FAILURES OCCURRED.           
C          ON A RETURN WITH KFLAG NEGATIVE, THE VALUES OF TN AND        
C          THE YH ARRAY ARE AS OF THE BEGINNING OF THE LAST             
C          STEP, AND H IS THE LAST STEP SIZE ATTEMPTED.                 
C MAXORD = THE MAXIMUM ORDER OF INTEGRATION METHOD TO BE ALLOWED.       
C METH/MITER = THE METHOD FLAGS.  SEE DESCRIPTION IN DRIVER.            
C N      = THE NUMBER OF FIRST-ORDER DIFFERENTIAL EQUATIONS.            
C-----------------------------------------------------------------------
      KFLAG = 0                                                         
      TOLD = TN                                                         
      NCF = 0                                                           
      IF (JSTART .GT. 0) GO TO 200                                      
      IF (JSTART .EQ. -1) GO TO 100                                     
      IF (JSTART .EQ. -2) GO TO 160                                     
C-----------------------------------------------------------------------
C ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER VARIABLES ARE     
C INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE INCREASED   
C IN A SINGLE STEP.  IT IS INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL   
C INITIAL H, BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE            
C OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT 2     
C FOR THE NEXT INCREASE.                                                
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1                                                 
      NQ = 1                                                            
      L = 2                                                             
      IALTH = 2                                                         
      RMAX = 10000.0D0                                                  
      RC = 0.0D0                                                        
      EL0 = 1.0D0                                                       
      CRATE = 0.7D0                                                     
      DELP = 0.0D0                                                      
      HOLD = H                                                          
      MEO = METH                                                        
      NSTEPJ = 0                                                        
      IRET = 3                                                          
      GO TO 140                                                         
C-----------------------------------------------------------------------
C THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN JSTART = -1.    
C IPUP IS SET TO MITER TO FORCE A MATRIX UPDATE.                        
C IF AN ORDER INCREASE IS ABOUT TO BE CONSIDERED (IALTH = 1),           
C IALTH IS RESET TO 2 TO POSTPONE CONSIDERATION ONE MORE STEP.          
C IF THE CALLER HAS CHANGED METH, CFODE IS CALLED TO RESET              
C THE COEFFICIENTS OF THE METHOD.                                       
C IF THE CALLER HAS CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT     
C ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN ACCORDINGLY.    
C IF H IS TO BE CHANGED, YH MUST BE RESCALED.                           
C IF H OR METH IS BEING CHANGED, IALTH IS RESET TO L = NQ + 1           
C TO PREVENT FURTHER CHANGES IN H FOR THAT MANY STEPS.                  
C-----------------------------------------------------------------------
 100  IPUP = MITER                                                      
      LMAX = MAXORD + 1                                                 
      IF (IALTH .EQ. 1) IALTH = 2                                       
      IF (METH .EQ. MEO) GO TO 110                                      
      CALL CFODE (METH, ELCO, TESCO)                                    
      MEO = METH                                                        
      IF (NQ .GT. MAXORD) GO TO 120                                     
      IALTH = L                                                         
      IRET = 1                                                          
      GO TO 150                                                         
 110  IF (NQ .LE. MAXORD) GO TO 160                                     
 120  NQ = MAXORD                                                       
      L = LMAX                                                          
      DO 125 I = 1,L                                                    
 125    EL(I) = ELCO(I,NQ)                                              
      NQNYH = NQ*NYH                                                    
      RC = RC*EL(1)/EL0                                                 
      EL0 = EL(1)                                                       
      CONIT = 0.5D0/FLOAT(NQ+2)                                        
C     DDN = VNORM (N, YH(1,L+1), EWT)/TESCO(1,L)                        
      DDN = VNORM (N, SAVF, EWT)/TESCO(1,L)                             
      EXDN = 1.0D0/FLOAT(L)                                            
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)                      
      RH = DMIN1(RHDN,1.0D0)                                            
      IREDO = 3                                                         
      IF (H .EQ. HOLD) GO TO 170                                        
      RH = DMIN1(RH,DABS(H/HOLD))                                       
      H = HOLD                                                          
      GO TO 175                                                         
C-----------------------------------------------------------------------
C CFODE IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS FOR THE       
C CURRENT METH.  THEN THE EL VECTOR AND RELATED CONSTANTS ARE RESET     
C WHENEVER THE ORDER NQ IS CHANGED, OR AT THE START OF THE PROBLEM.     
C-----------------------------------------------------------------------
 140  CALL CFODE (METH, ELCO, TESCO)                                    
 150  DO 155 I = 1,L                                                    
 155    EL(I) = ELCO(I,NQ)                                              
      NQNYH = NQ*NYH                                                    
      RC = RC*EL(1)/EL0                                                 
      EL0 = EL(1)                                                       
      CONIT = 0.5D0/FLOAT(NQ+2)                                        
      GO TO (160, 170, 200), IRET                                       
C-----------------------------------------------------------------------
C IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST              
C RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH IS SET TO     
C L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT MANY STEPS, UNLESS       
C FORCED BY A CONVERGENCE OR ERROR TEST FAILURE.                        
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200                                        
      RH = H/HOLD                                                       
      H = HOLD                                                          
      IREDO = 3                                                         
      GO TO 175                                                         
 170  RH = DMAX1(RH,HMIN/DABS(H))                                       
 175  RH = DMIN1(RH,RMAX)                                               
      RH = RH/DMAX1(1.0D0,DABS(H)*HMXI*RH)                              
      R = 1.0D0                                                         
      DO 180 J = 2,L                                                    
        R = R*RH                                                        
        DO 180 I = 1,N                                                  
 180      YH(I,J) = YH(I,J)*R                                           
      H = H*RH                                                          
      RC = RC*RH                                                        
      IALTH = L                                                         
      IF (IREDO .EQ. 0) GO TO 680                                       
C-----------------------------------------------------------------------
C THIS SECTION COMPUTES THE PREDICTED VALUES BY EFFECTIVELY             
C MULTIPLYING THE YH ARRAY BY THE PASCAL TRIANGLE MATRIX.               
C RC IS THE RATIO OF NEW TO OLD VALUES OF THE COEFFICIENT  H*EL(1).     
C WHEN RC DIFFERS FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER  
C TO FORCE PJAC TO BE CALLED, IF A JACOBIAN IS INVOLVED.                
C IN ANY CASE, PJAC IS CALLED AT LEAST EVERY 20-TH STEP.                
C-----------------------------------------------------------------------
 200  IF (DABS(RC-1.0D0) .GT. 0.3D0) IPUP = MITER                       
      IF (NST .GE. NSTEPJ+20) IPUP = MITER                              
      TN = TN + H                                                       
      I1 = NQNYH + 1                                                    
      DO 215 JB = 1,NQ                                                  
        I1 = I1 - NYH                                                   
        DO 210 I = I1,NQNYH                                             
 210      YH1(I) = YH1(I) + YH1(I+NYH)                                  
 215    CONTINUE                                                        
C-----------------------------------------------------------------------
C UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS        
C MADE ON THE R.M.S. NORM OF EACH CORRECTION, WEIGHTED BY THE ERROR     
C WEIGHT VECTOR EWT.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE  
C VECTOR ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE CORRECTOR LOOP.   
C-----------------------------------------------------------------------
 220  M = 0                                                             
      DO 230 I = 1,N                                                    
 230    Y(I) = YH(I,1)                                                  
      CALL F (NEQ, TN, Y, SAVF)                                         
      NFE = NFE + 1
      if(NFE.eq.617)then
      continue
      endif                                                     
      IF (IPUP .LE. 0) GO TO 250                                        
C-----------------------------------------------------------------------
C IF INDICATED, THE MATRIX P = I - H*EL(1)*J IS REEVALUATED AND         
C PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.  IPUP IS SET    
C TO 0 AS AN INDICATOR THAT THIS HAS BEEN DONE.                         
C-----------------------------------------------------------------------
      IPUP = 0                                                          
      RC = 1.0D0                                                        
      NSTEPJ = NST                                                      
      CRATE = 0.7D0                                                     
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)     
      IF (IER .NE. 0) GO TO 430                                         
 250  DO 260 I = 1,N                                                    
 260    ACOR(I) = 0.0D0                                                 
 270  IF (MITER .NE. 0) GO TO 350                                       
C-----------------------------------------------------------------------
C IN THE CASE OF FUNCTIONAL ITERATION, UPDATE Y DIRECTLY FROM           
C THE RESULT OF THE LAST FUNCTION EVALUATION.                           
C-----------------------------------------------------------------------
      DO 290 I = 1,N                                                    
        SAVF(I) = H*SAVF(I) - YH(I,2)                                   
 290    Y(I) = SAVF(I) - ACOR(I)                                        
      DEL = VNORM (N, Y, EWT)                                           
      DO 300 I = 1,N                                                    
        Y(I) = YH(I,1) + EL(1)*SAVF(I)                                  
 300    ACOR(I) = SAVF(I)                                               
      GO TO 400                                                         
C-----------------------------------------------------------------------
C IN THE CASE OF THE CHORD METHOD, COMPUTE THE CORRECTOR ERROR,         
C AND SOLVE THE LINEAR SYSTEM WITH THAT AS RIGHT-HAND SIDE AND          
C P AS COEFFICIENT MATRIX.                                              
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N                                                    
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))                          
      CALL SLVS (WM, IWM, Y, SAVF)                                      
      IF (IER .NE. 0) GO TO 410                                         
      DEL = VNORM (N, Y, EWT)                                           
      DO 380 I = 1,N                                                    
        ACOR(I) = ACOR(I) + Y(I)                                        
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)                                  
C-----------------------------------------------------------------------
C TEST FOR CONVERGENCE.  IF M.GT.0, AN ESTIMATE OF THE CONVERGENCE      
C RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.       
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = DMAX1(0.2D0*CRATE,DEL/DELP)                 
      DCON = DEL*DMIN1(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)           
      IF (DCON .LE. 1.0D0) GO TO 450                                    
      M = M + 1                                                         
      IF (M .EQ. 3) GO TO 410                                           
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410                 
      DELP = DEL                                                        
      CALL F (NEQ, TN, Y, SAVF)                                         
      NFE = NFE + 1                                                     
      GO TO 270                                                         
C-----------------------------------------------------------------------
C THE CORRECTOR ITERATION FAILED TO CONVERGE IN 3 TRIES.                
C IF MITER .NE. 0 AND THE JACOBIAN IS OUT OF DATE, PJAC IS CALLED FOR   
C THE NEXT TRY.  OTHERWISE THE YH ARRAY IS RETRACTED TO ITS VALUES      
C BEFORE PREDICTION, AND H IS REDUCED, IF POSSIBLE.  IF H CANNOT BE     
C REDUCED OR 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -2.           
C-----------------------------------------------------------------------
 410  IF (IPUP .EQ. 0) GO TO 430                                        
      IPUP = MITER                                                      
      GO TO 220                                                         
 430  TN = TOLD                                                         
      NCF = NCF + 1                                                     
      RMAX = 2.0D0                                                      
      I1 = NQNYH + 1                                                    
      DO 445 JB = 1,NQ                                                  
        I1 = I1 - NYH                                                   
        DO 440 I = I1,NQNYH                                             
 440      YH1(I) = YH1(I) - YH1(I+NYH)                                  
 445    CONTINUE                                                        
      IF (DABS(H) .LE. HMIN*1.00001D0) GO TO 670                        
      IF (NCF .EQ. 10) GO TO 670                                        
      RH = 0.25D0                                                       
      IPUP = MITER                                                      
      IREDO = 1                                                         
      GO TO 170                                                         
C-----------------------------------------------------------------------
C THE CORRECTOR HAS CONVERGED.  IPUP IS SET TO -1 IF MITER .NE. 0,      
C TO SIGNAL THAT THE JACOBIAN INVOLVED MAY NEED UPDATING LATER.         
C THE LOCAL ERROR TEST IS MADE AND CONTROL PASSES TO STATEMENT 500      
C IF IT FAILS.                                                          
C-----------------------------------------------------------------------
 450  IF (MITER .NE. 0) IPUP = -1                                       
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)                               
      IF (M .GT. 0) DSM = VNORM (N, ACOR, EWT)/TESCO(2,NQ)              
      IF (DSM .GT. 1.0D0) GO TO 500                                     
C-----------------------------------------------------------------------
C AFTER A SUCCESSFUL STEP, UPDATE THE YH ARRAY.                         
C CONSIDER CHANGING H IF IALTH = 1.  OTHERWISE DECREASE IALTH BY 1.     
C IF IALTH IS THEN 1 AND NQ .LT. MAXORD, THEN ACOR IS SAVED FOR         
C USE IN A POSSIBLE ORDER INCREASE ON THE NEXT STEP.                    
C IF A CHANGE IN H IS CONSIDERED, AN INCREASE OR DECREASE IN ORDER      
C BY ONE IS CONSIDERED ALSO.  A CHANGE IN H IS MADE ONLY IF IT IS BY A  
C FACTOR OF AT LEAST 1.1.  IF NOT, IALTH IS SET TO 3 TO PREVENT         
C TESTING FOR THAT MANY STEPS.                                          
C-----------------------------------------------------------------------
      KFLAG = 0                                                         
      IREDO = 0                                                         
      NST = NST + 1                                                     
      HU = H                                                            
      NQU = NQ                                                          
      DO 470 J = 1,L                                                    
        DO 470 I = 1,N                                                  
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)                             
      IALTH = IALTH - 1                                                 
      IF (IALTH .EQ. 0) GO TO 520                                       
      IF (IALTH .GT. 1) GO TO 690                                       
      IF (L .EQ. LMAX) GO TO 690                                        
      DO 490 I = 1,N                                                    
 490    YH(I,LMAX) = ACOR(I)                                            
      GO TO 690                                                         
C-----------------------------------------------------------------------
C THE ERROR TEST FAILED.  KFLAG KEEPS TRACK OF MULTIPLE FAILURES.       
C RESTORE TN AND THE YH ARRAY TO THEIR PREVIOUS VALUES, AND PREPARE     
C TO TRY THE STEP AGAIN.  COMPUTE THE OPTIMUM STEP SIZE FOR THIS OR     
C ONE LOWER ORDER.  AFTER 2 OR MORE FAILURES, H IS FORCED TO DECREASE   
C BY A FACTOR OF 0.2 OR LESS.                                           
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1                                                 
      TN = TOLD                                                         
      I1 = NQNYH + 1                                                    
      DO 515 JB = 1,NQ                                                  
        I1 = I1 - NYH                                                   
        DO 510 I = I1,NQNYH                                             
 510      YH1(I) = YH1(I) - YH1(I+NYH)                                  
 515    CONTINUE                                                        
      RMAX = 2.0D0                                                      
      IF (DABS(H) .LE. HMIN*1.00001D0) GO TO 660                        
      IF (KFLAG .LE. -3) GO TO 640                                      
      IREDO = 2                                                         
      RHUP = 0.0D0                                                      
      GO TO 540                                                         
C-----------------------------------------------------------------------
C REGARDLESS OF THE SUCCESS OR FAILURE OF THE STEP, FACTORS             
C RHDN, RHSM, AND RHUP ARE COMPUTED, BY WHICH H COULD BE MULTIPLIED     
C AT ORDER NQ - 1, ORDER NQ, OR ORDER NQ + 1, RESPECTIVELY.             
C IN THE CASE OF FAILURE, RHUP = 0.0 TO AVOID AN ORDER INCREASE.        
C THE LARGEST OF THESE IS DETERMINED AND THE NEW ORDER CHOSEN           
C ACCORDINGLY.  IF THE ORDER IS TO BE INCREASED, WE COMPUTE ONE         
C ADDITIONAL SCALED DERIVATIVE.                                         
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0                                                      
      IF (L .EQ. LMAX) GO TO 540                                        
      DO 530 I = 1,N                                                    
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)                                  
      DUP = VNORM (N, SAVF, EWT)/TESCO(3,NQ)                            
      EXUP = 1.0D0/FLOAT(L+1)                                          
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)                      
 540  EXSM = 1.0D0/FLOAT(L)                                            
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)                      
      RHDN = 0.0D0                                                      
      IF (NQ .EQ. 1) GO TO 560                                          
      DDN = VNORM (N, YH(1,L), EWT)/TESCO(1,NQ)                         
      EXDN = 1.0D0/FLOAT(NQ)                                           
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)                      
 560  IF (RHSM .GE. RHUP) GO TO 570                                     
      IF (RHUP .GT. RHDN) GO TO 590                                     
      GO TO 580                                                         
 570  IF (RHSM .LT. RHDN) GO TO 580                                     
      NEWQ = NQ                                                         
      RH = RHSM                                                         
      GO TO 620                                                         
 580  NEWQ = NQ - 1                                                     
      RH = RHDN                                                         
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0                  
      GO TO 620                                                         
 590  NEWQ = L                                                          
      RH = RHUP                                                         
      IF (RH .LT. 1.1D0) GO TO 610                                      
      R = EL(L)/FLOAT(L)                                               
      DO 600 I = 1,N                                                    
 600    YH(I,NEWQ+1) = ACOR(I)*R                                        
      GO TO 630                                                         
 610  IALTH = 3                                                         
      GO TO 690                                                         
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GO TO 610               
      IF (KFLAG .LE. -2) RH = DMIN1(RH,0.2D0)                           
C-----------------------------------------------------------------------
C IF THERE IS A CHANGE OF ORDER, RESET NQ, L, AND THE COEFFICIENTS.     
C IN ANY CASE H IS RESET ACCORDING TO RH AND THE YH ARRAY IS RESCALED.  
C THEN EXIT FROM 680 IF THE STEP WAS OK, OR REDO THE STEP OTHERWISE.    
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170                                       
 630  NQ = NEWQ                                                         
      L = NQ + 1                                                        
      IRET = 2                                                          
      GO TO 150                                                         
C-----------------------------------------------------------------------
C CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES HAVE OCCURED.      
C IF 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -1.                   
C IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE       
C YH ARRAY HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST             
C DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO 1.  THEN            
C H IS REDUCED BY A FACTOR OF 10, AND THE STEP IS RETRIED,              
C UNTIL IT SUCCEEDS OR H REACHES HMIN.                                  
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660                                     
      RH = 0.1D0                                                        
      RH = DMAX1(HMIN/DABS(H),RH)                                       
      H = H*RH                                                          
      DO 645 I = 1,N                                                    
 645    Y(I) = YH(I,1)                                                  
      CALL F (NEQ, TN, Y, SAVF)                                         
      NFE = NFE + 1                                                     
      DO 650 I = 1,N                                                    
 650    YH(I,2) = H*SAVF(I)                                             
      IPUP = MITER                                                      
      IALTH = 5                                                         
      IF (NQ .EQ. 1) GO TO 200                                          
      NQ = 1                                                            
      L = 2                                                             
      IRET = 3                                                          
      GO TO 150                                                         
C-----------------------------------------------------------------------
C ALL RETURNS ARE MADE THROUGH THIS SECTION.  H IS SAVED IN HOLD        
C TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP.                     
C-----------------------------------------------------------------------
 660  KFLAG = -1                                                        
      GO TO 700                                                         
 670  KFLAG = -2                                                        
      GO TO 700                                                         
 680  RMAX = 10.0D0                                                     
 690  R = 1.0D0/TESCO(2,NQU)                                            
      DO 695 I = 1,N                                                    
 695    ACOR(I) = ACOR(I)*R                                             
 700  HOLD = H                                                          
      JSTART = 1                                                        
      RETURN                                                            
C----------------------- END OF SUBROUTINE STODE -----------------------
      END                                                               
      SUBROUTINE CFODE (METH, ELCO, TESCO)                              
      INTEGER METH, I, IB, NQ, NQM1, NQP1                               
      DOUBLE PRECISION ELCO, TESCO, AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,  
     1   RQFAC, RQ1FAC, TSIGN, XPIN                                     
      DIMENSION ELCO(13,12), TESCO(3,12)                                
C-----------------------------------------------------------------------
C CFODE IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS         
C NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS            
C GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.         
C THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH = 2.   
C (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)               
C CFODE IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,                 
C AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.             
C                                                                       
C THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.                
C THE COEFFICIENTS EL(I), 1 .LE. I .LE. NQ+1, FOR THE METHOD OF         
C ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A GENETRATING   
C POLYNOMIAL, I.E.,                                                     
C     L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.                    
C FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY                      
C     DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) = 0.   
C FOR THE BDF METHODS, L(X) IS GIVEN BY                                 
C     L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,                                
C WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).               
C                                                                       
C THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE                  
C LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.         
C AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP            
C SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER     
C NQ + 1 IF K = 3.                                                      
C-----------------------------------------------------------------------
      DIMENSION PC(12)                                                  
C                                                                       
      GO TO (100, 200), METH                                            
C                                                                       
 100  ELCO(1,1) = 1.0D0                                                 
      ELCO(2,1) = 1.0D0                                                 
      TESCO(1,1) = 0.0D0                                                
      TESCO(2,1) = 2.0D0                                                
      TESCO(1,2) = 1.0D0                                                
      TESCO(3,12) = 0.0D0                                               
      PC(1) = 1.0D0                                                     
      RQFAC = 1.0D0                                                     
      DO 140 NQ = 2,12                                                  
C-----------------------------------------------------------------------
C THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL          
C     P(X) = (X+1)*(X+2)*...*(X+NQ-1).                                  
C INITIALLY, P(X) = 1.                                                  
C-----------------------------------------------------------------------
        RQ1FAC = RQFAC                                                  
        RQFAC = RQFAC/FLOAT(NQ)                                        
        NQM1 = NQ - 1                                                   
        FNQM1 = FLOAT(NQM1)                                            
        NQP1 = NQ + 1                                                   
C FORM COEFFICIENTS OF P(X)*(X+NQ-1). ----------------------------------
        PC(NQ) = 0.0D0                                                  
        DO 110 IB = 1,NQM1                                              
          I = NQP1 - IB                                                 
 110      PC(I) = PC(I-1) + FNQM1*PC(I)                                 
        PC(1) = FNQM1*PC(1)                                             
C COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X). -----------------------
        PINT = PC(1)                                                    
        XPIN = PC(1)/2.0D0                                              
        TSIGN = 1.0D0                                                   
        DO 120 I = 2,NQ                                                 
          TSIGN = -TSIGN                                                
          PINT = PINT + TSIGN*PC(I)/FLOAT(I)                           
 120      XPIN = XPIN + TSIGN*PC(I)/FLOAT(I+1)                         
C STORE COEFFICIENTS IN ELCO AND TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC                                        
        ELCO(2,NQ) = 1.0D0                                              
        DO 130 I = 2,NQ                                                 
 130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/FLOAT(I)                         
        AGAMQ = RQFAC*XPIN                                              
        RAGQ = 1.0D0/AGAMQ                                              
        TESCO(2,NQ) = RAGQ                                              
C       TESCO(1,NQP1) = RAGQ*RQFAC/FLOAT(NQP1)                         
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/FLOAT(NQP1)         
        TESCO(3,NQM1) = RAGQ                                            
 140    CONTINUE                                                        
      RETURN                                                            
C                                                                       
 200  PC(1) = 1.0D0                                                     
      RQ1FAC = 1.0D0                                                    
      DO 230 NQ = 1,5                                                   
C-----------------------------------------------------------------------
C THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL          
C     P(X) = (X+1)*(X+2)*...*(X+NQ).                                    
C INITIALLY, P(X) = 1.                                                  
C-----------------------------------------------------------------------
        FNQ = FLOAT(NQ)                                                
        NQP1 = NQ + 1                                                   
C FORM COEFFICIENTS OF P(X)*(X+NQ). ------------------------------------
        PC(NQP1) = 0.0D0                                                
        DO 210 IB = 1,NQ                                                
          I = NQ + 2 - IB                                               
 210      PC(I) = PC(I-1) + FNQ*PC(I)                                   
        PC(1) = FNQ*PC(1)                                               
C STORE COEFFICIENTS IN ELCO AND TESCO. --------------------------------
        DO 220 I = 1,NQP1                                               
 220      ELCO(I,NQ) = PC(I)/PC(2)                                      
        ELCO(2,NQ) = 1.0D0                                              
        TESCO(1,NQ) = RQ1FAC                                            
        TESCO(2,NQ) = FLOAT(NQP1)/ELCO(1,NQ)                           
        TESCO(3,NQ) = FLOAT(NQ+2)/ELCO(1,NQ)                           
        RQ1FAC = RQ1FAC/FNQ                                             
 230    CONTINUE                                                        
      RETURN                                                            
C----------------------- END OF SUBROUTINE CFODE -----------------------
      END                                                               
      SUBROUTINE PREPJ (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,      
     1   F, JAC)                                                        
      INTEGER NEQ, NYH, IWM, I, I1, I2, IER, II, IOWND, IOWNS, J, J1,   
     1   JJ, JSTART, KFLAG, L, LENP, MAXORD, MBA, MBAND, MEB1, MEBAND,  
     2   METH, MITER, ML, ML3, MU, N, NFE, NJE, NQ, NQU, NST            
      EXTERNAL F, JAC                                                   
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM,                      
     1   ROWND, ROWNS, EL0, H, HMIN, HMXI, HU, TN, UROUND,              
     2   CON, DI, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ, VNORM             
C     DIMENSION NEQ(1), Y(1), YH(NYH,1), EWT(1), FTEM(1), SAVF(1),      
C    1   WM(1), IWM(1) 
c     kearney changed dimensions below to solve array bound error
      DIMENSION NEQ(1), Y(1), YH(NYH,2), EWT(1), FTEM(1), SAVF(1),      
     1   WM(3), IWM(21)      
      COMMON /LS0001/ ROWND, ROWNS(210),                                
     1   EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(14), IOWNS(6),       
     4   IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE,   
     5   NJE, NQU                                                       
C-----------------------------------------------------------------------
C PREPJ IS CALLED BY STODE TO COMPUTE AND PROCESS THE MATRIX            
C P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN.      
C HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE JAC IF                
C MITER = 1 OR 4, OR BY FINITE DIFFERENCING IF MITER = 2, 3, OR 5.      
C IF MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED.                  
C J IS STORED IN WM AND REPLACED BY P.  IF MITER .NE. 3, P IS THEN      
C SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION       
C OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE          
C BY DGEFA IF MITER = 1 OR 2, AND BY DGBFA IF MITER = 4 OR 5.           
C                                                                       
C IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION          
C WITH PREPJ USES THE FOLLOWING..                                       
C Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.                    
C FTEM = WORK ARRAY OF LENGTH N (ACOR IN STODE).                        
C SAVF = ARRAY CONTAINING F EVALUATED AT PREDICTED Y.                   
C WM   = REAL WORK SPACE FOR MATRICES.  ON OUTPUT IT CONTAINS THE       
C        INVERSE DIAGONAL MATRIX IF MITER = 3 AND THE LU DECOMPOSITION  
C        OF P IF MITER IS 1, 2 , 4, OR 5.                               
C        STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).                    
C        WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..           
C        WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS.   
C        WM(2) = H*EL0, SAVED FOR LATER USE IF MITER = 3.               
C IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING AT   
C        IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS THE     
C        BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS 4 OR 5.
C EL0  = EL(1) (INPUT).                                                 
C IER  = OUTPUT ERROR FLAG,  = 0 IF NO TROUBLE, .NE. 0 IF               
C        P MATRIX FOUND TO BE SINGULAR.                                 
C THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND,       
C MITER, N, NFE, AND NJE.                                               
C-----------------------------------------------------------------------
      NJE = NJE + 1                                                     
      HL0 = H*EL0                                                       
      GO TO (100, 200, 300, 400, 500), MITER                            
C IF MITER = 1, CALL JAC AND MULTIPLY BY SCALAR. -----------------------
 100  LENP = N*N                                                        
      DO 110 I = 1,LENP                                                 
 110    WM(I+2) = 0.0D0                                                 
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)                             
      CON = -HL0                                                        
      DO 120 I = 1,LENP                                                 
 120    WM(I+2) = WM(I+2)*CON                                           
      GO TO 240                                                         
C IF MITER = 2, MAKE N CALLS TO F TO APPROXIMATE J. --------------------
 200  FAC = VNORM (N, SAVF, EWT)                                        
      R0 = 1000.0D0*DABS(H)*UROUND*FLOAT(N)*FAC                        
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0                                     
      SRUR = WM(1)                                                      
      J1 = 2                                                            
      DO 230 J = 1,N                                                    
        YJ = Y(J)                                                       
        R = DMAX1(SRUR*DABS(YJ),R0*EWT(J))                              
        Y(J) = Y(J) + R                                                 
        FAC = -HL0/R                                                    
        CALL F (NEQ, TN, Y, FTEM)                                       
        DO 220 I = 1,N                                                  
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC                            
        Y(J) = YJ                                                       
        J1 = J1 + N                                                     
 230    CONTINUE                                                        
      NFE = NFE + N                                                     
C ADD IDENTITY MATRIX. -------------------------------------------------
 240  J = 3                                                             
      DO 250 I = 1,N                                                    
        WM(J) = WM(J) + 1.0D0                                           
 250    J = J + (N + 1)                                                 
C DO LU DECOMPOSITION ON P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)                            
      RETURN                                                            
C IF MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND P. ---------
 300  WM(2) = HL0                                                       
      IER = 0                                                           
      R = EL0*0.1D0                                                     
      DO 310 I = 1,N                                                    
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))                           
      CALL F (NEQ, TN, Y, WM(3))                                        
      NFE = NFE + 1                                                     
      DO 320 I = 1,N                                                    
        R0 = H*SAVF(I) - YH(I,2)                                        
        DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))                           
        WM(I+2) = 1.0D0                                                 
        IF (DABS(R0) .LT. UROUND*EWT(I)) GO TO 320                      
        IF (DABS(DI) .EQ. 0.0D0) GO TO 330                              
        WM(I+2) = 0.1D0*R0/DI                                           
 320    CONTINUE                                                        
      RETURN                                                            
 330  IER = -1                                                          
      RETURN                                                            
C IF MITER = 4, CALL JAC AND MULTIPLY BY SCALAR. -----------------------
 400  ML = IWM(1)                                                       
      MU = IWM(2)                                                       
      ML3 = ML + 3                                                      
      MBAND = ML + MU + 1                                               
      MEBAND = MBAND + ML                                               
      LENP = MEBAND*N                                                   
      DO 410 I = 1,LENP                                                 
 410    WM(I+2) = 0.0D0                                                 
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)                    
      CON = -HL0                                                        
      DO 420 I = 1,LENP                                                 
 420    WM(I+2) = WM(I+2)*CON                                           
      GO TO 570                                                         
C IF MITER = 5, MAKE MBAND CALLS TO F TO APPROXIMATE J. ----------------
 500  ML = IWM(1)                                                       
      MU = IWM(2)                                                       
      MBAND = ML + MU + 1                                               
      MBA = MIN0(MBAND,N)                                               
      MEBAND = MBAND + ML                                               
      MEB1 = MEBAND - 1                                                 
      SRUR = WM(1)                                                      
      FAC = VNORM (N, SAVF, EWT)                                        
      R0 = 1000.0D0*DABS(H)*UROUND*FLOAT(N)*FAC                        
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0                                     
      DO 560 J = 1,MBA                                                  
        DO 530 I = J,N,MBAND                                            
          YI = Y(I)                                                     
          R = DMAX1(SRUR*DABS(YI),R0*EWT(I))                            
 530      Y(I) = Y(I) + R                                               
        CALL F (NEQ, TN, Y, FTEM)                                       
        DO 550 JJ = J,N,MBAND                                           
          Y(JJ) = YH(JJ,1)                                              
          YJJ = Y(JJ)                                                   
          R = DMAX1(SRUR*DABS(YJJ),R0*EWT(JJ))                          
          FAC = -HL0/R                                                  
          I1 = MAX0(JJ-MU,1)                                            
          I2 = MIN0(JJ+ML,N)                                            
          II = JJ*MEB1 - ML + 2                                         
          DO 540 I = I1,I2                                              
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC                          
 550      CONTINUE                                                      
 560    CONTINUE                                                        
      NFE = NFE + MBA                                                   
C ADD IDENTITY MATRIX. -------------------------------------------------
 570  II = MBAND + 2                                                    
      DO 580 I = 1,N                                                    
        WM(II) = WM(II) + 1.0D0                                         
 580    II = II + MEBAND                                                
C DO LU DECOMPOSITION OF P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)               
      RETURN                                                            
C----------------------- END OF SUBROUTINE PREPJ -----------------------
      END                                                               
      SUBROUTINE SOLSY (WM, IWM, X, TEM)                                
      INTEGER IWM, I, IER, IOWND, IOWNS, JSTART, KFLAG, L, MAXORD,      
     1   MEBAND, METH, MITER, ML, MU, N, NFE, NJE, NQ, NQU, NST         
      DOUBLE PRECISION WM, X, TEM,                                      
     1   ROWND, ROWNS, EL0, H, HMIN, HMXI, HU, TN, UROUND,              
     2   DI, HL0, PHL0, R                                               
C     DIMENSION WM(1), IWM(1), X(1), TEM(1)
c     kearney changed dimensions below to fix array bound error
      DIMENSION WM(3), IWM(21), X(1), TEM(1)                             
      COMMON /LS0001/ ROWND, ROWNS(210),                                
     1   EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(14), IOWNS(6),       
     4   IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE,   
     5   NJE, NQU                                                       
C-----------------------------------------------------------------------
C THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR SYSTEM ARISING FROM   
C A CHORD ITERATION.  IT IS CALLED BY STODE IF MITER .NE. 0.            
C IF MITER IS 1 OR 2, IT CALLS DGESL TO ACCOMPLISH THIS.                
C IF MITER = 3 IT UPDATES THE COEFFICIENT H*EL0 IN THE DIAGONAL         
C MATRIX, AND THEN COMPUTES THE SOLUTION.                               
C IF MITER IS 4 OR 5, IT CALLS DGBSL.                                   
C COMMUNICATION WITH SOLSY USES THE FOLLOWING VARIABLES..               
C WM  = REAL WORK SPACE CONTAINING THE INVERSE DIAGONAL MATRIX IF MITER 
C       IS 3 AND THE LU DECOMPOSITION OF THE MATRIX OTHERWISE.          
C       STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).                     
C       WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..            
C       WM(1) = SQRT(UROUND) (NOT USED HERE),                           
C       WM(2) = HL0, THE PREVIOUS VALUE OF H*EL0, USED IF MITER = 3.    
C IWM = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING AT    
C       IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS THE      
C       BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS 4 OR 5. 
C X   = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION VECTOR    
C       ON OUTPUT, OF LENGTH N.                                         
C TEM = VECTOR OF WORK SPACE OF LENGTH N, NOT USED IN THIS VERSION.     
C IER = OUTPUT FLAG (IN COMMON).  IER = 0 IF NO TROUBLE OCCURRED.       
C       IER = -1 IF A SINGULAR MATRIX AROSE WITH MITER = 3.             
C THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, MITER, AND N.     
C-----------------------------------------------------------------------
      IER = 0                                                           
      GO TO (100, 100, 300, 400, 400), MITER                            
 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)                           
      RETURN                                                            
C                                                                       
 300  PHL0 = WM(2)                                                      
      HL0 = H*EL0                                                       
      WM(2) = HL0                                                       
      IF (HL0 .EQ. PHL0) GO TO 330                                      
      R = HL0/PHL0                                                      
      DO 320 I = 1,N                                                    
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))                          
        IF (DABS(DI) .EQ. 0.0D0) GO TO 390                              
 320    WM(I+2) = 1.0D0/DI                                              
 330  DO 340 I = 1,N                                                    
 340    X(I) = WM(I+2)*X(I)                                             
      RETURN                                                            
 390  IER = -1                                                          
      RETURN                                                            
C                                                                       
 400  ML = IWM(1)                                                       
      MU = IWM(2)                                                       
      MEBAND = 2*ML + MU + 1                                            
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)              
      RETURN                                                            
C----------------------- END OF SUBROUTINE SOLSY -----------------------
      END                                                               
      SUBROUTINE EWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)                 
C-----------------------------------------------------------------------
C THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR EWT ACCORDING TO         
C     EWT(I) = RTOL(I)*DABS(YCUR(I)) + ATOL(I),  I = 1,...,N,           
C WITH THE SUBSCRIPT ON RTOL AND/OR ATOL POSSIBLY REPLACED BY 1 ABOVE,  
C DEPENDING ON THE VALUE OF ITOL.                                       
C-----------------------------------------------------------------------
      INTEGER N, ITOL, I                                                
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT, ATOLI, RTOLI              
      DIMENSION RTOL(1), ATOL(1), YCUR(N), EWT(N)                       
      RTOLI = RTOL(1)                                                   
      ATOLI = ATOL(1)                                                   
      DO 10 I = 1,N                                                     
        IF (ITOL .GE. 3) RTOLI = RTOL(I)                                
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)               
        EWT(I) = RTOLI*DABS(YCUR(I)) + ATOLI                            
 10     CONTINUE                                                        
      RETURN                                                            
C----------------------- END OF SUBROUTINE EWSET -----------------------
      END                                                               
      DOUBLE PRECISION FUNCTION VNORM (N, V, W)                         
C-----------------------------------------------------------------------
C THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM     
C OF THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS      
C CONTAINED IN THE ARRAY W OF LENGTH N..                                
C   VNORM = SQRT( (1/N) * SUM( V(I)/W(I) )**2 )                         
C-----------------------------------------------------------------------
      INTEGER N, I                                                      
      DOUBLE PRECISION V, W, SUM                                        
      DIMENSION V(N), W(N)                                              
      SUM = 0.0D0                                                       
      DO 10 I = 1,N                                                     
 10     SUM = SUM + (V(I)/W(I))**2                                      
      VNORM = DSQRT(SUM/FLOAT(N))                                      
      RETURN                                                            
C----------------------- END OF FUNCTION VNORM -------------------------
      END                                                               
      SUBROUTINE SVCOM (RSAV, ISAV)                                     
C-----------------------------------------------------------------------
C THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCKS    
C LS0001 AND EH0001, WHICH ARE USED INTERNALLY IN THE LSODE PACKAGE.    
C                                                                       
C RSAV = REAL ARRAY OF LENGTH 218 OR MORE.                              
C ISAV = INTEGER ARRAY OF LENGTH 35 OR MORE.                            
C-----------------------------------------------------------------------
      INTEGER ISAV, I, IEH, ILS, LENILS, LENRLS                         
      DOUBLE PRECISION RSAV, RLS                                        
      DIMENSION RSAV(1), ISAV(1)                                        
      COMMON /LS0001/ RLS(218), ILS(33)                                 
      COMMON /EH0001/ IEH(2)                                            
      DATA LENRLS/218/, LENILS/33/                                      
C                                                                       
      DO 10 I = 1,LENRLS                                                
 10     RSAV(I) = RLS(I)                                                
      DO 20 I = 1,LENILS                                                
 20     ISAV(I) = ILS(I)                                                
      ISAV(LENILS+1) = IEH(1)                                           
      ISAV(LENILS+2) = IEH(2)                                           
      RETURN                                                            
C----------------------- END OF SUBROUTINE SVCOM -----------------------
      END                                                               
      SUBROUTINE RSCOM (RSAV, ISAV)                                     
C-----------------------------------------------------------------------
C THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON       
C BLOCKS LS0001 AND EH0001, WHICH ARE USED INTERNALLY IN THE LSODE      
C PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS       
C OF SUBROUTINE SVCOM OR THE EQUIVALENT.                                
C-----------------------------------------------------------------------
      INTEGER ISAV, I, IEH, ILS, LENILS, LENRLS                         
      DOUBLE PRECISION RSAV, RLS                                        
      DIMENSION RSAV(1), ISAV(1)                                        
      COMMON /LS0001/ RLS(218), ILS(33)                                 
      COMMON /EH0001/ IEH(2)                                            
      DATA LENRLS/218/, LENILS/33/                                      
C                                                                       
      DO 10 I = 1,LENRLS                                                
 10     RLS(I) = RSAV(I)                                                
      DO 20 I = 1,LENILS                                                
 20     ILS(I) = ISAV(I)                                                
      IEH(1) = ISAV(LENILS+1)                                           
      IEH(2) = ISAV(LENILS+2)                                           
      RETURN                                                            
C----------------------- END OF SUBROUTINE RSCOM -----------------------
      END                                                               
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)                               
      INTEGER LDA,N,IPVT(1),INFO                                        
      DOUBLE PRECISION A(LDA,1)                                         
C                                                                       
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.  
C                                                                       
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED            
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .                   
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       DOUBLE PRECISION(LDA, N)                               
C                THE MATRIX TO BE FACTORED.                             
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         
C                WHICH WERE USED TO OBTAIN IT.                          
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       
C                                                                       
C        IPVT    INTEGER(N)                                             
C                AN INTEGER VECTOR OF PIVOT INDICES.                    
C                                                                       
C        INFO    INTEGER                                                
C                = 0  NORMAL VALUE.                                     
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO  
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE   
C                     INDICATION OF SINGULARITY.                        
C                                                                       
C     LINPACK. THIS VERSION DATED 08/14/78 .                            
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS DAXPY,DSCAL,IDAMAX                                           
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      DOUBLE PRECISION TT                                                
      INTEGER IDAMAX,J,K,KP1,L,NM1                                      
C                                                                       
C                                                                       
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        
C                                                                       
      INFO = 0                                                          
      NM1 = N - 1                                                       
      IF (NM1 .LT. 1) GO TO 70                                          
      DO 60 K = 1, NM1                                                  
         KP1 = K + 1                                                    
C                                                                       
C        FIND L = PIVOT INDEX                                           
C                                                                       
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1                             
         IPVT(K) = L                                                    
C                                                                       
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          
C                                                                       
         IF (A(L,K) .EQ. 0.0D0) GO TO 40                                
C                                                                       
C           INTERCHANGE IF NECESSARY                                    
C                                                                       
            IF (L .EQ. K) GO TO 10                                      
               TT = A(L,K)                                               
               A(L,K) = A(K,K)                                          
               A(K,K) = TT                                               
   10       CONTINUE                                                    
C                                                                       
C           COMPUTE MULTIPLIERS                                         
C                                                                       
            TT = -1.0D0/A(K,K)                                           
            CALL DSCAL(N-K,TT,A(K+1,K),1)                                
C                                                                       
C           ROW ELIMINATION WITH COLUMN INDEXING                        
C                                                                       
            DO 30 J = KP1, N                                            
               TT = A(L,J)                                               
               IF (L .EQ. K) GO TO 20                                   
                  A(L,J) = A(K,J)                                       
                  A(K,J) = TT                                            
   20          CONTINUE                                                 
               CALL DAXPY(N-K,TT,A(K+1,K),1,A(K+1,J),1)                  
   30       CONTINUE                                                    
         GO TO 50                                                       
   40    CONTINUE                                                       
            INFO = K                                                    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
   70 CONTINUE                                                          
      IPVT(N) = N                                                       
      IF (A(N,N) .EQ. 0.0D0) INFO = N                                   
      RETURN                                                            
      END                                                               
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)                              
      INTEGER LDA,N,IPVT(1),JOB                                         
      DOUBLE PRECISION A(LDA,1),B(1)                                    
C                                                                       
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM                          
C     A * X = B  OR  TRANS(A) * X = B                                   
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       DOUBLE PRECISION(LDA, N)                               
C                THE OUTPUT FROM DGECO OR DGEFA.                        
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C        IPVT    INTEGER(N)                                             
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  
C                                                                       
C        B       DOUBLE PRECISION(N)                                    
C                THE RIGHT HAND SIDE VECTOR.                            
C                                                                       
C        JOB     INTEGER                                                
C                = 0         TO SOLVE  A*X = B ,                        
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE            
C                            TRANS(A)  IS THE TRANSPOSE.                
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        B       THE SOLUTION VECTOR  X .                               
C                                                                       
C     ERROR CONDITION                                                   
C                                                                       
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0           
C        OR DGEFA HAS SET INFO .EQ. 0 .                                 
C                                                                       
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
C     WITH  P  COLUMNS                                                  
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                            
C           IF (RCOND IS TOO SMALL) GO TO ...                           
C           DO 10 J = 1, P                                              
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                        
C        10 CONTINUE                                                    
C                                                                       
C     LINPACK. THIS VERSION DATED 08/14/78 .                            
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS DAXPY,DDOT                                                   
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      DOUBLE PRECISION DDOT,TT                                           
      INTEGER K,KB,L,NM1                                                
C                                                                       
      NM1 = N - 1                                                       
      IF (JOB .NE. 0) GO TO 50                                          
C                                                                       
C        JOB = 0 , SOLVE  A * X = B                                     
C        FIRST SOLVE  L*Y = B                                           
C                                                                       
         IF (NM1 .LT. 1) GO TO 30                                       
         DO 20 K = 1, NM1                                               
            L = IPVT(K)                                                 
            TT = B(L)                                                    
            IF (L .EQ. K) GO TO 10                                      
               B(L) = B(K)                                              
               B(K) = TT                                                 
   10       CONTINUE                                                    
            CALL DAXPY(N-K,TT,A(K+1,K),1,B(K+1),1)                       
   20    CONTINUE                                                       
   30    CONTINUE                                                       
C                                                                       
C        NOW SOLVE  U*X = Y                                             
C                                                                       
         DO 40 KB = 1, N                                                
            K = N + 1 - KB                                              
            B(K) = B(K)/A(K,K)                                          
            TT = -B(K)                                                   
            CALL DAXPY(K-1,TT,A(1,K),1,B(1),1)                           
   40    CONTINUE                                                       
      GO TO 100                                                         
   50 CONTINUE                                                          
C                                                                       
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         
C        FIRST SOLVE  TRANS(U)*Y = B                                    
C                                                                       
         DO 60 K = 1, N                                                 
            TT = DDOT(K-1,A(1,K),1,B(1),1)                               
            B(K) = (B(K) - TT)/A(K,K)                                    
   60    CONTINUE                                                       
C                                                                       
C        NOW SOLVE TRANS(L)*X = Y                                       
C                                                                       
         IF (NM1 .LT. 1) GO TO 90                                       
         DO 80 KB = 1, NM1                                              
            K = N - KB                                                  
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)                 
            L = IPVT(K)                                                 
            IF (L .EQ. K) GO TO 70                                      
               TT = B(L)                                                 
               B(L) = B(K)                                              
               B(K) = TT                                                 
   70       CONTINUE                                                    
   80    CONTINUE                                                       
   90    CONTINUE                                                       
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                       
      INTEGER LDA,N,ML,MU,IPVT(1),INFO                                  
      DOUBLE PRECISION ABD(LDA,1)                                       
C                                                                       
C     DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION.      
C                                                                       
C     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED            
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        ABD     DOUBLE PRECISION(LDA, N)                               
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       
C                SEE THE COMMENTS BELOW FOR DETAILS.                    
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              
C                LDA MUST BE .GE. 2*ML + MU + 1 .                       
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE ORIGINAL MATRIX.                      
C                                                                       
C        ML      INTEGER                                                
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           
C                0 .LE. ML .LT. N .                                     
C                                                                       
C        MU      INTEGER                                                
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           
C                0 .LE. MU .LT. N .                                     
C                MORE EFFICIENT IF  ML .LE. MU .                        
C     ON RETURN                                                         
C                                                                       
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       
C                                                                       
C        IPVT    INTEGER(N)                                             
C                AN INTEGER VECTOR OF PIVOT INDICES.                    
C                                                                       
C        INFO    INTEGER                                                
C                = 0  NORMAL VALUE.                                     
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        
C                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF        
C                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE      
C                     INDICATION OF SINGULARITY.                        
C                                                                       
C     BAND STORAGE                                                      
C                                                                       
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      
C           WILL SET UP THE INPUT.                                      
C                                                                       
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)                
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                
C                   M = ML + MU + 1                                     
C                   DO 20 J = 1, N                                      
C                      I1 = MAX0(1, J-MU)                               
C                      I2 = MIN0(N, J+ML)                               
C                      DO 10 I = I1, I2                                 
C                         K = I - J + M                                 
C                         ABD(K,J) = A(I,J)                             
C                10    CONTINUE                                         
C                20 CONTINUE                                            
C                                                                       
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          
C                                                                       
C     LINPACK. THIS VERSION DATED 08/14/78 .                            
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS DAXPY,DSCAL,IDAMAX                                           
C     FORTRAN MAX0,MIN0                                                 
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      DOUBLE PRECISION TT                                                
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1             
C                                                                       
C                                                                       
      M = ML + MU + 1                                                   
      INFO = 0                                                          
C                                                                       
C     ZERO INITIAL FILL-IN COLUMNS                                      
C                                                                       
      J0 = MU + 2                                                       
      J1 = MIN0(N,M) - 1                                                
      IF (J1 .LT. J0) GO TO 30                                          
      DO 20 JZ = J0, J1                                                 
         I0 = M + 1 - JZ                                                
         DO 10 I = I0, ML                                               
            ABD(I,JZ) = 0.0D0                                           
   10    CONTINUE                                                       
   20 CONTINUE                                                          
   30 CONTINUE                                                          
      JZ = J1                                                           
      JU = 0                                                            
C                                                                       
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        
C                                                                       
      NM1 = N - 1                                                       
      IF (NM1 .LT. 1) GO TO 130                                         
      DO 120 K = 1, NM1                                                 
         KP1 = K + 1                                                    
C                                                                       
C        ZERO NEXT FILL-IN COLUMN                                       
C                                                                       
         JZ = JZ + 1                                                    
         IF (JZ .GT. N) GO TO 50                                        
         IF (ML .LT. 1) GO TO 50                                        
            DO 40 I = 1, ML                                             
               ABD(I,JZ) = 0.0D0                                        
   40       CONTINUE                                                    
   50    CONTINUE                                                       
C                                                                       
C        FIND L = PIVOT INDEX                                           
C                                                                       
         LM = MIN0(ML,N-K)                                              
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1                            
         IPVT(K) = L + K - M                                            
C                                                                       
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          
C                                                                       
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100                             
C                                                                       
C           INTERCHANGE IF NECESSARY                                    
C                                                                       
            IF (L .EQ. M) GO TO 60                                      
               TT = ABD(L,K)                                             
               ABD(L,K) = ABD(M,K)                                      
               ABD(M,K) = TT                                             
   60       CONTINUE                                                    
C                                                                       
C           COMPUTE MULTIPLIERS                                         
C                                                                       
            TT = -1.0D0/ABD(M,K)                                         
            CALL DSCAL(LM,TT,ABD(M+1,K),1)                               
C                                                                       
C           ROW ELIMINATION WITH COLUMN INDEXING                        
C                                                                       
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)                            
            MM = M                                                      
            IF (JU .LT. KP1) GO TO 90                                   
            DO 80 J = KP1, JU                                           
               L = L - 1                                                
               MM = MM - 1                                              
               TT = ABD(L,J)                                             
               IF (L .EQ. MM) GO TO 70                                  
                  ABD(L,J) = ABD(MM,J)                                  
                  ABD(MM,J) = TT                                         
   70          CONTINUE                                                 
               CALL DAXPY(LM,TT,ABD(M+1,K),1,ABD(MM+1,J),1)              
   80       CONTINUE                                                    
   90       CONTINUE                                                    
         GO TO 110                                                      
  100    CONTINUE                                                       
            INFO = K                                                    
  110    CONTINUE                                                       
  120 CONTINUE                                                          
  130 CONTINUE                                                          
      IPVT(N) = N                                                       
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N                                 
      RETURN                                                            
      END                                                               
      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)                      
      INTEGER LDA,N,ML,MU,IPVT(1),JOB                                   
      DOUBLE PRECISION ABD(LDA,1),B(1)                                  
C                                                                       
C     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM                     
C     A * X = B  OR  TRANS(A) * X = B                                   
C     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.                     
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        ABD     DOUBLE PRECISION(LDA, N)                               
C                THE OUTPUT FROM DGBCO OR DGBFA.                        
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE ORIGINAL MATRIX.                      
C                                                                       
C        ML      INTEGER                                                
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           
C                                                                       
C        MU      INTEGER                                                
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           
C                                                                       
C        IPVT    INTEGER(N)                                             
C                THE PIVOT VECTOR FROM DGBCO OR DGBFA.                  
C                                                                       
C        B       DOUBLE PRECISION(N)                                    
C                THE RIGHT HAND SIDE VECTOR.                            
C                                                                       
C        JOB     INTEGER                                                
C                = 0         TO SOLVE  A*X = B ,                        
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE           
C                            TRANS(A)  IS THE TRANSPOSE.                
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        B       THE SOLUTION VECTOR  X .                               
C                                                                       
C     ERROR CONDITION                                                   
C                                                                       
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     
C        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0           
C        OR DGBFA HAS SET INFO .EQ. 0 .                                 
C                                                                       
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
C     WITH  P  COLUMNS                                                  
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    
C           IF (RCOND IS TOO SMALL) GO TO ...                           
C           DO 10 J = 1, P                                              
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)                
C        10 CONTINUE                                                    
C                                                                       
C     LINPACK. THIS VERSION DATED 08/14/78 .                            
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS DAXPY,DDOT                                                   
C     FORTRAN MIN0                                                      
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      DOUBLE PRECISION DDOT,TT                                           
      INTEGER K,KB,L,LA,LB,LM,M,NM1                                     
C                                                                       
      M = MU + ML + 1                                                   
      NM1 = N - 1                                                       
      IF (JOB .NE. 0) GO TO 50                                          
C                                                                       
C        JOB = 0 , SOLVE  A * X = B                                     
C        FIRST SOLVE L*Y = B                                            
C                                                                       
         IF (ML .EQ. 0) GO TO 30                                        
         IF (NM1 .LT. 1) GO TO 30                                       
            DO 20 K = 1, NM1                                            
               LM = MIN0(ML,N-K)                                        
               L = IPVT(K)                                              
               TT = B(L)                                                 
               IF (L .EQ. K) GO TO 10                                   
                  B(L) = B(K)                                           
                  B(K) = TT                                              
   10          CONTINUE                                                 
               CALL DAXPY(LM,TT,ABD(M+1,K),1,B(K+1),1)                   
   20       CONTINUE                                                    
   30    CONTINUE                                                       
C                                                                       
C        NOW SOLVE  U*X = Y                                             
C                                                                       
         DO 40 KB = 1, N                                                
            K = N + 1 - KB                                              
            B(K) = B(K)/ABD(M,K)                                        
            LM = MIN0(K,M) - 1                                          
            LA = M - LM                                                 
            LB = K - LM                                                 
            TT = -B(K)                                                   
            CALL DAXPY(LM,TT,ABD(LA,K),1,B(LB),1)                        
   40    CONTINUE                                                       
      GO TO 100                                                         
   50 CONTINUE                                                          
C                                                                       
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         
C        FIRST SOLVE  TRANS(U)*Y = B                                    
C                                                                       
         DO 60 K = 1, N                                                 
            LM = MIN0(K,M) - 1                                          
            LA = M - LM                                                 
            LB = K - LM                                                 
            TT = DDOT(LM,ABD(LA,K),1,B(LB),1)                            
            B(K) = (B(K) - TT)/ABD(M,K)                                  
   60    CONTINUE                                                       
C                                                                       
C        NOW SOLVE TRANS(L)*X = Y                                       
C                                                                       
         IF (ML .EQ. 0) GO TO 90                                        
         IF (NM1 .LT. 1) GO TO 90                                       
            DO 80 KB = 1, NM1                                           
               K = N - KB                                               
               LM = MIN0(ML,N-K)                                        
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)             
               L = IPVT(K)                                              
               IF (L .EQ. K) GO TO 70                                   
                  TT = B(L)                                              
                  B(L) = B(K)                                           
                  B(K) = TT                                              
   70          CONTINUE                                                 
   80       CONTINUE                                                    
   90    CONTINUE                                                       
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)                            
C                                                                       
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      DOUBLE PRECISION DX(1),DY(1),DA                                   
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 
C                                                                       
      IF(N.LE.0)RETURN                                                  
      IF (DA .EQ. 0.0D0) RETURN                                         
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
C                                                                       
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
C          NOT EQUAL TO 1                                               
C                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        DY(IY) = DY(IY) + DA*DX(IX)                                     
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
C                                                                       
C                                                                       
C        CLEAN-UP LOOP                                                  
C                                                                       
   20 M = MOD(N,4)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        DY(I) = DY(I) + DA*DX(I)                                        
   30 CONTINUE                                                          
      IF( N .LT. 4 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,4                                                 
        DY(I) = DY(I) + DA*DX(I)                                        
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)                            
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)                            
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)                            
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSCAL(N,DA,DX,INCX)                                    
C                                                                       
C     SCALES A VECTOR BY A CONSTANT.                                    
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.                   
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      DOUBLE PRECISION DA,DX(1)                                         
      INTEGER I,INCX,M,MP1,N,NINCX                                      
C                                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        DX(I) = DA*DX(I)                                                
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
C                                                                       
C        CLEAN-UP LOOP                                                  
C                                                                       
   20 M = MOD(N,5)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        DX(I) = DA*DX(I)                                                
   30 CONTINUE                                                          
      IF( N .LT. 5 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,5                                                 
        DX(I) = DA*DX(I)                                                
        DX(I + 1) = DA*DX(I + 1)                                        
        DX(I + 2) = DA*DX(I + 2)                                        
        DX(I + 3) = DA*DX(I + 3)                                        
        DX(I + 4) = DA*DX(I + 4)                                        
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)                 
C                                                                       
C     FORMS THE DOT PRODUCT OF TWO VECTORS.                             
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      DOUBLE PRECISION DX(1),DY(1),DTEMP                                
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 
C                                                                       
      DDOT = 0.0D0                                                      
      DTEMP = 0.0D0                                                     
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
C                                                                       
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
C          NOT EQUAL TO 1                                               
C                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        DTEMP = DTEMP + DX(IX)*DY(IY)                                   
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      DDOT = DTEMP                                                      
      RETURN                                                            
C                                                                       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
C                                                                       
C                                                                       
C        CLEAN-UP LOOP                                                  
C                                                                       
   20 M = MOD(N,5)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        DTEMP = DTEMP + DX(I)*DY(I)                                     
   30 CONTINUE                                                          
      IF( N .LT. 5 ) GO TO 60                                           
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,5                                                 
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +             
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE                                                          
   60 DDOT = DTEMP                                                      
      RETURN                                                            
      END                                                               
      INTEGER FUNCTION IDAMAX(N,DX,INCX)                                
C                                                                       
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      DOUBLE PRECISION DX(1),DMAX                                       
      INTEGER I,INCX,IX,N                                               
C                                                                       
      IDAMAX = 0                                                        
      IF( N .LT. 1 ) RETURN                                             
      IDAMAX = 1                                                        
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      IX = 1                                                            
      DMAX = DABS(DX(1))                                                
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5                               
         IDAMAX = I                                                     
         DMAX = DABS(DX(IX))                                            
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 DMAX = DABS(DX(1))                                                
      DO 30 I = 2,N                                                     
         IF(DABS(DX(I)).LE.DMAX) GO TO 30                               
         IDAMAX = I                                                     
         DMAX = DABS(DX(I))                                             
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      DOUBLE PRECISION FUNCTION D1MACH (IDUM)                           
      INTEGER IDUM                                                      
C-----------------------------------------------------------------------
C THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE IN DOUBLE      
C PRECISION.  THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER   
C U SUCH THAT  1.0D0 + U .NE. 1.0D0 (IN DOUBLE PRECISION).              
C-----------------------------------------------------------------------
      DOUBLE PRECISION U, COMP                                          
      U = 1.0D0                                                         
 10   U = U*0.5D0                                                       
      COMP = 1.0D0 + U                                                  
      IF (COMP .NE. 1.0D0) GO TO 10                                     
      D1MACH = U*2.0D0                                                  
      RETURN                                                            
C----------------------- END OF FUNCTION D1MACH ------------------------
      END                                                               
      SUBROUTINE XERRWV (NERR, IERT, NI, I1, I2, NR, R1, R2) 
      INTEGER  NERR, IERT, NI, I1, I2, NR,                    
     1   I, LUN, LUNIT, MESFLG, NWDS  
      DOUBLE PRECISION R1, R2                                           
C-----------------------------------------------------------------------
C SUBROUTINES XERRWV, XSETF, AND XSETUN, AS GIVEN HERE, CONSTITUTE      
C A SIMPLIFIED VERSION OF THE SLATEC ERROR HANDLING PACKAGE.            
C WRITTEN BY A. C. HINDMARSH AT LLL.  VERSION OF JANUARY 23, 1980.      
C THIS VERSION IS IN DOUBLE PRECISION.                                  
C                                                                       
C ALL ARGUMENTS ARE INPUT ARGUMENTS.                                    
C                                                                       
C NERR   = THE ERROR NUMBER.                               
C IERT   = THE ERROR TYPE..                                             
C          1 MEANS RECOVERABLE (CONTROL RETURNS TO CALLER).             
C          2 MEANS FATAL (RUN IS ABORTED--SEE NOTE BELOW).              
C NI     = NUMBER OF INTEGERS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.  
C I1,I2  = INTEGERS TO BE PRINTED, DEPENDING ON NI.                     
C NR     = NUMBER OF REALS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.     
C R1,R2  = REALS TO BE PRINTED, DEPENDING ON NI.                        
C                                                                       
C NOTE..  THIS ROUTINE IS MACHINE-DEPENDENT AND SPECIALIZED FOR USE     
C IN LIMITED CONTEXT, IN THE FOLLOWING WAYS..                           
C 1. THE NUMBER OF HOLLERITH CHARACTERS STORED PER WORD IS ASSUMED      
C    TO BE 4.                                                           
C 2. THE VALUE OF NMES IS ASSUMED TO BE A MULTIPLE OF 4 AND AT MOST 60. 
C    (MULTI-LINE MESSAGES ARE GENERATED BY REPEATED CALLS.)             
C 3. IF IERT = 2, CONTROL PASSES TO THE STATEMENT   STOP                
C    TO ABORT THE RUN.  THIS STATEMENT MAY BE MACHINE-DEPENDENT.        
C 4. R1 AND R2 ARE ASSUMED TO BE IN DOUBLE PRECISION AND ARE PRINTED    
C    IN D21.13 FORMAT.                                                  
C 5. THE COMMON BLOCK /EH0001/ BELOW IS DATA-LOADED (A MACHINE-         
C    DEPENDENT FEATURE) WITH DEFAULT VALUES.                            
C    THIS BLOCK IS NEEDED FOR PROPER RETENTION OF PARAMETERS USED BY    
C    THIS ROUTINE WHICH THE USER CAN RESET BY CALLING XSETF OR XSETUN.  
C    THE VARIABLES IN THIS BLOCK ARE AS FOLLOWS..                       
C       MESFLG = PRINT CONTROL FLAG..                                   
C                1 MEANS PRINT ALL MESSAGES (THE DEFAULT).              
C                0 MEANS NO PRINTING.                                   
C       LUNIT  = LOGICAL UNIT NUMBER FOR MESSAGES.                      
C                                                                       
C-----------------------------------------------------------------------
C THE FOLLOWING ARE INSTRUCTIONS FOR INSTALLING THIS ROUTINE            
C IN DIFFERENT MACHINE ENVIRONMENTS.                                    
C                                                                       
C TO CHANGE THE DEFAULT OUTPUT UNIT, CHANGE THE DATA STATEMENT          
C IN THE BLOCK DATA SUBPROGRAM BELOW.                                   
C                                                                       
C FOR A DIFFERENT NUMBER OF CHARACTERS PER WORD, CHANGE THE             
C STATEMENTS SETTING NWDS BELOW, AND FORMAT 10.                         
C                                                                       
C FOR A DIFFERENT RUN-ABORT COMMAND, CHANGE THE STATEMENT FOLLOWING     
C STATEMENT 100 AT THE END.                                             
C-----------------------------------------------------------------------
      COMMON /EH0001/ MESFLG, LUNIT                                     
C                                                                       
      IF (MESFLG .EQ. 0) GO TO 100                                      
C GET LOGICAL UNIT NUMBER. ---------------------------------------------
      LUN = LUNIT                                                       
C WRITE THE MESSAGE. ---------------------------------------------------
      if(nerr .eq. 101)then
c        WRITE (LUN, *)'LSODE--  WARNING.. INTERNAL TT',R1,' AND H ',R2
c        WRITE (LUN, *)'ARE SUCH THAT IN THE MACHINE, TT + H = TT ON THE'
c	  WRITE (LUN, *)' NEXT STEP (H = STEP SIZE).' 
c	  WRITE (LUN, *)' SOLVER WILL CONTINUE ANYWAY'
      endif 
      if(nerr .eq. 102)then
c        WRITE (LUN, *)'LSODE--  ABOVE WARNING HAS BEEN ISSUED ',I1,   
c     &  'TIMES. IT WILL NOT BE ISSUED AGAIN FOR THIS PROBLEM'                                 
      endif
      if(nerr .eq. 301)then
c        WRITE (LUN, *)'LSODE--  REPEATED CALLS WITH ISTATE = 1',
c     &  ' AND TOUT = TT',R1
      endif
      if(nerr .eq. 201)then
c        WRITE (LUN, *)'LSODE--  AT CURRENT TT',R1,' MXSTEP ',I1,
c     &  ' STEPS TAKEN ON THIS CALL BEFORE REACHING TOUT'
      endif 
      if(nerr .eq. 202)then
c        WRITE (LUN, *)'LSODE--  AT TT ',R1,' EWT ',I1,' HAS BECOME',
c     &  R2,' .LE. 0.0' 
      endif
      if(nerr .eq. 203)then
c        WRITE (LUN, *)'LSODE--  AT TT ',R1,' TOO MUCH ACCURACY REQUESTED    
c     &  ',' FOR PRECISION OF MACHINE..  SEE ',TOLSF ,R2   
      endif
      if(nerr .eq. 204)then
c        WRITE (LUN, *)'LSODE--  AT TT ',R1,' AND STEP SIZE H ',R2,
c     &  ' THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H) = ',R2 
      endif
      if(nerr .eq. 205)then
c        WRITE (LUN, *)'LSODE--  AT TT ',R1,' AND STEP SIZE H ',R2,
c     &  ' THE CORRECTOR CONVERGENCE FAILED REPEATEDLY ',
c     &  ' OR WITH ABS(H) = ',R2                    
      endif
      if(nerr .eq. 1)then
c        WRITE (LUN, *)'LSODE--  ISTATE ',I1,' ILLEGAL INPUT'                 
      endif
      if(nerr .eq. 2)then
c        WRITE (LUN, *)'LSODE--  ITASK ',I1,' ILLEGAL.'   
      endif
      if(nerr .eq. 3)then
c        WRITE (LUN, *)'LSODE--  ISTATE .GT. 1 BUT LSODE NOT INITIALIZED'   
      endif
      if(nerr .eq. 4)then
c        WRITE (LUN, *)'LSODE--  NEQ = ',I1,' .LT. 0'                  
      endif
      if(nerr .eq. 5)then
c        WRITE (LUN, *)'LSODE--  ISTATE = 3 AND NEQ INCREASED (',I1,
c     &  ' TO ',I2,' )'
      endif
      if(nerr .eq. 6)then
c        WRITE (LUN, *)'LSODE--  Tolerance ITOL (',I1,') ILLEGAL.'    
      endif
      if(nerr .eq. 7)then
c        WRITE (LUN, *)'LSODE--  IOPT (',I1,' ) ILLEGAL.'    
      endif
      if(nerr .eq. 8)then
c        WRITE (LUN, *)'32HLSODE--  MF (',I1,') ILLEGAL.'    
      endif
      if(nerr .eq. 9)then
c        WRITE (LUN, *)'LSODE--  ML (',I1,') ILLEGAL.. .LT. 0 OR .GE. ',
c     &  'NEQ (',I2,').'
      endif
      if(nerr .eq. 10)then
c        WRITE (LUN, *)'LSODE--  MU (',I1,') ILLEGAL.. .LT. 0 OR .GE. ',
c     &  'NEQ (',I2,').'
      endif
      if(nerr .eq. 11)then
c        WRITE (LUN, *)'LSODE--  MAXORD (',I1,') .LT. 0'                 
      endif
      if(nerr .eq. 12)then
c        WRITE (LUN, *)'LSODE--  MXSTEP (',I1,') .LT. 0'                 
      endif
      if(nerr .eq. 13)then
c        WRITE (LUN, *)'LSODE--  MXHNIL (',I1,') .LT. 0'               
      endif
      if(nerr .eq. 14)then
c        WRITE (LUN, *)'40HLSODE--  TOUT (',R1,') BEHIND TT (',R2,')'        ,          
      endif
      if(nerr .eq. 15)then
c        WRITE (LUN, *)'LSODE--  HMAX (',R1,') .LT. 0.0'                  
      endif
      if(nerr .eq. 16)then
c        WRITE (LUN, *)'LSODE--  HMIN (',R1,') .LT. 0.0'                 
      endif
      if(nerr .eq. 17)then
c        WRITE (LUN, *)'LSODE--  RWORK LENGTH NEEDED, LENRW (',I1,'),',
c     &  'EXCEEDS LRW (',I2,')'
      endif
      if(nerr .eq. 18)then
c        WRITE (LUN, *)'LSODE--  IWORK LENGTH NEEDED, LENIW (',I1,')',
c     &  'EXCEEDS LIW (',I2,')'
      endif
      if(nerr .eq. 19)then
c        WRITE (LUN, *)'LSODE--  RTOL(',I1,') IS ',I2,' .LT. 0.0 '         
      endif
      if(nerr .eq. 20)then
c        WRITE (LUN, *)'LSODE--  ATOL(',I1,') IS ',R1,' .LT. 0.0'         
      endif
      if(nerr .eq. 21)then
c        WRITE (LUN, *)'LSODE--  EWT(',I1,') IS ',R1,' .LE. 0.0'        
      endif
      if(nerr .eq. 22)then
c        WRITE (LUN, *)'LSODE--  TOUT (',R1,') TOO CLOSE TO TT ',
c     &  '(',R2,') TO START INTEGRATION.'
      endif
      if(nerr .eq. 23)then
c        WRITE (LUN, *)'LSODE--  ITASK = ',I1,' AND TOUT (',R1,
c     &  ') BEHIND TCUR - HU (= ',R2,').'
      endif
      if(nerr .eq. 24)then
c        WRITE (LUN, *)'LSODE--  ITASK = 4 OR 5 AND TCRIT (',R1,
c     &  ') BEHIND TCUR (',R2,').'    
      endif
      if(nerr .eq. 25)then
c        WRITE (LUN, *)'LSODE--  ITASK = 4 OR 5 AND TCRIT (',R1,
c     &  ') BEHIND TOUT (',R2,')'
      endif
      if(nerr .eq. 26)then
c        WRITE (LUN, *)'LSODE--  AT START OF PROBLEM, TOO MUCH ACCURACY'    
      endif
      if(nerr .eq. 27)then
c        WRITE (LUN, *)'LSODE--  TROUBLE FROM INTDY. ITASK = ',I1, 
c     &  'TOUT = ',R1 
      endif
      if(nerr .eq. 302)then
c        WRITE (LUN, *)'LSODE--  REPEATED OCCURRENCES OF ILLEGAL INPUT'     
      endif
      if(nerr .eq. 303)then
c        WRITE (LUN, *)'LSODE--  RUN ABORTED.. APPARENT INFINITE LOOP'      
      endif
      if(nerr .eq. 51)then
c        WRITE (LUN, *)'INTDY--  K (',I1,') ILLEGAL'                   
      endif
      if(nerr .eq. 52)then
c        WRITE (LUN, *)'INTDY--  TT (',R1,') ILLEGAL; TT NOT IN INTERVAL
c     &  ','TCUR - HU (= ',R1,') TO TCUR (',R2,')'   
      endif
  
c      IF (NI .EQ. 1) WRITE (LUN, 20) I1                                 
c20    FORMAT(6X,23HIN ABOVE MESSAGE,  I1 =,I10)
c      IF (NI .EQ. 2) WRITE (LUN, 30) I1,I2                              
c 30   FORMAT(6X,23HIN ABOVE MESSAGE,  I1 =,I10,3X,4HI2 =,I10)           
c      IF (NR .EQ. 1) WRITE (LUN, 40) R1                                 
c 40   FORMAT(6X,23HIN ABOVE MESSAGE,  R1 =,D21.13)                      
c      IF (NR .EQ. 2) WRITE (LUN, 50) R1,R2                              
c 50   FORMAT(6X,15HIN ABOVE,  R1 =,D21.13,3X,4HR2 =,D21.13)             
C ABORT THE RUN IF IERT = 2. -------------------------------------------
 100  IF (IERT .NE. 2) RETURN                                           
      STOP                                                              
C----------------------- END OF SUBROUTINE XERRWV ----------------------
      END                                                               
      SUBROUTINE XSETF (MFLAG)                                          
C                                                                       
C THIS ROUTINE RESETS THE PRINT CONTROL FLAG MFLAG.                     
C                                                                       
      INTEGER MFLAG, MESFLG, LUNIT                                      
      COMMON /EH0001/ MESFLG, LUNIT                                     
C                                                                       
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) MESFLG = MFLAG                
      RETURN                                                            
C----------------------- END OF SUBROUTINE XSETF -----------------------
      END                                                               
      SUBROUTINE XSETUN (LUN)                                           
C                                                                       
C THIS ROUTINE RESETS THE LOGICAL UNIT NUMBER FOR MESSAGES.             
C                                                                       
      INTEGER LUN, MESFLG, LUNIT                                        
      COMMON /EH0001/ MESFLG, LUNIT                                     
C                                                                       
      IF (LUN .GT. 0) LUNIT = LUN                                       
      RETURN                                                            
C----------------------- END OF SUBROUTINE XSETUN ----------------------
      END                                                               
      BLOCK DATA INIT                                                       
C-----------------------------------------------------------------------
C THE FOLLOWING DATA STATEMENT LOADS VARIABLES INTO THE INTERNAL COMMON 
C BLOCKS USED BY LSODE.  THE VARIABLES ARE DEFINED AS FOLLOWS..         
C   ILLIN  = COUNTER FOR THE NUMBER OF CONSECUTIVE TIMES LSODE WAS      
C            CALLED WITH ILLEGAL INPUT.  THE RUN IS STOPPED WHEN        
C            ILLIN REACHES 5.                                           
C   NTREP  = COUNTER FOR THE NUMBER OF CONSECUTIVE TIMES LSODE WAS      
C            CALLED WITH ISTATE = 1 AND TOUT = TT.  THE RUN IS STOPPED   
C            WHEN NTREP REACHES 5.                                      
C   MESFLG = FLAG TO CONTROL PRINTING OF ERROR MESSAGES.  1 MEANS PRINT,
C            0 MEANS NO PRINTING.                                       
C   LUNIT  = DEFAULT VALUE OF LOGICAL UNIT NUMBER FOR PRINTING OF ERROR 
C            MESSAGES.                                                  
C-----------------------------------------------------------------------
      INTEGER ILLIN, IDUMA, NTREP, IDUMB, IOWNS, ICOMM, MESFLG, LUNIT   
      DOUBLE PRECISION ROWND, ROWNS, RCOMM                              
      COMMON /LS0001/ ROWND, ROWNS(210), RCOMM(7),                      
     1   ILLIN, IDUMA(10), NTREP, IDUMB(2), IOWNS(6), ICOMM(13)         
      COMMON /EH0001/ MESFLG, LUNIT                                     
      DATA ILLIN/0/, NTREP/0/                                           
C     DATA MESFLG/1/, LUNIT/3/                                          
      DATA MESFLG/1/, LUNIT/6/                                          
C                                                                       
C----------------------- END OF BLOCK DATA -----------------------------
      END