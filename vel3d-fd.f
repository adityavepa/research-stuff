      SUBROUTINE VEL3D_FD (RAD,ACELL,CA,RMU,TOL,DT,DTP,NTS,C,N2,NVFR,
     >IV,NNEI,IN,NVT,NVTP,PHI,FT2,V,VP,XNA,QW,QWP,QM)
      implicit real*8 (a-h,O-Z)
      INCLUDE 'nvfrmax.h'
c     sets NPCH and NVTCH:
      INCLUDE 'CH3D.h'
      INCLUDE 'NVTMAX.h'
      INCLUDE 'NTMAX.h'
      INCLUDE 'NEIBMAX.h'
      PARAMETER (NPCMAX=50, NVCMAX=10000,
     >NP0=NPCMAX-1, NVSMAX=NTMAX/2+2)
      REAL*8 L1CH,L1,L2,NTOMAX
      DIMENSION U(3),X0CH(3,NPCH),E1CH(3,NPCH),L1CH(NPCH),
     >XNWCH(3,NPCH+2),INBCH(NPCH+1),INFCH(NPCH+1),VCH(3,NVTCH),
     >PHICH(NVTCH),QWCH(3,NVTCH)

      DIMENSION V(3,*),VP(3,*),VOLD(3,NVTMAX),PHI(NVTMAX),XNA(3,*),
     >QW(3,*),QWP(3,*),QWOLD(3,NVTMAX),IM(NVTMAX,0:NP0+2),SX(3,NVTMAX),
     >SY(3,NVTMAX),XCW(3),INB(0:NP0+2),INF(0:NP0+2),INBH(0:NP0+2),
     >INFH(0:NP0+2),INV(NVTMAX),XNW(3,NP0+2),F(3,NVTMAX),DM(3,3)
      DIMENSION L1(NP0), X0(3,NP0+1), N1(NP0),
     >E1(3,NP0), E2(3),DX(3),SM(3),STR(3,3)

      DIMENSION R(3,NVTMAX),P(3,NVTMAX),A(3,NVTMAX),
     >RT(3,NVTMAX),PT(3,NVTMAX),AT(3,NVTMAX),
     >SXT(3,NVTMAX),SYT(3,NVTMAX)

      DIMENSION INB1(NP0),INF1(NP0),V1(2,NVCMAX),
     >DS1(NP0), XMC(2,NVFRMAX),IDT(3,20*NVFRMAX)

      DIMENSION URBP(3),OMRBP(3)

      DIMENSION XP(3),IV(3,*),NNEI(*), IN(NEIBMAX,*),
     >QM(*),FIT(3,NVSMAX),FT2(*),CUR(NVSMAX),DY(3),
     >CC(2,3), UD(3)

      LOGICAL SOL (NPCMAX)

      COMMON / XP / XP,RLENG,CURMAX,FSD
      COMMON /CH3D/DEPTH,X0CH,E1CH,L1CH,XNWCH,VCH,PHICH,QWCH,
     >INBCH,INFCH

      NVS=NTS/2+2
      NVSP=NVS
C--------------------------------------------------------------------
      CEST2D=0.8
      EXPRATE= 0.5D0
      QLIM=0.03D0
      NTOMAX= 6.D10
      SHRNK= 1.d0 - 1.d-6
      SHRNK= 1.d0
      pi=dacos(-1.d0)
C----------------------------------------------------------
      DO I=1,NVS
      PHI(I)=0.D0
      ENDDO

      FSD=0.D0
      DO JT=1,NTS
      DO K=1,2
      DO J=1,3
      JJ1=IV(K+1,JT)
      JJ2=IV(1,JT)
      CC(K,J)=V(J,JJ1)-V(J,JJ2)
      ENDDO
      ENDDO
      CX1=CC(1,2)*CC(2,3)-CC(2,2)*CC(1,3)
      CY1=-CC(1,1)*CC(2,3)+CC(2,1)*CC(1,3)
      CZ1=CC(1,1)*CC(2,2)-CC(2,1)*CC(1,2)
      DS= .5D0*DSQRT(CX1**2+CY1**2+CZ1**2)
      FSD=FSD+DS
      DO KJ=1,3
      I= IV(KJ,JT)
      PHI(I)=PHI(I)+DS/3.D0
      ENDDO
      ENDDO

      DO K=1,3
      XP(K)=0
      ENDDO

      DO I=1,NVS
      DO K=1,3
      XP(K)=XP(K)+V(K,I)*PHI(I)
      ENDDO
      ENDDO

      DO K=1,3
      XP(K)=XP(K)/FSD
      ENDDO

      PRINT *, 'XP=', XP

      DO I=1,NVS
C!
      V(3,I)=V(3,I)-XP(3)
      ENDDO
C!
      XP(3)=0.D0

321   FORMAT ('ZONE N=', I5, '  E=', I5, '  F=FEPOINT, ET=TRIANGLE')

866    FORMAT (1X, 3F13.7)
      OPEN (7, FILE='DROPMESH.DAT')
      WRITE (7,*) 'VARIABLES = "X","Y","Z"'
c      WRITE (7,*) 'VARIABLES = "X","Y"'
      WRITE (7,321)NVS, NTS
         DO I=1, NVS
         WRITE (7,866) V(1,I),V(2,I),V(3,I)
C         WRITE (7,86) V(1,I),V(2,I)
         ENDDO
         DO NT=1,NTS
         WRITE (7,*) IV(1,NT), IV(2,NT), IV(3,NT)
        ENDDO
      close (7)
      INB(0)=1
      INF(0)=NVS
      INBH(0)=1
      INFH(0)=NVS
      print *, ' start gim'
      CALL GIM (1,NVS,V,IN,NNEI,FIT,XNA)
      print *, ' end gim'
c     volume test:
      print *, ' XP=', XP
      vol=0.d0
      do i=1,nvs
      rnx=0.d0
      do k=1,3
      rnx=rnx + (V(K,I)-XP(K))*XNA(K,I)
      ENDDO
      vol=vol+rnx*phi(i)/3.d0
      enddo
      volt= 4.d0*pi*RAD**3/3.d0
c      print *, ' vol=', vol, ' volt=', volt
      print *, ' VOL=', vol
      print *, ' V0 =',volt

      factor2 = (volt/vol)**(1.d0/3.d0)

C      volume renormalization (lucas)
      do i=1,nvs
        phi(i) = phi(i)*factor2**2
        do k=1,3
          v(k,i) = xp(k)+factor2*(v(k,i)-xp(k))
          fit(k,i) = fit(k,i)/factor2
        enddo
      enddo

      CURMAX=0.d0
      DO I=1,NVS
      CUR(I)=-FIT(1,I)-FIT(3,I)
      FT2(I)=FIT(1,I)**2+FIT(3,I)**2+FIT(2,I)**2/2.D0
      curmax= max(curmax,CUR(I))
      ENDDO
C      rad= (3.d0*vol/(4.d0*pi))**(1.d0/3.d0)
      CURMAX=CURMAX*RAD
      print *, ' CURMAX=', CURMAX

      NVTOLD=NVT
      NVSOLD=NVS
      DO J=1,NVT
      DO K=1,3
      VOLD(K,J)=V(K,J)
      QWOLD(K,J)=QW(K,J)
      ENDDO
      ENDDO


      Print *, 'Call CONTOR'
c      Print *, 'X0=', x0
      RD=0.D0
      DO I=1,NVS
      RD1=((V(1,I)-XP(1))**2+(V(2,I)-XP(2))**2+
     >(V(3,I)-XP(3))**2)**(1.D0/2.D0)
      RD=MAX(RD1,RD)
      ENDDO
      print *, ' RD=', RD

      A0=MAX(ACELL*RAD,1.5*RD)
      print *, ' A0=', A0
C      CALL CONTOR (XP,A0,X0CH,E1CH,L1CH,NPC,X0,SOL)
C      infinitely long channel changes
      X0(1,1) = XP(1) - A0
C      for channel height of 1
      X0(2,1) = 1

      X0(1,2) = X0(1,1)
      X0(2,2) = 0

      X0(1,4) = XP(1) + A0
      X0(2,4) = 1

      X0(1,3) = X0(1,4)
      X0(2,3) = 0
C      set number of corner points + 1 to 5
      NPC = 5

      Print *, 'CONTOR called'

      OPEN (13, FILE='CELL.DAT')
      DO IP=1,NPC
      WRITE (13,*) X0(1,IP), X0(2,IP)
      ENDDO
      CLOSE(13)

      NP=NPC-1
      WRITE (4,*) ' NP=', NP
      PRINT *, ' NP=', NP, 'NPC= ', NPC

      E2(1)=0.D0
      E2(2)=0.D0
      E2(3)=1.D0
      TOTLENG=0.D0
      DO IP=1,NP
      JP=IP+1
      IF (JP.GT.NP) JP=JP-NP
      TMP=0.D0
      DO K=1,2
      TMP=TMP+ (X0(K,JP)-X0(K,IP))**2
      ENDDO
      L1(IP)=DSQRT(TMP)
c      PRINT *, 'L1', L1(IP)
      DO K=1,2
      E1(K,IP)= (X0(K,JP)-X0(K,IP))/L1(IP)
      ENDDO
      X0(3,IP)= -DEPTH/2.D0
      E1(3,IP)=0.D0
      TOTLENG=TOTLENG+L1(IP)
      N1(IP)= L1(IP)/C +1.d0
      ENDDO

      I= 0
      DO IP=1,NP
      Call PANEL2D(I,NVCMAX,N1(IP),L1(IP),X0(1,IP),
     >E1(1,IP),INB1(IP),INF1(IP),V1,XNW(1,IP),DS1(IP))
      ENDDO
C     SLIGHTLY CONCERVATIVE OVERESTIMATION
      NVT= NVS + 2*(I*N2+2*NVFR)
      IF (NVT.GT.NVTMAX) THEN
      PRINT *, 'NVS=', nvs
      PRINT *, 'NVFR=', nvfr
      PRINT *, 'NVT=', NVT
       PRINT 44, NVT
44     FORMAT ('INCREASE NVTMAX TO ', I7, 'RECOMPILE AND RERUN')
       STOP
      END IF
      L2=DEPTH
      DH= DEPTH/2
      I=NVS
      DO IP=1,NP
      CALL PNL3DCH(I,INB1(IP),INF1(IP),V1,DS1(IP),
     >L2,N2,INB(IP),INF(IP),INBH(IP),INFH(IP),INV,V,PHI)
      ENDDO
      NVT=I

C     OMIT COMMENTED OUT LINES BELOW?
C      DO K=1,2
C      XCW(K)=0.D0
C      ENDDO
C      DO IP=1,NP
C      JP=IP+1
C      IF (JP.GT.NP) JP=1
C      DO K=1,2
C      XCW(K)=XCW(K)+ 0.5D0*(X0(K,IP)+X0(K,JP))*L1(IP)
C      ENDDO
C      ENDDO

C      DO K=1,2
C      XCW(K)=XCW(K)/TOTLENG
C      ENDDO

C      DO IP=1,NP
C      DO K=1,2
C      X0(K,IP)= XCW(K)+SHRNK*(X0(K,IP)-XCW(K))
C      ENDDO
C      ENDDO
      CALL INIPOINT (NP,NVFR,X0,CEST2D,EXPRATE,NTOMAX,C2D,INE,
     >FRAREA,I0,XMC)
      PRINT *, ' C2D=', C2D, ' INE=', INE

C      OMIT LINES BELOW?
C      DO IP=1,NP
C      DO K=1,2
C      X0(K,IP)= XCW(K)+(X0(K,IP)-XCW(K))/SHRNK
C      ENDDO
C      ENDDO

      STOT=0.D0
      AF=   0.1255
      FACTOR=1.D20

      MTSEL=0
      MTB= NVFR-I0
      DO I=NVT+1,NVT+MTB
      PHI(I)=0.D0
      DO K=1,2
      V(K,I)= XMC(K,I-NVT+I0)
      ENDDO
      V(3,I)= -DH
      ENDDO

      PRINT *, 'CALL DTESSTV'
      DO I=1,NVFR-1
      CALL DTESSTV (I,NVFR,XMC,AF,QLIM,FACTOR,I0,V(1,NVT+1),
     >PHI(NVT+1),MTSEL,MTB,IDT,X0,NP)
      ENDDO
      print *, 'DTESSTV CALLED'

      NTTOT= MTSEL
      DO J=1,MTB
      IF (J.GT.NVFR-I0) V(3,NVT+J)= -DH
      STOT=STOT+PHI(NVT+J)
      ENDDO

c     necessary exclusion, if at all, of some possible near-corner nodes i
c     with phi (i)=0:
      MTBGOOD=0
      DO J=1,MTB
      IF (PHI(NVT+J).GT. 1.D-9*STOT/MTB) THEN
        MTBGOOD=MTBGOOD+1
        PHI(NVT+MTBGOOD)= PHI(NVT+J)
        DO K=1,3
        V(K,NVT+MTBGOOD)= V(K,NVT+J)
	ENDDO
      END IF
      ENDDO
      PRINT *, ' EXCL=', MTB-MTBGOOD
      MTB=MTBGOOD

      PRINT *, ' STOT=', STOT, ' NTTOT=', NTTOT, ' MTB=', MTB
      PRINT *, ' FRAREA=', FRAREA
      OPEN (2, FILE='meshFRBK.dat')
      write (2,*) 'variables = "x","y"'
      write (2,321) NVFR, NTTOT
      DO I=1,NVFR
      WRITE (2,86) XMC(1,I),XMC(2,I)
      ENDDO
      DO NT=1,NTTOT
      WRITE (2,*) IDT(1,NT),IDT(2,NT),IDT(3,NT)
      ENDDO
      CLOSE (2)
c      stop
86    FORMAT (1X, 2F13.7)

      INB(NP+1)= NVT+1
      INF(NP+1)= NVT+MTB
      inb(np+2)= NVT+MTB+1
      INF(NP+2)= NVT+2*MTB
      DO IP=NP+1,NP+2
      INBH(IP)=INB(IP)
      INFH(IP)=INF(IP)
      XNW(1,IP)=0.D0
      XNW(2,IP)=0.D0
      ENDDO
      XNW(3,NP+1)= -1.D0
      XNW(3,NP+2)=  1.D0
      DO I=1,MTB
C     CONTINUATION TO FRONT PANEL:
      V(1,NVT+MTB+I)= V(1,NVT+I)
      V(2,NVT+MTB+I)= V(2,NVT+I)
      V(3,NVT+MTB+I)= DH
      PHI(NVT+MTB+I)= PHI(NVT+I)
      INV(NVT+I)= NVT+MTB+I
      INV(NVT+MTB+I)= NVT+I
      ENDDO
      NVT=NVT+2*MTB
      PRINT *, ' NVT=', NVT

C     BUILDING INITIAL APPROXIMATION FOR QW:
      IF (NVTOLD.EQ.0) THEN
	DTP=1.D20
        NVTP=NVT
        NVSP=NVS
	DO I=1,NVT
	DO K=1,3
	QW(K,I)=0.D0
	QWP(K,I)=0.D0
	VP(K,I)=V(K,I)
	ENDDO
	ENDDO

      ELSE
      DO I=1,NVS
      DO K=1,3
c!      QW(K,I)=QW(K,I)+ (QW(K,I)-QWP(K,I))*DT/DTP
      ENDDO
      ENDDO
C     ----------------
	PMIN=1.D20
        DO IP=1,NP+1
        DO I=INBH(IP),INFH(IP)
	RRMIN=1.D20
	DO 66 J=NVSOLD+1,NVTOLD
	IF (VOLD(3,J).GT.0.D0) GOTO 66
	PP=0.D0
	DO K=1,3
	PP=PP+(VOLD(K,J)-V(K,I))**2
	ENDDO
	IF (PP.LT.RRMIN) THEN
	  RRMIN=PP
	  JMIN=J
	END IF
66      CONTINUE

        RRNEXT=1.D20
        DO 67 J=NVSOLD+1,NVTOLD
        IF (VOLD(3,J).GT.0.D0. OR. J.EQ.JMIN) GOTO 67
        PP=0.D0
        DO K=1,3
        PP=PP+(VOLD(K,J)-V(K,I))**2
        ENDDO
        IF (PP.LT.RRNEXT) THEN
          RRNEXT=PP
          JNEXT=J
        END IF
67      CONTINUE
        R2=0.D0
        PP=0.D0
        DO K=1,3
        R2=R2+(VOLD(K,JMIN)-VOLD(K,JNEXT))**2
        PP=PP+ (V(K,I)-VOLD(K,JNEXT))*(VOLD(K,JMIN)-VOLD(K,JNEXT))
        ENDDO
        PP=PP/R2
        PMIN=MIN(PMIN,PP)

	DO K=1,3
        QW(K,I)=PP*QWOLD(K,JMIN)+(1-PP)*QWOLD(K,JNEXT)
        ENDDO

	RRMIN=1.D20
	DO 660 J=NVSP+1,NVTP
	IF (VP(3,J).GT.0.D0) GOTO 660
	PP=0.D0
	DO K=1,3
	PP=PP+(VP(K,J)-V(K,I))**2
	ENDDO
	IF (PP.LT.RRMIN) THEN
	  RRMIN=PP
	  JMIN=J
	END IF
660     CONTINUE

        RRNEXT=1.D20

        DO 670 J=NVSP+1,NVTP
        IF (VP(3,J).GT.0.D0. OR. J.EQ.JMIN) GOTO 670
        PP=0.D0
        DO K=1,3
        PP=PP+(VP(K,J)-V(K,I))**2
        ENDDO
        IF (PP.LT.RRNEXT) THEN
          RRNEXT=PP
          JNEXT=J
        END IF
670     CONTINUE
        R2=0.D0
        PP=0.D0
        DO K=1,3
        R2=R2+(VP(K,JMIN)-VP(K,JNEXT))**2
        PP=PP+ (V(K,I)-VP(K,JNEXT))*(VP(K,JMIN)-VP(K,JNEXT))
        ENDDO
        PP=PP/R2
        PMIN=MIN(PMIN,PP)

	DO K=1,3
        TEMP=PP*QWP(K,JMIN)+(1-PP)*QWP(K,JNEXT)
c!	QW(K,I)=QW(K,I)+ (QW(K,I)-TEMP)*DT/DTP
        ENDDO
        I1=INV(I)
        QW(1,I1)=  QW(1,I)
        QW(2,I1)=  QW(2,I)
        QW(3,I1)= -QW(3,I)
        ENDDO
        ENDDO

        IF (DABS(DT).GT. 1.D-14) THEN
	  NVTP=NVTOLD
	  NVSP=NVSOLD
	  DO I=1, NVTP
	  DO K=1,3
	  VP(K,I)=VOLD(K,I)
	  QWP(K,I)=QWOLD(K,I)
	  ENDDO
	  ENDDO
	  DTP=DT
        END IF
      END IF
      PRINT  *, ' PMIN=', PMIN
C---------END OF INITIAL APPROXIMATON----------------------

      DO K=1,3
      XCW(K)=0.D0
      ENDDO
      SW=0.D0
      DO IP=1,NP+2
      DO I=INB(IP),INF(IP)
      DO K=1,3
      XNA(K,I)= XNW(K,IP)
      XCW(K)=XCW(K)+V(K,I)*PHI(I)
      ENDDO
      SW=SW+ PHI(I)
      ENDDO
      ENDDO

      DO K=1,3
      XCW(K)=XCW(K)/SW
      ENDDO

      CALL DEFLMAT (NVS+1,NVT,XCW,V,PHI,DM)

      PI=DACOS(-1.D0)
c     work with a different sign now (vs. try9cg.f)
      C2= -3.D0/(2.D0*PI)
      print *, ' start IM'
      DO IP=0,NP+1
      DO I=INBH(IP),INFH(IP)
      IF (IP.EQ.0) IM(I,0)=I
      DO JP=1,NP+1
      IF (ip.eq.jp) then
	IM(I,JP)=I
      ELSE
        RRMIN=1.D20
	INBTMP= INBH(JP)
	INFTMP= INFH(JP)
	IF (IP.EQ.0) THEN
	  INBTMP= INB(JP)
	  INFTMP= INF(JP)
	END IF

	DO J=INBTMP,INFTMP
        PP=0.D0
        DO K=1,3
        PP=PP+ (V(K,J)-V(K,I))**2
        ENDDO
        IF (PP.LT.RRMIN) THEN
	  RRMIN=PP
	  IM(I,JP)=J
        END IF
	ENDDO
      END IF
      ENDDO
      IM(I,NP+2)= IM(I,NP+1) + MTB
      ENDDO
      ENDDO

      DO IP=1,NP+1
      DO I=INBH(IP),INFH(IP)
      RRMIN=1.D20
      DO J=1,NVS
      PP=0.D0
      DO K=1,3
      PP=PP+ (V(K,J)-V(K,I))**2
      ENDDO
      IF (PP.LT.RRMIN) THEN
         RRMIN=PP
	  IM(I,0)=J
      END IF
      ENDDO
      ENDDO
      ENDDO

      PRINT *, ' END IM'
c-------------------------------
C      print *, ' start VELCH3D'
      print *, ' start BOUSVEL'

      call bousflux(1.d0,depth,1.d-8,qfact)

      do i=1,nvt
      do k=1,3
      F(k,i)=0.d0
      enddo

C      IF (I.LE.NVS) CALL VELCH3D (V(1,I),F(1,I))
C      infinite long straight channel velocity calcs (BOUS)
      if (i.le.nvs) then
        yb = v(2,i)
        zb = v(3,i)+depth/2.d0
        call bousvel(yb,zb,1.d0,depth,1.d-8,ubuss)
C      v_norm = velocity*flowrate/flux  (from binew3d.f flowrate=QF(1)=1.d0*DEPTH)
        F(1,i) = ubuss*(1.d0*depth)/qfact
      end if
C      F(1,inv(i)) = F(1,i)
      enddo
C      print *, ' VELCH3D finished'
      print *, ' BOUSVEL finished'


      FFACT=-1.D0/(4.D0*PI*CA)

!     F-CALCULATION: SELF-INTERACTION

       RLENG=0.d0
       DO I=1,NVS
       DO J=I+1,NVS
       RR=0.D0
       RNX=0.D0
       RNY=0.D0

       DO K=1,3
       DX(K)=V(K,J)-V(K,I)
       RR=RR+DX(K)**2
       RNX=RNX+DX(K)*XNA(K,J)
       RNY=RNY+DX(K)*XNA(K,I)
       ENDDO
       RLENG= max (RLENG,RR)
       R1= (CUR(J)-CUR(I))*FFACT/DSQRT(RR)
       R1X= R1*PHI(J)
       RNX=RNX*R1X/RR
       R1Y= R1*PHI(I)
       RNY=RNY*R1Y/RR
       DO K=1,3
       F(K,I)= F(K,I) + R1X*XNA(K,J) + RNX*DX(K)
       F(K,J)= F(K,J) - R1Y*XNA(K,I) - RNY*DX(K)
       ENDDO
       ENDDO
       ENDDO

       RLENG= dsqrt(RLENG)/RAD
       print *, ' LENG=', RLENG

!      F-CALCULATION: INTERACTION BTWN THE DROP AND THE MOVING FRAME

       DO IP=1,NP+1
       DO I=INBH(IP),INFH(IP)
       J1=IM(I,0)
       DO J=1,NVS
       RR=0.D0
       RNX=0.D0
       DO K=1,3
       DX(K)=V(K,J)-V(K,I)
       RR=RR+DX(K)**2
       RNX=RNX+DX(K)*XNA(K,J)
       ENDDO
       R1= PHI(J)*FFACT*(CUR(J)-CUR(J1))/DSQRT(RR)
       RNX=RNX*R1/RR
       DO K=1,3
       F(K,I)=F(K,I) + R1*XNA(K,J) + RNX*DX(K)
       ENDDO
       ENDDO
       ENDDO
       ENDDO

      ITTOT=0

775   DO i=1,nvt
      do k=1,3
      R(k,I)= -F(K,I)
      IF (I.LE.NVS) R(K,I)= 2.D0*F(K,I)/(RMU+1.D0)
      enddo
      enddo

C     MODIFIED FOR DROP PRESENCE
C     ASSUMING
      C4= 3.D0/(PI*(RMU+1.D0))

C     HERE, QW AND IM ASSUMED TO BE CONTINUED TO JP=NP+2:
      DO IP=0,NP+1
      DO 330 JP=1,NP+2
      IF (JP.EQ.IP) GOTO 330
      IF (INFH(JP).LT.INBH(JP)) GOTO 330
      ID= 2*(JP-NP) -3
      DO I=INBH(IP),INFH(IP)
      IF (JP.LE.NP) THEN
c        CALL RECSTREX (V(1,I),X0(1,JP),E1(1,JP),L1(JP),E2,L2,
        CALL NEWRECST (V(1,I),X0(1,JP),E1(1,JP),L1(JP),E2,L2,
     >  XNW(1,JP),STR)
      ELSE
	CALL BKFRSTR(V(1,I),ID,NP,L2,X0,E1,L1,STR)
      END IF
      CTEMP= C2
      IF (IP.EQ.0) CTEMP=C4
      DO K=1,3
      DO L=1,3
      R(K,I)=R(K,I)+CTEMP*STR(K,L)*QW(L,IM(I,JP))
      ENDDO
      ENDDO
      ENDDO
330   CONTINUE
      ENDDO

      CALL DEFL(NVS+1,NVT,V,PHI,XNA,XCW,DM,SW,QW,0,-1,URBP,OMRBP,R)
C     INTERACTION BETWEEN SIDE WALLS
34     DO IP=1,NP
      DO JP=IP+1,NP
      DO I=INBH(IP), INFH(IP)
      DO K=1,3
      SY(K,I)=0.D0
      ENDDO
      ENDDO

      DO J=INBH(JP), INFH(JP)
      DO K=1,3
      SX(K,J)=0.D0
      ENDDO
      ENDDO

      DO I=INBH(IP), INFH(IP)
      IMI=IM(I,JP)
      DO J=INBH(JP), INFH(JP)
      IMJ=IM(J,IP)
      PP=0.D0
      PX=0.D0
      PY=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (QW(K,J)-QW(K,IMI))*DX(K)
      PY=PY+ (QW(K,I)-QW(K,IMJ))*DX(K)
      ENDDO
      R3M=V(3,J)-V(3,I)
      R3P=V(3,J)+V(3,I)
      RR=PP+R3M**2
      R5= 1.D0/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= 1.D0/(RR**2*DSQRT(RR))

      RQX=  (PX+R3M*(QW(3,J)-QW(3,IMI)))*R5*PHI(J)
      RQXST=  (PX+R3P*(QW(3,J)+QW(3,IMI)))*R5ST*PHI(J)

      RQY= -(PY+R3M*(QW(3,I)-QW(3,IMJ)))*R5*PHI(I)
      RQYST= (-PY+R3P*(QW(3,I)+QW(3,IMJ)))*R5ST*PHI(I)

      RX=RQX+RQXST
      RY=RQY+RQYST
      DO K=1,2
      SY(K,I)=SY(K,I)+RX*DX(K)
      SX(K,J)=SX(K,J)+RY*DX(K)
      ENDDO
      SY(3,I)=SY(3,I)+RQX*R3M-RQXST*R3P
      SX(3,J)=SX(3,J)+RQY*R3M+RQYST*R3P
      ENDDO
      ENDDO

      IY=INBH(IP)
      IX=INBH(JP)

      DO I=INBH(IP),INFH(IP)
      RNX=0.D0
      DO K=1,3
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*C2
      DO K=1,3
      R(K,I)=R(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO

      DO J=INBH(JP),INFH(JP)
      RNY=0.D0
      DO K=1,3
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*C2
      DO K=1,3
      R(K,J)=R(K,J)+ RNY*SX(K,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

c     NEW:
C     INTERACTION BETWEEN SIDE (IP from 1 to NP) AND TOP/BOTTOM (jp=NP+1,np+2)
C     PANELS  (ASSUMING THAT PANEL np+1 IS BACK (V(3,J)= -DH), AND PANEL NP+2
C     IS FRONT (V(3,J)=DH); CONTRIBUTIONS FROM FRONT AND BACK PANELS ARE JOINED;
C     DEPTH IS FULL CHANEL DEPTH
      JP=NP+1
      DH=DEPTH/2

      DO IP=1,NP

      DO J=INB(JP), INF(JP)
      DO K=1,3
      SX(K,J)=0.D0
      ENDDO
      ENDDO

      DO I=INBH(IP), INFH(IP)
      IMI=IM(I,JP)
      R3M= -DH-V(3,I)
      R3P= V(3,I)-DH
      RNX= -R3M*C2
      RNXST= -R3P*C2
      DO J=INB(JP), INF(JP)
      IMJ=IM(J,IP)
      PP=0.D0
      PX=0.D0
      PY=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (QW(K,J)-QW(K,IMI))*DX(K)
      PY=PY+ (QW(K,I)-QW(K,IMJ))*DX(K)
      ENDDO
C     R3M=V(3,J)-V(3,I)
C     R3P=V(3,J)+V(3,I)
      RR=PP+R3M**2
      R5= 1.D0/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= 1.D0/(RR**2*DSQRT(RR))

      RQX=  (PX+R3M*(QW(3,J)-QW(3,IMI)))*RNX*R5*PHI(J)
      RQXST=  (PX+R3P*(QW(3,J)-QW(3,IMI)))*RNXST*R5ST*PHI(J)

      RQY= -(PY+R3M*(QW(3,I)-QW(3,IMJ)))*R5*PHI(I)
      RQYST= (-PY+R3P*(QW(3,I)+QW(3,IMJ)))*R5ST*PHI(I)

      RX=RQX+RQXST
      RY=RQY+RQYST
      DO K=1,2
      R(K,I)=R(K,I)+RX*DX(K)
      SX(K,J)=SX(K,J)+RY*DX(K)
      ENDDO
      R(3,I)=R(3,I)+RQX*R3M-RQXST*R3P
      SX(3,J)=SX(3,J)+RQY*R3M+RQYST*R3P
      ENDDO
      ENDDO

      IY=INBH(IP)

      DO J=INB(JP),INF(JP)
C      SIGN CORRECT
      RNY=0.D0
      DO K=1,3
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*C2
      DO K=1,3
      R(K,J)=R(K,J)+ RNY*SX(K,J)
      ENDDO
      ENDDO
      ENDDO

c     INTERACTION BETWEEN BACK AND FRONT PANELS:
      D2=DEPTH**2
      IP=NP+1
      C22= C2*DEPTH
      DO I=INB(IP), INF(IP)
      DO J=I+1, INF(IP)
      RR=D2
      RQ=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      RR=RR+ DX(K)**2
      RQ=RQ+ (QW(K,J)-QW(K,I))*DX(K)
      ENDDO
      R5= C22/(RR**2*DSQRT(RR))
      TMP= DEPTH*(QW(3,J)-QW(3,I))
      RX= (RQ-TMP)*R5*PHI(J)
      RY= (RQ+TMP)*R5*PHI(I)
      DO K=1,2
      R(K,I)=R(K,I)+RX*DX(K)
      R(K,J)=R(K,J)-RY*DX(K)
      ENDDO
      R(3,I)=R(3,I) +RX*DEPTH
      R(3,J)=R(3,J) +RY*DEPTH
      ENDDO
      ENDDO

c     new: INSERT #1R: DROP TO ALL SOLID WALLS CONTRIBUTION:
C     ASSUMING
      C3= -(RMU-1.D0)*3.D0/(4.D0*PI)
      DO IP=1,NP+1
      DO I=INBH(IP), INFH(IP)
      J1= IM(I,0)
      DO J=1,NVS
      PP=0.D0
      PX=0.D0
      RNX=0.D0
      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      RNX=RNX+DX(K)*XNA(K,J)
      PX=PX+ (QW(K,J)-QW(K,J1))*DX(K)
      ENDDO
      RNX=C3*RNX*PHI(J)/(PP**2*DSQRT(PP))
      RNX=RNX*PX
      DO K=1,3
      R(K,I)=R(K,I)+RNX*DX(K)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C------------------------------------------------
C     INSERT #2R: CONTRIBUTION FROM SIDE PANELS TO DROP
C     ASSUMING
      C4= 3.D0/((RMU+1.D0)*PI)
      DO JP=1,NP
      IX=INBH(JP)
      DO I= 1,NVS
      RNX=0.D0
      DO K=1,3
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*C4
      DO K=1,3
      SY(K,I)=0.D0
      ENDDO
      J1=IM(I,JP)
      DO J=INBH(JP), INFH(JP)
      PP=0.D0
      PX=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (QW(K,J)-QW(K,J1))*DX(K)
      ENDDO
      R3M=V(3,J)-V(3,I)
      R3P=V(3,J)+V(3,I)
      RR=PP+R3M**2
      R5= PHI(J)/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= PHI(J)/(RR**2*DSQRT(RR))

      RQX=  (PX+R3M*(QW(3,J)-QW(3,J1)))*R5
      RQXST=  (PX+R3P*(QW(3,J)+QW(3,J1)))*R5ST

      RX=RQX+RQXST
      DO K=1,2
      SY(K,I)=SY(K,I)+RX*DX(K)
      ENDDO
      SY(3,I)=SY(3,I)+RQX*R3M-RQXST*R3P
      ENDDO

      DO K=1,3
      R(K,I)=R(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO
      ENDDO
C-------------------------------------------------------------
cc      insert #3R: front and back walls to drop contribution
c     assuming
      C4= 3.D0/((RMU+1.D0)*PI)
      DO I=1,NVS
      J1=IM(I,NP+1)
      R3M= -DH-V(3,I)
      R3P= V(3,I)-DH
      RNX= -R3M*C4
      RNXST= -R3P*C4

      DO J=INB(NP+1), INF(NP+1)
      PP=0.D0
      PX=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (QW(K,J)-QW(K,J1))*DX(K)
      ENDDO
      RR=PP+R3M**2
      R5= PHI(J)/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= PHI(J)/(RR**2*DSQRT(RR))

      RQX=  (PX+R3M*(QW(3,J)-QW(3,J1)))*RNX*R5
      RQXST=  (PX+R3P*(QW(3,J)-QW(3,J1)))*RNXST*R5ST
      RX=RQX+RQXST

      DO K=1,2
      R(K,I)=R(K,I)+RX*DX(K)
      ENDDO

      R(3,I)=R(3,I)+RQX*R3M-RQXST*R3P
      ENDDO
      ENDDO
C---------------------------------------------------------
cc      insert #4R: drop self-interaction
c     assuming
      C5= (RMU-1.D0)/(RMU+1.D0)
      C6= C5*3.D0/(2.D0*PI)
      DO I=1,NVS
C     ADDED-BACK TERM (DUE TO SINGYLARITY SUBTRACTION)
      DO K=1,3
      R(K,I)=R(K,I)+C5*QW(K,I)
      ENDDO
      DO J=I+1,NVS
      RR=0.D0
      RQ=0.D0
      RNX=0.D0
      RNY=0.D0
      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      RR=RR+DX(K)**2
      RQ=RQ+DX(K)*(QW(K,J)-QW(K,I))
      RNX=RNX+ DX(K)*XNA(K,J)
      RNY=RNY+ DX(K)*XNA(K,I)
      ENDDO
      R5=C6*RQ/(RR**2*DSQRT(RR))
      RX=R5*RNX*PHI(J)
      RY=R5*RNY*PHI(I)
      DO K=1,3
      R(K,I)=R(K,I) +RX*DX(K)
      R(K,J)=R(K,J) +RY*DX(K)
      ENDDO
      ENDDO
      ENDDO
C------------------------------------------------------------
      do IP=1,NP+1
      DO I=INBH(IP),INFH(IP)
      I1=INV(I)
      R(1,I1)=  R(1,I)
      R(2,I1)=  R(2,I)
      R(3,I1)= -R(3,I)
      ENDDO
      ENDDO

      DO I=1,NVT
      DO K=1,3
      QWNEW=R(K,I)
      R(K,I)=QW(K,I)-R(K,I)
      P(K,I)=   R(K,I)
      RT(K,I)=  R(K,I)
      PT(K,I)=  RT(K,I)
      ENDDO
      ENDDO

      IT=0

      erra=1.d10

777   IT=IT+1
      ITTOT=ITTOT+1

C     CONVERGENCE CRITERION NOT MODIFIED YET
      ERRMAX=0.D0
      RR=0.D0
      ERRA=0.D0
      QW2=0.D0
      SPHI=0.D0
      DO IP=0,NP+1
      DO I= INBH(IP),INFH(IP)
      ERR=0.D0
      TEMP=0.D0
      DO K=1,3
      RR=RR+R(K,I)*RT(K,I)*PHI(I)
      ERR=ERR+R(K,I)**2
      TEMP=TEMP+QW(K,I)**2
      ENDDO
      ERRA=ERRA+ERR*PHI(I)
      QW2= QW2+TEMP*PHI(I)
      SPHI=SPHI+PHI(I)
      ERRMAX=MAX(ERRMAX,ERR)
      ENDDO
      ENDDO

      QW2= DSQRT(QW2/SPHI)
      IF (ITTOT.GT.1) THEN
        ERRMAX=DSQRT(ERRMAX)/QW2
        ERRA= DSQRT(ERRA/SPHI)/QW2

c	PRINT *, ' ITTOT=', ITTOT
	PRINT *, 'IT=', IT, ' ERRMAX=', ERRMAX, ' ERRA=', ERRA
	PRINT *, ' QW2=', QW2
	IF (ERRA.LT.TOL) GOTO 778
      END IF

      DO I=1,NVT
      DO K=1,3
      A(K,I)=0.D0
      AT(K,I)=0.D0
      ENDDO
      ENDDO

!      ..............CALCULATE A AND AT
C     MODIFIED FOR DROP PRESENCE
C     ASSUMING
      C4= 3.D0/(PI*(RMU+1.D0))

      DO IP=0,NP+1
      DO 430 JP=1,NP+2
      IF (JP.EQ.IP) GOTO 430
      IF (INFH(JP).LT.INBH(JP)) GOTO 430
      ID= 2*(JP-NP) -3
      DO I=INBH(IP),INFH(IP)
      IF (JP.LE.NP) THEN
c        CALL RECSTREX (V(1,I),X0(1,JP),E1(1,JP),L1(JP),E2,L2,
        CALL NEWRECST (V(1,I),X0(1,JP),E1(1,JP),L1(JP),E2,L2,
     >  XNW(1,JP),STR)
      ELSE
	CALL BKFRSTR(V(1,I),ID,NP,L2,X0,E1,L1,STR)
      END IF
      J=IM(I,JP)
      J2=J
      SIG=1.D0
C      THIS WAY, EASIER TO GENERALIZE FOR THE PRESENCE OF DROP/PARTICLE
      IF (J.GT. INFH(JP).OR. JP.EQ.NP+2) THEN
        J2=INV(J)
        SIG= -1.D0
      END IF

C     P MUST BE CONTINUED
      CTEMP=C2
      IF (IP.EQ.0) CTEMP=C4
      DO K=1,3
      DO L=1,3
      A(K,I)=A(K,I)+CTEMP*STR(K,L)*P(L,J)
      TEMP=CTEMP*STR(K,L)*PT(K,I)*PHI(I)
      IF (L.EQ.3) TEMP=SIG*TEMP
      AT(L,J2)=AT(L,J2)+ TEMP
      ENDDO
      ENDDO
      ENDDO
430   CONTINUE
      ENDDO

!     INTERACTION BETWEEN side WALLS
      DO IP=1,NP
      DO JP=IP+1,NP

      IY=INBH(IP)
      IX=INBH(JP)

      DO I=INBH(IP),INFH(IP)
      RNX=0.D0
      DO K=1,3
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*C2
      DO K=1,3
      SYT(k,i)= RNX*PT(K,I)*PHI(I)
      SY(K,I)=0.D0
      ENDDO
      ENDDO

      DO J=INBH(JP),INFH(JP)
      RNY=0.D0
      DO K=1,3
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*C2
      DO K=1,3
      SXT(K,J)= RNY*PT(K,J)*PHI(J)
      SX(K,J)=0.D0
      ENDDO
      ENDDO

      DO I=INBH(IP), INFH(IP)
      DO J=INBH(JP), INFH(JP)
      I1=IM(J,IP)
      J1=IM(I,JP)

      PP=0.D0
      PX=0.D0
      PY=0.D0
      SYTR=0.D0
      SXTR=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (P(K,J)-P(K,J1))*DX(K)
      PY=PY+ (P(K,I)-P(K,I1))*DX(K)
      SYTR=SYTR+SYT(K,I)*DX(K)
      SXTR=SXTR+SXT(K,J)*DX(K)
      ENDDO
      R3M=V(3,J)-V(3,I)
      R3P=V(3,J)+V(3,I)
      RR=PP+R3M**2
      R5= 1.D0/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= 1.D0/(RR**2*DSQRT(RR))

      R7=  R5*R3M-R5ST*R3P
      R8=  R5*R3M+R5ST*R3P
      R9=  R5*R3M**2-R5ST*R3P**2
      R10= R5*R3M**2+R5ST*R3P**2

      R11= ( (R5+R5ST)*SYTR+ R7*SYT(3,I) )*PHI(J)
      R12= ( (R5+R5ST)*SXTR+ R8*SXT(3,J) )*PHI(I)
      RQX=  (PX+R3M*(P(3,J)-P(3,J1)))*R5*PHI(J)
      RQXST=  (PX+R3P*(P(3,J)+P(3,J1)))*R5ST*PHI(J)

      RQY= -(PY+R3M*(P(3,I)-P(3,I1)))*R5*PHI(I)
      RQYST= (-PY+R3P*(P(3,I)+P(3,I1)))*R5ST*PHI(I)

      RX=RQX+RQXST
      RY=RQY+RQYST
      DO K=1,2
      SY(K,I)=SY(K,I)+RX*DX(K)
      SX(K,J)=SX(K,J)+RY*DX(K)
      TEMP=R11*DX(K)
      AT(K,J)=AT(K,J)+TEMP
      AT(K,J1)=AT(K,J1)-TEMP
      TEMP= R12*DX(K)
      AT(K,I)=AT(K,I)-TEMP
      AT(K,I1)=AT(K,I1)+TEMP
      ENDDO
      SY(3,I)=SY(3,I)+RQX*R3M-RQXST*R3P
      SX(3,J)=SX(3,J)+RQY*R3M+RQYST*R3P
      AT(3,J)= AT(3,J) + (R8*SYTR+R9*SYT(3,I))*PHI(J)
      AT(3,J1)=AT(3,J1)- (R7*SYTR+R10*SYT(3,I))*PHI(J)
      AT(3,I)= AT(3,I) - (R7*SXTR+R9*SXT(3,J))*PHI(I)
      AT(3,I1)=AT(3,I1)+ (R8*SXTR+R10*SXT(3,J))*PHI(I)

      ENDDO
      ENDDO

      IY=INBH(IP)
      IX=INBH(JP)

      DO I=INBH(IP),INFH(IP)
      RNX=0.D0
      DO K=1,3
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*C2
      DO K=1,3
      A(K,I)=A(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO

      DO J=INBH(JP),INFH(JP)
      RNY=0.D0
      DO K=1,3
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*C2
      DO K=1,3
      A(K,J)=A(K,J)+ RNY*SX(K,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C     INTERACTION BETWEEN SIDE AND FRONT/BACK PANELS:

      DO IP=1,NP
      JP=NP+1

      IY=INBH(IP)

      DO I=INBH(IP),INFH(IP)
      DO K=1,3
      SYT(k,i)= PT(K,I)*PHI(I)
      ENDDO
      ENDDO

      DO J=INBH(JP),INFH(JP)
      RNY=0.D0
      DO K=1,3
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*C2
      DO K=1,3
      SXT(K,J)= RNY*PT(K,J)*PHI(J)
      SX(K,J)=0.D0
      ENDDO
      ENDDO

      DO I=INBH(IP), INFH(IP)
      J1=IM(I,JP)
      R3M= -DH-V(3,I)
      R3P= V(3,I)-DH
      RNX= -R3M*C2
      RNXST= -R3P*C2

      DO J=INB(JP), INF(JP)
      I1=IM(J,IP)
      PP=0.D0
      PX=0.D0
      PY=0.D0
      SYTR=0.D0
      SXTR=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (P(K,J)-P(K,J1))*DX(K)
      PY=PY+ (P(K,I)-P(K,I1))*DX(K)
      SYTR=SYTR+ SYT(K,I)*DX(K)
      SXTR=SXTR+ SXT(K,J)*DX(K)
      ENDDO
C     R3M=V(3,J)-V(3,I)
C     R3P=V(3,J)+V(3,I)
      RR=PP+R3M**2
      R5= 1.D0/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= 1.D0/(RR**2*DSQRT(RR))

      RQX=  (PX+R3M*(P(3,J)-P(3,J1)))*RNX*R5*PHI(J)
      RQXST=  (PX+R3P*(P(3,J)-P(3,J1)))*RNXST*R5ST*PHI(J)
      RQY= -(PY+R3M*(P(3,I)-P(3,I1)))*R5*PHI(I)
      RQYST= (-PY+R3P*(P(3,I)+P(3,I1)))*R5ST*PHI(I)

      RX=RQX+RQXST
      RY=RQY+RQYST

      R11= ( (RNX*R5+RNXST*R5ST)*SYTR+ SYT(3,I)*
     >(RNX*R5*R3M -RNXST*R5ST*R3P) )*PHI(J)

      R7=  R5*R3M-R5ST*R3P
      R8=  R5*R3M+R5ST*R3P
      R9=  R5*R3M**2-R5ST*R3P**2
      R10= R5*R3M**2+R5ST*R3P**2

      R12= ( (R5+R5ST)*SXTR+ R8*SXT(3,J) )*PHI(I)

      DO K=1,2
      A(K,I)=A(K,I)+RX*DX(K)
      SX(K,J)=SX(K,J)+RY*DX(K)
      TEMP=R11*DX(K)
      AT(K,J)=AT(K,J)+TEMP
      AT(K,J1)=AT(K,J1)-TEMP
      TEMP= R12*DX(K)
      AT(K,I)=AT(K,I)-TEMP
      AT(K,I1)=AT(K,I1)+TEMP
      ENDDO

      A(3,I)=A(3,I)+RQX*R3M-RQXST*R3P
      SX(3,J)=SX(3,J)+RQY*R3M+RQYST*R3P
      TEMP= (SYTR*(R3M*RNX*R5+R3P*RNXST*R5ST)+
     >SYT(3,I)*(R3M**2*RNX*R5-R3P**2*RNXST*R5ST))*PHI(J)
      AT(3,J)=AT(3,J)+TEMP
      AT(3,J1)=AT(3,J1)-TEMP
      AT(3,I)= AT(3,I) - (R7*SXTR+R9*SXT(3,J))*PHI(I)
      AT(3,I1)=AT(3,I1)+ (R8*SXTR+R10*SXT(3,J))*PHI(I)
      ENDDO
      ENDDO

      IY=INBH(IP)

      DO J=INB(JP),INF(JP)
C      SIGN CORRECT
      RNY=0.D0
      DO K=1,3
      RNY=RNY-(V(K,IY)-V(K,J))*XNW(K,IP)
      ENDDO
      RNY=RNY*C2
      DO K=1,3
      A(K,J)=A(K,J)+ RNY*SX(K,J)
      ENDDO
      ENDDO
      ENDDO

c     new: INSERT #1: DROP TO ALL SOLID WALLS CONTRIBUTION:
C     ASSUMING
      C3= -(RMU-1.D0)*3.D0/(4.D0*PI)
      DO IP=1,NP+1
      DO I=INBH(IP), INFH(IP)
      J1= IM(I,0)
      DO J=1,NVS
      PP=0.D0
      PX=0.D0
      PYT=0.D0
      RNX=0.D0
      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      RNX=RNX+DX(K)*XNA(K,J)
      PX=PX+ (P(K,J)-P(K,J1))*DX(K)
      PYT=PYT+ PT(K,I)*DX(K)
      ENDDO
      RNX=C3*RNX*PHI(J)/(PP**2*DSQRT(PP))
      RNXT=RNX*PYT*PHI(I)
      RNX=RNX*PX
      DO K=1,3
      A(K,I)=A(K,I)+RNX*DX(K)
      AT(K,J)=  AT(K,J)  + RNXT*DX(K)
      AT(K,J1)= AT(K,J1) - RNXT*DX(K)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C------------------------------------------------
C     INSERT #2: CONTRIBUTION FROM SIDE PANELS TO DROP
C     ASSUMING
      C4= 3.D0/((RMU+1.D0)*PI)
      DO JP=1,NP
      IX=INBH(JP)

      DO I= 1,NVS
      RNX=0.D0
      DO K=1,3
      RNX=RNX+(V(K,IX)-V(K,I))*XNW(K,JP)
      ENDDO
      RNX=RNX*C4
      DO K=1,3
      SYT(k,i)= RNX*PT(K,I)*PHI(I)
      SY(K,I)=0.D0
      ENDDO
      J1=IM(I,JP)
      SIG=1.D0
      J2=J1
      IF (J1.GT. INFH(JP)) THEN
	J2=INV(J1)
	SIG= -1.D0
      END IF

      DO J=INBH(JP), INFH(JP)

      PP=0.D0
      PX=0.D0
      SYTR=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (P(K,J)-P(K,J1))*DX(K)
      SYTR=SYTR+SYT(K,I)*DX(K)
      ENDDO
      R3M=V(3,J)-V(3,I)
      R3P=V(3,J)+V(3,I)
      RR=PP+R3M**2
      R5= PHI(J)/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= PHI(J)/(RR**2*DSQRT(RR))

      R7=  R5*R3M-R5ST*R3P
      R8=  R5*R3M+R5ST*R3P
      R9=  R5*R3M**2-R5ST*R3P**2
      R10= R5*R3M**2+R5ST*R3P**2

      R11=  (R5+R5ST)*SYTR+ R7*SYT(3,I)
      RQX=  (PX+R3M*(P(3,J)-P(3,J1)))*R5
      RQXST=  (PX+R3P*(P(3,J)+P(3,J1)))*R5ST

      RX=RQX+RQXST
      DO K=1,2
      SY(K,I)=SY(K,I)+RX*DX(K)
      TEMP=R11*DX(K)
      AT(K,J)=AT(K,J)+TEMP
      AT(K,J2)=AT(K,J2)-TEMP
      ENDDO
      SY(3,I)=SY(3,I)+RQX*R3M-RQXST*R3P
      AT(3,J)= AT(3,J) + R8*SYTR+R9*SYT(3,I)
      AT(3,J2)=AT(3,J2)- SIG*(R7*SYTR+R10*SYT(3,I))
      ENDDO

      DO K=1,3
      A(K,I)=A(K,I)+ RNX*SY(K,I)
      ENDDO
      ENDDO
      ENDDO
C-------------------------------------------------------------
cc      INSERT #3: front and back walls to drop contribution
c     assuming
      C4= 3.D0/((RMU+1.D0)*PI)

      DO I=1,NVS
      DO K=1,3
      SYT(k,i)= PT(K,I)*PHI(I)
      ENDDO

      J1=IM(I,NP+1)
      R3M= -DH-V(3,I)
      R3P= V(3,I)-DH
      RNX= -R3M*C4
      RNXST= -R3P*C4

      DO J=INB(NP+1), INF(NP+1)
      PP=0.D0
      PX=0.D0
      SYTR=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      PP=PP+DX(K)**2
      PX=PX+ (P(K,J)-P(K,J1))*DX(K)
      SYTR=SYTR+ SYT(K,I)*DX(K)
      ENDDO
      RR=PP+R3M**2
      R5= PHI(J)/(RR**2*DSQRT(RR))
      RR=PP+R3P**2
      R5ST= PHI(J)/(RR**2*DSQRT(RR))

      RQX=  (PX+R3M*(P(3,J)-P(3,J1)))*RNX*R5
      RQXST=  (PX+R3P*(P(3,J)-P(3,J1)))*RNXST*R5ST
      RX=RQX+RQXST

      R11= (RNX*R5+RNXST*R5ST)*SYTR+ SYT(3,I)*
     >(RNX*R5*R3M -RNXST*R5ST*R3P)

      DO K=1,2
      A(K,I)=A(K,I)+RX*DX(K)
      TEMP=R11*DX(K)
      AT(K,J)=AT(K,J)+TEMP
      AT(K,J1)=AT(K,J1)-TEMP
      ENDDO

      A(3,I)=A(3,I)+RQX*R3M-RQXST*R3P
      TEMP= SYTR*(R3M*RNX*R5+R3P*RNXST*R5ST)+
     >SYT(3,I)*(R3M**2*RNX*R5-R3P**2*RNXST*R5ST)
      AT(3,J)=AT(3,J)+TEMP
      AT(3,J1)=AT(3,J1)-TEMP
      ENDDO
      ENDDO
C---------------------------------------------------------
      do ip=0,np+1
      DO I= inbh(ip),infh(ip)
      DO K=1,3
      AT(K,I)=AT(K,I)/PHI(I)
      ENDDO
      ENDDO
      ENDDO

C     INTERACTION BETWEEN BACK AND FRONT PANELS:
      IP=NP+1
      D=DEPTH
      DO I=INB(IP), INF(IP)
      DO J=I+1, INF(IP)
      RR=D2
      RQ=0.D0
      RPTM=0.D0
      RPTP=0.D0
      DO K=1,2
      DX(K)=V(K,J)-V(K,I)
      RR=RR+ DX(K)**2
      RQ=RQ+ (P(K,J)-P(K,I))*DX(K)
      RPTM=RPTM+(PT(K,J)-PT(K,I))*DX(K)
      RPTP=RPTP+(PT(K,J)+PT(K,I))*DX(K)
      ENDDO
      R5= C22/(RR**2*DSQRT(RR))
      TMP= DEPTH*(P(3,J)-P(3,I))
      RX= (RQ-TMP)*R5*PHI(J)
      RY= (RQ+TMP)*R5*PHI(I)
      TMP=(D*(PT(3,I)+PT(3,J))-RPTM)*R5
      TMPX= -TMP*PHI(J)
      TMPY=  TMP*PHI(I)

      DO K=1,2
      A(K,I)=A(K,I)+RX*DX(K)
      A(K,J)=A(K,J)-RY*DX(K)
      AT(K,I)=AT(K,I)+TMPX*DX(K)
      AT(K,J)=AT(K,J)+TMPY*DX(K)
      ENDDO
      A(3,I)=A(3,I) +RX*D
      A(3,J)=A(3,J) +RY*D
      TMP= (D*(PT(3,I)-PT(3,J))+RPTP)*D*R5
      AT(3,I)=AT(3,I)+TMP*PHI(J)
      AT(3,J)=AT(3,J)-TMP*PHI(I)
      ENDDO
      ENDDO
c----------------------------------------------------------------
cc      insert #4: drop self-interaction
c     assuming
      C5= (RMU-1.D0)/(RMU+1.D0)
      C6= C5*3.D0/(2.D0*PI)
      DO I=1,NVS
C     ADDED-BACK TERM (DUE TO SINGYLARITY SUBTRACTION)
      DO K=1,3
      A(K,I)=A(K,I)+C5*P(K,I)
      AT(K,I)=AT(K,I)+C5*PT(K,I)
      ENDDO
      DO J=I+1,NVS
      RR=0.D0
      RQ=0.D0
      RPTX=0.D0
      RPTY=0.D0
      RNX=0.D0
      RNY=0.D0
      DO K=1,3
      DX(K)=V(K,J)-V(K,I)
      RR=RR+DX(K)**2
      RQ=RQ+DX(K)*(P(K,J)-P(K,I))
      RPTX=RPTX+DX(K)*PT(K,J)
      RPTY=RPTY+DX(K)*PT(K,I)
      RNX=RNX+ DX(K)*XNA(K,J)
      RNY=RNY+ DX(K)*XNA(K,I)
      ENDDO
      R5=C6/(RR**2*DSQRT(RR))
      RX=R5*RQ*RNX*PHI(J)
      RY=R5*RQ*RNY*PHI(I)
      TMP= (RNX-RNY)*(RPTX-RPTY) -(RNX+RNY)*(RPTX+RPTY)
      TMP= TMP*R5/2.D0
      TMPX=TMP*PHI(J)
      TMPY=TMP*PHI(I)
      DO K=1,3
      A(K,I)=A(K,I) +RX*DX(K)
      A(K,J)=A(K,J) +RY*DX(K)
      AT(K,I)=AT(K,I)+ TMPX*DX(K)
      AT(K,J)=AT(K,J)- TMPY*DX(K)
      ENDDO
      ENDDO
      ENDDO
C------------------------------------------------------------------------------
C!    MAKE SURE P AND PT ARE CONTINUED TO 2ND HALF
      CALL DEFL(NVS+1,NVT,V,PHI,XNA,XCW,DM,SW,P,0,-1,URBP,OMRBP,A)

      CALL DEFL(NVS+1,NVT,V,PHI,XNA,XCW,DM,SW,PT,0,-1,URBP,OMRBP,AT)

      DO IP=0,NP+1
      DO I=INBH(IP),INFH(IP)
      DO K=1,3
      A(K,I)=P(K,I) - A(K,I)
      AT(K,I)=PT(K,I) -AT(K,I)
      ENDDO
      ENDDO
      ENDDO

      RR=0.D0
      AL=0.D0
      DO IP=0,NP+1
      DO I= INBH(IP),INFH(IP)
      DO K=1,3
      RR=RR+R(K,I)*RT(K,I)*PHI(I)
      AL=AL+PT(K,I)*A(K,I)*PHI(I)
      ENDDO
      ENDDO
      ENDDO

      TE=0.D0
      DO IP= 0,NP+1
      DO I= INBH(IP),INFH(IP)
      DO K=1,3
      TE=TE+P(K,I)*AT(K,I)*PHI(I)
      ENDDO
      ENDDO
      ENDDO

      print *, ' PHI=',PHI(I)
      print *, ' PT(K,I)=',PT(K,I)
      print *, ' A(K,I)=',A(K,I)

      PRINT 323,  AL, TE
323   format (' al=', d25.17,' te=', d25.17)
      PRINT *, ' RR=', RR

      AL= RR/AL
      RRNEW=0.D0
      DO IP=0,NP+1
      DO I=INBH(IP),INFH(IP)
      DO K=1,3
      R(K,I)= R(K,I) - AL*A(K,I)
      RT(K,I)= RT(K,I) - AL*AT(K,I)
      RRNEW=RRNEW+R(K,I)*RT(K,I)*PHI(I)
      ENDDO
      ENDDO
      ENDDO

      BET=RRNEW/RR
      DO IP=0,NP+1
      DO I=INBH(IP),INFH(IP)
      DO K=1,3
      QW(K,I)= QW(K,I) - AL*P(K,I)
      P(K,I)= R(K,I) + BET*P(K,I)
      PT(K,I)= RT(K,I) + BET*PT(K,I)
      ENDDO
      IF (IP.NE.0) THEN
        I1=INV(I)
        QW(1,I1)=  QW(1,I)
        QW(2,I1)=  QW(2,I)
        QW(3,I1)= -QW(3,I)

	P(1,I1)=  P(1,I)
        P(2,I1)=  P(2,I)
        P(3,I1)= -P(3,I)

	R(1,I1)=  R(1,I)
        R(2,I1)=  R(2,I)
        R(3,I1)= -R(3,I)

	PT(1,I1)=  PT(1,I)
        PT(2,I1)=  PT(2,I)
        PT(3,I1)= -PT(3,I)

	RT(1,I1)=  RT(1,I)
        RT(2,I1)=  RT(2,I)
        RT(3,I1)= -RT(3,I)

      END IF
      ENDDO
      ENDDO

c     DROP velocity and RMS:
      do k=1,3
      ud(k)=0.d0
      enddo
      QMSQ=0.d0
      DO I=1,NVS
      QM(I)=0.D0
      DO K=1,3
      QM(I)=QM(I)+QW(K,I)*XNA(K,I)
      ENDDO
      do k=1,3
      UD(k)=UD(k)+ (V(K,I)-XP(K))*QM(I)*PHI(i)/VOL
      ENDDO
      QMSQ=QMSQ+ qm(i)**2*PHI(i)
      ENDDO
      QMSQ= dsqrt(QMSQ/FSD)
      print *, ' IT=', IT, ' QMSQ=', QMSQ
      print *, 'UD=', UD

      IF (IT.NE.90) THEN
        GOTO 777
      ELSE
	GOTO 775
      END IF
778   continue
      PRINT *, 'XP=', XP(1), XP(2), XP(3)
c      STOP
      RETURN
      END

      SUBROUTINE DEFL(INB,INF,V,PHI,XNA,XC,DM,FS,QW,MRBM,MNORM,UAV,
     >OM,DL)
      IMPLICIT REAL *8(A-H,O-Z)
      DIMENSION QW(3,*),DL(3,*),XC(3),DM(3,3),PHI(*),XNA(3,*),
     >V(3,*),SM(3),DX(3),OM(3),UAV(3)

      DO K=1,3
      SM(K)=0.D0
      UAV(K)=0.D0
      ENDDO
      QN=0.D0
      DO 502 I=INB,INF
      P=0.D0
      DO K=1,3
      UAV(K)=UAV(K)+QW(K,I)*PHI(I)
      DX(K)=V(K,I)-XC(K)
      P=P+QW(K,I)*XNA(K,I)
      ENDDO
      QN=QN+P*PHI(I)
      SM(1)=SM(1)+(DX(2)*QW(3,I)-DX(3)*QW(2,I))*PHI(I)
      SM(2)=SM(2)+(DX(3)*QW(1,I)-DX(1)*QW(3,I))*PHI(I)
      SM(3)=SM(3)+(DX(1)*QW(2,I)-DX(2)*QW(1,I))*PHI(I)
502   CONTINUE
      QN=MNORM*QN/FS
      DO K=1,3
      UAV(K)=UAV(K)/FS
      OM(K)=0.D0
      DO L=1,3
      OM(K)=OM(K)+DM(K,L)*SM(L)
      ENDDO
      ENDDO

      DO 503 I=INB,INF
      DO K=1,3
      DX(K)=(V(K,I)-XC(K))
      ENDDO
      DL(1,I)= DL(1,I) -(UAV(1)+OM(2)*DX(3)-OM(3)*DX(2))*MRBM +
     >QN*XNA(1,I)
      DL(2,I)= DL(2,I) -(UAV(2)+OM(3)*DX(1)-OM(1)*DX(3))*MRBM +
     >QN*XNA(2,I)
      DL(3,I)= DL(3,I) -(UAV(3)+OM(1)*DX(2)-OM(2)*DX(1))*MRBM +
     >QN*XNA(3,I)
503   CONTINUE
      RETURN
      END

      SUBROUTINE DEFLMAT (INB,INF,XC,V,PHI,DM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XC(3),V(3,*),PHI(*),T(3,6),DM(3,3)
      DO I=1,3
      DO J=1,6
       T(I,J)=0.D0
      ENDDO
       T(I,I+3)=1.D0
      ENDDO
      do I=INB,INF
      R2=(V(1,I)-XC(1))**2 + (V(2,I)-XC(2))**2 + (V(3,I)-XC(3))**2
      DO K=1,3
      DO L=1,3
      DEL=0.D0
      IF(K.EQ.L) DEL=1.D0
      T(K,L)=T(K,L)+(DEL*R2-(V(K,I)-XC(K))*(V(L,I)-XC(L)))*PHI(I)
      ENDDO
      ENDDO
      ENDDO
      CALL GM36(3,6,T)
      DO K=1,3
      DO L=1,3
      DM(K,L)=T(K,L+3)
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE GM36(N,M,T)
      real*8 T(3,6),R,U
      INTEGER S
      DO 1 I=1,N
      S=I
      R=T(I,I)
      DO 2 K=I+1,N
      IF(DABS(T(K,I)).GT.DABS(R)) THEN
         S=K
	 R=T(K,I)
      END IF
2     CONTINUE
      T(S,I)=T(I,I)
      DO 3 J=I+1,M
      U=T(S,J)/R
      T(S,J)=T(I,J)
      T(I,J)=U
      DO 4 K=I+1,N
4     T(K,J)=T(K,J)-T(K,I)*U
3     CONTINUE
1     CONTINUE
      DO 5 I=N-1,1,-1
      DO 5 J=N+1,M
      DO 5 K=I+1,N
      R=T(I,K)
      T(I,J)=T(I,J)-R*T(K,J)
5     CONTINUE
      RETURN
      END

      SUBROUTINE GMRES(N,M,T)
      PARAMETER (NRESMAX=20)
      real*8 T(NRESMAX,NRESMAX+1),R,U
      INTEGER S
      DO 1 I=1,N
      S=I
      R=T(I,I)
      DO 2 K=I+1,N
      IF(DABS(T(K,I)).GT.DABS(R)) THEN
         S=K
	 R=T(K,I)
      END IF
2     CONTINUE
      T(S,I)=T(I,I)
      DO 3 J=I+1,M
      U=T(S,J)/R
      T(S,J)=T(I,J)
      T(I,J)=U
      DO 4 K=I+1,N
4     T(K,J)=T(K,J)-T(K,I)*U
3     CONTINUE
1     CONTINUE
      DO 5 I=N-1,1,-1
      DO 5 J=N+1,M
      DO 5 K=I+1,N
      R=T(I,K)
      T(I,J)=T(I,J)-R*T(K,J)
5     CONTINUE
      RETURN
      END

      SUBROUTINE SPHAREA(Y,X0,DSSPH)
!     CALCULATES THE AREA OF A SPHERICAL TRIANGLE PROJECTED ON UNIT SPHERE
!     Y(j, j=1,2,3)IS THE SPHERE CENTER, AND TRIANGLE VERTICES
!     BEFORE PROJECTION ARE X0(K,J
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION XNA(3,3),SIDE1(3),SIDE2(3),Y(3),X0(3,*)
!     FOR THIS APPLICATION, X0(3,*) INSTEAD OF X0(3,3) IS IMPORTANT
      DO K=1,3
      P=0.D0
      DO J=1,3
      XNA(J,K)=X0(J,K)-Y(J)
      P=P+XNA(J,K)**2
      ENDDO
      P=1.D0/DSQRT(P)
      DO J=1,3
      XNA(J,K)=P*XNA(J,K)
      ENDDO
      ENDDO

      PI=DACOS(-1.D0)
      SUM=0.D0

      DO 446 K=1,3
      K1=K+1
      IF (K1.GT.3) K1=K1-3
      K2=K+2
      IF (K2.GT.3) K2=K2-3
      P1=0.D0
      P2=0.D0
      DO J=1,3
      SIDE1(J)=XNA(J,K1)-XNA(J,K)
      SIDE2(J)=XNA(J,K2)-XNA(J,K)
      P1=P1+ SIDE1(J)*XNA(J,K)
      P2=P2+ SIDE2(J)*XNA(J,K)
      ENDDO
      A1=0.D0
      A2=0.D0
      P=0.D0
      DO J=1,3
      SIDE1(J)=SIDE1(J)-P1*XNA(J,K)
      A1=A1+SIDE1(J)**2
      SIDE2(J)=SIDE2(J)-P2*XNA(J,K)
      A2=A2+SIDE2(J)**2
      P=P+SIDE1(J)*SIDE2(J)
      ENDDO
      P=P/DSQRT(A1*A2)
      p=p*(1.d0-1.d-15)
      SUM=SUM+DACOS(P)
446   CONTINUE
      DSSPH= SUM-PI
      RETURN
      END

      SUBROUTINE RECSTREX (Y,X0,E1,L1,E2,L2,XN,STR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION Y(3),X0(3),E1(3),E2(3),XN(3),STR(3,3),E(3,4),X(3,4),
     >L(4)
      REAL *8 L1,L2,L

      DO I=1,3
      DO J=1,3
      STR(I,J)=0.D0
      ENDDO
      ENDDO

      DO J=1,3
      X(J,1)=X0(J)
      E(J,1)=E1(J)

      X(J,2)=X0(J)+L1*E1(J)
      E(J,2)=E2(J)

      X(J,3)=X(J,2)+L2*E2(J)
      E(J,3)= -E1(J)

      X(J,4)=X0(J)+L2*E2(J)
      E(J,4)= -E2(J)
      ENDDO

      L(1)=L1
      L(2)=L2
      L(3)=L1
      L(4)=L2

      DO K=1,4
      CALL SEGCONT (Y,X(1,K),E(1,K),L(K),STR)
      ENDDO

      P=XN(1)*(E1(2)*E2(3)-E1(3)*E2(2))+XN(2)*
     >(E1(3)*E2(1)-E1(1)*E2(3))+ XN(3)*(E1(1)*E2(2)-E1(2)*E2(1))
      IF (P.LT.0.D0) THEN
        DO I=1,3
        DO J=1,3
	STR(I,J)=-STR(I,J)
	ENDDO
	ENDDO
      END IF

      RNX=0.D0
      DO J=1,3
      RNX=RNX+XN(J)*(X0(J)-Y(J))
      ENDDO

      CALL SPHAREA(Y,X,S1)
      DO J=1,3
      X(J,2)=X(J,3)
      X(J,3)=X(J,4)
      ENDDO
      CALL SPHAREA(Y,X,S2)
      S=(S1+S2)/3.D0
      IF(RNX.LT.0.D0) S= -S
      DO I=1,3
      STR(I,I)=STR(I,I)+S
      ENDDO
c     to make it perfectly symmetrical:
      do i=1,2
      do j=i+1,3
      str(i,j)= str(j,i)
      enddo
      enddo
      RETURN
      END

C     ANALYTICAL INTEGRATION OVER BACK OR FRONT PANEL TO OBTAIN STR-TENSOR
C     Y IS OBSERVATION  POINT (ASSUMED TO BE BETWEEN FRONT AND BACK PANELS
C     IN X3-DIRECTION)
C     id = -1 for back panel ibegration, id=1 for front panel integration
c     NP IS NUMBER OF *SIDE* PANELS (ON THE WALL FOR BINEW3D, OR ON THE MF
C     FOR FUTURE VEL3D WITH FINITE DEPTH
C     depth IS FULL CHANNEL DEPTH
C     X0 IS ARRAY OF CORNER POINTS FOR CHANNEL 2d PROFILE (CHANNEL DOMAIN
C     REMAINS ON THE LEFT, AS BEFORE), E1 ARE UNIT DIRECTORS FROM
C     X0(., JP) TO THE NEXT CORNER POINT, L1(JP) IS EDGE LENGTH TO THAT
C     NEXT CORNER POINT; THIRD COORDINATES OF X0 AND E1 NOT USED HERE.

      SUBROUTINE BKFRSTR(Y,ID,NP,DEPTH,X0,E1,L1,STR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION STR(3,3),Y(3),X0(3,*),E1(3,*),L1(*),X(3),DIR(3),
     >XN(2),XX(3,3)
      REAL *8 L1
      DH=DEPTH/2
      DO I=1,3
      DO J=1,3
      STR(I,J)=0.D0
      ENDDO
      ENDDO

      S=0.D0
      DO 1 JP=1,NP

      DO K=1,2
      X(K)=X0(K,JP)
      DIR(K)=E1(K,JP)
      ENDDO
      X(3)=ID*DH
      DIR(3)=0.D0
      CALL SEGCONT (Y,X,DIR,L1(JP),STR)
      XN(1)= DIR(2)
      XN(2)= -DIR(1)
      RNX=0.D0
      DO K=1,2
      RNX=RNX+ (X(K)-Y(K))*XN(K)
      ENDDO
      DO K=1,3
      XX(K,1)= X(K)
      XX(K,2)= X(K)+DIR(K)*L1(JP)
      XX(K,3)=Y(K)
      ENDDO
      XX(3,3)= XX(3,3)+ID
      CALL SPHAREA(Y,XX,OM)
      IF (RNX.LT.0.D0) OM= -OM
      S=S+OM
1     CONTINUE

      DO I=1,3
      DO J=1,3
      STR(I,J)=STR(I,J)*ID
      ENDDO
      STR(I,I)=STR(I,I)+ S/3.D0
      ENDDO
c     to make it perfectly symmetrical:
      do i=1,2
      do j=i+1,3
      str(i,j)= str(j,i)
      enddo
      enddo
      RETURN
      END

      SUBROUTINE SEGCONT (Y,X0,DIR,L,STR)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 L
      DIMENSION Y(3),X0(3),DIR(3),STR(3,3),DX(3),E(3,3)
      RR=0.D0
      P=0.D0
      DO K=1,3
      DX(K)=Y(K)-X0(K)
      RR=RR+DX(K)**2
      P=P+DX(K)*DIR(K)
      ENDDO
      H=DSQRT(RR-P**2)
      T1=-P/H
      T2= (L-P)/H
      RI1=T2/DSQRT(1.D0+T2**2) - T1/DSQRT(1.D0+T1**2)
      RI2=1.D0/DSQRT(1.D0+T1**2) -1.D0/DSQRT(1.D0+T2**2)
      RI1=RI1/3.D0
      RI2=RI2/3.D0
      E(1,3)=DIR(2)*DX(3)-DIR(3)*DX(2)
      E(2,3)=DIR(3)*DX(1)-DIR(1)*DX(3)
      E(3,3)=DIR(1)*DX(2)-DIR(2)*DX(1)
      TEMP=0.D0
      DO K=1,3
      TEMP=TEMP+E(K,3)**2
      E(K,2)=DX(K)-P*DIR(K)
      E(K,1)=DIR(K)
      ENDDO
      TEMP=1.D0/DSQRT(TEMP)
      DO K=1,3
      E(K,2)=E(K,2)/H
      E(K,3)=E(K,3)*TEMP
      ENDDO
      DO K=1,3
      DO J=1,3
      STR(K,J)=STR(K,J)+(RI2*E(J,1)-RI1*E(J,2))*E(K,3)
      ENDDO
      ENDDO
      RETURN
      END

	SUBROUTINE PANEL2D(I,NVCMAX,N1,L1,X0,E1,INB,INF,V,XNW,DSW)
c     change real to match original
	IMPLICIT REAL*8 (A-H,O-Z)
	REAL *8 L1
	DIMENSION X0(3),E1(3),V(2,*),XNW(3)
	H1=L1/N1
	INB=I + 1
	INF = I + N1
	DSW = H1
	XNW(1)=E1(2)
	XNW(2)=-E1(1)
	xnw(3)=0.d0
	DO I1 = 1,N1
	I = I + 1
	IF (I.GT.NVCMAX) THEN
	PRINT *,'panel2d: Out of Bounds '
	STOP
	ENDIF
	DO K=1,2
	V(K,I)=X0(K)+(I1-0.5D0)*H1*E1(K)
	ENDDO
	ENDDO
	RETURN
	END

      SUBROUTINE PNL3DCH(I,INB1,INF1,V1,DS1,
     >L2,N2,INB,INF,INBH,INFH,INV,V,PHI)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 L2,L2H
      DIMENSION V(3,*),INV(*),PHI(*),V1(2,*)
!     SYMMETRY assumed ABOUT Z=0 PLANE, L2 is FULL channel depth
      N1= INF1-INB1+1
      L2H=0.5D0*L2
      H= L2H/N2
      DS2=H
      I0=I
      I=0
      DO 44 I1 = INB1,INF1
        DO J=1,N2
        DO K=1,2
        V(K,I0+I+J)=V1(K,I1)
        ENDDO
        V(3,I0+I+J)= -L2H+ (J-0.5D0)*H
        PHI(I0+I+J)= DS1*H
        ENDDO
      I=I+N2
44    CONTINUE
      INB=I0+1
      INF=I0+2*I
      INBH=INB
      INFH=I0+I

      DO J=I0+1,I0+I
      INV(J)=J+I
      INV(J+I)=J
      DO K=1,2
      V(K,INV(J))=V(K,J)
      ENDDO
      V(3,INV(J)) = -V(3,J)
      PHI(INV(J))=PHI(J)
      ENDDO
      I=INF
      RETURN
      END

      SUBROUTINE BOUSFLUX(H,L,TOL,Q)
      IMPLICIT REAL*8 (A-H,L,N-Z)
	  PI=DACOS(-1.D0)

	  N=0.D0
	  SUM1=0.D0
	  G=1.D0

	  A=G*h**3.d0*L/12.d0
	  B=16.d0*h**4.d0/PI**5.d0

1     N=N+1

      BETA=(2.D0*N-1.D0)*PI/h
      C=(1.D0-2.D0*EXP(-BETA*L)+EXP(-2.D0*BETA*L))
      D=(1.D0-EXP(-2.D0*BETA*L))
	  E=(2.d0*N-1.d0)**5.d0
	  SUM1=SUM1+C/(D*E)

	  Q=A-B*SUM1
	  RN = 1.d0/(N**4.d0)

      IF (RN.GT.TOL)GOTO 1

	  PRINT *, 'Q=', Q

	  RETURN
	  END

      subroutine NEWRECST (Y,X0,E1,L1,E2,L2,XN,STR)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION Y(3),X0(3),E1(3),E2(3),XN(3),STR(3,3),
     >X(3,4),E(3,4),RA(3),RB(3),HN(3),HNXE(3)
      REAL *8 L1,L2, LA,LB
      DO I=1,3
      DO J=1,3
      STR(I,J)=0.D0
      ENDDO
      ENDDO

      DO J=1,3
      X(J,1)=X0(J)
      E(J,1)=E1(J)

      X(J,2)=X0(J)+L1*E1(J)
      E(J,2)=E2(J)

      X(J,3)=X(J,2)+L2*E2(J)
      E(J,3)= -E1(J)

      X(J,4)=X0(J)+L2*E2(J)
      E(J,4)= -E2(J)
      ENDDO

      DO 1 K=1,4
      K2=K+1
      IF (K2.GT.4) K2=K2-4
      LA=0.D0
      LB=0.D0
      XA=0.D0
      XB=0.D0
      DO J=1,3
      RA(J)= X(J,K)-Y(J)
      LA=LA+RA(J)**2
      RB(J)= X(J,K2)-Y(J)
      LB=LB+ RB(J)**2
      XA=XA+RA(J)*E(J,K)
      XB=XB+RB(J)*E(J,K)
      ENDDO

      LA=DSQRT(LA)
      LB=DSQRT(LB)

      IF (DABS(XA/LA).GT. DABS(XB/LB)) THEN
        H=0.D0
        DO J=1,3
        HN(J)=RA(J)-XA*E(J,K)
        H=H+HN(J)**2
        ENDDO
      ELSE
        H=0.D0
        DO J=1,3
        HN(J)=RB(J)-XB*E(J,K)
        H=H+HN(J)**2
        ENDDO
      END IF
      H=DSQRT(H)
      DO J=1,3
      HN(J)=HN(J)/H
      ENDDO

      HNXE(1)=HN(2)*E(3,K)-HN(3)*E(2,K)
      HNXE(2)=HN(3)*E(1,K)-HN(1)*E(3,K)
      HNXE(3)=HN(1)*E(2,K)-HN(2)*E(1,K)
      FE= (H/LA-H/LB)/3.D0
      FH= (XB/LB-XA/LA)/3.D0
      DO I=1,3
      DO J=1,I
      STR(I,J)=STR(I,J)+(FE*E(I,K)+FH*HN(I))*HNXE(J)
      ENDDO
      ENDDO
1     CONTINUE

      P=XN(1)*(E1(2)*E2(3)-E1(3)*E2(2))+XN(2)*
     >(E1(3)*E2(1)-E1(1)*E2(3))+ XN(3)*(E1(1)*E2(2)-E1(2)*E2(1))
      IF (P.LT.0.D0) THEN
        DO I=1,3
        DO J=1,I
	STR(I,J)= -STR(I,J)
	ENDDO
	ENDDO
      END IF

      RNX=0.D0
      DO J=1,3
      RNX=RNX+XN(J)*(X0(J)-Y(J))
      ENDDO

      CALL SPHAREA(Y,X,S1)
      DO J=1,3
      X(J,2)=X(J,3)
      X(J,3)=X(J,4)
      ENDDO
      CALL SPHAREA(Y,X,S2)
      S=(S1+S2)/3.D0
      IF(RNX.LT.0.D0) S= -S
      DO I=1,3
      STR(I,I)=STR(I,I)+S
      ENDDO
c     to make it perfectly symmetrical:
      do i=1,2
      do j=i+1,3
      str(i,j)= str(j,i)
      enddo
      enddo
      RETURN
      END

      SUBROUTINE CONTOR (XP,A,X0,E1,L1,NPC,X0C,SOL)
C     A FIXED VERSION OF CONTOR, TO ALLOW FOR DOUBLE-CONNECTED INTERSECTION
C     AREA OF THE SQUARE AND CHANNEL;IF THIS HAPPENS, CHOOSES THE PIECE
C     CONTAINING PARTICLE CENTER.
C       DID NOT TRY YET TO FIX SOL; NOT USED SO FAR, ANYWAY.
!     ************************************************************************
!     * XP(1),XP(2) PARTICLE CENTER LOCATION                                 *
!     * A IS HALFWIDTH OF THE SQUARE BOX (NOT TRIMMED YET) AROUND XP(1),XP(2)*
!     * X0, E1, L1 ARE ARRAYS SPECIFYING THE CHANNEL (EXACTLY LIKE IN BINEW  *
!     * AND TRAJNEW.                                                         *
!     *                                                                      *
!     * OUTPUTS:                                                             *
!     *      NPC IS THE NUMBER OF CORNER POINTS OF THE CELL (BUT THE         *
!     *      FIRST AND LAST POINT COUNTED TWICE)                             *
!     *                                                                      *
!     *     X0C IS ARRAYS OF CORNER POINTS FOR THE CELL                      *
!     *     (LISTED COUNTERCLOCKWISE)                                        *
!     *     SOL(IP) IS .TRUE. IF THE CELL EDGE ORIGINATING                   *
!     *     FROM X0(1,IP),X0(2,IP)                                           *
!     *     BELONG TO THE CHANNEL BOUNDARY, AND SOL(IP)=.FALSE. OTHERWISE.   *
!     ************************************************************************

      IMPLICIT REAL *8 (A-H,O-Z)
      INCLUDE 'CH3D.h'
      PARAMETER (NPCMAX=50, NVT=50000)
      REAL *8 L1
      LOGICAL SOL
!     XP(*) ESSENTIAL;
C     1ST DIMENSION OF X0, E1 AND X0C IS 3, UNLIKE IN INFINITE DEPTH PACKAGE
      DIMENSION XP(3),L1(*), X0(3,*), E1(3,*),V(2,NVT),X0SQ(2,4),
     >INB(NPCH),INF(NPCH),X0INT(2),X0C(3,*),SOL(*), X(2,4),Y(2)
C     INDF CALLED HERE IS FROM PACKAGE INIPOINT.F, NOT FROM INFINITE DEPTH PAC.
C     (MUCH FASTER, WITH X0(3,*) ASSUMED)
      INTEGER COLOR(4)
      IF (INDF(XP,NPCH,X0).EQ.0) THEN
	PRINT *, ' PARTICLE CENTER OUTSIDE CHANNEL'
	STOP
      END IF
      EPS=1.D-14*A**2
      DO IP=1,NPCMAX
      SOL(IP)=.FALSE.
      ENDDO

      S=0.D0
      DO IP=1,NPCH
      S=S+L1(IP)
      ENDDO

      I=0
      DO IP=1,NPCH
      N1=L1(IP)*NVT/S
      H1=L1(IP)/N1
      INB(IP)=I + 1
      INF(IP) = I + N1
      DO I1 = 1,N1
      I = I + 1
      DO K=1,2
      V(K,I)=X0(K,IP)+(I1-1)*H1*E1(K,IP)
      ENDDO
      ENDDO
      ENDDO
!     ************************************************************************
C      OPEN (2,FILE='V.DAT')
C      DO IP=1,NPCH
C      DO I=INB(IP),INF(IP)
C      WRITE (2,*) V(1,I),V(2,I)
C      ENDDO
C      ENDDO
C      CLOSE(2)
!      V OK
!     ************************************************************************

      X0SQ(1,1)=XP(1)-A
      X0SQ(2,1)=XP(2)-A

      X0SQ(1,2)=XP(1)+A
      X0SQ(2,2)=XP(2)-A

      X0SQ(1,3)=XP(1)+A
      X0SQ(2,3)=XP(2)+A

      X0SQ(1,4)=XP(1)-A
      X0SQ(2,4)=XP(2)+A

      OPEN (2,FILE='SQ.DAT')
      DO IP=1,4
      WRITE (2,*) X0SQ(1,IP),X0SQ(2,IP)
      ENDDO
      WRITE (2,*) X0SQ(1,1),X0SQ(2,1)
      CLOSE(2)

      IPSTART=1

51    IP=IPSTART

      NPC=0
      IPSQ0=0
      DO JP=1,4
      COLOR(JP)=0
      ENDDO
      IND=1

2     DO 1 I=INB(IP),INF(IP)
      I1=I+1
      IF (I1.GT.INF(NPCH))I1=1
      SIGOLD= -1.D20
      SIGNEW= -1.D20
      DO K=1,2
      SIGOLD=MAX(SIGOLD, DABS(V(K,I)-XP(K)) -A)
      SIGNEW=MAX(SIGNEW, DABS(V(K,I1)-XP(K)) -A)
      ENDDO
      IF (SIGNEW*SIGOLD.GT.0.D0) GOTO 1
      IF( (DABS(V(1,I)-XP(1))-A)*(DABS(V(1,I1)-XP(1))-A).LT.0.D0 ) THEN
	  K=1
      ELSE
        K=2
      END IF

      IF ( (V(K,I)-XP(K)-A)*(V(K,I1)-XP(K)-A).LT.0.D0 ) THEN
	  X0INT(K)=XP(K)+A
	  IF (K.EQ.1) THEN
	    IPSQ=2
	  ELSE
	    IPSQ=3
	  END IF
      ELSE
        X0INT(K)=XP(K)-A
	  IF (K.EQ.1) THEN
	    IPSQ=4
	  ELSE
	    IPSQ=1
	  END IF
      END IF

      M=3-K
      X0INT(M)=V(M,I)+(V(M,I1)-V(M,I))*(X0INT(K)-V(K,I))/
     >(V(K,I1)-V(K,I))
C      PRINT *, ' IP=', IP, ' K=', K, ' IPSQ=', IPSQ, ' SIGNEW=',SIGNEW
C      PRINT *, ' SIGOLD=', SIGOLD
C      PRINT *, ' IND=', IND, ' IPSQ0=', IPSQ0
C      PRINT *, ' NPC=', NPC
      IF (SIGNEW.GT.0.D0 .AND. IND.EQ.0) GOTO 1
      IF (SIGNEW.GT.0.D0 .AND. IND.EQ.1) THEN
        IPSQ0=IPSQ
	IP0=IP
	NPC=NPC+1
        IF (NPC.GT.NPCMAX) THEN
	    PRINT *, ' NPCMAX TOO SMALL'
	    STOP
        END IF
	DIFF=0.D0
        DO K=1,2
        X0C(K,NPC)=X0INT(K)
        DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
	ENDDO
        IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
      ELSE

        IF(IPSQ0.NE.0) THEN
	    DO K=1,2
	    Y(K)=0.5D0*(X0C(K,NPC)+X0INT(K))
	    ENDDO
	    IND= INDF(Y,NPCH,X0)
	    IF(IND.EQ.0  .AND.  MOD(IP-IP0-1,NPCH).EQ.0) GOTO 1
	    IF(IPSQ.EQ.IPSQ0  .AND. IND.EQ.0) IPSQ=IPSQ+4
	    IF(IPSQ.LT.IPSQ0) IPSQ=IPSQ+4
            IND=1
	    DO 44 J=IPSQ0,IPSQ
	    J1=MOD(J-1,4)+1

          IF (J.EQ.IPSQ0) GOTO 45
	  NPC=NPC+1
          IF (NPC.GT.NPCMAX) THEN
	      PRINT *, ' NPCMAX TOO SMALL'
	      STOP
          END IF
	    DIFF=0.D0
	    DO K=1,2
	    X0C(K,NPC)=X0SQ(K,J1)
          DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
	    ENDDO
            IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
45          IF (COLOR(J1).EQ.1) THEN
             DO K=1,2
             Y(K)=0.5D0*(X0C(K,NPC)+X(K,J1))
             ENDDO
             IF(INDF(Y,NPCH,X0).EQ.1)  THEN
              NPC=NPC+1
              IF (NPC.GT.NPCMAX) THEN
                PRINT *, ' NPCMAX TOO SMALL'
                STOP
              END IF

              DO K=1,2
              X0C(K,NPC)=X(K,J1)
              ENDDO
              RETURN
             END IF
            END IF
44          CONTINUE

	END IF

	NPC=NPC+1
        IF (NPC.GT.NPCMAX) THEN
          PRINT *, ' NPCMAX TOO SMALL'
          STOP
        END IF
	IF (IPSQ.GT.4) IPSQ=IPSQ-4
	DIFF=0.D0
        DO K=1,2
        X0C(K,NPC)=X0INT(K)
        DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
	X(K,IPSQ)=X0C(K,NPC)
	ENDDO
	COLOR(IPSQ)=1
        SOL(NPC)=.TRUE.
        IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
      END IF
1     CONTINUE

      IF (SIGNEW.LT.0.D0) THEN
        NPC=NPC+1
        IF (NPC.GT.NPCMAX) THEN
          PRINT *, ' NPCMAX TOO SMALL'
          STOP
        END IF
        DIFF=0.D0
        DO K=1,2
        X0C(K,NPC)=V(K,I1)
        DIFF=DIFF+(X0C(K,NPC)-X0C(K,1))**2
        ENDDO
        SOL(NPC)=.TRUE.
        IF (NPC.NE.1 .AND. DIFF.LT.EPS) GOTO 50
      END IF

      IP=IP+1
      IF (IP.GT.NPCH) IP=1
      IF (IP.NE.IPSTART. OR. NPC.NE.0) GOTO 2

      DO K=1,2
      DO IPC=1,4
      X0C(K,IPC)=X0SQ(K,IPC)
      ENDDO
      X0C(K,5)=X0SQ(K,1)
      ENDDO
      NPC=5
      RETURN

50    IF (INDF(XP,NPC-1,X0C).EQ.0) THEN
        IPSTART=IPSTART+1
        IF (IPSTART.GT.NPCH) IPSTART=1
	GOTO 51
      END IF
      RETURN
      END

C      add BOUSVEL subroutine here in order to do infinite channel calcs
      SUBROUTINE BOUSVEL(Y0,Z0,H0,L0,TOL,U)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL *8 L,MU,L0
        G=1.D0
        MU=1.D0
        PI=DACOS(-1.D0)
        Y=Y0
        Z=Z0
        H=H0
        L=L0

        IF(MIN(((L0-Z0)/H0),Z0/H0).LT.MIN(((H0-Y0)/L0),Y0/L0)) THEN
          Y=Z0
          Z=Y0
          H=L0
          L=H0
        END IF

        N=0
        SUM=0.D0

        A=G*Y*(H-Y)/(2.d0*MU)
        B=4.D0*G*H**2.D0/(MU*PI**3.D0)

2     N=N+1

      BETA=(2.D0*N-1.D0)*PI/h
      C=EXP(BETA*(Z-L))*(1-EXP(-2.d0*BETA*Z))
      D=EXP(-BETA*Z)*(1-EXP(-2.D0*BETA*(L-Z)))
      E=(1.d0-EXP(-2.D0*BETA*L))
      F=(2.d0*N-1.d0)**3.d0
      SUM=SUM+(C+D)/(E*F)*SIN(BETA*Y)

      U=A-B*SUM
      RN = 1.d0/(N**4)

      IF (RN.GT.TOL)GOTO 2

C     PRINT *, 'U=', U
      RETURN
      END
