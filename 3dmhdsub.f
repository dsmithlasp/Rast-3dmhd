C**********************************************************************
	SUBROUTINE SLENGTH(STRNG,JLEN)
C
	CHARACTER*80 STRNG
C
	ILEN = LEN(STRNG)
	DO 10 I=1,ILEN
		IF (STRNG(I:I).EQ.'@') THEN
			JLEN = I-1
			GOTO 20
	        ENDIF
10	CONTINUE
C
20	RETURN
	END
C**********************************************************************
        SUBROUTINE MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,ICOORD)
C----------------------------------------------------------------------
C  If ICOORD=1 then: gets XCODE and returns XPHYS, DXXDX and D2XXDX2
C  If ICOORD=2 then: gets YCODE and returns YPHYS, DYYDY and D2YYDY2
C  If ICOORD=3 then: gets ZCODE and returns ZPHYS, DZZDX and D2ZZDZ2
C----------------------------------------------------------------------
        INCLUDE '3dmhdparam.f'
C
        COMMON/BOUNDS/XMAX,YMAX,ZMAX
C----------------------------------------------------------------------	
C  Grid selection
C----------------------------------------------------------------------	
        SELECT CASE (NGRID)
C----------------------------------------------------------------------	
C  Uniform grid. 
C----------------------------------------------------------------------
	CASE(0)
	  SELECT CASE (ICOORD)
	  CASE(1)
	    SMAX=XMAX
	  CASE(2)
	    SMAX=YMAX
	  CASE(3)
	    SMAX=ZMAX
          CASE DEFAULT
	    WRITE(6,*)'MKGRID: Invalid ICOORD number'
            CALL MPI_FINALIZE(IERR)
            STOP
          END SELECT	  
C
	  SPHYS   = SCODE*SMAX
          DSSDS   = 1.0E00/SMAX
          D2SSDS2 = 0.0E00
C----------------------------------------------------------------------	
C  NOTE:  The nonuniform grid NGRID > 0 has some problems that must be
C  fixed before using.  These are:
C	1.  Horizontal grid stretching must be symmetric about domain 
C           center so that the hard wired periodicity of the domain 
C	    Makes sense.  The domain is periodic over XMAX+DX.
C	2.  Horizontal inegration for mean values (SPLINEX and SPLINEY)
C           must be extended to include ghost point on right (I=2,NX).
C       3.  Normalization of integration for mean values must be changed
C           from XMAX or YMAX to XMAX+DELTAX or YMAX+DELTAY where DELTAX
C           and DELTAY are the physical grid distances to the ghost point.
C       4.  In general the derviative metrics should be calculated 
C           numerically, using the same difference scheme and the code,
C           rather than analytically to take advatage of truncation error 
C	    cancelation.
C----------------------------------------------------------------------	
C  Original grid (by Mark Rast)	
C----------------------------------------------------------------------	
        CASE(1)
          SELECT CASE (ICOORD)
          CASE(1)
            S1   = XX1
            S2   = XX2
            SMAX = XMAX
          CASE(2)
            S1   = YY1
            S2   = YY2
            SMAX = YMAX
          CASE(3)
            S1   = ZZ1
            S2   = ZZ2
            SMAX = ZMAX
          CASE DEFAULT
	    WRITE(6,*)'MKGRID: Invalid ICOORD number'
            CALL MPI_FINALIZE(IERR)
            STOP
          END SELECT	  
C            
	  A1=(S2-S1)/SMAX
	  A3=ATAN(S1)
	  A2=ATAN(S2)-A3
          IF (SCODE.EQ.0.0E00) THEN
            SPHYS=0.0E00
          ELSE IF (SCODE.EQ.1.0E00) THEN
            SPHYS=SMAX
          ELSE
            SPHYS=(TAN(A2*SCODE+A3)-S1)/A1
          ENDIF
	  DSSDS=A1/A2/(1.0E00+(S1+A1*SPHYS)**2)
	  D2SSDS2=-2.0E00*A1**2*(S1+A1*SPHYS)/A2/
     1            (1.0E00+(S1+A1*SPHYS)**2)**2
C----------------------------------------------------------------------	
C  Second grid (by Thierry Emonet 98/05/30)	
C----------------------------------------------------------------------	
	CASE(2)
	  PI=2.0E00*ASIN(1.0E00)
          SELECT CASE (ICOORD)
          CASE(1)
            A=XA
            B=XB
            C=XC
            D=XD
            SMAX = XMAX
          CASE(2)
            A=YA
            B=YB
            C=YC
            D=YD
            SMAX = YMAX
          CASE(3)
            A=ZA
            B=ZB
            C=ZC
            D=ZD
            SMAX = ZMAX
          CASE DEFAULT
	    WRITE(6,*)'MKGRID: Invalid ICOORD number'
            CALL MPI_FINALIZE(IERR)
            STOP
          END SELECT	  
C            
	  IF (C.LE.1.0E-09) THEN
            SPHYS   = SCODE*SMAX
	    DSSDS   = 1.0E00/SMAX
	    D2SSDS2 = 0.0E00
	  ELSE
	    SH=(1.0E00-D)*PI / ( D*ATAN((A-B)*C) + 2.0E00*ATAN(B*C) - D
     1        *ATAN((A+B)*C) )
	    SK=2.0E00*C*PI*SMAX / (2.0E00*C*PI + SH*(2.0E00*C*
     1        (ATAN((A-B)*C)*(A-B) - ATAN((A+B)*C)*(A+B) +
     2         ATAN((A-B-1.0E00)*C) +
     3         ATAN((A+B-1.0E00)*C)*(A+B-1.0E00) +
     4         ATAN((1.0E00-A+B)*C)*(A-B) ) -
     5         LOG(1.0E00+(A-B)**2*C**2) +
     6         LOG(1.0E00+(A+B)**2*C**2) -
     7         LOG(1.0E00+(A+B-1.0E00)**2*C**2) +
     8         LOG(1.0E00+(1.0E00-A+B)**2*C**2) ))
	    DSSDS=1.0E00 / (SK * (1.0E00 + SH/PI *
     1            (ATAN(C*(SCODE-A-B))+ATAN(C*(A-SCODE-B))) ))
	    D2SSDS2=-DSSDS**3 * C*SH*SK/PI *
     1              (1.0E00/(1.0E00+C**2*(A+B-SCODE)**2) - 
     2               1.0E00/(1.0E00+C**2*(A-B-SCODE)**2))
	    SPHYS=SK/(2.0E00*C*PI) * (2.0E00*C*PI*SCODE + SH*(2.0E00*C*
     1            (ATAN((A-B)*C)*(A-B) -
     2             ATAN((A+B)*C)*(A+B) +
     3             ATAN((A-B-SCODE)*C)*SCODE +
     4             ATAN((A+B-SCODE)*C)*(A+B-SCODE) +
     5             ATAN((SCODE-A+B)*C)*(A-B) ) -
     6             LOG(1.0E00+(A-B)**2*C**2) +
     7             LOG(1.0E00+(A+B)**2*C**2) -
     8             LOG(1.0E00+(A+B-SCODE)**2*C**2) +
     9             LOG(1.0E00+(SCODE-A+B)**2*C**2) ))
	  ENDIF
	  IF (SCODE.EQ.0.0E00) SPHYS=0.0E00
	  IF (SCODE.EQ.1.0E00) SPHYS=SMAX
C----------------------------------------------------------------------	
C  If error then... and calculation of dxx
C----------------------------------------------------------------------	
        CASE DEFAULT
	    WRITE(6,*)'MKGRID: Invalid NGRID number'
            CALL MPI_FINALIZE(IERR)
            STOP
        END SELECT	  
C
	RETURN
	END
C**********************************************************************
	SUBROUTINE XJACOBI
C----------------------------------------------------------------------
C  Returns the x coordinate Jacobian transformation elements.
C----------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
        DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
	DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
     	DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
C
	COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
C----------------------------------------------------------------------
C  Computational grid has evenly spaced x values between zero and one.
C----------------------------------------------------------------------
	ORX=1.0E00/FLOAT(NX-IX-1)
	DO 10 I=0,NX-IX-1
	  EXX(I+2)=FLOAT(I)*ORX
10	CONTINUE
	EXX(2)=0.0E00
	EXX(NX-IX+1)=1.0E00
C
	DO 20 I=2,NX-IX+1
	  SCODE = EXX(I)
	  CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,1)
	  EXX(I)     = SPHYS
	  DXXDX(I)   = DSSDS
	  D2XXDX2(I) = D2SSDS2
	  DDX(I)     = ORX*(1.0E00/DSSDS)
20	CONTINUE
C
	RETURN
	END
C**********************************************************************
	SUBROUTINE YJACOBI
C----------------------------------------------------------------------
C  Returns the y coordinate Jacobian transformation elements.
C  Modified to look like ZJACOBI (M. Rempel).
C----------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
	DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
	DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
	DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
	DIMENSION YY(NPY)
C
	COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
	COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C----------------------------------------------------------------------
C  Full computational grid has evenly spaced y values between zero and one.
C----------------------------------------------------------------------
        ORY=1.0E00/FLOAT(NPY-1)
        DO 10 J=0,NPY-1
                YY(J+1)=FLOAT(J)*ORY
10      CONTINUE
        YY(1)   = 0.0E00
        YY(NPY) = 1.0E00
C----------------------------------------------------------------------
C  Domain split between processors and jacobian calculation
C----------------------------------------------------------------------
        IF (MYPEY.EQ.0) THEN
          DO 20 J=1,IY/2
            SCODE = 0.0E00
            CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,2)
            WYY(J)     = SPHYS
            DYYDY(J)   = DSSDS
            D2YYDY2(J) = D2SSDS2
20        CONTINUE
          DO 30 J=IY/2+1,NRY+IY
            SCODE = YY(J-IY/2)
            CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,2)
            WYY(J)     = SPHYS
            DYYDY(J)   = DSSDS
            D2YYDY2(J) = D2SSDS2
30        CONTINUE
        ELSE IF (MYPEY.EQ.NPEY-1) THEN
          DO 40 J=1,NRY+IY/2
            SCODE = YY(J+NPY-NRY-IY/2)
            CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,2)
            WYY(J)     = SPHYS
            DYYDY(J)   = DSSDS
            D2YYDY2(J) = D2SSDS2
40        CONTINUE
          DO 50 J=NRY+IY/2+1,NRY+IY
            SCODE = 1.0E00
            CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,2)
            WYY(J)     = SPHYS
            DYYDY(J)   = DSSDS
            D2YYDY2(J) = D2SSDS2
50        CONTINUE
        ELSE
          DO 60 J=1,NY
            SCODE = YY(J+MYPEY*NRY-IY/2)
            CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,2)
            WYY(J)     = SPHYS
            DYYDY(J)   = DSSDS
            D2YYDY2(J) = D2SSDS2
60        CONTINUE
        ENDIF
C
        DDY=ORY*(1.0E00/DYYDY)
C
	RETURN
	END
C**********************************************************************
        SUBROUTINE ZJACOBI
C----------------------------------------------------------------------
C  Returns the z coordinate Jacobian transformation elements.
C  Modified (M. Rempel).
C----------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
	DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
	DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
     	DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
        DIMENSION ZZ(NPZ)
C
	COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
        COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C----------------------------------------------------------------------
C  Full computational grid has evenly space z values between zero and one.
C----------------------------------------------------------------------
        ORZ=1.0E00/FLOAT(NPZ-1)
        DO 10 K=0,NPZ-1
                ZZ(K+1)=FLOAT(K)*ORZ
10      CONTINUE
	ZZ(1)   = 0.0E00
        ZZ(NPZ) = 1.0E00
C----------------------------------------------------------------------
C  Domain split between processors and jacobian calculation
C----------------------------------------------------------------------
        IF (MYPEZ.EQ.0) THEN
	  DO 20 K=1,ILAP/2
	    SCODE = 0.0E00
	    CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,3)
	    ZEE(K)     = SPHYS
	    DZZDZ(K)   = DSSDS
	    D2ZZDZ2(K) = D2SSDS2
20	  CONTINUE
	  DO 30 K=ILAP/2+1,NRZ+ILAP
	    SCODE = ZZ(K-ILAP/2)
	    CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,3)
	    ZEE(K)     = SPHYS
	    DZZDZ(K)   = DSSDS
	    D2ZZDZ2(K) = D2SSDS2
30	  CONTINUE
	ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
	  DO 40 K=1,NRZ+ILAP/2
	    SCODE = ZZ(K+NPZ-NRZ-ILAP/2)
	    CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,3)
	    ZEE(K)     = SPHYS
	    DZZDZ(K)   = DSSDS
	    D2ZZDZ2(K) = D2SSDS2
40	  CONTINUE
	  DO 50 K=NRZ+ILAP/2+1,NRZ+ILAP
	    SCODE = 1.0E00
	    CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,3)
	    ZEE(K)     = SPHYS
	    DZZDZ(K)   = DSSDS
	    D2ZZDZ2(K) = D2SSDS2
50	  CONTINUE
	ELSE
	  DO 60 K=1,NZ
	    SCODE = ZZ(K+MYPEZ*NRZ-ILAP/2)
	    CALL MKGRID(SCODE,SPHYS,DSSDS,D2SSDS2,3)
	    ZEE(K)     = SPHYS
	    DZZDZ(K)   = DSSDS
	    D2ZZDZ2(K) = D2SSDS2
60	  CONTINUE
        ENDIF
C
        DDZ=ORZ*(1.0E00/DZZDZ)
C
        RETURN
        END
C**********************************************************************
        SUBROUTINE STATIC
C----------------------------------------------------------------------
C  Generates initial profile.  Upper linear adiabatic temperature
C  gradient, lower linear stable gradient, linked across inversion layer
C  It calls also TUBE which creates the BX, BY, BZ and modifies the
C  thermodynamic variables adequatelly. The calling process is as
C  follows: IF LMAG=.FALSE.: no magnetic field is created.
C           IF LMAG=.TRUE. : if AMPB=/=0: a horizontal layer of magnetic
C                                         field is constructed
C                      	     if AMPB=0: the routine TUBE is called
C				if NTUBE=0 no layer, no tube, but MHD.
C----------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
        include 'mpif.h'
C
        PARAMETER(NV=2)
C
        DIMENSION Y(NV),DYDX(NV)
        DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2		 ,TT(NX,NY,NZ)
        DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
        DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2		 ,FT(NX,NY,NZ)
        DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
     2		 ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
        DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
C
	DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
        DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
C
	DIMENSION RND(NX,NY)
	DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
        DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
        DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
        DIMENSION RKAPA(NZ),DKAPA(NZ)
        DIMENSION T(NPZ),R(NPZ),RK(NPZ),DRK(NPZ)
        DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
     2           ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
     3           ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
     4           ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
C
     5           ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
     6           ,SP26(IPAD)
C
        DIMENSION IIR(97)
C
        DIMENSION ISTATUS(MPI_STATUS_SIZE)
C
        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4		  ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C
	COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
	COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
	COMMON/CMAG/ORM,RM,OBETA,AMPB,BFH,BZP
        COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
	COMMON/GRID/DD,HX,H2X,HY,H2Y,HZ,H2Z,C13,C23,C43
        COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
	COMMON/BOUNDS/XMAX,YMAX,ZMAX
C----------------------------------------------------------------------
C  If LREM: File to store RKAPPA in and rescaled boundary temperatures.
C----------------------------------------------------------------------
	CHARACTER*50 FNAME
	COMMON/SPECIALBOUND/TU,DZTB,DZTU
C
        DZZ=1.0E00/FLOAT(NPZ-1)
        EPS=1.0E-12
c        EPS=1.0E-06
C----------------------------------------------------------------------
C  Zero velocity fields.
C----------------------------------------------------------------------
        RU = 0.0E00
        RV = 0.0E00
        RW = 0.0E00
C
	IF (.NOT.LSHR) THEN
C
	IF (LREM) THEN
C----------------------------------------------------------------------
C  New (M. Rempel).  Factor 8.07*REPR scales RK so that 
C  THETA=ETA=F/rho/Cp/T/Cs.
C----------------------------------------------------------------------
	   DO K=1,NPZ
		SZZ=FLOAT(K-1)/FLOAT(NPZ-1)
		CALL MKGRID(SZZ,SZ,SDZZDZ,SD2ZZDZ2,3)
           	CALL GETBACKGROUND(SZ,TZ,RZ)
           	CALL KAPPA(SZ,RZ,TZ,RKZ)
           	T(K) =TZ
           	R(K) =RZ
           	RK(K)=RKZ
           END DO
C  Derivative of kappa
           RK=RK*8.07*THETA*REPR
           DO K=2,NPZ-1
		SZZ=FLOAT(K-1)/FLOAT(NPZ-1)
		CALL MKGRID(SZZ,SZ,SDZZDZ,SD2ZZDZ2,3)
           	DRK(K)=(RK(K+1)-RK(K-1))*HZ*SDZZDZ
           END DO
	   SZZ=0
	   CALL MKGRID(SZZ,SZ,SDZZDZ,SD2ZZDZ2,3)
           DRK(1)=(-3.0E00*RK(1)+4.0E00*RK(2)-RK(3))*HZ*SDZZDZ
	   SZZ=1
	   CALL MKGRID(SZZ,SZ,SDZZDZ,SD2ZZDZ2,3)
           DRK(NPZ)=(3.0E00*RK(NPZ)-4.0E00*RK(NPZ-1)+RK(NPZ-2))
     2						     *HZ*SDZZDZ
C
           IF (MYPE.EQ.0) THEN
           	CALL SETUP(FINP,FOUT,IPAR,PAR)
		I1    = INDEX(FOUT,'  ')-1
		FNAME = FOUT(1:I1)//'.strat'
		OPEN(999,FILE=FNAME,STATUS='UNKNOWN')
		WRITE(999,*)T,R,RK/REPR,DRK/REPR
		CLOSE(999)
	   ENDIF
C
	   TU   = T(1)
C----------------------------------------------------------------------
C  DZTB and DZTU equal -2/3*DZ*dT/dz.
C----------------------------------------------------------------------
           DZTB = T(NPZ)-4.0/3.0*T(NPZ-1)+1.0/3.0*T(NPZ-2)
           DZTU = T(1)-4.0/3.0*T(2)+1.0/3.0*T(3)
	ELSE
C----------------------------------------------------------------------
C   Assign upper boundary values.
C----------------------------------------------------------------------
           T(1) = 1.0E00
           R(1) = 1.0E00
C----------------------------------------------------------------------
C   Extrapolate for remaining values.
C----------------------------------------------------------------------
           ZZ=0.0E00
	   FFZ=0.0E00
           PLN=0.0E00
	   Y(1)=FFZ
	   Y(2)=PLN
	   DO 10 K=2,NPZ
		HTRY=DZZ
               	CALL DERIVS(ZZ,Y,NV,DYDX)
               	CALL BSSTEP(Y,DYDX,NV,ZZ,HTRY,EPS,HDID,HNEXT)
C
               	IF (HDID.NE.HTRY) THEN
                       	WRITE(6,*)'STATIC: Static structure error.'
			CALL MPI_FINALIZE(IERR)
                       	STOP
               	ENDIF
C
		T(K) = 1.0E00+Y(1)
		R(K) = EXP(Y(2))/(1.0E00+Y(1))
10         CONTINUE
	ENDIF
C
	ELSE
C
	DO 11 K=1,NPZ
		T(K)=1.0E00
		R(K)=1.0E00
11	CONTINUE
C
	ENDIF
C
	TB=T(NPZ)
	RB=R(NPZ)
C----------------------------------------------------------------------
C   Divide among processors.
C----------------------------------------------------------------------
	IF (MYPEZ.EQ.0) THEN
		TT(2:NX-IX+1,2:NY-IY+1,1:ILAP/2)=T(1)
		TT(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ)
     2		       =SPREAD(SPREAD(T(1:NRZ+ILAP/2),1,NX-IX),2,NY-IY)
		RO(2:NX-IX+1,2:NY-IY+1,1:ILAP/2)=R(1)
		RO(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ)
     2		       =SPREAD(SPREAD(R(1:NRZ+ILAP/2),1,NX-IX),2,NY-IY)
        	IF (LREM) THEN
			RKAPA(1:ILAP/2)    = RK(1)
                	RKAPA(ILAP/2+1:NZ) = RK(1:NRZ+ILAP/2)
                	DKAPA(1:ILAP/2)    = DRK(1)
                	DKAPA(ILAP/2+1:NZ) = DRK(1:NRZ+ILAP/2)
		ENDIF
	ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
		TT(2:NX-IX+1,2:NY-IY+1,1:NRZ+ILAP/2)
     2	               =SPREAD(SPREAD(T(NPZ-NRZ-ILAP/2+1:NPZ),1,NX-IX)
     3							      ,2,NY-IY)
		TT(2:NX-IX+1,2:NY-IY+1,NRZ+ILAP/2+1:NRZ+ILAP)=TB
		RO(2:NX-IX+1,2:NY-IY+1,1:NRZ+ILAP/2)
     2		       =SPREAD(SPREAD(R(NPZ-NRZ-ILAP/2+1:NPZ),1,NX-IX)
     3							      ,2,NY-IY)
		RO(2:NX-IX+1,2:NY-IY+1,NRZ+ILAP/2+1:NRZ+ILAP)=RB
		IF (LREM) THEN
			RKAPA(1:NRZ+ILAP/2) = RK(NPZ-NRZ-ILAP/2+1:NPZ)
                	RKAPA(NRZ+ILAP/2+1:NRZ+ILAP) = RK(NPZ)
                	DKAPA(1:NRZ+ILAP/2) = DRK(NPZ-NRZ-ILAP/2+1:NPZ)
                	DKAPA(NRZ+ILAP/2+1:NRZ+ILAP) = DRK(NPZ)
		ENDIF
	ELSE
		I1=MYPEZ*NRZ-ILAP/2+1
		I2=(MYPEZ+1)*NRZ+ILAP/2
		TT(2:NX-IX+1,2:NY-IY+1,:)
     2		       	     =SPREAD(SPREAD(T(I1:I2),1,NX-IX),2,NY-IY)
		RO(2:NX-IX+1,2:NY-IY+1,:)
     2			     =SPREAD(SPREAD(R(I1:I2),1,NX-IX),2,NY-IY)
		IF (LREM) THEN
			RKAPA(:) = RK(I1:I2)
			DKAPA(:) = DRK(I1:I2)
		ENDIF
	ENDIF
C----------------------------------------------------------------------
C   Initiate magnetic field layer of full-width at half-maximum equal
C   to BFH and of maximum amplitude AMPB at depth BZP.
C----------------------------------------------------------------------
        IF (LMAG) THEN
                IF (AMPB.NE.0.0E00) THEN
                        CLN=-4.0E00*LOG(2.0E00)/BFH/BFH
                        DO 15 K=ILAP/2+1,NZ-ILAP/2
                                BX(:,:,K)=AMPB*EXP(CLN*(ZEE(K)-BZP)**2)
15                      CONTINUE
                        BY=0.0E00
                        BZ=0.0E00
C
                        DO 16 K=ILAP/2+1,NZ-ILAP/2
                        DO 16 J=2,NY-IY+1
                        DO 16 I=2,NX-IX+1
                                RO(I,J,K)=RO(I,J,K)
     2                            -BX(I,J,K)*BX(I,J,K)/TT(I,J,K)*OBETA
16                      CONTINUE
                ELSE
		  CALL TUBE
                ENDIF
        ENDIF
C----------------------------------------------------------------------
C  Initiate temperature perturbations.  Switch ISW to set either random
C  or sinusoidal perturbations (M. Rempel).  Random  perturbations 
C  changed to yield identical initial state independent of processor 
C  configuration (Rast 3/20/02).
C----------------------------------------------------------------------
        IF (AMPT.NE.0.0E00) THEN
C
	ISW=1
	IF (ISW.EQ.1) THEN
C
	IDUM=-062659
	DO 18 NNZ=1,NPEZ
	   DO 19 K=ILAP/2+1,NZ-ILAP/2
	      ITAG=K
	      DO 20 NNY=1,NPEY
	      	 IF (MYPE.EQ.0) THEN
	   	    DO 21 J=IY/2+1,NY-IY/2
	   	    DO 21 I=IX/2+1,NX-IX/2
		       RND(I,J)=1.0E00+AMPT*(RAN2(IDUM,IIY,IIR)-0.5E00)
21	   	    CONTINUE
		    IF ((NNY+NPEY*(NNZ-1)-1).NE.0) THEN
	      	     CALL MPI_SEND(RND,NX*NY,MPISIZE,NNY+NPEY*(NNZ-1)-1
     2				             ,ITAG,MPI_COMM_WORLD,IERR)
		    ELSE
		     WW1(:,:,K)=RND
		    ENDIF
		 ELSE	
		    IF (MYPE.EQ.(NNY+NPEY*(NNZ-1)-1)) THEN
		       CALL MPI_RECV(WW1(:,:,K),NX*NY,MPISIZE,0,ITAG
     2					  ,MPI_COMM_WORLD,ISTATUS,IERR)
		    ENDIF
		 ENDIF
20	     CONTINUE
19         CONTINUE
18      CONTINUE
C
	DO 40 K=ILAP/2+1,NZ-ILAP/2
        DO 40 J=IY/2+1,NY-IY/2
        DO 40 I=IX/2+1,NX-IX/2
		IF (.NOT.LSHR) THEN
			TT(I,J,K)=TT(I,J,K)*WW1(I,J,K)
		ELSE
			RV(I,J,K)=RV(I,J,K)+WW1(I,J,K)-1.0E00
		ENDIF
40	CONTINUE
C
	ELSE
C
	RPI=2.0E00*ASIN(1.0E00)
	DO K=ILAP/2+1,NZ-ILAP/2
	   DO J=2,NY-IY+1
	      DO I=2,NX-IX+1
                    RKY=WYY(J)/YMAX*FLOAT(NPY)/FLOAT(NPY+1)*2.0*RPI*8.0
                    RKZ=ZEE(K)/ZMAX*2.0*RPI
		    IF (.NOT.LSHR) THEN
                    TT(I,J,K)=TT(I,J,K)+AMPT*SIN(RKY)*SIN(RKZ)/RO(I,J,K)
		    ELSE 
		    RV(I,J,K)=RV(I,J,K)+AMPT*SIN(RKY)*SIN(RKZ)/RO(I,J,K)
		    ENDIF
                 END DO
              END DO
           END DO
C
        ENDIF
C
        ENDIF
C
	RETURN                                   
	END        
C**********************************************************************
        SUBROUTINE DERIVS(ZZ,Y,NV,DYDX)
C
        INCLUDE '3dmhdparam.f'
C
        DIMENSION Y(NV),DYDX(NV)
	DIMENSION RKAPA(NZ),DKAPA(NZ)
C
	COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
        COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
C-----------------------------------------------------------------------
C  Calculate Z AND DZZDZ
C-----------------------------------------------------------------------
	CALL MKGRID(ZZ,Z,DZZDZ,D2ZZDZ2,3)
C----------------------------------------------------------------------
C  Calculate dF/dZZ and dlnP/dZZ (DYDX).
C----------------------------------------------------------------------
        IF (PZP.EQ.0.0E00) THEN
                RKAP=1.0E00
        ELSE
                RKAP=1.0E00
     2              +(RKAPST-1.0E00)/2.0E00*(1.0E00+TANH((Z-PZP)/SIGMA))
        ENDIF
        DYDX(1)=THETA/RKAP/DZZDZ
	DYDX(2)=GRAV/(1.0E00+Y(1))/DZZDZ
C
        RETURN
        END
C**********************************************************************
	SUBROUTINE STEP
C
	INCLUDE '3dmhdparam.f'
C
	include 'mpif.h'
C
	DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2           ,TT(NX,NY,NZ)
	DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
	DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2           ,FT(NX,NY,NZ)
        DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
     2           ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
        DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
C
	DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
        DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
C
	DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
	DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
	DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
	DIMENSION RKAPA(NZ),DKAPA(NZ)
C
        DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
     2           ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
     3           ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
     4           ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
C
     5           ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
     6           ,SP26(IPAD)
C
        DIMENSION WMIN(4),WMOUT(4)
C
	COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4            ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C
	COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
	COMMON/CMAG/ORM,RM,OBETA,AMPB,BFH,BZP
	COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
	COMMON/GRID/DD,HX,H2X,HY,H2Y,HZ,H2Z,C13,C23,C43
        COMMON/TRACE/UMACH
        COMMON/CTIM/DT,TIMT,TIMC,TIMI
        COMMON/ITER/NTOTAL,NSTEP0,NIT
	COMMON/RUNGKU/GAM1,GAM2,GAM3,ZETA1,ZETA2
	COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C----------------------------------------------------------------------
C  Calculate the time-step and maximum Mach number.
C----------------------------------------------------------------------
        UMACH=0.0E00
	RMIN=1.0E09
	VMAX=0.0E00
	OGAMMA=1.0E00/GAMMA
C
c        UU=SQRT((RU**2+RV**2+RW**2)/RO**2)
        DO 10 K=ILAP/2+1,NZ-ILAP/2
        DO 10 J=2,NY-IY+1
        DO 10 I=2,NX-IX+1
                UU(I,J,K)=(RU(I,J,K)*RU(I,J,K)+RV(I,J,K)*RV(I,J,K)
     2                   +RW(I,J,K)*RW(I,J,K))*(1.0E00/RO(I,J,K))**2
10      CONTINUE
C
c        VV=SQRT(GAMMA*TT)
c        UMACH=MAXVAL(UU/VV)
        DO 20 K=ILAP/2+1,NZ-ILAP/2
        DO 20 J=2,NY-IY+1
        DO 20 I=2,NX-IX+1
		UMACH=MAX(UMACH,OGAMMA*UU(I,J,K)/TT(I,J,K))
20      CONTINUE
	UMACH=SQRT(UMACH)
C
c        VMAX=MAXVAL(UU+VV)
	IF (LMAG) THEN
	     DO 30 K=ILAP/2+1,NZ-ILAP/2
             DO 30 J=2,NY-IY+1
             DO 30 I=2,NX-IX+1
               	UU(I,J,K)=UU(I,J,K)+GAMMA*TT(I,J,K)
     2			       	     +(BX(I,J,K)*BX(I,J,K)
     3				      +BY(I,J,K)*BY(I,J,K)
     4                  	      +BZ(I,J,K)*BZ(I,J,K))
     5					   *OBETA/RO(I,J,K)
     6                               +2.0E00*SQRT(UU(I,J,K)
     7					     *(GAMMA*TT(I,J,K)
     8					       +(BX(I,J,K)*BX(I,J,K)
     9                          		+BY(I,J,K)*BY(I,J,K)
     1                                          +BZ(I,J,K)*BZ(I,J,K))
     2						    *OBETA/RO(I,J,K)))
30           CONTINUE
	ELSE
       	     DO 40 K=ILAP/2+1,NZ-ILAP/2
       	     DO 40 J=2,NY-IY+1
       	     DO 40 I=2,NX-IX+1
	    	UU(I,J,K)=UU(I,J,K)+GAMMA*TT(I,J,K)
     2		   	     	     +2.0E00*SQRT(UU(I,J,K)
     3					     *GAMMA*TT(I,J,K))
40	     CONTINUE
	ENDIF
        DO 50 K=ILAP/2+1,NZ-ILAP/2
        DO 50 J=2,NY-IY+1
        DO 50 I=2,NX-IX+1
                VMAX=MAX(VMAX,UU(I,J,K))
		RMIN=MIN(RMIN,RO(I,J,K))
50      CONTINUE
        VMAX=SQRT(VMAX)
C
	CALL  MPI_ALLREDUCE(UMACH,WMOUT,1,MPISIZE,MPI_MAX,
     2			    MPI_COMM_WORLD,IERR)
        UMACH=WMOUT(1)
C
	ISW=1
	IF (ISW.EQ.0) THEN
C----------------------------------------------------------------------
C  Find global minimum timestep. Note that this approach yields a 
C  timestep possibly smaller than the minimum pointwise value.
C  Values of DD and RKAPM are unchanged from initial time step evaluation.
C----------------------------------------------------------------------
	WMIN(1)=RMIN
	WMIN(2)=-1.0E00*VMAX
	CALL MPI_ALLREDUCE(WMIN,WMOUT,2,MPISIZE,MPI_MIN,
     2			   MPI_COMM_WORLD,IERR)
	RMIN=WMOUT(1)
	VMAX=-1.0E00*WMOUT(2)
C
        DT1=DD/VMAX
        IF (LREM) THEN
C----------------------------------------------------------------------
C  Thermal conductivity split into a constant and depth dependant part.
C----------------------------------------------------------------------
                DT2=0.5*DD*DD*REPR*CV*RMIN/(1.0+RKAPM)
        ELSE
                DT2=0.5E00*DD*DD*REPR*CV*RMIN/RKAPM
        ENDIF
        DT3=0.375E00*DD*DD*RE*RMIN
	IF (LMAG) THEN
                DT4=0.5E00*DD*DD*RM
                DT=SF*MIN(DT1,DT2,DT3,DT4)
        ELSE
                DT=SF*MIN(DT1,DT2,DT3)
        ENDIF
	ELSE
C----------------------------------------------------------------------
C  Find the pointwise minimum timestep.
C----------------------------------------------------------------------
	DO 2 K=ILAP/2+1,NZ-ILAP/2
        DO 2 J=2,NY-IY+1
        DO 2 I=2,NX-IX+1
		IF (NX.GT.IX+1) THEN
                	DD=MIN(DDX(I),DDY(J),DDZ(K))
		ELSE
			DD=MIN(DDY(J),DDZ(K))
		ENDIF
                WW1(I,J,K)=DD/SQRT(UU(I,J,K))
        IF (LREM) THEN
C----------------------------------------------------------------------
C  Thermal conductivity split into a constant and depth dependant part.
C----------------------------------------------------------------------
                WW2(I,J,K)=0.5*DD*DD*REPR*CV*RO(I,J,K)/(1.0+RKAPA(K))
        ELSE
                WW2(I,J,K)=0.5E00*DD*DD*REPR*CV*RO(I,J,K)/RKAPA(K)
        ENDIF
                WW3(I,J,K)=0.375E00*DD*DD*RE*RO(I,J,K)
                IF (LMAG) THEN
                        VV(I,J,K)=0.5E00*DD*DD*RM
                ENDIF
2       CONTINUE
C
        WMIN(1)=MINVAL(WW1(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
        WMIN(2)=MINVAL(WW2(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
        WMIN(3)=MINVAL(WW3(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	MINCNT=3
	IF (LMAG) THEN
	     WMIN(4)=MINVAL(VV(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	     MINCNT=4
	ENDIF
	CALL MPI_ALLREDUCE(WMIN,WMOUT,MINCNT,MPISIZE,MPI_MIN,
     2			   MPI_COMM_WORLD,IERR)
        IF (LMAG) THEN
             DT=SF*MIN(WMOUT(1),WMOUT(2),WMOUT(3),WMOUT(4))
        ELSE
             DT=SF*MIN(WMOUT(1),WMOUT(2),WMOUT(3))
        ENDIF
	ENDIF
C----------------------------------------------------------------------
C  Third-order Runge-Kutta timestepping scheme of Wray (Spalart, P.R.,
C  Moser, R.D., & Rogers M.M., Spectral Methods for the Navier-Stokes
C  Equations with One Infinite and Two Periodic Directions,
C  J. Comp. Phys., 96 297-324 1991).
C  Calculate the first Runge-Kutta substep.
C----------------------------------------------------------------------	
        DO 71 K=ILAP/2+1,NZ-ILAP/2
        DO 71 J=2,NY-IY+1
        DO 71 I=2,NX-IX+1
                ZRU(I,J,K)=RU(I,J,K)
71      CONTINUE
        DO 81 K=ILAP/2+1,NZ-ILAP/2
        DO 81 J=2,NY-IY+1
        DO 81 I=2,NX-IX+1
                ZRV(I,J,K)=RV(I,J,K)
81      CONTINUE
        DO 91 K=ILAP/2+1,NZ-ILAP/2
        DO 91 J=2,NY-IY+1
        DO 91 I=2,NX-IX+1
                ZRW(I,J,K)=RW(I,J,K)
91     CONTINUE
        DO 101 K=ILAP/2+1,NZ-ILAP/2
        DO 101 J=2,NY-IY+1
        DO 101 I=2,NX-IX+1
                ZRO(I,J,K)=RO(I,J,K)
101     CONTINUE
        DO 111 K=ILAP/2+1,NZ-ILAP/2
        DO 111 J=2,NY-IY+1
        DO 111 I=2,NX-IX+1
                ZTT(I,J,K)=TT(I,J,K)
111     CONTINUE
        IF (LMAG) THEN
                DO 121 K=ILAP/2+1,NZ-ILAP/2
                DO 121 J=2,NY-IY+1
                DO 121 I=2,NX-IX+1
                        ZBX(I,J,K)=BX(I,J,K)
121             CONTINUE
                DO 131 K=ILAP/2+1,NZ-ILAP/2
                DO 131 J=2,NY-IY+1
                DO 131 I=2,NX-IX+1
                        ZBY(I,J,K)=BY(I,J,K)
131             CONTINUE
                DO 141 K=ILAP/2+1,NZ-ILAP/2
                DO 141 J=2,NY-IY+1
                DO 141 I=2,NX-IX+1
                        ZBZ(I,J,K)=BZ(I,J,K)
141             CONTINUE
        ENDIF
C
	CALL FLUXES
C
	COEF=GAM1*DT
C
	DO 70 K=ILAP/2+1,NZ-ILAP/2
        DO 70 J=2,NY-IY+1
        DO 70 I=2,NX-IX+1
                RU(I,J,K)=ZRU(I,J,K)+COEF*FU(I,J,K)
70      CONTINUE
        DO 80 K=ILAP/2+1,NZ-ILAP/2
        DO 80 J=2,NY-IY+1
        DO 80 I=2,NX-IX+1
                RV(I,J,K)=ZRV(I,J,K)+COEF*FV(I,J,K)
80      CONTINUE
        DO 90 K=ILAP/2+1,NZ-ILAP/2
        DO 90 J=2,NY-IY+1
        DO 90 I=2,NX-IX+1
                RW(I,J,K)=ZRW(I,J,K)+COEF*FW(I,J,K)
90      CONTINUE
        DO 100 K=ILAP/2+1,NZ-ILAP/2
        DO 100 J=2,NY-IY+1
        DO 100 I=2,NX-IX+1
                RO(I,J,K)=ZRO(I,J,K)+COEF*FR(I,J,K)
100     CONTINUE
        DO 110 K=ILAP/2+1,NZ-ILAP/2
        DO 110 J=2,NY-IY+1
        DO 110 I=2,NX-IX+1
                TT(I,J,K)=ZTT(I,J,K)+COEF*FT(I,J,K)
110     CONTINUE
	IF (LMAG) THEN
        	DO 120 K=ILAP/2+1,NZ-ILAP/2
        	DO 120 J=2,NY-IY+1
        	DO 120 I=2,NX-IX+1
                	BX(I,J,K)=ZBX(I,J,K)+COEF*WW1(I,J,K)
120     	CONTINUE
        	DO 130 K=ILAP/2+1,NZ-ILAP/2
        	DO 130 J=2,NY-IY+1
        	DO 130 I=2,NX-IX+1
                	BY(I,J,K)=ZBY(I,J,K)+COEF*WW2(I,J,K)
130     	CONTINUE
        	DO 140 K=ILAP/2+1,NZ-ILAP/2
        	DO 140 J=2,NY-IY+1
        	DO 140 I=2,NX-IX+1
               		BZ(I,J,K)=ZBZ(I,J,K)+COEF*WW3(I,J,K)
140     	CONTINUE
	ENDIF
C
        CALL BCON
C
        CALL COMMUNICATE
C----------------------------------------------------------------------
C   Calculate second Runge-Kutta substep.
C----------------------------------------------------------------------
	COEF=ZETA1*DT
C
	DO 190 K=ILAP/2+1,NZ-ILAP/2
        DO 190 J=2,NY-IY+1
        DO 190 I=2,NX-IX+1
                ZRU(I,J,K)=RU(I,J,K)+COEF*FU(I,J,K)
190      CONTINUE
        DO 200 K=ILAP/2+1,NZ-ILAP/2
        DO 200 J=2,NY-IY+1
        DO 200 I=2,NX-IX+1
                ZRV(I,J,K)=RV(I,J,K)+COEF*FV(I,J,K)
200      CONTINUE
        DO 210 K=ILAP/2+1,NZ-ILAP/2
        DO 210 J=2,NY-IY+1
        DO 210 I=2,NX-IX+1
                ZRW(I,J,K)=RW(I,J,K)+COEF*FW(I,J,K)
210      CONTINUE
        DO 220 K=ILAP/2+1,NZ-ILAP/2
        DO 220 J=2,NY-IY+1
        DO 220 I=2,NX-IX+1
                ZRO(I,J,K)=RO(I,J,K)+COEF*FR(I,J,K)
220     CONTINUE
        DO 230 K=ILAP/2+1,NZ-ILAP/2
        DO 230 J=2,NY-IY+1
        DO 230 I=2,NX-IX+1
                ZTT(I,J,K)=TT(I,J,K)+COEF*FT(I,J,K)
230     CONTINUE
        IF (LMAG) THEN
                DO 240 K=ILAP/2+1,NZ-ILAP/2
                DO 240 J=2,NY-IY+1
                DO 240 I=2,NX-IX+1
                        ZBX(I,J,K)=BX(I,J,K)+COEF*WW1(I,J,K)
240             CONTINUE
                DO 250 K=ILAP/2+1,NZ-ILAP/2
                DO 250 J=2,NY-IY+1
                DO 250 I=2,NX-IX+1
                        ZBY(I,J,K)=BY(I,J,K)+COEF*WW2(I,J,K)
250             CONTINUE
                DO 260 K=ILAP/2+1,NZ-ILAP/2
                DO 260 J=2,NY-IY+1
                DO 260 I=2,NX-IX+1
                        ZBZ(I,J,K)=BZ(I,J,K)+COEF*WW3(I,J,K)
260             CONTINUE
        ENDIF
C	
	CALL FLUXES
C	
	COEF=GAM2*DT
C
        DO 270 K=ILAP/2+1,NZ-ILAP/2
        DO 270 J=2,NY-IY+1
        DO 270 I=2,NX-IX+1
                RU(I,J,K)=ZRU(I,J,K)+COEF*FU(I,J,K)
270      CONTINUE
        DO 280 K=ILAP/2+1,NZ-ILAP/2
        DO 280 J=2,NY-IY+1
        DO 280 I=2,NX-IX+1
                RV(I,J,K)=ZRV(I,J,K)+COEF*FV(I,J,K)
280      CONTINUE
        DO 290 K=ILAP/2+1,NZ-ILAP/2
        DO 290 J=2,NY-IY+1
        DO 290 I=2,NX-IX+1
                RW(I,J,K)=ZRW(I,J,K)+COEF*FW(I,J,K)
290      CONTINUE
        DO 300 K=ILAP/2+1,NZ-ILAP/2
        DO 300 J=2,NY-IY+1
        DO 300 I=2,NX-IX+1
                RO(I,J,K)=ZRO(I,J,K)+COEF*FR(I,J,K)
300      CONTINUE
        DO 310 K=ILAP/2+1,NZ-ILAP/2
        DO 310 J=2,NY-IY+1
        DO 310 I=2,NX-IX+1
                TT(I,J,K)=ZTT(I,J,K)+COEF*FT(I,J,K)
310      CONTINUE
	IF (LMAG) THEN
		DO 320 K=ILAP/2+1,NZ-ILAP/2
        	DO 320 J=2,NY-IY+1
        	DO 320 I=2,NX-IX+1
                	BX(I,J,K)=ZBX(I,J,K)+COEF*WW1(I,J,K)
320      	CONTINUE
		DO 330 K=ILAP/2+1,NZ-ILAP/2
                DO 330 J=2,NY-IY+1
                DO 330 I=2,NX-IX+1
                        BY(I,J,K)=ZBY(I,J,K)+COEF*WW2(I,J,K)
330      	CONTINUE
		DO 340 K=ILAP/2+1,NZ-ILAP/2
                DO 340 J=2,NY-IY+1
                DO 340 I=2,NX-IX+1
                        BZ(I,J,K)=ZBZ(I,J,K)+COEF*WW3(I,J,K)
340      	CONTINUE
	ENDIF
C
	CALL BCON
C
        CALL COMMUNICATE
C----------------------------------------------------------------------
C   Calculate third Runge-Kutta substep.
C----------------------------------------------------------------------
        COEF=ZETA2*DT
C
        DO 390 K=ILAP/2+1,NZ-ILAP/2
        DO 390 J=2,NY-IY+1
        DO 390 I=2,NX-IX+1
                ZRU(I,J,K)=RU(I,J,K)+COEF*FU(I,J,K)
390      CONTINUE
        DO 400 K=ILAP/2+1,NZ-ILAP/2
        DO 400 J=2,NY-IY+1
        DO 400 I=2,NX-IX+1
                ZRV(I,J,K)=RV(I,J,K)+COEF*FV(I,J,K)
400      CONTINUE
        DO 410 K=ILAP/2+1,NZ-ILAP/2
        DO 410 J=2,NY-IY+1
        DO 410 I=2,NX-IX+1
                ZRW(I,J,K)=RW(I,J,K)+COEF*FW(I,J,K)
410      CONTINUE
        DO 420 K=ILAP/2+1,NZ-ILAP/2
        DO 420 J=2,NY-IY+1
        DO 420 I=2,NX-IX+1
                ZRO(I,J,K)=RO(I,J,K)+COEF*FR(I,J,K)
420     CONTINUE
        DO 430 K=ILAP/2+1,NZ-ILAP/2
        DO 430 J=2,NY-IY+1
        DO 430 I=2,NX-IX+1
                ZTT(I,J,K)=TT(I,J,K)+COEF*FT(I,J,K)
430     CONTINUE
        IF (LMAG) THEN
                DO 440 K=ILAP/2+1,NZ-ILAP/2
                DO 440 J=2,NY-IY+1
                DO 440 I=2,NX-IX+1
                        ZBX(I,J,K)=BX(I,J,K)+COEF*WW1(I,J,K)
440             CONTINUE
                DO 450 K=ILAP/2+1,NZ-ILAP/2
                DO 450 J=2,NY-IY+1
                DO 450 I=2,NX-IX+1
                        ZBY(I,J,K)=BY(I,J,K)+COEF*WW2(I,J,K)
450             CONTINUE
                DO 460 K=ILAP/2+1,NZ-ILAP/2
                DO 460 J=2,NY-IY+1
                DO 460 I=2,NX-IX+1
                        ZBZ(I,J,K)=BZ(I,J,K)+COEF*WW3(I,J,K)
460             CONTINUE
        ENDIF
C
	CALL FLUXES
C
        COEF=GAM3*DT
C
        DO 470 K=ILAP/2+1,NZ-ILAP/2
        DO 470 J=2,NY-IY+1
        DO 470 I=2,NX-IX+1
                RU(I,J,K)=ZRU(I,J,K)+COEF*FU(I,J,K)
470      CONTINUE
        DO 480 K=ILAP/2+1,NZ-ILAP/2
        DO 480 J=2,NY-IY+1
        DO 480 I=2,NX-IX+1
                RV(I,J,K)=ZRV(I,J,K)+COEF*FV(I,J,K)
480      CONTINUE
        DO 490 K=ILAP/2+1,NZ-ILAP/2
        DO 490 J=2,NY-IY+1
        DO 490 I=2,NX-IX+1
                RW(I,J,K)=ZRW(I,J,K)+COEF*FW(I,J,K)
490      CONTINUE
        DO 500 K=ILAP/2+1,NZ-ILAP/2
        DO 500 J=2,NY-IY+1
        DO 500 I=2,NX-IX+1
                RO(I,J,K)=ZRO(I,J,K)+COEF*FR(I,J,K)
500      CONTINUE
        DO 510 K=ILAP/2+1,NZ-ILAP/2
        DO 510 J=2,NY-IY+1
        DO 510 I=2,NX-IX+1
                TT(I,J,K)=ZTT(I,J,K)+COEF*FT(I,J,K)
510      CONTINUE
        IF (LMAG) THEN
                DO 520 K=ILAP/2+1,NZ-ILAP/2
                DO 520 J=2,NY-IY+1
                DO 520 I=2,NX-IX+1
                        BX(I,J,K)=ZBX(I,J,K)+COEF*WW1(I,J,K)
520             CONTINUE
                DO 530 K=ILAP/2+1,NZ-ILAP/2
                DO 530 J=2,NY-IY+1
                DO 530 I=2,NX-IX+1
                        BY(I,J,K)=ZBY(I,J,K)+COEF*WW2(I,J,K)
530             CONTINUE
                DO 540 K=ILAP/2+1,NZ-ILAP/2
                DO 540 J=2,NY-IY+1
                DO 540 I=2,NX-IX+1
                        BZ(I,J,K)=ZBZ(I,J,K)+COEF*WW3(I,J,K)
540             CONTINUE
        ENDIF
C
        CALL BCON
C
        CALL COMMUNICATE
C
	NIT=NIT+1
	TIMC=TIMC+DT
	TIMT=TIMI+TIMC
C
	RETURN                                            
	END	  	                    
C**********************************************************************
	SUBROUTINE FLUXES
C
	INCLUDE '3dmhdparam.f'
C
        include 'mpif.h'
C
	DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2           ,TT(NX,NY,NZ)
        DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
        DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2           ,FT(NX,NY,NZ)
        DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
     2           ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
        DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
C
        DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
        DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
C
        DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
        DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
        DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
	DIMENSION WWY(NY),WWZ(NZ)
        DIMENSION RKAPA(NZ),DKAPA(NZ)
	DIMENSION TTM(NZ),ROM(NZ),HRAD(NZ),FCONM(NZ),FRADM(NZ),
     2		  ALPHA(NZ),CORRECT(NZ),ADDSUM(NPE),ENDVAL(NPE)
        DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
     2           ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
     3           ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
     4           ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
C
     5           ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
     6           ,SP26(IPAD)
C
        DIMENSION ISTATUS(MPI_STATUS_SIZE)
c sleak: fix common block size mismatch for splinex and spliney
 	DIMENSION KLOX(NPX+3),KHIX(NPX+3),HHX(NPX+3),SIGX(NPX),
     2		  AAX(NPX+3),BBX(NPX+3),XPINV(NPX)
        DIMENSION KLOY(NRY+3),KHIY(NRY+3),HHY(NRY+3),SIGY(NRY),
     2            AAY(NRY+3),BBY(NRY+3),YPINV(NRY)
C
        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4            ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C
        COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
	COMMON/GRID/DD,HX,H2X,HY,H2Y,HZ,H2Z,C13,C23,C43
	COMMON/CPER/TP,XP,YP,ZP,TC,QFH,HH
	COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
	COMMON/CROT/OMX,OMZ
	COMMON/CMAG/ORM,RM,OBETA,AMPB,BFH,BZP
	COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
	COMMON/BOUNDS/XMAX,YMAX,ZMAX
        COMMON/BCT/IXC,IYC,IZC,ITC,IBC
	COMMON/SPLINEX/KLOX,KHIX,HHX,SIGX,AAX,BBX,XPINV,XHH,ISEGX
        COMMON/SPLINEY/KLOY,KHIY,HHY,SIGY,AAY,BBY,YPINV,YHH,ISEGY
        COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
	COMMON/CTIM/DT,TIMT,TIMC,TIMI
	COMMON/RELAX/TSTART,TOFF,RLAX
C
	DATA ICALL/0/
	DATA TTM_MIN/1D99/
        DATA TTM_MAX/-1D99/
        DATA TFAC/1D0/
        SAVE ICALL,TTM_MAX,TTM_MIN,TFAC
C----------------------------------------------------------------------
C  Take care of FR.
C----------------------------------------------------------------------
c       FR=(CSHIFT(RU,1,1)-CSHIFT(RU,-1,1))*HX*DXXDX
c       FR=-FR
c	WW1=(CSHIFT(RV,1,2)-CSHIFT(RV,-1,2))*HY*DYYDY
        IF ((IXC.EQ.0).AND.(IYC.EQ.0)) THEN
		DO 10 K=ILAP/2+1,NZ-ILAP/2
		DO 10 J=2,NY-IY+1
			TMPY=HY*DYYDY(J)
		DO 10 I=2,NX-IX+1
			FR(I,J,K)=(RU(I-1,J,K)-RU(I+1,J,K))
     2							*HX*DXXDX(I)
			FR(I,J,K)=FR(I,J,K)-(RV(I,J+1,K)-RV(I,J-1,K))
     2								*TMPY
10		CONTINUE
	ELSE
		WRITE(6,*)'FLUXES: Non-periodic horizontal boundaries'
		CALL MPI_FINALIZE(IERR)
		STOP
	ENDIF
c       WW1=(CSHIFT(RW,1,3)-CSHIFT(RW,-1,3))*HZ*DZZDZ
        DO 40 K=ILAP/2+1,NZ-ILAP/2
		TMPZ=HZ*DZZDZ(K)
        DO 40 J=2,NY-IY+1
        DO 40 I=2,NX-IX+1
                WW1(I,J,K)=(RW(I,J,K+1)-RW(I,J,K-1))*TMPZ
40      CONTINUE
        IF (MYPEZ.EQ.0) THEN
		TMPZ=HZ*DZZDZ(ILAP/2+1)
        	DO 50 J=2,NY-IY+1
        	DO 50 I=2,NX-IX+1
                	WW1(I,J,ILAP/2+1)=(4.0E00*RW(I,J,ILAP/2+2)
     2			      		   -RW(I,J,ILAP/2+3))*TMPZ
50              CONTINUE
	ENDIF
        IF (MYPEZ.EQ.NPEZ-1) THEN
		TMPZ=HZ*DZZDZ(NZ-ILAP/2)
		DO 60 J=2,NY-IY+1
		DO 60 I=2,NX-IX+1
                	WW1(I,J,NZ-ILAP/2)=(RW(I,J,NZ-ILAP/2-2)
     2                  	-4.0E00*RW(I,J,NZ-ILAP/2-1))*TMPZ
60		CONTINUE
        ENDIF
c       FR=FR-WW1
	DO 70 K=ILAP/2+1,NZ-ILAP/2
	DO 70 J=2,NY-IY+1
	DO 70 I=2,NX-IX+1
		FR(I,J,K)=FR(I,J,K)-WW1(I,J,K)
70	CONTINUE
C----------------------------------------------------------------------
C  Diffusion of density perturbations. Density ghost points loaded by 
C  mirror points so that diffusion on boundary is centered to avoid
C  forward/backward difference instability (ITP Santa Barbara).
C----------------------------------------------------------------------
        IF ((ID.NE.0).OR.LREM) THEN
	     IF (MYPEZ.EQ.0) THEN
            	RO(:,:,1:ILAP/2)=RO(:,:,ILAP/2+2:ILAP+1)
	     ENDIF
	     IF (MYPEZ.EQ.NPEZ-1) THEN
	    	RO(:,:,NZ-ILAP/2+1:NZ)=RO(:,:,NZ-ILAP:NZ-ILAP/2-1)
	     ENDIF
C
	     CALL HORIZONTAL_MEAN(ROM,RO)
C
             DO 9120 K=1,NZ
             DO 9120 J=1,NY
             DO 9120 I=1,NX
                   WW1(I,J,K)=RO(I,J,K)-ROM(K)
9120         CONTINUE
	ENDIF
        IF (ID.NE.0) THEN
C-----------------------------------------------------------------------
C  Specify vertical variation in density diffusion (=1, goes as 1/rhobar,
C  =2 exponentially confined to upper and lower boundaries).
C-----------------------------------------------------------------------
             IF (ID.EQ.1) THEN
                WWZ=1.0E00/ROM
             ENDIF
             IF (ID.EQ.2) THEN
                SGM=0.1E00
                CLN=-4.0E00*LOG(2.0E00)/SGM/SGM
                WWZ=EXP(CLN*ZEE**2)+EXP(CLN*(ZEE-ZMAX)**2)
             ENDIF
        IF (DH.NE.0.0E00) THEN
             DO 9130 K=ILAP/2+1,NZ-ILAP/2
             DO 9130 J=2,NY-IY+1
             DO 9130 I=2,NX-IX+1
                   WW2(I,J,K)=(WW1(I+1,J,K)-WW1(I-1,J,K))*HX*D2XXDX2(I)
     2                   +(WW1(I+1,J,K)-2.0E00*WW1(I,J,K)+WW1(I-1,J,K))
     3                              *H2X*DXXDX(I)*DXXDX(I)
9130    CONTINUE
             DO 9135 K=ILAP/2+1,NZ-ILAP/2
             DO 9135 J=2,NY-IY+1
		   TMPY1=HY*D2YYDY2(J)
		   TMPY2=H2Y*DYYDY(J)*DYYDY(J)
             DO 9135 I=2,NX-IX+1
                   WW3(I,J,K)=(WW1(I,J+1,K)-WW1(I,J-1,K))*TMPY1
     2                   +(WW1(I,J+1,K)-2.0E00*WW1(I,J,K)+WW1(I,J-1,K))
     3                              *TMPY2
9135         CONTINUE
             DO 9140 K=ILAP/2+1,NZ-ILAP/2
                TMP=DH*ORE*WWZ(K)
             DO 9140 J=2,NY-IY+1
             DO 9140 I=2,NX-IX+1
                FR(I,J,K)=FR(I,J,K)+(WW2(I,J,K)+WW3(I,J,K))*TMP
9140         CONTINUE
        ENDIF
        IF (DV.NE.0.0E00) THEN
               DO 9145 K=ILAP/2+1,NZ-ILAP/2
		   TMPZ1=HZ*D2ZZDZ2(K)
		   TMPZ2=H2Z*DZZDZ(K)*DZZDZ(K)
               DO 9145 J=2,NY-IY+1
               DO 9145 I=2,NX-IX+1
                   WW2(I,J,K)=(WW1(I,J,K+1)-WW1(I,J,K-1))*TMPZ1
     2                   +(WW1(I,J,K+1)-2.0E00*WW1(I,J,K)+WW1(I,J,K-1))
     3                              *TMPZ2
9145         CONTINUE
             DO 9150 K=ILAP/2+1,NZ-ILAP/2
                TMP=DV*ORE*WWZ(K)
             DO 9150 J=2,NY-IY+1
             DO 9150 I=2,NX-IX+1
                FR(I,J,K)=FR(I,J,K)+WW2(I,J,K)*TMP
9150         CONTINUE
        ENDIF
        ENDIF
C----------------------------------------------------------------------
C  Calculate forces.  Buoyancy ...
C----------------------------------------------------------------------
c	FW=GRAV*RO
	DO 80 K=ILAP/2+1,NZ-ILAP/2
	DO 80 J=2,NY-IY+1
	DO 80 I=2,NX-IX+1
		FW(I,J,K)=GRAV*RO(I,J,K)
80	CONTINUE
C----------------------------------------------------------------------
C  Calculate pressure and 1/rho.
C----------------------------------------------------------------------
c	WW1=RO*TT
c	RO=1.0E00/RO
	DO 90 K=1,NZ
        DO 90 J=1,NY
        DO 90 I=1,NX
                WW1(I,J,K)=RO(I,J,K)*TT(I,J,K)
		RO(I,J,K)=1.0E00/RO(I,J,K)
90      CONTINUE
C----------------------------------------------------------------------
C  Grad P ...
C----------------------------------------------------------------------
c       FW=FW-(CSHIFT(WW1,1,3)-CSHIFT(WW1,-1,3))*HZ*DZZDZ
        DO 120 K=ILAP/2+1,NZ-ILAP/2
		TMPZ=HZ*DZZDZ(K)
        DO 120 J=2,NY-IY+1
		TMPY=HY*DYYDY(J)
        DO 120 I=2,NX-IX+1
		FW(I,J,K)=FW(I,J,K)-(WW1(I,J,K+1)-WW1(I,J,K-1))*TMPZ
		FV(I,J,K)=(WW1(I,J-1,K)-WW1(I,J+1,K))*TMPY
		FU(I,J,K)=(WW1(I-1,J,K)-WW1(I+1,J,K))*HX*DXXDX(I)
120     CONTINUE
C----------------------------------------------------------------------
C  Calculate velocities.
C----------------------------------------------------------------------
c	UU=RU*RO
c	VV=RV*RO
c	WW=RW*RO
        DO 160 K=1,NZ
        DO 160 J=1,NY
        DO 160 I=1,NX
                UU(I,J,K)=RU(I,J,K)*RO(I,J,K)
160      CONTINUE
        DO 170 K=1,NZ
        DO 170 J=1,NY
        DO 170 I=1,NX
                VV(I,J,K)=RV(I,J,K)*RO(I,J,K)
170      CONTINUE
        DO 180 K=1,NZ
        DO 180 J=1,NY
        DO 180 I=1,NX
                WW(I,J,K)=RW(I,J,K)*RO(I,J,K)
180      CONTINUE
C----------------------------------------------------------------------
C  Advection ...
C----------------------------------------------------------------------
c       WW1=RU*UU
c       FU=FU-(CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
c	WW1=RU*VV
c       FU=FU-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
c	WW1=RU*WW
c       FU=FU-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HZ*DZZDZ
        DO 210 K=ILAP/2+1,NZ-ILAP/2
		TMPZ=HZ*DZZDZ(K)
        DO 210 J=2,NY-IY+1
		TMPY=HY*DYYDY(J)
        DO 210 I=2,NX-IX+1
                FU(I,J,K)=FU(I,J,K)-(RU(I+1,J,K)*UU(I+1,J,K)
     2				    -RU(I-1,J,K)*UU(I-1,J,K))
     3							*HX*DXXDX(I)
		FU(I,J,K)=FU(I,J,K)-(RU(I,J+1,K)*VV(I,J+1,K)
     2				    -RU(I,J-1,K)*VV(I,J-1,K))*TMPY
		FU(I,J,K)=FU(I,J,K)-(RU(I,J,K+1)*WW(I,J,K+1)
     2                              -RU(I,J,K-1)*WW(I,J,K-1))*TMPZ
210      CONTINUE
c       FV=FV-(CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
c	WW1=RV*VV
c	FV=FV-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
c	WW1=RV*WW
c       FV=FV-(CSHIFT(WW1,1,3)-CSHIFT(WW1,-1,3))*HZ*DZZDZ
	DO 280 K=ILAP/2+1,NZ-ILAP/2
	 	TMPZ=HZ*DZZDZ(K)
	DO 280 J=2,NY-IY+1
		TMPY=HY*DYYDY(J)
	DO 280 I=2,NX-IX+1
		FV(I,J,K)=FV(I,J,K)-(RV(I+1,J,K)*UU(I+1,J,K)
     2                              -RV(I-1,J,K)*UU(I-1,J,K))
     3                                                  *HX*DXXDX(I)
                FV(I,J,K)=FV(I,J,K)-(RV(I,J+1,K)*VV(I,J+1,K)
     2                              -RV(I,J-1,K)*VV(I,J-1,K))*TMPY
                FV(I,J,K)=FV(I,J,K)-(RV(I,J,K+1)*WW(I,J,K+1)
     2                              -RV(I,J,K-1)*WW(I,J,K-1))*TMPZ
280      CONTINUE
c       FW=FW-(CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
c 	FW=FW-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
c	WW1=RW*WW
c       FW=FW-(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HZ*DZZDZ
        DO 330 K=ILAP/2+1,NZ-ILAP/2
		TMPZ=HZ*DZZDZ(K)
        DO 330 J=2,NY-IY+1
		TMPY=HY*DYYDY(J)
        DO 330 I=2,NX-IX+1
		FW(I,J,K)=FW(I,J,K)-(RW(I+1,J,K)*UU(I+1,J,K)
     2                              -RW(I-1,J,K)*UU(I-1,J,K))
     3                                                  *HX*DXXDX(I)
                FW(I,J,K)=FW(I,J,K)-(RW(I,J+1,K)*VV(I,J+1,K)
     2                              -RW(I,J-1,K)*VV(I,J-1,K))*TMPY
                FW(I,J,K)=FW(I,J,K)-(RW(I,J,K+1)*WW(I,J,K+1)
     2                              -RW(I,J,K-1)*WW(I,J,K-1))*TMPZ
330      CONTINUE
C----------------------------------------------------------------------
C  Calculate energy fluxes.  Advection ...
C----------------------------------------------------------------------
c       WW1=(CSHIFT(TT,1,3)-CSHIFT(TT,-1,3))*HZ*DZZDZ
c       WW2=(CSHIFT(TT,1,1)-CSHIFT(TT,-1,1))*HX*DXXDX
c	WW3=(CSHIFT(TT,1,2)-CSHIFT(TT,-1,2))*HY*DYYDY
c	FT=-WW*WW1
c	FT=FT-UU*WW2
c	FT=FT-VV*WW3
C
	IF (.NOT.LSHR) THEN
C
        DO 450 K=ILAP/2+1,NZ-ILAP/2
		TMPZ=HZ*DZZDZ(K)
        DO 450 J=2,NY-IY+1
		TMPY=HY*DYYDY(J)
        DO 450 I=2,NX-IX+1
		FT(I,J,K)=-1.0E00*UU(I,J,K)*(TT(I+1,J,K)-TT(I-1,J,K))
     2							 *HX*DXXDX(I)
		FT(I,J,K)=FT(I,J,K)-VV(I,J,K)*(TT(I,J+1,K)-TT(I,J-1,K))
     2								  *TMPY
		FT(I,J,K)=FT(I,J,K)-WW(I,J,K)*(TT(I,J,K+1)-TT(I,J,K-1))
     2								  *TMPZ
450      CONTINUE
C----------------------------------------------------------------------
C  Diffusion ...
C----------------------------------------------------------------------
c       WW3=WW3*(1.0E00/DYYDY)*D2YYDY2
c     2    +(CSHIFT(TT,1,2)-2.0E00*TT+CSHIFT(TT,-1,2))*H2Y*DYYDY*DYYDY
c     	WW2=WW2*(1.0E00/DXXDX)*D2XXDX2
c     2	   +(CSHIFT(TT,1,1)-2.0E00*TT+CSHIFT(TT,-1,1))*H2X*DXXDX*DXXDX
c     	WW1=WW1*(1.0E00/DZZDZ)*D2ZZDZ2
c     2	   +(CSHIFT(TT,1,3)-2.0E00*TT+CSHIFT(TT,-1,3))*H2Z*DZZDZ*DZZDZ
c	WW1=WW1+WW2+WW3
c	FT=FT+WW1*RO*OCV*(1.0E00/REPR)
C
	IF ((RLAX.GT.0.0E00).OR.LREM) THEN
C----------------------------------------------------------------------
C  Calculate horizontal mean temperature and temperature perturbations.
C----------------------------------------------------------------------
		CALL HORIZONTAL_MEAN(TTM,TT)
C
           DO K=1,NZ
              DO J=1,NY
                 DO I=1,NX
                    WW1(I,J,K)=TT(I,J,K)-TTM(K)
                 END DO
              END DO
           END DO
	ENDIF
C
	ENDIF
C
	IF (LREM) THEN
C----------------------------------------------------------------------
C  Diffuse only temeprature perturbations (M. Rempel)
C----------------------------------------------------------------------
		TMP=OCV/REPR
                DO 518 K=ILAP/2+1,NZ-ILAP/2
                        TMPZ1=HZ*D2ZZDZ2(K)
                        TMPZ2=H2Z*DZZDZ(K)*DZZDZ(K)
                DO 518 J=2,NY-IY+1
                        TMPY1=HY*D2YYDY2(J)
                        TMPY2=H2Y*DYYDY(J)*DYYDY(J)
                DO 518 I=2,NX-IX+1
                        FT(I,J,K)=FT(I,J,K)
     2                          +((WW1(I+1,J,K)-WW1(I-1,J,K))
     3                                          *HX*D2XXDX2(I)
     4                   +(WW1(I+1,J,K)-2.0E00*WW1(I,J,K)+WW1(I-1,J,K))
     5                                          *H2X*DXXDX(I)*DXXDX(I))
     6                                          	*RO(I,J,K)*TMP
                        FT(I,J,K)=FT(I,J,K)
     2                     +((WW1(I,J+1,K)-WW1(I,J-1,K))*TMPY1
     3                   +(WW1(I,J+1,K)-2.0E00*WW1(I,J,K)+WW1(I,J-1,K))
     4                                          *TMPY2)
     5                                          	*RO(I,J,K)*TMP
                        FT(I,J,K)=FT(I,J,K)
     2                     +((WW1(I,J,K+1)-WW1(I,J,K-1))*TMPZ1
     3                   +(WW1(I,J,K+1)-2.0E00*WW1(I,J,K)+WW1(I,J,K-1))
     4                                          *TMPZ2)
     5                                          	*RO(I,J,K)*TMP
518             CONTINUE
C----------------------------------------------------------------------
C  Radiative heating (M. Rempel)
C----------------------------------------------------------------------
        	DO K=ILAP/2+1,NZ-ILAP/2
           		HRAD(K)=((TTM(K+1)-2.0E00*TTM(K)+TTM(K-1))
     2						*H2Z*DZZDZ(K)*DZZDZ(K)
     3        		       +(TTM(K+1)-TTM(K-1))
     4					      *HZ*D2ZZDZ2(K))*RKAPA(K)
     5        		       +(TTM(K+1)-TTM(K-1))*HZ*DZZDZ(K)*DKAPA(K)
        	END DO
		DO K=ILAP/2+1,NZ-ILAP/2
                DO J=2,NY-IY+1
                DO I=2,NX-IX+1
			FT(I,J,K)=FT(I,J,K)+RO(I,J,K)*TMP*HRAD(K)
		END DO
		END DO
		END DO
C----------------------------------------------------------------------
	ELSE
		TMP=OCV/REPR
		DO 519 K=ILAP/2+1,NZ-ILAP/2
			TMPZ1=HZ*D2ZZDZ2(K)
			TMPZ2=H2Z*DZZDZ(K)*DZZDZ(K)
       		DO 519 J=2,NY-IY+1
			TMPY1=HY*D2YYDY2(J)
			TMPY2=H2Y*DYYDY(J)*DYYDY(J)
        	DO 519 I=2,NX-IX+1
                	FT(I,J,K)=FT(I,J,K)
     2			  	+((TT(I+1,J,K)-TT(I-1,J,K))
     3						*HX*D2XXDX2(I)
     4			    +(TT(I+1,J,K)-2.0E00*TT(I,J,K)+TT(I-1,J,K))
     5						*H2X*DXXDX(I)*DXXDX(I))
     6						*RO(I,J,K)*TMP*RKAPA(K)
     			FT(I,J,K)=FT(I,J,K)
     2                     +((TT(I,J+1,K)-TT(I,J-1,K))*TMPY1
     3                      +(TT(I,J+1,K)-2.0E00*TT(I,J,K)+TT(I,J-1,K))
     4                                          *TMPY2)
     5                                          *RO(I,J,K)*TMP*RKAPA(K)
     			FT(I,J,K)=FT(I,J,K)
     2                     +((TT(I,J,K+1)-TT(I,J,K-1))*TMPZ1
     3                      +(TT(I,J,K+1)-2.0E00*TT(I,J,K)+TT(I,J,K-1))
     4                                          *TMPZ2)
     5                                          *RO(I,J,K)*TMP*RKAPA(K)
     			FT(I,J,K)=FT(I,J,K)+(TT(I,J,K+1)-TT(I,J,K-1))
     2				    *HZ*DZZDZ(K)*RO(I,J,K)*TMP*DKAPA(K)
519      	CONTINUE
	ENDIF
C----------------------------------------------------------------------
C  Flux relaxation (M. Rempel).
C----------------------------------------------------------------------
	IF (RLAX.GT.0.0E00) THEN
C
	   IF ((TIMT.GT.TSTART).AND.(TFAC.GT.1D-4)) THEN
C
	      IF((ITC.EQ.0).AND.(IZC.EQ.1)) THEN
		 IF (MYPE.EQ.NPE-1) THEN
		    IF(TTM(NZ-ILAP/2).LE.TTM_MAX) THEN
		       TFAC    = TFAC-DT/(3.0*TOFF)
		    ELSE
		       TTM_MAX = TTM(NZ-ILAP/2)
                       TFAC    = TFAC+DT/(3.0*TOFF)
                       TFAC    = MIN(TFAC,2.0)
		    ENDIF   
		 ENDIF
		 CALL MPI_BCAST(TFAC,1,MPISIZE,NPE-1,
     2						 MPI_COMM_WORLD,IERR)
	      ELSE IF((ITC.EQ.1).AND.(IZC.EQ.0)) THEN
                 IF (MYPE.EQ.0) THEN
                    IF(TTM(ILAP/2+1).GE.TTM_MIN) THEN
                       TFAC    = TFAC-DT/(3.0*TOFF)
                    ELSE
                       TTM_MIN = TTM(ILAP/2+1)
                       TFAC    = TFAC+DT/(3.0*TOFF)
                       TFAC    = MIN(TFAC,2.0)
                    ENDIF
                 ENDIF
		 CALL MPI_BCAST(TFAC,1,MPISIZE,0,MPI_COMM_WORLD,IERR)
	      ELSE
		 WRITE(6,*)'Relax: check boundary temperature'
		 CALL MPI_FINALIZE(IERR)
                 STOP
              ENDIF
C----------------------------------------------------------------------
C  Compute mean convective and radiative flux.
C----------------------------------------------------------------------
              DO K=ILAP/2+1,NZ
                 DO J=1+IY/2,NY-IY+1
                    DO I=1+IX/2,NX-IX/2
                       WW2(I,J,K)=-RW(I,J,K)*(2.5*WW1(I,J,K)
     2                                +0.5*(UU(I,J,K)*UU(I,J,K)
     3                                     +VV(I,J,K)*VV(I,J,K)
     4                                     +WW(I,J,K)*WW(I,J,K)))
                    END DO
                 END DO
              END DO
C    
              CALL HORIZONTAL_MEAN(FCONM,WW2)
C    
              DO K=ILAP/2+1,NZ-ILAP/2
                 FRADM(K)=(TTM(K+1)-TTM(K-1))
     2                *HZ*DZZDZ(K)*RKAPA(K)/REPR
              END DO
              IF (MYPEZ.EQ.0) THEN
                 FRADM(ILAP/2+1)=(-3.0*TTM(ILAP/2+1)+4.0*TTM(ILAP/2+2)
     2                -TTM(ILAP/2+3))*HZ*DZZDZ(ILAP/2+1)
     3                *RKAPA(ILAP/2+1)/REPR
              ENDIF
              IF (MYPEZ.EQ.NPEZ-1) THEN
                 FRADM(NZ-ILAP/2)=(3.0*TTM(NZ-ILAP/2)
     2                -4.0*TTM(NZ-ILAP/2-1)
     3                +TTM(NZ-ILAP/2-2))*HZ*DZZDZ(NZ-ILAP/2)
     4                *RKAPA(NZ-ILAP/2)/REPR
              ENDIF
C----------------------------------------------------------------------
C  Communication to get FRADM(NZ,MYPE)=FRADM(ILAP/2+1,MYPE+NPEY)
C----------------------------------------------------------------------
              ITAG = 100
              IF (MYPEZ.EQ.0) THEN
                 CALL MPI_RECV(FRADM(NZ),1,MPISIZE,MYPE+NPEY,ITAG,
     2                    MPI_COMM_WORLD,ISTATUS,IERR)
              ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
                 CALL MPI_SEND(FRADM(ILAP/2+1),1,MPISIZE,
     2                    MYPE-NPEY,ITAG,MPI_COMM_WORLD,IERR)
              ELSE
C                print*, "pt 1 mype, npey", mype, npey
                 CALL MPI_SENDRECV(FRADM(ILAP/2+1),1,MPISIZE,
     2                    MYPE-NPEY,ITAG,FRADM(NZ),1,MPISIZE,
     3                    MYPE+NPEY,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
              ENDIF
C
              IF (LREM) THEN
                 FTOT=8.07*0.4*THETA
              ELSE
                 FTOT=THETA/REPR
              ENDIF
C
              DO K=ILAP/2+1,NZ
                 ALPHA(K) = 1.0-(FCONM(K)+FRADM(K))/FTOT
              END DO
              DO K=ILAP/2+1,NZ
                 ALPHA(K)=ALPHA(K)/(1.0+4.0*ABS(ALPHA(K)))
              END DO
              CORRECT(ILAP/2+1)=0.0
              DO K=ILAP/2+2,NZ
                 CORRECT(K)=CORRECT(K-1)+0.5*(ALPHA(K)+ALPHA(K-1))
     2                     *(ZEE(K)-ZEE(K-1))
              END DO
C
              IF(MYPEY.EQ.0) THEN
                 ENDVAL(MYPEZ+1)=CORRECT(NZ)
C
                 ITAG = 200
                 IF (MYPEZ.GT.0) THEN
                    CALL MPI_SEND(ENDVAL(MYPEZ+1),1,MPISIZE,0,
     2                           ITAG,MPI_COMM_WORLD,IERR)
                 ENDIF
                 IF (MYPEZ.EQ.0) THEN
                    ADDSUM(1)=ENDVAL(1)
                    DO IPE=1,NPEZ-1
                       CALL MPI_RECV(ADDSUM(IPE+1),1,MPISIZE,IPE*NPEY,
     2                              ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
                    END DO
                 ENDIF
              ENDIF
              CALL MPI_BCAST(ADDSUM,NPEZ,MPISIZE,0,MPI_COMM_WORLD,IERR)
C
              SUM=0.0
              IF((ITC.EQ.0).AND.(IZC.EQ.1)) THEN
                 IF(MYPEZ.GT.0) THEN
                    DO IPE=1,MYPEZ
                       SUM=SUM+ADDSUM(IPE)
                    END DO
                 END IF
              ELSE IF((ITC.EQ.1).AND.(IZC.EQ.0)) THEN
                 DO IPE=NPEZ-1,MYPEZ,-1
                    SUM=SUM-ADDSUM(IPE+1)
                    SUM=SUM-ADDSUM(IPE+1)
                 END DO
              ELSE
                 WRITE(6,*)'Relax: check boundary temperature.'
                 CALL MPI_FINALIZE(IERR)
                 STOP
              END IF
C
              DO K=ILAP/2+1,NZ-ILAP/2
                 CORRECT(K)=CORRECT(K)+SUM
              END DO
           ELSE
             DO K=ILAP/2+1,NZ-ILAP/2
                CORRECT(K) = 0.0
             END DO
          ENDIF
C
          DO K=ILAP/2+1,NZ-ILAP/2
             DO J=1+IY/2,NY-IY+1
                DO I=1+IX/2,NX-IX/2
                   FT(I,J,K)=FT(I,J,K)+RLAX*TFAC*CORRECT(K)
                END DO
             END DO
          END DO
C----------------------------------------------------------------------
C  Write progress of relaxation in output file every 25 time steps.
C----------------------------------------------------------------------
          IF((MOD(TIMT,25.0).LE.DT).AND.(MOD(ICALL,3).EQ.0)) THEN
             IF((ITC.EQ.0).AND.(IZC.EQ.1).AND.(MYPE.EQ.NPE-1)) THEN
                WRITE(*,'(A7,5(D15.6))') 'Relax: ',TIMT,TFAC,
     &               CORRECT(NZ-ILAP/2),TTM(NZ-ILAP/2),TTM_MAX
             ENDIF
             IF((ITC.EQ.1).AND.(IZC.EQ.0).AND.(MYPE.EQ.0)) THEN
                WRITE(*,'(A7,5(D15.6))') 'Relax: ',TIMT,TFAC,
     &               CORRECT(ILAP/2+1),TTM(ILAP/2+1),TTM_MIN
             ENDIF
          ENDIF
          ICALL=ICALL+1
        ENDIF
C----------------------------------------------------------------------
C  Localized embedded heat loss, ramps up to a constant value TP.
C----------------------------------------------------------------------
        IF (ITC.EQ.2) THEN
            CLN=-4.0E00*LOG(2.0E00)/HH/HH
            IF (TC.NE.0.0E00) THEN
                DO 521 K=ILAP/2+1,NZ-ILAP/2
                DO 521 J=2,NY-IY+1
                DO 521 I=2,NX-IX+1
                    FT(I,J,K)=FT(I,J,K)
     2                 -TP*OCV*0.5E00*(1.0E00+TANH(LOG(3.0E00)/QFH
     2                                          *(TIMT-TC)))*RO(I,J,K)
     3                         		*EXP(CLN*(EXX(I)-XP)**2)/HH
     4                         		*EXP(CLN*(WYY(J)-YP)**2)/HH
     5                         		*EXP(CLN*(ZEE(K)-ZP)**2)/HH
521             CONTINUE
            ELSE
                DO 522 K=ILAP/2+1,NZ-ILAP/2
                DO 522 J=2,NY-IY+1
                DO 522 I=2,NX-IX+1
                        FT(I,J,K)=FT(I,J,K)-TP*OCV*RO(I,J,K)
     2                                  *EXP(CLN*(EXX(I)-XP)**2)/HH
     3                                  *EXP(CLN*(WYY(J)-YP)**2)/HH
     4                                  *EXP(CLN*(ZEE(K)-ZP)**2)/HH
522             CONTINUE
            ENDIF
        ENDIF
C----------------------------------------------------------------------
C  Calculate the viscous terms and add them to FU, FV, and FW.
C----------------------------------------------------------------------
c	WW1=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*D2XXDX2
c     2	   +(CSHIFT(UU,1,1)-2.0E00*UU+CSHIFT(UU,-1,1))*H2X*DXXDX*DXXDX
c	WW1=C43*WW1
c	WW1=WW1+(CSHIFT(UU,1,2)-CSHIFT(UU,-1,2))*HY*D2YYDY2
c     2    +(CSHIFT(UU,1,2)-2.0E00*UU+CSHIFT(UU,-1,2))*H2Y*DYYDY*DYYDY
c	WW1=WW1+(CSHIFT(UU,1,3)-CSHIFT(UU,-1,3))*HZ*D2ZZDZ2
c     2     +(CSHIFT(UU,1,3)-2.0E00*UU+CSHIFT(UU,-1,3))*H2Z*DZZDZ*DZZDZ
c	FU=FU+WW1*ORE
        DO 530 K=ILAP/2+1,NZ-ILAP/2
		TMPZ1=HZ*D2ZZDZ2(K)
		TMPZ2=H2Z*DZZDZ(K)*DZZDZ(K)
        DO 530 J=2,NY-IY+1
		TMPY1=HY*D2YYDY2(J)
		TMPY2=H2Y*DYYDY(J)*DYYDY(J)
        DO 530 I=2,NX-IX+1
		FU(I,J,K)=FU(I,J,K)+ORE*(
     2		  	C43*((UU(I+1,J,K)-UU(I-1,J,K))*HX*D2XXDX2(I)
     3		  	 +(UU(I+1,J,K)-2.0E00*UU(I,J,K)+UU(I-1,J,K))
     4		 		      	      *H2X*DXXDX(I)*DXXDX(I))
     5		 		     +(UU(I,J+1,K)-UU(I,J-1,K))*TMPY1
     6		 	  +(UU(I,J+1,K)-2.0E00*UU(I,J,K)+UU(I,J-1,K))
     7						       	       *TMPY2
     8		 		     +(UU(I,J,K+1)-UU(I,J,K-1))*TMPZ1
     9		          +(UU(I,J,K+1)-2.0E00*UU(I,J,K)+UU(I,J,K-1))
     1						      	       *TMPZ2)
530	CONTINUE
c	WW1=(CSHIFT(VV,1,1)-CSHIFT(VV,-1,1))*HX*D2XXDX2
c     2     +(CSHIFT(VV,1,1)-2.0E00*VV+CSHIFT(VV,-1,1))*H2X*DXXDX*DXXDX
c       WW1=WW1+C43*(CSHIFT(VV,1,2)-CSHIFT(VV,-1,2))*HY*D2YYDY2
c     2    +(CSHIFT(VV,1,2)-2.0E00*VV+CSHIFT(VV,-1,2))*H2Y*DYYDY*DYYDY)
c	WW1=WW1+(CSHIFT(VV,1,3)-CSHIFT(VV,-1,3))*HZ*D2ZZDZ2
c     2     +(CSHIFT(VV,1,3)-2.0E00*VV+CSHIFT(VV,-1,3))*H2Z*DZZDZ*DZZDZ
c	FV=FV+WW1*ORE
        DO 540 K=ILAP/2+1,NZ-ILAP/2
                TMPZ1=HZ*D2ZZDZ2(K)
                TMPZ2=H2Z*DZZDZ(K)*DZZDZ(K)
        DO 540 J=2,NY-IY+1
                TMPY1=HY*D2YYDY2(J)
                TMPY2=H2Y*DYYDY(J)*DYYDY(J)
        DO 540 I=2,NX-IX+1
                FV(I,J,K)=FV(I,J,K)+ORE*(
     2			  (VV(I+1,J,K)-VV(I-1,J,K))*HX*D2XXDX2(I)
     3                    +(VV(I+1,J,K)-2.0E00*VV(I,J,K)+VV(I-1,J,K))
     4                                         *H2X*DXXDX(I)*DXXDX(I)
     5                   +C43*((VV(I,J+1,K)-VV(I,J-1,K))*TMPY1
     6                    +(VV(I,J+1,K)-2.0E00*VV(I,J,K)+VV(I,J-1,K))
     7                                                        *TMPY2)
     8                   +(VV(I,J,K+1)-VV(I,J,K-1))*TMPZ1
     9                    +(VV(I,J,K+1)-2.0E00*VV(I,J,K)+VV(I,J,K-1))
     1                                                         *TMPZ2)
540     CONTINUE
c	WW1=(CSHIFT(WW,1,3)-CSHIFT(WW,-1,3))*HZ*D2ZZDZ2
c     2     +(CSHIFT(WW,1,3)-2.0E00*WW+CSHIFT(WW,-1,3))*H2Z*DZZDZ*DZZDZ
c	WW1=C43*WW1
c	WW1=WW1+(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*D2XXDX2
c     2     +(CSHIFT(WW,1,1)-2.0E00*WW+CSHIFT(WW,-1,1))*H2X*DXXDX*DXXDX
c       WW1=WW1+(CSHIFT(WW,1,2)-CSHIFT(WW,-1,2))*HY*D2YYDY2
c     2     +(CSHIFT(WW,1,2)-2.0E00*WW+CSHIFT(WW,-1,2))*H2Y*DYYDY*DYYDY
c	FW=FW+WW1*ORE
	DO 550 K=ILAP/2+1,NZ-ILAP/2
                TMPZ1=HZ*D2ZZDZ2(K)
                TMPZ2=H2Z*DZZDZ(K)*DZZDZ(K)
        DO 550 J=2,NY-IY+1
                TMPY1=HY*D2YYDY2(J)
                TMPY2=H2Y*DYYDY(J)*DYYDY(J)
        DO 550 I=2,NX-IX+1
                FW(I,J,K)=FW(I,J,K)+ORE*(
     2                    (WW(I+1,J,K)-WW(I-1,J,K))*HX*D2XXDX2(I)
     3                    +(WW(I+1,J,K)-2.0E00*WW(I,J,K)+WW(I-1,J,K))
     4                                         *H2X*DXXDX(I)*DXXDX(I)
     5                    +(WW(I,J+1,K)-WW(I,J-1,K))*TMPY1
     6                    +(WW(I,J+1,K)-2.0E00*WW(I,J,K)+WW(I,J-1,K))
     7                                                        *TMPY2
     8                   +C43*((WW(I,J,K+1)-WW(I,J,K-1))*TMPZ1
     9                    +(WW(I,J,K+1)-2.0E00*WW(I,J,K)+WW(I,J,K-1))
     1                                                        *TMPZ2))
550     CONTINUE
C----------------------------------------------------------------------
C  Viscous heating, compressional heating, and viscous dissipation ...
C----------------------------------------------------------------------
c       WW1=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*DXXDX
	DO 650 K=ILAP/2+1,NZ-ILAP/2
        DO 650 J=1,NY
        DO 650 I=2,NX-IX+1
                WW1(I,J,K)=(UU(I+1,J,K)-UU(I-1,J,K))*HX*DXXDX(I)
650     CONTINUE
c       WW2=(CSHIFT(VV,1,1)-CSHIFT(VV,-1,1))*HX*DXXDX
        DO 670 K=ILAP/2+1,NZ-ILAP/2
        DO 670 J=1,NY
        DO 670 I=2,NX-IX+1
                WW2(I,J,K)=(VV(I+1,J,K)-VV(I-1,J,K))*HX*DXXDX(I)
670     CONTINUE
C
	IF (.NOT.LSHR) THEN
C
c       WW3=(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*DXXDX
        DO 690 K=ILAP/2+1,NZ-ILAP/2
        DO 690 J=2,NY-IY+1
        DO 690 I=2,NX-IX+1
                WW3(I,J,K)=(WW(I+1,J,K)-WW(I-1,J,K))*HX*DXXDX(I)
690     CONTINUE
c       FT=FT+(C43*WW1*WW1+WW2*WW2+WW3*WW3)*RO*OCV*ORE
        TMP=OCV*ORE
        DO 700 K=ILAP/2+1,NZ-ILAP/2
        DO 700 J=2,NY-IY+1
        DO 700 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)+(C43*WW1(I,J,K)*WW1(I,J,K)
     2                          	+WW2(I,J,K)*WW2(I,J,K)
     3                              	+WW3(I,J,K)*WW3(I,J,K))
     4                                           *RO(I,J,K)*TMP
700     CONTINUE
c	FT=FT-TT*OCV*WW1
	DO 710 K=ILAP/2+1,NZ-ILAP/2
        DO 710 J=2,NY-IY+1
        DO 710 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)-OCV*WW1(I,J,K)*TT(I,J,K)
710     CONTINUE
C
	ENDIF
C
c	FU=FU+C13*ORE*(CSHIFT(WW2,1,2)-CSHIFT(WW2,-1,2))*HY*DYYDY
	TMP=C13*ORE
	DO 720 K=ILAP/2+1,NZ-ILAP/2
	DO 720 J=2,NY-IY+1
	DO 720 I=2,NX-IX+1
		FU(I,J,K)=FU(I,J,K)+TMP*(WW2(I,J+1,K)-WW2(I,J-1,K))
     2							*HY*DYYDY(J)
720	CONTINUE
c       FV=FV+C13*ORE*(CSHIFT(WW1,1,2)-CSHIFT(WW1,-1,2))*HY*DYYDY
	DO 730 K=ILAP/2+1,NZ-ILAP/2
        DO 730 J=2,NY-IY+1
        DO 730 I=2,NX-IX+1
                FV(I,J,K)=FV(I,J,K)+TMP*(WW1(I,J+1,K)-WW1(I,J-1,K))
     2                                                  *HY*DYYDY(J)
730     CONTINUE
C
	IF (.NOT.LSHR) THEN
C
c       WW1=(CSHIFT(UU,1,2)-CSHIFT(UU,-1,2))*HY*DYYDY
        DO 740 K=ILAP/2+1,NZ-ILAP/2
        DO 740 J=2,NY-IY+1
        DO 740 I=2,NX-IX+1
                WW1(I,J,K)=(UU(I,J+1,K)-UU(I,J-1,K))*HY*DYYDY(J)
740     CONTINUE
c       FT=FT+(WW1*WW1+2.0E00*WW1*WW2)*RO*OCV*ORE
        TMP=OCV*ORE
        DO 750 K=ILAP/2+1,NZ-ILAP/2
        DO 750 J=2,NY-IY+1
        DO 750 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)+(WW1(I,J,K)*WW1(I,J,K)
     2                              +2.0E00*WW1(I,J,K)*WW2(I,J,K))
     3                                         	    *RO(I,J,K)*TMP
750     CONTINUE
c       WW2=(CSHIFT(VV,1,2)-CSHIFT(VV,-1,2))*HY*DYYDY
        DO 755 K=ILAP/2+1,NZ-ILAP/2
        DO 755 J=2,NY-IY+1
        DO 755 I=2,NX-IX+1
                WW2(I,J,K)=(VV(I,J+1,K)-VV(I,J-1,K))*HY*DYYDY(J)
755     CONTINUE
c       FT=FT+C43*WW2*WW2*RO*OCV*ORE
        TMP=C43*OCV*ORE
        DO 760 K=ILAP/2+1,NZ-ILAP/2
        DO 760 J=2,NY-IY+1
        DO 760 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)+WW2(I,J,K)*WW2(I,J,K)
     2                                          *RO(I,J,K)*TMP
760    CONTINUE
c       FT=FT-TT*OCV*WW2
        DO 770 K=ILAP/2+1,NZ-ILAP/2
        DO 770 J=2,NY-IY+1
        DO 770 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)-OCV*WW2(I,J,K)*TT(I,J,K)
770     CONTINUE
c       WW1=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*DXXDX
        DO 780 K=ILAP/2+1,NZ-ILAP/2
        DO 780 J=2,NY-IY+1
        DO 780 I=2,NX-IX+1
                WW1(I,J,K)=(UU(I+1,J,K)-UU(I-1,J,K))*HX*DXXDX(I)
780     CONTINUE
C
	ENDIF
C
c       WW3=(CSHIFT(WW,1,3)-CSHIFT(WW,-1,3))*HZ*DZZDZ
        DO 790 K=ILAP/2+1,NZ-ILAP/2
        DO 790 J=1,NY
        DO 790 I=1,NX
                WW3(I,J,K)=(WW(I,J,K+1)-WW(I,J,K-1))*HZ*DZZDZ(K)
790     CONTINUE
C
	IF (.NOT.LSHR) THEN
C
c       FT=FT+C43*(WW3*WW3-WW1*WW3-WW1*WW2-WW2*WW3)*RO*OCV*ORE
        TMP=C43*OCV*ORE
	DO 820 K=ILAP/2+1,NZ-ILAP/2
        DO 820 J=2,NY-IY+1
        DO 820 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)+(WW3(I,J,K)*WW3(I,J,K)
     2                              -WW1(I,J,K)*WW3(I,J,K)
     3                              -WW1(I,J,K)*WW2(I,J,K)
     4                              -WW2(I,J,K)*WW3(I,J,K))
     5                                       *RO(I,J,K)*TMP
820     CONTINUE
c       FT=FT-TT*OCV*WW3
        DO 830 K=ILAP/2+1,NZ-ILAP/2
        DO 830 J=2,NY-IY+1
        DO 830 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)-OCV*WW3(I,J,K)*TT(I,J,K)
830     CONTINUE
C
	ENDIF
C
c       FU=FU+C13*ORE*(CSHIFT(WW3,1,1)-CSHIFT(WW3,-1,1))*HX*DXXDX
        TMP=C13*ORE
        DO 840 K=ILAP/2+1,NZ-ILAP/2
        DO 840 J=2,NY-IY+1
        DO 840 I=2,NX-IX+1
                FU(I,J,K)=FU(I,J,K)+TMP*(WW3(I+1,J,K)-WW3(I-1,J,K))
     2                                                 *HX*DXXDX(I)
840     CONTINUE
c       FV=FV+C13*ORE*(CSHIFT(WW3,1,2)-CSHIFT(WW3,-1,2))*HY*DYYDY
        TMP=C13*ORE
        DO 850 K=ILAP/2+1,NZ-ILAP/2
        DO 850 J=2,NY-IY+1
        DO 850 I=2,NX-IX+1
                FV(I,J,K)=FV(I,J,K)+TMP*(WW3(I,J+1,K)-WW3(I,J-1,K))
     2                                                 *HY*DYYDY(J)
850     CONTINUE
c       WW1=(CSHIFT(UU,1,3)-CSHIFT(UU,-1,3))*HZ*DZZDZ
        DO 860 K=ILAP/2+1,NZ-ILAP/2
        DO 860 J=2,NY-IY+1
        DO 860 I=1,NX
                WW1(I,J,K)=(UU(I,J,K+1)-UU(I,J,K-1))*HZ*DZZDZ(K)
860     CONTINUE
c       WW2=(CSHIFT(VV,1,3)-CSHIFT(VV,-1,3))*HZ*DZZDZ
        DO 880 K=ILAP/2+1,NZ-ILAP/2
        DO 880 J=1,NY
        DO 880 I=2,NX-IX+1
                WW2(I,J,K)=(VV(I,J,K+1)-VV(I,J,K-1))*HZ*DZZDZ(K)
880     CONTINUE
C
	IF (.NOT.LSHR) THEN
C
c       WW3=(CSHIFT(WW,1,2)-CSHIFT(WW,-1,2))*HY*DYYDY
        DO 900 K=ILAP/2+1,NZ-ILAP/2
        DO 900 J=2,NY-IY+1
        DO 900 I=2,NX-IX+1
                WW3(I,J,K)=(WW(I,J+1,K)-WW(I,J-1,K))*HY*DYYDY(J)
900     CONTINUE
c       FT=FT+(WW1*WW1+WW2*WW2+2.0E00*WW2*WW3+WW3*WW3)*RO*OCV*ORE
        TMP=OCV*ORE
        DO 910 K=ILAP/2+1,NZ-ILAP/2
        DO 910 J=2,NY-IY+1
        DO 910 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)+(WW1(I,J,K)*WW1(I,J,K)
     2				    +WW2(I,J,K)*WW2(I,J,K)
     3				    +WW3(I,J,K)*WW3(I,J,K)
     4				    +2.0E00*WW2(I,J,K)*WW3(I,J,K))
     5                                              *RO(I,J,K)*TMP
910     CONTINUE
C
	ENDIF
C
c       FW=FW+C13*ORE*((CSHIFT(WW1,1,1)-CSHIFT(WW1,-1,1))*HX*DXXDX
c	     +(CSHIFT(WW2,1,2)-CSHIFT(WW2,-1,2))*HY*DYYDY)
	TMP=C13*ORE
        DO 920 K=ILAP/2+1,NZ-ILAP/2
        DO 920 J=2,NY-IY+1
        DO 920 I=2,NX-IX+1
                FW(I,J,K)=FW(I,J,K)+TMP*((WW1(I+1,J,K)-WW1(I-1,J,K))
     2                                                  *HX*DXXDX(I)
     3					+(WW2(I,J+1,K)-WW2(I,J-1,K))
     4                                                  *HY*DYYDY(J))
920     CONTINUE
C
	IF (.NOT.LSHR) THEN
C
c       WW3=(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*DXXDX
        DO 930 K=ILAP/2+1,NZ-ILAP/2
        DO 930 J=2,NY-IY+1
        DO 930 I=2,NX-IX+1
                WW3(I,J,K)=(WW(I+1,J,K)-WW(I-1,J,K))*HX*DXXDX(I)
930     CONTINUE
c       FT=FT+2.0E00*WW1*WW3*RO*OCV*ORE
        TMP=OCV*ORE
        DO 940 K=ILAP/2+1,NZ-ILAP/2
        DO 940 J=2,NY-IY+1
        DO 940 I=2,NX-IX+1
                FT(I,J,K)=FT(I,J,K)+2.0E00*WW1(I,J,K)*WW3(I,J,K)
     2                                            *RO(I,J,K)*TMP
940     CONTINUE
C
	ENDIF
C
C----------------------------------------------------------------------
C  Add rotation terms.
C----------------------------------------------------------------------
	IF (LROT) THEN
c		FU=FU+RV*OMZ
        	DO 950 K=ILAP/2+1,NZ-ILAP/2
        	DO 950 J=2,NY-IY+1
        	DO 950 I=2,NX-IX+1
                	FU(I,J,K)=FU(I,J,K)+OMZ*RV(I,J,K)
950     	CONTINUE
c		FV=FV-RU*OMZ+RW*OMX
		DO 960 K=ILAP/2+1,NZ-ILAP/2
        	DO 960 J=2,NY-IY+1
        	DO 960 I=2,NX-IX+1
                	FV(I,J,K)=FV(I,J,K)-OMZ*RU(I,J,K)+OMX*RW(I,J,K)
960     	CONTINUE
c        	FW=FW-RV*OMX
        	DO 970 K=ILAP/2+1,NZ-ILAP/2
        	DO 970 J=2,NY-IY+1
        	DO 970 I=2,NX-IX+1
                	FW(I,J,K)=FW(I,J,K)-OMX*RV(I,J,K)
970     	CONTINUE
	ENDIF
C----------------------------------------------------------------------
C  Magnetic fields ...
C----------------------------------------------------------------------
	IF (LMAG) THEN
c       	RU=(CSHIFT(UU,1,1)-CSHIFT(UU,-1,1))*HX*DXXDX
        	DO 1000 K=ILAP/2+1,NZ-ILAP/2
        	DO 1000 J=2,NY-IY+1
        	DO 1000 I=2,NX-IX+1
                	RU(I,J,K)=(UU(I+1,J,K)-UU(I-1,J,K))*HX*DXXDX(I)
1000    	CONTINUE
c		RV=(CSHIFT(VV,1,1)-CSHIFT(VV,-1,1))*HX*DXXDX
                DO 1010 K=ILAP/2+1,NZ-ILAP/2
                DO 1010 J=2,NY-IY+1
                DO 1010 I=2,NX-IX+1
                        RV(I,J,K)=(VV(I+1,J,K)-VV(I-1,J,K))*HX*DXXDX(I)
1010            CONTINUE
c		RW=(CSHIFT(WW,1,1)-CSHIFT(WW,-1,1))*HX*DXXDX
                DO 1020 K=ILAP/2+1,NZ-ILAP/2
                DO 1020 J=2,NY-IY+1
                DO 1020 I=2,NX-IX+1
                        RW(I,J,K)=(WW(I+1,J,K)-WW(I-1,J,K))*HX*DXXDX(I)
1020            CONTINUE
c		WW2=BX*RV-BY*RU
		DO 1030 K=ILAP/2+1,NZ-ILAP/2
                DO 1030 J=2,NY-IY+1
                DO 1030 I=2,NX-IX+1
                        WW2(I,J,K)=BX(I,J,K)*RV(I,J,K)
     2						-BY(I,J,K)*RU(I,J,K)
1030            CONTINUE
c		WW3=BX*RW-BZ*RU
		DO 1040 K=ILAP/2+1,NZ-ILAP/2
                DO 1040 J=2,NY-IY+1
                DO 1040 I=2,NX-IX+1
                        WW3(I,J,K)=BX(I,J,K)*RW(I,J,K)
     2                                          -BZ(I,J,K)*RU(I,J,K)
1040            CONTINUE
c		RU=(CSHIFT(UU,1,2)-CSHIFT(UU,-1,2))*HY*DYYDY
                DO 1050 K=ILAP/2+1,NZ-ILAP/2
                DO 1050 J=2,NY-IY+1
                DO 1050 I=2,NX-IX+1
                        RU(I,J,K)=(UU(I,J+1,K)-UU(I,J-1,K))*HY*DYYDY(J)
1050            CONTINUE
c               RV=(CSHIFT(VV,1,2)-CSHIFT(VV,-1,2))*HY*DYYDY
                DO 1060 K=ILAP/2+1,NZ-ILAP/2
                DO 1060 J=2,NY-IY+1
                DO 1060 I=2,NX-IX+1
                        RV(I,J,K)=(VV(I,J+1,K)-VV(I,J-1,K))*HY*DYYDY(J)
1060            CONTINUE
c               RW=(CSHIFT(WW,1,2)-CSHIFT(WW,-1,2))*HY*DYYDY
                DO 1070 K=ILAP/2+1,NZ-ILAP/2
                DO 1070 J=2,NY-IY+1
                DO 1070 I=2,NX-IX+1
                        RW(I,J,K)=(WW(I,J+1,K)-WW(I,J-1,K))*HY*DYYDY(J)
1070             CONTINUE
c               WW1=BY*RU-BX*RV
                DO 1080 K=ILAP/2+1,NZ-ILAP/2
                DO 1080 J=2,NY-IY+1
                DO 1080 I=2,NX-IX+1
                        WW1(I,J,K)=BY(I,J,K)*RU(I,J,K)
     2                                          -BX(I,J,K)*RV(I,J,K)
1080            CONTINUE
c               WW3=WW3+BY*RW-BZ*RV
                DO 1090 K=ILAP/2+1,NZ-ILAP/2
                DO 1090 J=2,NY-IY+1
                DO 1090 I=2,NX-IX+1
                        WW3(I,J,K)=WW3(I,J,K)+BY(I,J,K)*RW(I,J,K)
     2                                          -BZ(I,J,K)*RV(I,J,K)
1090            CONTINUE
c               RU=(CSHIFT(UU,1,3)-CSHIFT(UU,-1,3))*HZ*DZZDZ
                DO 1100 K=ILAP/2+1,NZ-ILAP/2
                DO 1100 J=2,NY-IY+1
                DO 1100 I=2,NX-IX+1
                        RU(I,J,K)=(UU(I,J,K+1)-UU(I,J,K-1))*HZ*DZZDZ(K)
1100            CONTINUE
c               RV=(CSHIFT(VV,1,3)-CSHIFT(VV,-1,3))*HZ*DZZDZ
                DO 1110 K=ILAP/2+1,NZ-ILAP/2
                DO 1110 J=2,NY-IY+1
                DO 1110 I=2,NX-IX+1
                        RV(I,J,K)=(VV(I,J,K+1)-VV(I,J,K-1))*HZ*DZZDZ(K)
1110            CONTINUE
c               RW=(CSHIFT(WW,1,3)-CSHIFT(WW,-1,3))*HZ*DZZDZ
                DO 1120 K=ILAP/2+1,NZ-ILAP/2
                DO 1120 J=2,NY-IY+1
                DO 1120 I=2,NX-IX+1
                        RW(I,J,K)=(WW(I,J,K+1)-WW(I,J,K-1))*HZ*DZZDZ(K)
1120            CONTINUE
		IF (MYPEZ.EQ.0) THEN
                    DO 1125 J=2,NY-IY+1
                    DO 1125 I=2,NX-IX+1
                       	RW(I,J,ILAP/2+1)=
     2			     (4.0E00*WW(I,J,ILAP/2+2)-WW(I,J,ILAP/2+3))
     3                        			    *HZ*DZZDZ(ILAP/2+1)
1125                CONTINUE
        	ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
                    DO 1126 J=2,NY-IY+1
                    DO 1126 I=2,NX-IX+1
                        RW(I,J,NZ-ILAP/2)=
     2			(WW(I,J,NZ-ILAP/2-2)-4.0E00*WW(I,J,NZ-ILAP/2-1))
     3                                              *HZ*DZZDZ(NZ-ILAP/2)
1126                CONTINUE
        	ENDIF
c               WW1=WW1+BZ*RU-BX*RW
                DO 1130 K=ILAP/2+1,NZ-ILAP/2
                DO 1130 J=2,NY-IY+1
                DO 1130 I=2,NX-IX+1
                        WW1(I,J,K)=WW1(I,J,K)+BZ(I,J,K)*RU(I,J,K)
     2                                          -BX(I,J,K)*RW(I,J,K)
1130             CONTINUE
c               WW2=WW2+BZ*RV-BY*RW
                DO 1140 K=ILAP/2+1,NZ-ILAP/2
                DO 1140 J=2,NY-IY+1
                DO 1140 I=2,NX-IX+1
                        WW2(I,J,K)=WW2(I,J,K)+BZ(I,J,K)*RV(I,J,K)
     2                                          -BY(I,J,K)*RW(I,J,K)
1140             CONTINUE
c		RU=(CSHIFT(BX,1,1)-CSHIFT(BX,-1,1))*HX*DXXDX
        	DO 1150 K=ILAP/2+1,NZ-ILAP/2
        	DO 1150 J=2,NY-IY+1
        	DO 1150 I=2,NX-IX+1
                	RU(I,J,K)=(BX(I+1,J,K)-BX(I-1,J,K))
     2                                                  *HX*DXXDX(I)
1150     	CONTINUE
c               RV=(CSHIFT(BY,1,1)-CSHIFT(BY,-1,1))*HX*DXXDX
                DO 1160 K=ILAP/2+1,NZ-ILAP/2
                DO 1160 J=2,NY-IY+1
                DO 1160 I=2,NX-IX+1
                        RV(I,J,K)=(BY(I+1,J,K)-BY(I-1,J,K))
     2                                                  *HX*DXXDX(I)
1160     	CONTINUE
c               RW=(CSHIFT(BZ,1,1)-CSHIFT(BZ,-1,1))*HX*DXXDX
                DO 1170 K=ILAP/2+1,NZ-ILAP/2
                DO 1170 J=2,NY-IY+1
                DO 1170 I=2,NX-IX+1
                        RW(I,J,K)=(BZ(I+1,J,K)-BZ(I-1,J,K))
     2                                                  *HX*DXXDX(I)
1170     	CONTINUE
c               TT=(CSHIFT(BX,1,2)-CSHIFT(BX,-1,2))*HY*DYYDY
                DO 1180 K=ILAP/2+1,NZ-ILAP/2
                DO 1180 J=2,NY-IY+1
                DO 1180 I=2,NX-IX+1
                        TT(I,J,K)=(BX(I,J+1,K)-BX(I,J-1,K))
     2                                                  *HY*DYYDY(J)
1180             CONTINUE
c       	FU=FU+OBETA*BY*(TT-RV)-OBETA*BZ*RW
        	DO 1190 K=ILAP/2+1,NZ-ILAP/2
        	DO 1190 J=2,NY-IY+1
        	DO 1190 I=2,NX-IX+1
                	FU(I,J,K)=FU(I,J,K)+OBETA*BY(I,J,K)
     2					      *(TT(I,J,K)-RV(I,J,K))
     3					   -OBETA*BZ(I,J,K)*RW(I,J,K)
1190      	CONTINUE
c               FV=FV+OBETA*BX*(RV-TT)
                DO 1200 K=ILAP/2+1,NZ-ILAP/2
                DO 1200 J=2,NY-IY+1
                DO 1200 I=2,NX-IX+1
                        FV(I,J,K)=FV(I,J,K)+OBETA*BX(I,J,K)
     2                                        *(RV(I,J,K)-TT(I,J,K))
1200            CONTINUE
c               FW=FW+OBETA*BX*RW
                DO 1210 K=ILAP/2+1,NZ-ILAP/2
                DO 1210 J=2,NY-IY+1
                DO 1210 I=2,NX-IX+1
                        FW(I,J,K)=FW(I,J,K)+OBETA*BX(I,J,K)*RW(I,J,K)
1210             CONTINUE
c       	FT=FT+OBETA*ORM*RO*OCV*(TT*TT+RV*RV+RW*RW-2.0E00*RV*TT)
        	TMP=OBETA*ORM*OCV
        	DO 1220 K=ILAP/2+1,NZ-ILAP/2
        	DO 1220 J=2,NY-IY+1
        	DO 1220 I=2,NX-IX+1
                	FT(I,J,K)=FT(I,J,K)+(TT(I,J,K)*TT(I,J,K)
     2				    	 +RV(I,J,K)*RV(I,J,K)
     3				         +RW(I,J,K)*RW(I,J,K)
     4				         -2.0E00*RV(I,J,K)*TT(I,J,K))
     5                                                 *RO(I,J,K)*TMP
1220     CONTINUE
c       	WW1=WW1-UU*RU-VV*TT
        	DO 1230 K=ILAP/2+1,NZ-ILAP/2
                DO 1230 J=2,NY-IY+1
                DO 1230 I=2,NX-IX+1
                        WW1(I,J,K)=WW1(I,J,K)-UU(I,J,K)*RU(I,J,K)
     2                                       -VV(I,J,K)*TT(I,J,K)
1230             CONTINUE
c               WW2=WW2-UU*RV
                DO 1240 K=ILAP/2+1,NZ-ILAP/2
                DO 1240 J=2,NY-IY+1
                DO 1240 I=2,NX-IX+1
                        WW2(I,J,K)=WW2(I,J,K)-UU(I,J,K)*RV(I,J,K)
1240             CONTINUE
c               WW3=WW3-UU*RW
                DO 1250 K=ILAP/2+1,NZ-ILAP/2
                DO 1250 J=2,NY-IY+1
                DO 1250 I=2,NX-IX+1
                        WW3(I,J,K)=WW3(I,J,K)-UU(I,J,K)*RW(I,J,K)
1250             CONTINUE
c		RU=RU*(1.0E00/DXXDX)*D2XXDX2
c     2    		+(CSHIFT(BX,1,1)-2.0E00*BX+CSHIFT(BX,-1,1))
c     3						       *H2X*DXXDX*DXXDX
	        DO 1260 K=ILAP/2+1,NZ-ILAP/2
        	DO 1260 J=2,NY-IY+1
        	DO 1260 I=2,NX-IX+1
                	RU(I,J,K)=(BX(I+1,J,K)-BX(I-1,J,K))
     2							*HX*D2XXDX2(I)
     3                     +(BX(I+1,J,K)-2.0E00*BX(I,J,K)+BX(I-1,J,K))
     4                              		*H2X*DXXDX(I)*DXXDX(I)
1260     	CONTINUE
c               RV=RV*(1.0E00/DXXDX)*D2XXDX2
c     2                 +(CSHIFT(BY,1,1)-2.0E00*BY+CSHIFT(BY,-1,1))
c     3                                                *H2X*DXXDX*DXXDX
                DO 1265 K=ILAP/2+1,NZ-ILAP/2
                DO 1265 J=2,NY-IY+1
                DO 1265 I=2,NX-IX+1
                        RV(I,J,K)=(BY(I+1,J,K)-BY(I-1,J,K))
     2                                                  *HX*D2XXDX2(I)
     3                     +(BY(I+1,J,K)-2.0E00*BY(I,J,K)+BY(I-1,J,K))
     4                                          *H2X*DXXDX(I)*DXXDX(I)
1265     	CONTINUE
c               RW=RW*(1.0E00/DXXDX)*D2XXDX2
c     2                 +(CSHIFT(BY,1,1)-2.0E00*BZ+CSHIFT(BZ,-1,1))
c     3                                                *H2X*DXXDX*DXXDX
c                DO 329 K=ILAP/2+1,NZ-ILAP/2
c                DO 329 J=2,NY-IY+1
c                DO 329 I=2,NX-IX+1
c                        RW(I,J,K)=(BZ(I+1,J,K)-BZ(I-1,J,K))
c     2                                                  *HX*D2XXDX2(I)
c     3                     +(BZ(I+1,J,K)-2.0E00*BZ(I,J,K)+BZ(I-1,J,K))
c     4                                          *H2X*DXXDX(I)*DXXDX(I)
c329     	CONTINUE
c       	TT=TT*(1.0E00/DYYDY)*D2YYDY2
c     2    		+(CSHIFT(BX,1,2)-2.0E00*BX+CSHIFT(BX,-1,2))
c     3						      *H2Y*DYYDY*DYYDY
        	DO 1270 K=ILAP/2+1,NZ-ILAP/2
        	DO 1270 J=2,NY-IY+1
        	DO 1270 I=2,NX-IX+1
                	TT(I,J,K)=(BX(I,J+1,K)-BX(I,J-1,K))
     2							*HY*D2YYDY2(J)
     2                     +(BX(I,J+1,K)-2.0E00*BX(I,J,K)+BX(I,J-1,K))
     3                              		*H2Y*DYYDY(J)*DYYDY(J)
1270     CONTINUE
c               WW1=WW1+ORM*(RU+TT)
                DO 1280 K=ILAP/2+1,NZ-ILAP/2
                DO 1280 J=2,NY-IY+1
                DO 1280 I=2,NX-IX+1
                        WW1(I,J,K)=WW1(I,J,K)+ORM*(RU(I,J,K)+TT(I,J,K))
1280             CONTINUE
c               WW2=WW2+ORM*RV
                DO 1290 K=ILAP/2+1,NZ-ILAP/2
                DO 1290 J=2,NY-IY+1
                DO 1290 I=2,NX-IX+1
                        WW2(I,J,K)=WW2(I,J,K)+ORM*RV(I,J,K)
1290             CONTINUE
c               TT=RW*(1.0E00/DXXDX)*D2XXDX2
c     2                 +(CSHIFT(BY,1,1)-2.0E00*BZ+CSHIFT(BZ,-1,1))
c     3                                                *H2X*DXXDX*DXXDX
                DO 1300 K=ILAP/2+1,NZ-ILAP/2
                DO 1300 J=2,NY-IY+1
                DO 1300 I=2,NX-IX+1
                        TT(I,J,K)=(BZ(I+1,J,K)-BZ(I-1,J,K))
     2                                                  *HX*D2XXDX2(I)
     3                     +(BZ(I+1,J,K)-2.0E00*BZ(I,J,K)+BZ(I-1,J,K))
     4                                          *H2X*DXXDX(I)*DXXDX(I)
1300             CONTINUE
c               WW3=WW3+ORM*TT
                DO 1310 K=ILAP/2+1,NZ-ILAP/2
                DO 1310 J=2,NY-IY+1
                DO 1310 I=2,NX-IX+1
                        WW3(I,J,K)=WW3(I,J,K)+ORM*TT(I,J,K)
1310             CONTINUE
c               RU=(CSHIFT(BX,1,3)-CSHIFT(BX,-1,3))*HZ*DZZDZ
                DO 1320 K=ILAP/2+1,NZ-ILAP/2
                DO 1320 J=2,NY-IY+1
                DO 1320 I=2,NX-IX+1
                        RU(I,J,K)=(BX(I,J,K+1)-BX(I,J,K-1))
     2                                                  *HZ*DZZDZ(K)
1320             CONTINUE
c               RV=(CSHIFT(BY,1,2)-CSHIFT(BY,-1,2))*HY*DYYDY
                DO 1330 K=ILAP/2+1,NZ-ILAP/2
                DO 1330 J=2,NY-IY+1
                DO 1330 I=2,NX-IX+1
                        RV(I,J,K)=(BY(I,J+1,K)-BY(I,J-1,K))
     2                                                  *HY*DYYDY(J)
1330            CONTINUE
c               FU=FU+OBETA*BZ*RU
                DO 1340 K=ILAP/2+1,NZ-ILAP/2
                DO 1340 J=2,NY-IY+1
                DO 1340 I=2,NX-IX+1
                        FU(I,J,K)=FU(I,J,K)+OBETA*BZ(I,J,K)*RU(I,J,K)
1340             CONTINUE
c               FW=FW-OBETA*BX*RU
                DO 1350 K=ILAP/2+1,NZ-ILAP/2
                DO 1350 J=2,NY-IY+1
                DO 1350 I=2,NX-IX+1
                        FW(I,J,K)=FW(I,J,K)-OBETA*BX(I,J,K)*RU(I,J,K)
1350             CONTINUE
c               FT=FT+OBETA*ORM*RO*OCV*(RU*RU-2.0E00*RU*RW)
                TMP=OBETA*ORM*OCV
                DO 1360 K=ILAP/2+1,NZ-ILAP/2
                DO 1360 J=2,NY-IY+1
                DO 1360 I=2,NX-IX+1
                        FT(I,J,K)=FT(I,J,K)+(RU(I,J,K)*RU(I,J,K)
     2                                   -2.0E00*RU(I,J,K)*RW(I,J,K))
     3                                                 *RO(I,J,K)*TMP
1360     	CONTINUE
c               WW1=WW1-WW*RU
                DO 1370 K=ILAP/2+1,NZ-ILAP/2
                DO 1370 J=2,NY-IY+1
                DO 1370 I=2,NX-IX+1
                        WW1(I,J,K)=WW1(I,J,K)-WW(I,J,K)*RU(I,J,K)
1370             CONTINUE
c               WW2=WW2-VV*RV
                DO 1380 K=ILAP/2+1,NZ-ILAP/2
                DO 1380 J=2,NY-IY+1
                DO 1380 I=2,NX-IX+1
                        WW2(I,J,K)=WW2(I,J,K)-VV(I,J,K)*RV(I,J,K)
1380             CONTINUE
c               RU=RU*(1.0E00/DZZDZ)*D2ZZDZ2
c     2                 +(CSHIFT(BX,1,3)-2.0E00*BX+CSHIFT(BX,-1,3))
c     3                                                *H2Z*DZZDZ*DZZDZ
                DO 1390 K=ILAP/2+1,NZ-ILAP/2
                DO 1390 J=2,NY-IY+1
                DO 1390 I=2,NX-IX+1
                        RU(I,J,K)=(BX(I,J,K+1)-BX(I,J,K-1))
     2                                                  *HZ*D2ZZDZ2(K)
     3                     +(BX(I,J,K+1)-2.0E00*BX(I,J,K)+BX(I,J,K-1))
     4                                          *H2Z*DZZDZ(K)*DZZDZ(K)
1390            CONTINUE
		IF (MYPEZ.EQ.0) THEN
                    DO 1391 J=2,NY-IY+1
                    DO 1391 I=2,NX-IX+1
                        RU(I,J,ILAP/2+1)=(-3.0E00*BX(I,J,ILAP/2+1)
     2                        +4.0E00*BX(I,J,ILAP/2+2)-BX(I,J,ILAP/2+3))
     3                                             *HZ*D2ZZDZ2(ILAP/2+1)
     4		       +(2.0E00*BX(I,J,ILAP/2+1)-5.0E00*BX(I,J,ILAP/2+2)
     5			 +4.0E00*BX(I,J,ILAP/2+3)-BX(I,J,ILAP/2+4))
     6			       *H2Z*DZZDZ(ILAP/2+1)*DZZDZ(ILAP/2+1)
1391                CONTINUE
                ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
                    DO 1392 J=2,NY-IY+1
                    DO 1392 I=2,NX-IX+1
                        RU(I,J,NZ-ILAP/2)=(3.0E00*BX(I,J,NZ-ILAP/2)
     2                  -4.0E00*BX(I,J,NZ-ILAP/2-1)+BX(I,J,NZ-ILAP/2-2))
     3                                            *HZ*D2ZZDZ2(NZ-ILAP/2)
     4		   +(2.0E00*BX(I,J,NZ-ILAP/2)-5.0E00*BX(I,J,NZ-ILAP/2-1)
     5		    +4.0E00*BX(I,J,NZ-ILAP/2-2)-BX(I,J,NZ-ILAP/2-3))
     6			      *H2Z*DZZDZ(NZ-ILAP/2)*DZZDZ(NZ-ILAP/2)
1392                CONTINUE
                ENDIF
c               RV=RV*(1.0E00/DYYDY)*D2YYDY2
c     2                 +(CSHIFT(BY,1,2)-2.0E00*BY+CSHIFT(BY,-1,2))
c     3                                                *H2Y*DYYDY*DYYDY
                DO 1400 K=ILAP/2+1,NZ-ILAP/2
                DO 1400 J=2,NY-IY+1
                DO 1400 I=2,NX-IX+1
                        RV(I,J,K)=(BY(I,J+1,K)-BY(I,J-1,K))
     2                                                  *HY*D2YYDY2(J)
     3                     +(BY(I,J+1,K)-2.0E00*BY(I,J,K)+BY(I,J-1,K))
     4                                          *H2Y*DYYDY(J)*DYYDY(J)
1400            CONTINUE
c               WW1=WW1+ORM*RU
                DO 1410 K=ILAP/2+1,NZ-ILAP/2
                DO 1410 J=2,NY-IY+1
                DO 1410 I=2,NX-IX+1
                        WW1(I,J,K)=WW1(I,J,K)+ORM*RU(I,J,K)
1410            CONTINUE
c               WW2=WW2+ORM*RV
                DO 1420 K=ILAP/2+1,NZ-ILAP/2
                DO 1420 J=2,NY-IY+1
                DO 1420 I=2,NX-IX+1
                        WW2(I,J,K)=WW2(I,J,K)+ORM*RV(I,J,K)
1420             CONTINUE
c               RV=(CSHIFT(BY,1,3)-CSHIFT(BY,-1,3))*HZ*DZZDZ
                DO 1430 K=ILAP/2+1,NZ-ILAP/2
                DO 1430 J=2,NY-IY+1
                DO 1430 I=2,NX-IX+1
                        RV(I,J,K)=(BY(I,J,K+1)-BY(I,J,K-1))
     2                                                  *HZ*DZZDZ(K)
1430             CONTINUE
c               RW=(CSHIFT(BZ,1,2)-CSHIFT(BZ,-1,2))*HY*DYYDY
                DO 1440 K=ILAP/2+1,NZ-ILAP/2
                DO 1440 J=2,NY-IY+1
                DO 1440 I=2,NX-IX+1
                        RW(I,J,K)=(BZ(I,J+1,K)-BZ(I,J-1,K))
     2                                                  *HY*DYYDY(J)
1440            CONTINUE
c		RU=RV-RW
		DO 1450 K=ILAP/2+1,NZ-ILAP/2
                DO 1450 J=2,NY-IY+1
                DO 1450 I=2,NX-IX+1
			RU(I,J,K)=RV(I,J,K)-RW(I,J,K)
1450		CONTINUE
c               FV=FV+OBETA*BZ*RU
                DO 1460 K=ILAP/2+1,NZ-ILAP/2
                DO 1460 J=2,NY-IY+1
                DO 1460 I=2,NX-IX+1
                        FV(I,J,K)=FV(I,J,K)+OBETA*BZ(I,J,K)*RU(I,J,K)
1460            CONTINUE
c               FW=FW-OBETA*BY*RU
                DO 1470 K=ILAP/2+1,NZ-ILAP/2
                DO 1470 J=2,NY-IY+1
                DO 1470 I=2,NX-IX+1
                        FW(I,J,K)=FW(I,J,K)-OBETA*BY(I,J,K)*RU(I,J,K)
1470            CONTINUE
c               FT=FT+OBETA*ORM*RO*OCV*(RV*RV+RW*RW-2.0E00*RV*RW)
                TMP=OBETA*ORM*OCV
                DO 1480 K=ILAP/2+1,NZ-ILAP/2
                DO 1480 J=2,NY-IY+1
                DO 1480 I=2,NX-IX+1
                        FT(I,J,K)=FT(I,J,K)+(RV(I,J,K)*RV(I,J,K)
     2					    +RW(I,J,K)*RW(I,J,K)
     3                                   -2.0E00*RV(I,J,K)*RW(I,J,K))
     4                                                 *RO(I,J,K)*TMP
1480            CONTINUE
c               WW2=WW2-WW*RV
                DO 1490 K=ILAP/2+1,NZ-ILAP/2
                DO 1490 J=2,NY-IY+1
                DO 1490 I=2,NX-IX+1
                        WW2(I,J,K)=WW2(I,J,K)-WW(I,J,K)*RV(I,J,K)
1490            CONTINUE
c               WW3=WW3-VV*RW
                DO 1500 K=ILAP/2+1,NZ-ILAP/2
                DO 1500 J=2,NY-IY+1
                DO 1500 I=2,NX-IX+1
                        WW3(I,J,K)=WW3(I,J,K)-VV(I,J,K)*RW(I,J,K)
1500            CONTINUE
c               RV=RV*(1.0E00/DZZDZ)*D2ZZDZ2
c     2                 +(CSHIFT(BY,1,3)-2.0E00*BY+CSHIFT(BY,-1,3))
c     3                                                *H2Z*DZZDZ*DZZDZ
                DO 1510 K=ILAP/2+1,NZ-ILAP/2
                DO 1510 J=2,NY-IY+1
                DO 1510 I=2,NX-IX+1
                        RV(I,J,K)=(BY(I,J,K+1)-BY(I,J,K-1))
     2                                                  *HZ*D2ZZDZ2(K)
     3                     +(BY(I,J,K+1)-2.0E00*BY(I,J,K)+BY(I,J,K-1))
     4                                          *H2Z*DZZDZ(K)*DZZDZ(K)
1510            CONTINUE
		IF (MYPEZ.EQ.0) THEN
                    DO 1511 J=2,NY-IY+1
                    DO 1511 I=2,NX-IX+1
                        RV(I,J,ILAP/2+1)=(-3.0E00*BY(I,J,ILAP/2+1)
     2                        +4.0E00*BY(I,J,ILAP/2+2)-BY(I,J,ILAP/2+3))
     3                                             *HZ*D2ZZDZ2(ILAP/2+1)
     4                 +(2.0E00*BY(I,J,ILAP/2+1)-5.0E00*BY(I,J,ILAP/2+2)
     5                   +4.0E00*BY(I,J,ILAP/2+3)-BY(I,J,ILAP/2+4))
     6                         *H2Z*DZZDZ(ILAP/2+1)*DZZDZ(ILAP/2+1)
1511                CONTINUE
                ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
                    DO 1513 J=2,NY-IY+1
                    DO 1513 I=2,NX-IX+1
                        RV(I,J,NZ-ILAP/2)=(3.0E00*BY(I,J,NZ-ILAP/2)
     2                  -4.0E00*BY(I,J,NZ-ILAP/2-1)+BY(I,J,NZ-ILAP/2-2))
     3                                            *HZ*D2ZZDZ2(NZ-ILAP/2)
     4             +(2.0E00*BY(I,J,NZ-ILAP/2)-5.0E00*BY(I,J,NZ-ILAP/2-1)
     5              +4.0E00*BY(I,J,NZ-ILAP/2-2)-BY(I,J,NZ-ILAP/2-3))
     6                        *H2Z*DZZDZ(NZ-ILAP/2)*DZZDZ(NZ-ILAP/2)
1513                CONTINUE
                ENDIF
c               RW=RW*(1.0E00/DYYDY)*D2YYDY2
c     2                 +(CSHIFT(BZ,1,2)-2.0E00*BZ+CSHIFT(BZ,-1,2))
c     3                                                *H2Y*DYYDY*DYYDY
                DO 1520 K=ILAP/2+1,NZ-ILAP/2
                DO 1520 J=2,NY-IY+1
                DO 1520 I=2,NX-IX+1
                        RW(I,J,K)=(BZ(I,J+1,K)-BZ(I,J-1,K))
     2                                                  *HY*D2YYDY2(J)
     3                     +(BZ(I,J+1,K)-2.0E00*BZ(I,J,K)+BZ(I,J-1,K))
     4                                          *H2Y*DYYDY(J)*DYYDY(J)
1520            CONTINUE
c               WW2=WW2+ORM*RV
                DO 1530 K=ILAP/2+1,NZ-ILAP/2
                DO 1530 J=2,NY-IY+1
                DO 1530 I=2,NX-IX+1
                        WW2(I,J,K)=WW2(I,J,K)+ORM*RV(I,J,K)
1530            CONTINUE
c               WW3=WW3+ORM*RW
                DO 1540 K=ILAP/2+1,NZ-ILAP/2
                DO 1540 J=2,NY-IY+1
                DO 1540 I=2,NX-IX+1
                        WW3(I,J,K)=WW3(I,J,K)+ORM*RW(I,J,K)
1540            CONTINUE
c               RW=(CSHIFT(BZ,1,3)-CSHIFT(BZ,-1,3))*HZ*DZZDZ
                DO 1550 K=ILAP/2+1,NZ-ILAP/2
                DO 1550 J=2,NY-IY+1
                DO 1550 I=2,NX-IX+1
                        RW(I,J,K)=(BZ(I,J,K+1)-BZ(I,J,K-1))
     2                                                  *HZ*DZZDZ(K)
1550            CONTINUE
c               WW3=WW3-WW*RW
                DO 1560 K=ILAP/2+1,NZ-ILAP/2
                DO 1560 J=2,NY-IY+1
                DO 1560 I=2,NX-IX+1
                        WW3(I,J,K)=WW3(I,J,K)-WW(I,J,K)*RW(I,J,K)
1560            CONTINUE
c               RW=RW*(1.0E00/DZZDZ)*D2ZZDZ2
c     2                 +(CSHIFT(BZ,1,3)-2.0E00*BZ+CSHIFT(BZ,-1,3))
c     3                                                *H2Z*DZZDZ*DZZDZ
                DO 1570 K=ILAP/2+1,NZ-ILAP/2
                DO 1570 J=2,NY-IY+1
                DO 1570 I=2,NX-IX+1
                        RW(I,J,K)=(BZ(I,J,K+1)-BZ(I,J,K-1))
     2                                                  *HZ*D2ZZDZ2(K)
     3                     +(BZ(I,J,K+1)-2.0E00*BZ(I,J,K)+BZ(I,J,K-1))
     4                                          *H2Z*DZZDZ(K)*DZZDZ(K)
1570            CONTINUE
c               WW3=WW3+ORM*RW
                DO 1580 K=ILAP/2+1,NZ-ILAP/2
                DO 1580 J=2,NY-IY+1
                DO 1580 I=2,NX-IX+1
                        WW3(I,J,K)=WW3(I,J,K)+ORM*RW(I,J,K)
1580             CONTINUE
	ENDIF
C
	RETURN
	END
C**********************************************************************
	SUBROUTINE BCON
C
	INCLUDE '3dmhdparam.f'
C
        DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2           ,TT(NX,NY,NZ)
        DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
        DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2           ,FT(NX,NY,NZ)
        DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
     2           ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
        DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
C
        DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
        DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
C
        DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
        DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
        DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
	DIMENSION RKAPA(NZ),DKAPA(NZ)
        DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
     2           ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
     3           ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
     4           ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
C
     5           ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
     6           ,SP26(IPAD)
C
        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4            ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C
        COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
	COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
        COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
	COMMON/GRID/DD,HX,H2X,HY,H2Y,HZ,H2Z,C13,C23,C43
        COMMON/CPER/TP,XP,YP,ZP,TC,QFH,HH
	COMMON/CTIM/DT,TIMT,TIMC,TIMI
        COMMON/BCT/IXC,IYC,IZC,ITC,IBC
	COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C----------------------------------------------------------------------
C  Special boundary for T (M. Rempel)
C----------------------------------------------------------------------
	COMMON/SPECIALBOUND/TU,DZTB,DZTU
C----------------------------------------------------------------------
C  Ensure that div(B)=0
C----------------------------------------------------------------------
c        IF ((LMAG).AND.(LPOT)) THEN
c                DO 250 K=ILAP/2+1,NZ-ILAP/2
c                DO 250 J=2,NY-IY+1
c                        BX(1,J,K)=BX(NX-IX+1,J,K)
c                        BX(NX-IX+2,J,K)=BX(2,J,K)
c                        BY(1,J,K)=BY(NX-IX+1,J,K)
c                        BY(NX-IX+2,J,K)=BY(2,J,K)
c                        BZ(1,J,K)=BZ(NX-IX+1,J,K)
c                        BZ(NX-IX+2,J,K)=BZ(2,J,K)
c250             CONTINUE
c                DO 260 K=ILAP/2+1,NZ-ILAP/2
c                DO 260 I=1,NX
c                        BX(I,1,K)=BX(I,NY-IY+1,K)
c                        BX(I,NY-IY+2,K)=BX(I,2,K)
c                        BY(I,1,K)=BY(I,NY-IY+1,K)
c                        BY(I,NY-IY+2,K)=BY(I,2,K)
c                        BZ(I,1,K)=BZ(I,NY-IY+1,K)
c                        BZ(I,NY-IY+2,K)=BZ(I,2,K)
c260             CONTINUE
c                CALL BARRIER()
c                IF (MYPE.EQ.0) THEN
c                DO 15 K=1,ILAP/2
c                BX(:,:,K)=BX(:,:,ILAP/2+1)
c                CALL SHMEM_GGG(BX(1,1,NZ-ILAP/2+K),BX(1,1,ILAP/2+K)
c     2                                                    ,NX*NY,1)
c                BY(:,:,K)=BY(:,:,ILAP/2+1)
c                CALL SHMEM_GGG(BY(1,1,NZ-ILAP/2+K),BY(1,1,ILAP/2+K)
c     2                                                    ,NX*NY,1)
c                BZ(:,:,K)=BZ(:,:,ILAP/2+1)
c                CALL SHMEM_GGG(BZ(1,1,NZ-ILAP/2+K),BZ(1,1,ILAP/2+K)
c     2                                                    ,NX*NY,1)
c15              CONTINUE
c                ELSE IF (MYPE.EQ.NPE-1) THEN
c                DO 25 K=1,ILAP/2
c                CALL SHMEM_GGG(BX(1,1,K),BX(1,1,NZ-ILAP+K),NX*NY,NPE-2)
c                BX(:,:,NZ-ILAP/2+K)=BX(:,:,NZ-ILAP/2)
c                CALL SHMEM_GGG(BY(1,1,K),BY(1,1,NZ-ILAP+K),NX*NY,NPE-2)
c                BY(:,:,NZ-ILAP/2+K)=BY(:,:,NZ-ILAP/2)
c                CALL SHMEM_GGG(BZ(1,1,K),BZ(1,1,NZ-ILAP+K),NX*NY,NPE-2)
c                BZ(:,:,NZ-ILAP/2+K)=BZ(:,:,NZ-ILAP/2)
c25              CONTINUE
c                ELSE
c                DO 35 K=1,ILAP/2
c                CALL SHMEM_GGG(BX(1,1,K),BX(1,1,NZ-ILAP+K),NX*NY,MYPE-1)
c                CALL SHMEM_GGG(BX(1,1,NZ-ILAP/2+K),BX(1,1,ILAP/2+K)
c     2                                                    ,NX*NY,MYPE+1)
c                CALL SHMEM_GGG(BY(1,1,K),BY(1,1,NZ-ILAP+K),NX*NY,MYPE-1)
c                CALL SHMEM_GGG(BY(1,1,NZ-ILAP/2+K),BY(1,1,ILAP/2+K)
c     2                                                    ,NX*NY,MYPE+1)
c                CALL SHMEM_GGG(BZ(1,1,K),BZ(1,1,NZ-ILAP+K),NX*NY,MYPE-1)
c                CALL SHMEM_GGG(BZ(1,1,NZ-ILAP/2+K),BZ(1,1,ILAP/2+K)
c     2                                                    ,NX*NY,MYPE+1)
c35              CONTINUE
c                ENDIF
c                CALL BARRIER()
cC
c                CALL POTENTIAL
c        ENDIF
C----------------------------------------------------------------------
C  Satisfy boundary conditions.
C----------------------------------------------------------------------
	IF (LSHR) THEN
                RRT=1.0E00
                RRB=-1.0E00
		RPP=0.0E00
c		RPP=2.0E00*ASIN(1.0E00)/5.0E00
        ENDIF
C----------------------------------------------------------------------
C  Top.
C----------------------------------------------------------------------
	IF (MYPEZ.EQ.0) THEN
C
		IF (LSHR) THEN
			DO 10 J=2,NY-IY+1
                        DO 10 I=2,NX-IX+1
                                RU(I,J,ILAP/2+1)=0.0E00
				RV(I,J,ILAP/2+1)=RRT*COS(RPP*TIMT)
     2						    *RO(I,J,ILAP/2+1)
10			CONTINUE
		ELSE
C
		DO 20 J=2,NY-IY+1
		DO 20 I=2,NX-IX+1
			RU(I,J,ILAP/2+1)=RO(I,J,ILAP/2+1)
     2				*(C43*RU(I,J,ILAP/2+2)/RO(I,J,ILAP/2+2)
     3	     	        	-C13*RU(I,J,ILAP/2+3)/RO(I,J,ILAP/2+3))
20		CONTINUE
		DO 30 J=2,NY-IY+1
		DO 30 I=2,NX-IX+1
			RV(I,J,ILAP/2+1)=RO(I,J,ILAP/2+1)
     2				*(C43*RV(I,J,ILAP/2+2)/RO(I,J,ILAP/2+2)
     3		        	-C13*RV(I,J,ILAP/2+3)/RO(I,J,ILAP/2+3))
30		CONTINUE
C
		ENDIF
C
		DO 40 J=2,NY-IY+1
		DO 40 I=2,NX-IX+1
			RW(I,J,ILAP/2+1)=0.0E00
40		CONTINUE
                IF ((TP.NE.0.0E00).AND.((ITC.EQ.0).OR.(ITC.EQ.1))) THEN
C----------------------------------------------------------------------
C  Temperature or flux peturbation held on upper boundary.
C----------------------------------------------------------------------
                    IF (TC.NE.0.0E00) THEN
                        TPR=TP+(1.0E00-TP)*EXP(-TIMT/TC)
                    ELSE
                        TPR=TP
                    ENDIF
		    IF (ITC.EQ.0) THEN
C----------------------------------------------------------------------
C  Gaussian temperature profile with full width at half maximum
C  equal to HH and minimum equal to TP.
C----------------------------------------------------------------------
			CLN=-4.0E00*LOG(2.0E00)/HH/HH
			DO 50 J=2,NY-IY+1
			DO 50 I=2,NX-IX+1
				TT(I,J,ILAP/2+1)=1.0E00-(1.0E00-TPR)
     2					   *EXP(CLN*(EXX(I)-XP)**2)/HH
     3					   *EXP(CLN*(WYY(J)-YP)**2)/HH
50			CONTINUE
		    ELSE
C----------------------------------------------------------------------
C  Gaussian flux profile with full width at half maximum equal to 1.0
C  and maximum equal to TP.
C----------------------------------------------------------------------
			CLN=-4.0E00*LOG(2.0E00)/HH/HH
                	DZ=1.0E00/FLOAT(NPZ-1)/DZZDZ(ILAP/2+1)
			DO 60 J=2,NY-IY+1
			DO 60 I=2,NX-IX+1
                		TT(I,J,ILAP/2+1)=C43*TT(I,J,ILAP/2+2)
     2                                	    -C13*TT(I,J,ILAP/2+3)
     3                                	    -C23*DZ*(THETA-(THETA-TPR)
     4                                	    *EXP(CLN*(EXX(I)-XP)**2)/HH
     5                                	    *EXP(CLN*(WYY(J)-YP)**2)/HH)
60			CONTINUE
        	    ENDIF  
		ELSE
                    IF ((ITC.EQ.0).OR.(ITC.EQ.2).OR.(ITC.EQ.4)) THEN
C----------------------------------------------------------------------
C  Constant temperature upper boundary, no plume perturbation. Embedded 
C  heat loss plume hard wired for constant temp upper boundary codition.
C----------------------------------------------------------------------
			DO 61 J=2,NY-IY+1
                    	DO 61 I=2,NX-IX+1
				IF (LREM) THEN
					TT(I,J,ILAP/2+1)=TU
				ELSE
                     			TT(I,J,ILAP/2+1)=1.0E00
					 IF (ITC.EQ.4)THEN
C----------------------------------------------------------------------
C Half domain TP cooler.
C----------------------------------------------------------------------
					IF (I.LE.NX/2) THEN
						TT(I,J,ILAP/2+1)=
     2							1.0E00-TP
					ENDIF
					ENDIF
C----------------------------------------------------------------------
				ENDIF
61		    	CONTINUE
		    ENDIF 
		    IF (ITC.EQ.1) THEN
C----------------------------------------------------------------------
C  Constant flux upper boundary, no plume perturbation.
C----------------------------------------------------------------------
			DZ=1.0E00/FLOAT(NPZ-1)/DZZDZ(ILAP/2+1)
			DO 62 J=2,NY-IY+1
                        DO 62 I=2,NX-IX+1
		    	IF (LREM) THEN
			    TT(I,J,ILAP/2+1)=C43*TT(I,J,ILAP/2+2)
     2                                      -C13*TT(I,J,ILAP/2+3)
     3                                      +DZTU
			ELSE
                            TT(I,J,ILAP/2+1)=C43*TT(I,J,ILAP/2+2)
     2					    -C13*TT(I,J,ILAP/2+3)
     3					    -C23*DZ*THETA
			ENDIF
62                      CONTINUE
	       	    ENDIF
		ENDIF
		IF (LMAG) THEN
			IF (IBC.EQ.0) THEN
				DO 65 J=2,NY-IY+1
				DO 65 I=2,NX-IX+1
					BZ(I,J,ILAP/2+1)=0.0E00
65				CONTINUE
			ENDIF
		ENDIF
        ENDIF
C----------------------------------------------------------------------
C  Bottom.
C----------------------------------------------------------------------
	IF (MYPEZ.EQ.NPEZ-1) THEN
C
		IF (LSHR) THEN
                        DO 69 J=2,NY-IY+1
                        DO 69 I=2,NX-IX+1
                                RU(I,J,NZ-ILAP/2)=0.0E00
                                RV(I,J,NZ-ILAP/2)=RRB*COS(RPP*TIMT)
     2						     *RO(I,J,NZ-ILAP/2)
69                      CONTINUE
        	ELSE
C
		DO 70 J=2,NY-IY+1
		DO 70 I=2,NX-IX+1
    			RU(I,J,NZ-ILAP/2)=RO(I,J,NZ-ILAP/2)
     2			      	      *(C43*RU(I,J,NZ-ILAP/2-1)
     3				       /RO(I,J,NZ-ILAP/2-1)
     4			      	      -C13*RU(I,J,NZ-ILAP/2-2)
     5			      	       /RO(I,J,NZ-ILAP/2-2))
70		CONTINUE
		DO 80 J=2,NY-IY+1
		DO 80 I=2,NX-IX+1
    			RV(I,J,NZ-ILAP/2)=RO(I,J,NZ-ILAP/2)
     1			      	      *(C43*RV(I,J,NZ-ILAP/2-1)
     2			      		/RO(I,J,NZ-ILAP/2-1)
     3			      	      -C13*RV(I,J,NZ-ILAP/2-2)
     4		  	      		/RO(I,J,NZ-ILAP/2-2))
80		CONTINUE
C
		ENDIF
C
		DO 90 J=2,NY-IY+1
		DO 90 I=2,NX-IX+1
			RW(I,J,NZ-ILAP/2)=0.0E00
90		CONTINUE
		IF (IZC.EQ.0) THEN
C----------------------------------------------------------------------
C  Constant temperature lower boundary.
C----------------------------------------------------------------------
			DO 100 J=2,NY-IY+1
			DO 100 I=2,NX-IX+1
				TT(I,J,NZ-ILAP/2)=TB
100			CONTINUE
		ELSE IF (IZC.EQ.1) THEN
C----------------------------------------------------------------------
C  Constant flux lower boundary.
C----------------------------------------------------------------------
                        DZ=1.0E00/FLOAT(NPZ-1)/DZZDZ(NZ-ILAP/2)
                        DO 105 J=2,NY-IY+1
                        DO 105 I=2,NX-IX+1
			IF (LREM) THEN
			    TT(I,J,NZ-ILAP/2)=C43*TT(I,J,NZ-ILAP/2-1)
     2                                   -C13*TT(I,J,NZ-ILAP/2-2)
     3                                   +DZTB
			ELSE
                            TT(I,J,NZ-ILAP/2)=C43*TT(I,J,NZ-ILAP/2-1)
     2                                   -C13*TT(I,J,NZ-ILAP/2-2)
     3                                   +C23*DZ*THETA/RKAPA(NZ-ILAP/2)
			ENDIF
105			CONTINUE
		ELSE
			WRITE(6,*)'BCON:  Invalid IZCON'
			CALL MPI_FINALIZE(IERR)
                	STOP
		ENDIF
		IF (LMAG) THEN
			IF (IBC.EQ.0) THEN
				DO 110 J=2,NY-IY+1
				DO 110 I=2,NX-IX+1
					BZ(I,J,NZ-ILAP/2)=0.0E00
110				CONTINUE
			ENDIF
		ENDIF
	ENDIF
C
	RETURN
	END
C**********************************************************************
        SUBROUTINE PEEXTN(PESTRNG)
C----------------------------------------------------------------------
C   Creates character string from processor number with padded zeros.
C   Written by Joe Werne.
C----------------------------------------------------------------------
	COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C
	CHARACTER*4 PEZERO,PESTRNG
C
	PEZERO='0000'
	WRITE(PESTRNG,'(I4)')MYPE
	DO 10 K=1,4
		I0=INDEX(PESTRNG,' ')
		IF (I0.EQ.0) GOTO 20
		PESTRNG(1:I0)=PEZERO(1:I0)
10	CONTINUE
C
20	RETURN
	END
C**********************************************************************
        SUBROUTINE COMMUNICATE
C----------------------------------------------------------------------
C  Does communications between processors.
C----------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
	include 'mpif.h'
C
        DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2           ,TT(NX,NY,NZ)
        DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
        DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2           ,FT(NX,NY,NZ)
        DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
     2           ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
        DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
C
        DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
        DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
C
        DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
     2           ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
     3           ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
     4           ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
C
     5           ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
     6           ,SP26(IPAD)
C
        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4            ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C
	COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C
	CALL COMM_MPI(RU)
        CALL COMM_MPI(RV)
        CALL COMM_MPI(RW)
        CALL COMM_MPI(TT)
        CALL COMM_MPI(RO)
        IF(LMAG) THEN
           CALL COMM_MPI(BX)
           CALL COMM_MPI(BY)
           CALL COMM_MPI(BZ)
        ENDIF
C
	RETURN
	END
C**********************************************************************
	SUBROUTINE COMM_MPI(VAR)
C----------------------------------------------------------------------
C  Subroutine COMM_MPI now also enforces horizontal periodicity (M. Rempel).
C----------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
        include 'mpif.h'
C
	DIMENSION VAR(NX,NY,NZ),ISTATUS(MPI_STATUS_SIZE)
C
	COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C----------------------------------------------------------------------
C  Vertical communication.
C----------------------------------------------------------------------
        IF(NPEZ.GT.1) THEN
           ITAG = 100
           IF (MYPEZ.EQ.0) THEN
              CALL MPI_SEND(VAR(1,1,NZ-ILAP+1),NX*NY*(ILAP/2),MPISIZE,
     2                     MYPE+NPEY,ITAG,MPI_COMM_WORLD,IERR)
           ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
              CALL MPI_RECV(VAR(1,1,1),NX*NY*(ILAP/2),MPISIZE,MYPE-NPEY,
     2                     ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
           ELSE
C              print*, "pt 2 mype, npey", mype, npey
              CALL MPI_SENDRECV(VAR(1,1,NZ-ILAP+1),NX*NY*(ILAP/2),
     2                    MPISIZE,MYPE+NPEY,ITAG,VAR(1,1,1),
     3                    NX*NY*(ILAP/2),MPISIZE,MYPE-NPEY,ITAG,
     4                    MPI_COMM_WORLD,ISTATUS,IERR)
           ENDIF
C
           ITAG = 200
           IF (MYPEZ.EQ.0) THEN
              CALL MPI_RECV(VAR(1,1,NZ-ILAP/2+1),NX*NY*(ILAP/2),MPISIZE,
     2                    MYPE+NPEY,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
           ELSE IF (MYPEZ.EQ.NPEZ-1) THEN
              CALL MPI_SEND(VAR(1,1,ILAP/2+1),NX*NY*(ILAP/2),MPISIZE,
     2                    MYPE-NPEY,ITAG,MPI_COMM_WORLD,IERR)
           ELSE
C                print*, "pt 3 mype, npey", mype, npey
              CALL MPI_SENDRECV(VAR(1,1,ILAP/2+1),NX*NY*(ILAP/2),
     2                    MPISIZE,MYPE-NPEY,ITAG,VAR(1,1,NZ-ILAP/2+1),
     3                    NX*NY*(ILAP/2),MPISIZE,MYPE+NPEY,ITAG,
     4                    MPI_COMM_WORLD,ISTATUS,IERR)
           ENDIF
        ENDIF
C----------------------------------------------------------------------
C  Horizontal communication in y-direction.
C----------------------------------------------------------------------
        IF(NPEY.GT.1) THEN
           ITAG = 300
           IF (MYPEY.EQ.0) THEN
              CALL MPI_SEND(VAR(:,NY-IY+1:NY-IY/2,:),NX*NZ*(IY/2),
     2                    MPISIZE,MYPE+1,ITAG,MPI_COMM_WORLD,IERR)
           ELSE IF (MYPEY.EQ.NPEY-1) THEN
              CALL MPI_RECV(VAR(:,1:IY/2,:),NX*NZ*(IY/2),MPISIZE,
     2                    MYPE-1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
           ELSE
C                print*, "pt 4 mype, npey", mype, npey
              CALL MPI_SENDRECV(VAR(:,NY-IY+1:NY-IY/2,:),NX*NZ*(IY/2),
     2                    MPISIZE,MYPE+1,ITAG,VAR(:,1:IY/2,:),
     3                    NX*NZ*(IY/2),MPISIZE,MYPE-1,ITAG,
     4                    MPI_COMM_WORLD,ISTATUS,IERR)
           ENDIF
C
           ITAG = 400
           IF (MYPEY.EQ.0) THEN
              CALL MPI_RECV(VAR(:,NY-IY/2+1:NY,:),NX*NZ*(IY/2),MPISIZE,
     2                    MYPE+1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
           ELSE IF (MYPEY.EQ.NPEY-1) THEN
              CALL MPI_SEND(VAR(:,IY/2+1:IY,:),NX*NZ*(IY/2),MPISIZE,
     2                    MYPE-1,ITAG,MPI_COMM_WORLD,IERR)
           ELSE
C                print*, "pt 5 mype, npey", mype, npey
              CALL MPI_SENDRECV(VAR(:,IY/2+1:IY,:),NX*NZ*(IY/2),
     2                    MPISIZE,MYPE-1,ITAG,VAR(:,NY-IY/2+1:NY,:),
     3                    NX*NZ*(IY/2),MPISIZE,MYPE+1,ITAG,
     4                    MPI_COMM_WORLD,ISTATUS,IERR)
           ENDIF
C----------------------------------------------------------------------
C  Communication for periodicity in y-direction.
C----------------------------------------------------------------------
           ITAG = 500
           IF (MYPEY.EQ.0) THEN
              CALL MPI_SEND(VAR(:,IY/2+1:IY,:),NX*NZ*(IY/2),MPISIZE,
     2                     MYPE+NPEY-1,ITAG,MPI_COMM_WORLD,IERR)
           ELSE IF (MYPEY.EQ.NPEY-1) THEN
              CALL MPI_RECV(VAR(:,NY-IY/2+1:NY,:),NX*NZ*(IY/2),
     2        MPISIZE,MYPE-NPEY+1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
           ENDIF
C
           ITAG = 600
           IF (MYPEY.EQ.0) THEN
              CALL MPI_RECV(VAR(:,1:IY/2,:),NX*NZ*(IY/2),MPISIZE,
     2                MYPE+NPEY-1,ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
           ELSE IF (MYPEY.EQ.NPEY-1) THEN
              CALL MPI_SEND(VAR(:,NY-IY+1:NY-IY/2,:),NX*NZ*(IY/2),
     2              MPISIZE,MYPE-NPEY+1,ITAG,MPI_COMM_WORLD,IERR)
           ENDIF

        ELSE
           VAR(:,1:IY/2,:)       = VAR(:,NY-IY+1:NY-IY/2,:)
           VAR(:,NY-IY/2+1:NY,:) = VAR(:,IY/2+1:IY,:)
        ENDIF
C----------------------------------------------------------------------
C  Periodicity in x-direction.
C----------------------------------------------------------------------
	VAR(1:IX/2,:,:)       = VAR(NX-IX+1:NX-IX/2,:,:)
        VAR(NX-IX/2+1:NX,:,:) = VAR(IX/2+1:IX,:,:)
C
        RETURN
        END
C**********************************************************************
        SUBROUTINE HORIZONTAL_MEAN(VARM,VAR)
C
        INCLUDE '3dmhdparam.f'
C
        include 'mpif.h'
C
        DIMENSION VAR(NX,NY,NZ),VARM(NZ),WWY(NY),WWZ(NZ),
     &            ISTATUS(MPI_STATUS_SIZE)
        DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
        DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
        DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
C
        COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     &                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
        COMMON/BOUNDS/XMAX,YMAX,ZMAX
        COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
C
        IF(NGRID.EQ.0) THEN
           DO K=1,NZ
              WWZ(K)=0.0
              DO J=IY/2+1,NY-IY/2
                 DO I=IX/2+1,NX-IX/2
                    WWZ(K)=WWZ(K)+VAR(I,J,K)
                 END DO
              END DO
              WWZ(K)=WWZ(K)/(FLOAT(NPX*NPY))
           END DO
        ELSE
           WRITE(*,*) 'Update spline interpolation'
           CALL MPI_FINALIZE(IERR)
           STOP
C----------------------------------------------------------------
C  The spline integration has to be extended to the ghost point
C  at the end of the domain in order to get the correct average.
C  For the uniform grid we use the total instead (M. Rempel)
C
c       IF (NX.EQ.IX+1) THEN
c          DO K=1,NZ
c             WWZ(K)=FSPLINEY(WYY,VAR(IX/2+1,:,K))/YMAX
c          END DO
c       ELSE
c           DO K=1,NZ
c              DO J=1+IY/2,NY-IY+1
c                WWY(J)=FSPLINEX(EXX,VAR(:,J,K))/XMAX
c             END DO
c              WWZ(K)=FSPLINEY(WYY,WWY)/YMAX
c          END DO
c        ENDIF
C----------------------------------------------------------------
        ENDIF
C
        IF(NPEY.EQ.1) THEN
           VARM=WWZ
        ELSE
           ITAG=100
           IF(MYPEY.EQ.0) THEN
              VARM=WWZ
              DO IPE=1,NPEY-1
                 CALL MPI_RECV(WWZ,NZ,MPISIZE,MYPE+IPE,
     &                    ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
                 VARM=VARM+WWZ
              END DO
           ELSE
              CALL MPI_SEND(WWZ,NZ,MPISIZE,MYPEZ*NPEY,
     &                    ITAG,MPI_COMM_WORLD,IERR)
           ENDIF
C
           ITAG=200
           IF(MYPEY.EQ.0) THEN
              DO IPE=1,NPEY-1
                 CALL MPI_SEND(VARM,NZ,MPISIZE,MYPE+IPE,
     &                    ITAG,MPI_COMM_WORLD,IERR)
              END DO
           ELSE
              CALL MPI_RECV(VARM,NZ,MPISIZE,MYPEZ*NPEY,
     &                    ITAG,MPI_COMM_WORLD,ISTATUS,IERR)
           ENDIF
        ENDIF
C
        RETURN
        END
C**********************************************************************
        SUBROUTINE BSSTEP (Y,DYDX,NV,X,HTRY,EPS,HDID,HNEXT)
C----------------------------------------------------------------------
C  Adapted from Numerical Recipes--Press et al., p. 563ff.
C  Burlisch-Stoer step with monitoring of local truncation error to
C  ensure accuracy and adjust stepsize. YSCAL=Y for constant fractional
C  error.
C----------------------------------------------------------------------
        IMPLICIT REAL*8(A-H,O-Z)
c        IMPLICIT REAL(A-H,O-Z)
C
        PARAMETER(NMAX=10,IMAX=11,NUSE=7,ONE=1.0E00,SHRINK=0.95E00,
     2          GROW=1.2E00)
C
        DIMENSION Y(NV),DYDX(NV),YERR(NMAX),YSAV(NMAX),DYSAV(NMAX),
     *          YSEQ(NMAX),NSEQ(IMAX),D(NMAX,NUSE)
C
        DATA NSEQ /2,4,6,8,12,16,24,32,48,64,96/
C
        H=HTRY   
        XSAV=X
C
        DO 10 I=1,NV
                YSAV(I)=Y(I)
                DYSAV(I)=DYDX(I)
10      CONTINUE   
1       DO 20 I=1,IMAX
                CALL MMID(YSAV,DYSAV,NV,XSAV,H,NSEQ(I),YSEQ)
                XEST=(H/NSEQ(I))**2
                CALL RZEXTR(I,XEST,YSEQ,Y,D,YERR,NV,NUSE)
                ERRMAX=0.E00
                DO 30 J=1,NV
                        ERRMAX=MAX(ERRMAX,ABS(YERR(J)/Y(J)))
30              CONTINUE
                ERRMAX=ERRMAX/EPS
                IF (ERRMAX.LT.ONE) THEN
                        X=X+H
                        HDID=H
                        IF (I.EQ.NUSE) THEN
                                HNEXT=H*SHRINK
                        ELSE IF (I.EQ.NUSE-1) THEN
                                HNEXT=H*GROW
                        ELSE
                                HNEXT=(H*NSEQ(NUSE-1))/NSEQ(I)
                        ENDIF
                        RETURN
                ENDIF
20      CONTINUE   
        H=0.25E00*H/2.0E00**((IMAX-NUSE)/2)
        IF (X+H.EQ.X) THEN
                WRITE(6,*)'BSSTEP:  Step size underflow'
		CALL MPI_FINALIZE(IERR)
                STOP
        ENDIF      
        GOTO 1
C
        END
C**********************************************************************
        SUBROUTINE MMID(Y,DYDX,NVAR,XS,HTOT,NSTP,YOUT)
C----------------------------------------------------------------------
C  Adapted from Numerical Recipes--Press et al., p. 560ff.
C  Modified midpoint method.
C----------------------------------------------------------------------
        IMPLICIT REAL*8(A-H,O-Z)
c        IMPLICIT REAL(A-H,O-Z)
C
        PARAMETER(NMAX=10)
C
        DIMENSION Y(NVAR),DYDX(NVAR),YOUT(NVAR),YM(NMAX),YN(NMAX)
C
        H=HTOT/NSTP
        DO 10 I=1,NVAR
                YM(I)=Y(I)
                YN(I)=Y(I)+H*DYDX(I)
10      CONTINUE
        X=XS+H
C
        CALL DERIVS(X,YN,NVAR,YOUT)
C
        H2=2.0E00*H
        DO 20 N=2,NSTP
                DO 30 I=1,NVAR
                        SWAP=YM(I)+H2*YOUT(I)
                        YM(I)=YN(I)
                        YN(I)=SWAP
30              CONTINUE
        X=X+H
C
        CALL DERIVS(X,YN,NVAR,YOUT)
C
20      CONTINUE
        DO 40 I=1,NVAR
                YOUT(I)=0.5E00*(YM(I)+YN(I)+H*YOUT(I))
40      CONTINUE
C
        RETURN
        END
C**********************************************************************
        SUBROUTINE RZEXTR(IEST,XEST,YEST,YZ,D,DY,NV,NUSE)
C----------------------------------------------------------------------
C  Adapted from Numerical Recipes--Press et al., p. 566f.
C  Diagonal rational function extrapolation to evaluate NV functions at
C  X=0 by fitting a diagonal rational function to a sequence of
C  estimates with progressively smaller values X=XEST, and
C  corresponding function vectors YEST.
C----------------------------------------------------------------------
        IMPLICIT REAL*8(A-H,O-Z)
c        IMPLICIT REAL(A-H,O-Z)
C
        PARAMETER(IMAX=11,NMAX=10,NCOL=7)
C
        DIMENSION X(IMAX),YEST(NV),YZ(NV),DY(NV),D(NMAX,NCOL),FX(NCOL)
	SAVE X
C
        X(IEST)=XEST
        IF (IEST.EQ.1) THEN
                DO 10 J=1,NV
                        YZ(J)=YEST(J)
                        D(J,1)=YEST(J)
                        DY(J)=YEST(J)
10              CONTINUE
        ELSE
                M1=MIN(IEST,NUSE)
                DO 20 K=1,M1-1
                        FX(K+1)=X(IEST-K)/XEST
20              CONTINUE
                DO 30 J=1,NV
                        YY=YEST(J)
                        V=D(J,1)
                        C=YY
                        D(J,1)=YY
                        DO 40 K=2,M1
                                B1=FX(K)*V
                                B=B1-C
                                IF (B.NE.0) THEN
                                        B=(C-V)/B
                                        DDY=C*B
                                        C=B1*B
                                ELSE
                                        DDY=V
                                ENDIF
                                V=D(J,K)
                                D(J,K)=DDY
                                YY=YY+DDY
40                      CONTINUE
                        DY(J)=DDY
                        YZ(J)=YY
30              CONTINUE
        ENDIF
C
        RETURN
        END
C**********************************************************************
	FUNCTION RAN2(IDUM,IIY,IIR)
C----------------------------------------------------------------------
C  Adapted from Numerical Recipes--Press et al., p. 191f. 
C  Returns a uniform random deviate between 0.0 and 1.0.
C----------------------------------------------------------------------
        IMPLICIT REAL*8(A-H,O-Z)
c        IMPLICIT REAL(A-H,O-Z)
C
	PARAMETER(M=714025,IA=1366,IC=150889,RM=1.0E00/714025.0E00)
C
	DIMENSION IIR(97)
C
	IF (IDUM.LT.0) THEN
		IDUM=MOD(IC-IDUM,M)
		DO 10 J=1,97
			IDUM=MOD(IA*IDUM+IC,M)
			IIR(J)=IDUM
10		CONTINUE
		IDUM=MOD(IA*IDUM+IC,M)
		IIY=IDUM
	ENDIF
	J=1+(97*IIY)/M
	IF ((J.GT.97).OR.(J.LT.1)) THEN
		WRITE(6,*)'RAN2:  Error.'
		CALL MPI_FINALIZE(IERR)
		STOP
	ENDIF
	IIY=IIR(J)
	RAN2=IIY*RM
	IDUM=MOD(IA*IDUM+IC,M)
	IIR(J)=IDUM
C
	RETURN
	END
C**********************************************************************
c        SUBROUTINE POTENTIAL
cC
c	INCLUDE '3dmhdparam.f'
cC
cc       INCLUDE '/opt/ctl/craylibs/craylibs/include/mpp/shmem.fh'
c        INCLUDE '/opt/ctl/craylibs_m/craylibs_m/include/mpp/mpp/shmem.fh'
cC
c        DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
c     2           ,TT(NX,NY,NZ)
c        DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
c        DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
c     2           ,FT(NX,NY,NZ)
c        DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
c     2           ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
c        DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
cC
c        DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
c        DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
cC
c        DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
c        DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
c        DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
c        DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
c     2           ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
c     3           ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
c     4           ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
cC
c     5           ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
c     6           ,SP26(IPAD)
cC
c        DIMENSION PWORK1(SHMEM_REDUCE_MIN_WRKDATA_SIZE)
c        DIMENSION ISYNC1(SHMEM_COLLECT_SYNC_SIZE)
c        DATA  ISYNC1 /SHMEM_COLLECT_SYNC_SIZE*SHMEM_SYNC_VALUE/
cC
c        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
c     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
c     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
c     4            ,SP18,WW1,SP19,WW2,SP20,WW3
cC
c     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
cC
c        COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
c     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
c        COMMON/GRID/DD,HX,H2X,HY,H2Y,HZ,H2Z,C13,C23,C43
c        COMMON/BCT/IXC,IYC,IZC,ITC,IBC
c        COMMON/COMMUN/MYPE,MPISIZE
cC----------------------------------------------------------------------
cC  Calculate div(B)
cC----------------------------------------------------------------------
c        DO 10 K=ILAP/2+1,NZ-ILAP/2
c        DO 10 J=2,NY-IY+1
c        DO 10 I=2,NX-IX+1
c                WW1(I,J,K)=(BZ(I,J,K+1)-BZ(I,J,K-1))*HZ*DZZDZ(K)
c10      CONTINUE
c        IF (IBC.EQ.0) THEN
c                IF (MYPE.EQ.0) THEN
c                        DO 20 J=2,NY-IY+1
c                        DO 20 I=2,NX-IX+1
c                                WW1(I,J,ILAP/2+1)=
c     2                                  (4.0E00*BZ(I,J,ILAP/2+2)
c     3                                         -BZ(I,J,ILAP/2+3))
c     4                                              *HZ*DZZDZ(ILAP/2+1)
c20                      CONTINUE
c                ELSE IF (MYPE.EQ.NPE-1) THEN
c                        DO 30 J=2,NY-IY+1
c                        DO 30 I=2,NX-IX+1
c                                WW1(I,J,NZ-ILAP/2)=
c     2                                  (BZ(I,J,NZ-ILAP/2-2)
c     3                                     -4.0E00*BZ(I,J,NZ-ILAP/2-1))
c     4                                             *HZ*DZZDZ(NZ-ILAP/2)
c30                      CONTINUE
c                ENDIF
c        ENDIF
c        DO 40 K=ILAP/2+1,NZ-ILAP/2
c        DO 40 J=2,NY-IY+1
c        DO 40 I=2,NX-IX+1
c                WW1(I,J,K)=WW1(I,J,K)
c     2                          +(BX(I+1,J,K)-BX(I-1,J,K))*HX*DXXDX(I)
c40      CONTINUE
c        DO 50 K=ILAP/2+1,NZ-ILAP/2
c        DO 50 J=2,NY-IY+1
c        DO 50 I=2,NX-IX+1
c                WW1(I,J,K)=WW1(I,J,K)
c     2                          +(BY(I,J+1,K)-BY(I,J-1,K))*HY*DYYDY(J)
c50      CONTINUE
cC----------------------------------------------------------------------
cC  Initial guess for potential
cC----------------------------------------------------------------------
c        DO 60 K=1,NZ
c        DO 60 J=1,NY
c        DO 60 I=1,NX
c                WW2(I,J,K)=0.0E00
c                WW3(I,J,K)=0.0E00
c60      CONTINUE
cC---------------------------------------------------------------------
cC  Solve for potential.
cC---------------------------------------------------------------------
c        ITMAX=1000
c        DTT=DD*DD/6.0E00
c        ERR=1.0E-10
c	RADC=0.99E00
c	SOR=1.0E00
c        DO 500 N=1,ITMAX
c            DO 501 K=ILAP/2+1,NZ-ILAP/2
c            DO 501 J=2,NY-IY+1
c            DO 501 I=2,NX-IX+1
c	    IF (MOD(I+J+K,2).EQ.MOD(N,2)) THEN
c                WW3(I,J,K)=(WW2(I,J,K+1)-WW2(I,J,K-1))*HZ*D2ZZDZ2(K)
c     2                  +(WW2(I,J,K+1)-2.0E00*WW2(I,J,K)+WW2(I,J,K-1))
c     3                                          *H2Z*DZZDZ(K)*DZZDZ(K)
c	    ENDIF
c501         CONTINUE
c            IF (MYPE.EQ.0) THEN
c                DO 502 J=2,NY-IY+1
c                DO 502 I=2,NX-IX+1
c		IF (MOD(I+J+ILAP/2+1,2).EQ.MOD(N,2)) THEN
c                    WW3(I,J,ILAP/2+1)=(2.0E00*WW2(I,J,ILAP/2+1)
c     2                                -5.0E00*WW2(I,J,ILAP/2+2)
c     3                                +4.0E00*WW2(I,J,ILAP/2+3)
c     4                                -WW2(I,J,ILAP/2+4))
c     5                                          *H2Z*DZZDZ(K)*DZZDZ(K)
c		ENDIF
c502             CONTINUE
c            ELSE IF (MYPE.EQ.NPE-1) THEN
c                DO 503 J=2,NY-IY+1
c                DO 503 I=2,NX-IX+1
c		IF (MOD(I+J+NZ-ILAP/2,2).EQ.MOD(N,2)) THEN
c                    WW3(I,J,NZ-ILAP/2)=(2.0E00*WW2(I,J,NZ-ILAP/2)
c     2                                 -5.0E00*WW2(I,J,NZ-ILAP/2-1)
c     3                                 +4.0E00*WW2(I,J,NZ-ILAP/2-2)
c     4                                 -WW2(I,J,NZ-ILAP/2-3))
c     5                                          *H2Z*DZZDZ(K)*DZZDZ(K)
c		ENDIF
c503             CONTINUE
c            ENDIF
c            DO 510 K=ILAP/2+1,NZ-ILAP/2
c            DO 510 J=2,NY-IY+1
c            DO 510 I=2,NX-IX+1
c	    IF (MOD(I+J+K,2).EQ.MOD(N,2)) THEN
c                WW3(I,J,K)=WW3(I,J,K)
c     2                  +(WW2(I+1,J,K)-WW2(I-1,J,K))*HX*D2XXDX2(I)
c     3                  +(WW2(I+1,J,K)-2.0E00*WW2(I,J,K)+WW2(I-1,J,K))
c     4                                          *H2X*DXXDX(I)*DXXDX(I)
c	    ENDIF
c510         CONTINUE
c            DO 520 K=ILAP/2+1,NZ-ILAP/2
c            DO 520 J=2,NY-IY+1
c            DO 520 I=2,NX-IX+1
c	    IF (MOD(I+J+K,2).EQ.MOD(N,2)) THEN
c                WW3(I,J,K)=WW3(I,J,K)
c     2                  +(WW2(I,J+1,K)-WW2(I,J-1,K))*HY*D2YYDY2(J)
c     3                  +(WW2(I,J+1,K)-2.0E00*WW2(I,J,K)+WW2(I,J-1,K))
c     4                                          *H2Y*DYYDY(J)*DYYDY(J)
c	    ENDIF
c520         CONTINUE
c            DO 530 K=ILAP/2+1,NZ-ILAP/2
c            DO 530 J=2,NY-IY+1
c            DO 530 I=2,NX-IX+1
c	    IF (MOD(I+J+K,2).EQ.MOD(N,2)) THEN
c                WW3(I,J,K)=WW3(I,J,K)+WW1(I,J,K)
c	    ENDIF
c530         CONTINUE
c            DO 540 K=ILAP/2+1,NZ-ILAP/2
c            DO 540 J=2,NY-IY+1
c            DO 540 I=2,NX-IX+1
c	    IF (MOD(I+J+K,2).EQ.MOD(N,2)) THEN
c                WW2(I,J,K)=WW2(I,J,K)+SOR*DTT*WW3(I,J,K)
c	    ENDIF
c540         CONTINUE
cC
c	    IF (N.EQ.1) THEN
c		SOR=1.0E00/(1.0E00-0.5E00*RADC**2)
c	    ELSE
c		SOR=1.0E00/(1.0E00-0.25E00*RADC**2*SOR)
c	    ENDIF
cC----------------------------------------------------------------------
cC  Set dphidz=0 on upper and lower boundary.
cC----------------------------------------------------------------------
c            IF (MYPE.EQ.0) THEN
c                DO 541 J=2,NY-IY+1
c                DO 541 I=2,NX-IX+1
c		IF (MOD(I+J+ILAP/2+1,2).EQ.MOD(N,2)) THEN
c                        WW2(I,J,ILAP/2+1)=C43*WW2(I,J,ILAP/2+2)
c     2                                          -C13*WW2(I,J,ILAP/2+3)
c		ENDIF
c541             CONTINUE
c            ELSE IF (MYPE.EQ.NPE-1) THEN
c                DO 542 J=2,NY-IY+1
c                DO 542 I=2,NX-IX+1
c		IF (MOD(I+J+NZ-ILAP/2,2).EQ.MOD(N,2)) THEN
c                        WW2(I,J,NZ-ILAP/2)=C43*WW2(I,J,NZ-ILAP/2-1)
c     2                                        -C13*WW2(I,J,NZ-ILAP/2-2)
c		ENDIF
c542             CONTINUE
c            ENDIF
cC----------------------------------------------------------------------
cC  Communicate.
cC----------------------------------------------------------------------
c            DO 550 K=ILAP/2+1,NZ-ILAP/2
c            DO 550 J=2,NY-IY+1
c                WW2(1,J,K)=WW2(NX-IX+1,J,K)
c                WW2(NX-IX+2,J,K)=WW2(2,J,K)
c550         CONTINUE
c            DO 560 K=ILAP/2+1,NZ-ILAP/2
c            DO 560 I=1,NX
c                WW2(I,1,K)=WW2(I,NY-IY+1,K)
c                WW2(I,NY-IY+2,K)=WW2(I,2,K)
c560         CONTINUE
c            CALL BARRIER()
c            IF (MYPE.EQ.0) THEN
c                DO 15 K=1,ILAP/2
c                WW2(:,:,K)=WW2(:,:,ILAP/2+1)
c                CALL SHMEM_GGG(WW2(1,1,NZ-ILAP/2+K),WW2(1,1,ILAP/2+K)
c     2                                                      ,NX*NY,1)
c15              CONTINUE
c            ELSE IF (MYPE.EQ.NPE-1) THEN
c                DO 25 K=1,ILAP/2
c                CALL SHMEM_GGG(WW2(1,1,K),WW2(1,1,NZ-ILAP+K)
c     2							,NX*NY,NPE-2)
c                WW2(:,:,NZ-ILAP/2+K)=WW2(:,:,NZ-ILAP/2)
c25              CONTINUE
c            ELSE
c                DO 35 K=1,ILAP/2
c                CALL SHMEM_GGG(WW2(1,1,K),WW2(1,1,NZ-ILAP+K)
c     2							,NX*NY,MYPE-1)
c                CALL SHMEM_GGG(WW2(1,1,NZ-ILAP/2+K),WW2(1,1,ILAP/2+K)
c     2                                                  ,NX*NY,MYPE+1)
c35              CONTINUE
c            ENDIF
c            CALL BARRIER()
cC
c            RSID=0.0E00
c            DO 535 K=ILAP/2+1,NZ-ILAP/2
c            DO 535 J=2,NY-IY+1
c            DO 535 I=2,NX-IX+1
c                RSID=MAX(RSID,ABS(WW3(I,J,K)))
c535         CONTINUE
c            CALL BARRIER()
c            CALL SHMEM_REAL8_MAX_TO_ALL(RSID,RSID,1,0,0,NPE,PWORK1,ISYNC1)
c            CALL BARRIER()
c            IF (RSID.LT.ERR) GOTO 1000
c	    if (mype.eq.0) then
c            	write(6,*)n,radc,sor,rsid
c	    endif
cC
c500     CONTINUE
cC
c1000    CONTINUE
cC
c        RETURN
c        END
C**********************************************************************
C integration routine grabbed from the Numerical Recipies
C I have modified them: changed from REAL to REAL*8	
C I changed also EPS=1.e-6 into EPS=1.e-9	
C**********************************************************************	
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
c      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-9, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
c      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.0d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software );.
      
      
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
c      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
c      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software );.
      
      
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
c      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
c      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software );.
C**********************************************************************
      SUBROUTINE INIT_SPLINEX(EX)
C----------------------------------------------------------------------
C  Adapted from the IDL int_tabulated function.  Integrates a tabulated
C  set of data WY(EX) on interval min(EX) to max(EX).
C----------------------------------------------------------------------
      	INCLUDE '3dmhdparam.f'
C
      	DIMENSION EX(NPX),Y2(NPX),EX2(NPX+3)
 	DIMENSION KLOX(NPX+3),KHIX(NPX+3),HHX(NPX+3),SIGX(NPX),
     2		  AAX(NPX+3),BBX(NPX+3),XPINV(NPX)
C
	COMMON/SPLINEX/KLOX,KHIX,HHX,SIGX,AAX,BBX,XPINV,XHH,ISEGX
C
      	ISEGX=NPX-1
      	DO WHILE (MOD(ISEGX,4).NE.0)
         	ISEGX=ISEGX+1
      	END DO
      	XMIN=MINVAL(EX)
      	XMAX=MAXVAL(EX)
   	XHH=(XMAX-XMIN)/FLOAT(ISEGX)
C
      	DO 10 I=0,ISEGX
         	EX2(I+1)=XHH*FLOAT(I)+XMIN
10	CONTINUE
C----------------------------------------------------------------------
C  The following from numerical recipes SPLINE (natural boundaries).
C----------------------------------------------------------------------
      	Y2(1)=0.0E00
      	DO 20 I=2,NPX-1
         	SIGX(I)=(EX(I)-EX(I-1))/(EX(I+1)-EX(I-1))
         	XPINV(I)=1.0E00/(SIGX(I)*Y2(I-1)+2.0E00)
         	Y2(I)=(SIGX(I)-1.0E00)*XPINV(I)
20	CONTINUE
      	DO 30 I=1,ISEGX+1
         	X=EX2(I)
         	KLOX(I)=1
         	KHIX(I)=NPX
         	DO WHILE (KHIX(I)-KLOX(I).GT.1)
            		K=(KHIX(I)+KLOX(I))/2
            		IF (EX(K).GT.X) THEN
               			KHIX(I)=K
            		ELSE
               			KLOX(I)=K
            		ENDIF
         	ENDDO
         	HHX(I)=EX(KHIX(I))-EX(KLOX(I))
         	IF (HHX(I).EQ.0) THEN
            		WRITE(6,*)'INIT_SPLINEX: 
     2					EX values must be distinct.',I
			CALL MPI_FINALIZE(IERR)
            		STOP
         	ENDIF
         	AAX(I)=(EX(KHIX(I))-X)/HHX(I)
         	BBX(I)=(X-EX(KLOX(I)))/HHX(I)
30	CONTINUE
C
      	RETURN
      	END
C**********************************************************************
      SUBROUTINE INIT_SPLINEY(EY)
C----------------------------------------------------------------------
C  Adapted from the IDL int_tabulated function.  Integrates a tabulated
C  set of data WY(EY) on interval min(EY) to max(EY).
C----------------------------------------------------------------------
        INCLUDE '3dmhdparam.f'
C
        DIMENSION EY(NRY),Y2(NRY),EY2(NRY+3)
        DIMENSION KLOY(NRY+3),KHIY(NRY+3),HHY(NRY+3),SIGY(NRY),
     2            AAY(NRY+3),BBY(NRY+3),YPINV(NRY)
C
        COMMON/SPLINEY/KLOY,KHIY,HHY,SIGY,AAY,BBY,YPINV,YHH,ISEGY
C
        ISEGY=NRY-1
        DO WHILE (MOD(ISEGY,4).NE.0)
                ISEGY=ISEGY+1
        END DO
        YMIN=MINVAL(EY)
        YMAX=MAXVAL(EY)
        YHH=(YMAX-YMIN)/FLOAT(ISEGY)
C
        DO 10 I=0,ISEGY
                EY2(I+1)=YHH*FLOAT(I)+YMIN
10      CONTINUE
C----------------------------------------------------------------------
C  The following from numerical recipes SPLINE (natural boundaries).
C----------------------------------------------------------------------
        Y2(1)=0.0E00
        DO 20 I=2,NRY-1
                SIGY(I)=(EY(I)-EY(I-1))/(EY(I+1)-EY(I-1))
                YPINV(I)=1.0E00/(SIGY(I)*Y2(I-1)+2.0E00)
                Y2(I)=(SIGY(I)-1.0E00)*YPINV(I)
20      CONTINUE
        DO 30 I=1,ISEGY+1
                Y=EY2(I)
                KLOY(I)=1
                KHIY(I)=NRY
                DO WHILE (KHIY(I)-KLOY(I).GT.1)
                        K=(KHIY(I)+KLOY(I))/2
                        IF (EY(K).GT.Y) THEN
                                KHIY(I)=K
                        ELSE
                                KLOY(I)=K
                        ENDIF
                ENDDO
                HHY(I)=EY(KHIY(I))-EY(KLOY(I))
                IF (HHY(I).EQ.0) THEN
                        WRITE(6,*)'INIT_SPLINEY: 
     2					EY values must be distinct.',I
                        CALL MPI_FINALIZE(IERR)
                        STOP
                ENDIF
                AAY(I)=(EY(KHIY(I))-Y)/HHY(I)
                BBY(I)=(Y-EY(KLOY(I)))/HHY(I)
30      CONTINUE
C
        RETURN
        END
C**********************************************************************
      	FUNCTION FSPLINEX(EX,WY)
C----------------------------------------------------------------------
C  Adapted from the IDL int_tabulated function.  Integrates a tabulated
C  set of data WY(EX) on interval min(EX) to max(EX).
C----------------------------------------------------------------------
        INCLUDE '3dmhdparam.f'
C
        DIMENSION EX(NX),Y2(NPX),WY(NX),WY2(NPX),U(NPX)
C----------------------------------------------------------------------
C  Whole EX and WY passed in to save copy, so either dimension EX and WY
C  down by IX-1 (as Jim Edwards (IBM) did REAL*8 EX(0:NPX-1), WY(0:NPX-1)),
C  or add IX-1 to indexing of EX and WY.
C----------------------------------------------------------------------
        DIMENSION KLOX(NPX+3),KHIX(NPX+3),HHX(NPX+3),SIGX(NPX),
     2            AAX(NPX+3),BBX(NPX+3),XPINV(NPX)
C
        COMMON/SPLINEX/KLOX,KHIX,HHX,SIGX,AAX,BBX,XPINV,XHH,ISEGX
C----------------------------------------------------------------------
C  The following from numerical recipes SPLINE (natural boundaries).
C----------------------------------------------------------------------
      	Y2(1)=0.0E00
      	U(1)=0.0E00
      	DO 10 I=2,NPX-1
        	Y2(I)=(SIGX(I)-1.0E00)*XPINV(I)
         	U(I)=(6.0E00*((WY(I+1+IX-1)-WY(I+IX-1))
     2					     /(EX(I+1+IX-1)-EX(I+IX-1))
     3        		     -(WY(I+IX-1)-WY(I-1+IX-1))
     4					     /(EX(I+IX-1)-EX(I-1+IX-1)))
     5        	    /(EX(I+1+IX-1)-EX(I-1+IX-1))-SIGX(I)*U(I-1))*XPINV(I)
10	CONTINUE
      	Y2(NPX)=0.0
      	DO 20 K=NPX-1,1,-1
         	Y2(K)=Y2(K)*Y2(K+1)+U(K)
20	CONTINUE
      	DO 30 I=1,ISEGX+1
             WY2(I)=AAX(I)*WY(KLOX(I)+IX-1)+BBX(I)*WY(KHIX(I)+IX-1)
     2             +((AAX(I)**3-AAX(I))*Y2(KLOX(I))
     3		    +(BBX(I)**3-BBX(I))*Y2(KHIX(I)))*(HHX(I)**2)/6.0E00
30	CONTINUE
C
      	FSPLINEX=0.0E00
      	DO 40 I=5,ISEGX+1,4
      		FSPLINEX=FSPLINEX+2.0E00*XHH*(7.0E00*(WY2(I-4)+WY2(I))
     2        		       		   +32.0E00*(WY2(I-3)+WY2(I-1))
     3        		       		   +12.0E00*WY2(I-2))/45.0E00
40	CONTINUE
C
      RETURN
      END
C**********************************************************************
        FUNCTION FSPLINEY(EY,WY)
C----------------------------------------------------------------------
C  Adapted from the IDL int_tabulated function.  Integrates a tabulated
C  set of data WY(EY) on interval min(EY) to max(EY).
C----------------------------------------------------------------------
        INCLUDE '3dmhdparam.f'
C
        DIMENSION EY(NY),Y2(NRY),WY(NY),WY2(NRY),U(NRY)
C----------------------------------------------------------------------
C  Whole EY and WY passed in to save copy, so either dimension EY and WY
C  down by IY-1 (as Jim Edwards (IBM) did REAL*8 EY(0:NRY-1), WY(0:NRY-1)),
C  or add IY-1 to indexing of EY and WY.
C----------------------------------------------------------------------
        DIMENSION KLOY(NRY+3),KHIY(NRY+3),HHY(NRY+3),SIGY(NRY),
     2            AAY(NRY+3),BBY(NRY+3),YPINV(NRY)
C
        COMMON/SPLINEY/KLOY,KHIY,HHY,SIGY,AAY,BBY,YPINV,YHH,ISEGY
C----------------------------------------------------------------------
C  The following from numerical recipes SPLINE (natural boundaries).
C----------------------------------------------------------------------
        Y2(1)=0.0E00
        U(1)=0.0E00
        DO 10 I=2,NRY-1
                Y2(I)=(SIGY(I)-1.0E00)*YPINV(I)
                U(I)=(6.0E00*((WY(I+1+IY-1)-WY(I+IY-1))
     2                                       /(EY(I+1+IY-1)-EY(I+IY-1))
     3                       -(WY(I+IY-1)-WY(I-1+IY-1))
     4                                       /(EY(I+IY-1)-EY(I-1+IY-1)))
     5             /(EY(I+1+IY-1)-EY(I-1+IY-1))-SIGY(I)*U(I-1))*YPINV(I)
10      CONTINUE
        Y2(NRY)=0.0
        DO 20 K=NRY-1,1,-1
                Y2(K)=Y2(K)*Y2(K+1)+U(K)
20      CONTINUE
        DO 30 I=1,ISEGY+1
             WY2(I)=AAY(I)*WY(KLOY(I)+IY-1)+BBY(I)*WY(KHIY(I)+IY-1)
     2             +((AAY(I)**3-AAY(I))*Y2(KLOY(I))
     3              +(BBY(I)**3-BBY(I))*Y2(KHIY(I)))*(HHY(I)**2)/6.0E00
30      CONTINUE
C
        FSPLINEY=0.0E00
        DO 40 I=5,ISEGY+1,4
                FSPLINEY=FSPLINEY+2.0E00*YHH*(7.0E00*(WY2(I-4)+WY2(I))
     2                                     +32.0E00*(WY2(I-3)+WY2(I-1))
     3                                     +12.0E00*WY2(I-2))/45.0E00
40      CONTINUE
C
      RETURN
      END
C**********************************************************************
        SUBROUTINE TUBE
C----------------------------------------------------------------------
C It creates BX, BY, BZ and modifies the thermodynamqic variables
C adequatelly. It is called from STATIC. The calling process is as 
C follows: IF LMAG=.FALSE.: no magnetic field is created.
C          IF LMAG=.TRUE. : if AMPB=/=0: a horizontal layer of magnetic
C                                  		  field is constructed
C                     	    if AMPB= 0: this routine TUBE is called
C----------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
C *** dimensions declarations for /BIG/ ***
        DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2		 ,TT(NX,NY,NZ)
        DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
        DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2		 ,FT(NX,NY,NZ)
        DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ)
     2		 ,ZRO(NX,NY,NZ),ZTT(NX,NY,NZ)
        DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
C
	DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
        DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
C *** dimension declarations for /AJACOBI/ *** 
	DIMENSION EXX(NX),DXXDX(NX),D2XXDX2(NX),DDX(NX)
        DIMENSION WYY(NY),DYYDY(NY),D2YYDY2(NY),DDY(NY)
        DIMENSION ZEE(NZ),DZZDZ(NZ),D2ZZDZ2(NZ),DDZ(NZ)
C *** dimension declarations for the padding in /BIG/ ***        
        DIMENSION SP1(IPAD),SP2(IPAD),SP3(IPAD),SP4(IPAD),SP5(IPAD)
     2           ,SP6(IPAD),SP7(IPAD),SP8(IPAD),SP9(IPAD),SP10(IPAD)
     3           ,SP11(IPAD),SP12(IPAD),SP13(IPAD),SP14(IPAD),SP15(IPAD)
     4           ,SP16(IPAD),SP17(IPAD),SP18(IPAD),SP19(IPAD),SP20(IPAD)
C
     5           ,SP21(IPAD),SP22(IPAD),SP23(IPAD),SP24(IPAD),SP25(IPAD)
     6           ,SP26(IPAD)
c sleak try to fix mismatch common block size warning:
	DIMENSION RKAPA(NZ),DKAPA(NZ)
C
C *** common block /BIG/ contains everything ***
        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4		  ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C *** common block /AJACOBI/ contains everything ***
	COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2                ,ZEE,DZZDZ,D2ZZDZ2,DDZ
	COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
	COMMON/CMAG/ORM,RM,OBETA,AMPB,BFH,BZP
        COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
        COMMON/BOUNDS/XMAX,YMAX,ZMAX
C *** tube related stuff ***       
        COMMON /FFCOM/ A,C  
        COMMON /NUMBERS/ ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX
        EXTERNAL FF
        REAL*8 LAMBDA
c        REAL LAMBDA
C *** CASE=2,4, and 5  related stuff
        REAL*8, DIMENSION(NTUBES) :: XN,ZN,RRN,EXPI,HN,QN,BYN,AN
c        REAL, DIMENSION(NTUBES) :: XN,ZN,RRN,EXPI,HN,QN,BYN,AN
C----------------------------------------------------------------------
C Constants 
C----------------------------------------------------------------------
        ZERO  = 0.0E00
        ONE   = 1.0E00
        TWO   = ONE+ONE
        THREE = ONE+TWO
        FOUR  = ONE+THREE
        FIVE  = ONE+FOUR
        SIX   = ONE+FIVE
        PI = 2.00E00*ASIN(1.00E00)
C----------------------------------------------------------------------
C INITIALIZE VARIABLES 
C----------------------------------------------------------------------
	IF (LMAG) THEN
        	BX=ZERO
        	BY=ZERO
        	BZ=ZERO
	ENDIF
        OGAMMA=ONE/GAMMA
        CALL SETUP(FINP,FOUT,IPAR,PAR)
        NTUBE = IPAR(16)
C**********************************************************************
C Selection of different models of tubes
C**********************************************************************
        SELECT CASE (NTUBE)
C**********************************************************************
C MODEL 0 : Model 0:  No tube, no layer.
C**********************************************************************
        CASE (0)
	  IF (LMAG) THEN
	  	BX = ZERO
          	BY = ZERO
          	BZ = ZERO
	  ENDIF
          RU = ZERO
          RV = ZERO
          RW = ZERO
C**********************************************************************
C MODEL 1 : 2D tube in equilibrium apart for the buoyancy force
C**********************************************************************
        CASE (1) 
C------------------
C Tube's Parameters
C------------------
          XCENT = PAR(34)
          ZCENT = PAR(35)
          R_MAX = PAR(36)
          C_MT  = PAR(37)
          A     = PAR(38)
          C     = EXP(-R_MAX*R_MAX)
C-----------------------------------------------------------
C  NOTE:  The following line assumes an adiabatically
C  stratified domain (entropy constant throughout).
C-----------------------------------------------------------
          ROGAPR= RO(3,2,3)**GAMMA/(RO(3,2,3)*TT(3,2,3))
C------------------------------------------------------------
C we calculate  2d tube in the plan (x,z), cut at JCUT=IY/2+1
C------------------------------------------------------------
          JCUT=IY/2+1
          I1 = 2
          I2 = NX-IX+1
          J1 = 2 
          J2 = NY-IY+1
          K1 = ILAP/2+1
          K2 = NZ-ILAP/2
          DO K=K1,K2
            DO I=I1,I2
              XNEW = EXX(I)-XCENT
              ZNEW = ZEE(K)-ZCENT
              RNEW = SQRT(XNEW**TWO + ZNEW**TWO) 
C-----------------------------
C Calculate the magnetic field.  NOTE: OBETA = TWO / BETA
C-----------------------------
              IF (RNEW.LT.R_MAX) THEN
                BY(I,JCUT,K) =  (EXP(-RNEW**TWO)-C) / (ONE-C)
                BPHI =  BY(I,JCUT,K) * C_MT * A*RNEW**THREE /
     &               (A*RNEW**THREE + ONE) 
                IF (RNEW.NE.ZERO) THEN
                  BX(I,JCUT,K) =  BPHI*ZNEW/RNEW
                  BZ(I,JCUT,K) = -BPHI*XNEW/RNEW
                ENDIF
                CALL QROMB(FF,RNEW,R_MAX,DPR_TOT)
                DPR = OBETA/(ONE-C)**TWO*DPR_TOT -
     &               (BX(I,JCUT,K)**TWO + BY(I,JCUT,K)**TWO
     &               + BZ(I,JCUT,K)**TWO)*OBETA/TWO 
                PRE = RO(I,JCUT,K)*TT(I,JCUT,K)+DPR
                RO(I,JCUT,K) = (ROGAPR * PRE)**OGAMMA
                TT(I,JCUT,K) = PRE/RO(I,JCUT,K)
              ENDIF
            ENDDO
          ENDDO
          RO(I1:I2,J1:J2,K1:K2)=SPREAD(RO(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          TT(I1:I2,J1:J2,K1:K2)=SPREAD(TT(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BX(I1:I2,J1:J2,K1:K2)=SPREAD(BX(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BY(I1:I2,J1:J2,K1:K2)=SPREAD(BY(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BZ(I1:I2,J1:J2,K1:K2)=SPREAD(BZ(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          RU = ZERO
          RV = ZERO
          RW = ZERO
C**********************************************************************
C MODEL 2 : 3D tube in equilibrium with drho=0
C**********************************************************************
        CASE (2) 
C------------------
C Tube's Parameters
C------------------
          XCENT = PAR(34)
          ZCENT = PAR(35)
          R_MAX = PAR(36)
          C_MT  = PAR(37)
          A     = PAR(38)
          C     = EXP(-R_MAX*R_MAX)
          LAMBDA= PAR(51)
          VZ0   = PAR(52)
C        
C------------------------------------------------------------
C we calculate  2d tube in the plan (x,z), cut at JCUT=IY/2+1
C------------------------------------------------------------
          RW=ZERO
C          
          JCUT=IY/2+1
          I1 = 2
          I2 = NX-IX+1
          J1 = 2 
          J2 = NY-IY+1
          K1 = ILAP/2+1
          K2 = NZ-ILAP/2
          DO K=K1,K2
            DO I=I1,I2
              XNEW = EXX(I)-XCENT
              ZNEW = ZEE(K)-ZCENT
              RNEW = SQRT(XNEW**TWO + ZNEW**TWO) 
C-----------------------------
C Calculate the magnetic field.  NOTE: OBETA = TWO / BETA
C-----------------------------
              IF (RNEW.LT.R_MAX) THEN
                BY(I,JCUT,K) =  (EXP(-RNEW**TWO)-C) / (ONE-C)
                BPHI =  BY(I,JCUT,K) * C_MT * A*RNEW**THREE /
     &               (A*RNEW**THREE + ONE) 
                IF (RNEW.NE.ZERO) THEN
                  BX(I,JCUT,K) =  BPHI*ZNEW/RNEW
                  BZ(I,JCUT,K) = -BPHI*XNEW/RNEW
                ENDIF
                CALL QROMB(FF,RNEW,R_MAX,DPR_TOT)
                DPR = OBETA/(ONE-C)**TWO*DPR_TOT -
     &               (BX(I,JCUT,K)**TWO + BY(I,JCUT,K)**TWO
     &               + BZ(I,JCUT,K)**TWO)*OBETA/TWO 
                PRE = RO(I,JCUT,K)*TT(I,JCUT,K)+DPR
                TT(I,JCUT,K) = PRE/RO(I,JCUT,K)
C-----note: RW is set to a negative!!! value because z is the depth!!!
                RW(I,JCUT,K) = -VZ0*RO(I,JCUT,K)
              ENDIF
            ENDDO
          ENDDO
          TT(I1:I2,J1:J2,K1:K2)=SPREAD(TT(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BX(I1:I2,J1:J2,K1:K2)=SPREAD(BX(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BY(I1:I2,J1:J2,K1:K2)=SPREAD(BY(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BZ(I1:I2,J1:J2,K1:K2)=SPREAD(BZ(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          RU = ZERO
          RV = ZERO
          RW(I1:I2,J1:J2,K1:K2)=SPREAD(RW(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          DO J=J1,J2
            RW(I1:I2,J,K1:K2)=RW(I1:I2,J,K1:K2)*
     &                        SIN(2.0E00*PI*WYY(J)/LAMBDA-PI/2.0E00)
          ENDDO
C**********************************************************************
C MODEL 3 : Same as MODEL 2 but 2D
C**********************************************************************
        CASE (3) 
C------------------
C Tube's Parameters
C------------------
          XCENT = PAR(34)
          ZCENT = PAR(35)
          R_MAX = PAR(36)
          C_MT  = PAR(37)
          A     = PAR(38)
          C     = EXP(-R_MAX*R_MAX)
          LAMBDA= PAR(51)
          VZ0   = PAR(52)
C        
C------------------------------------------------------------
C we calculate  2d tube in the plan (x,z), cut at JCUT=IY/2+1
C------------------------------------------------------------
          RW=ZERO
C          
          JCUT=IY/2+1
          I1 = 2
          I2 = NX-IX+1
          J1 = 2 
          J2 = NY-IY+1
          K1 = ILAP/2+1
          K2 = NZ-ILAP/2
          DO K=K1,K2
            DO I=I1,I2
              XNEW = EXX(I)-XCENT
              ZNEW = ZEE(K)-ZCENT
              RNEW = SQRT(XNEW**TWO + ZNEW**TWO) 
C-----------------------------
C Calculate the magnetic field.  NOTE: OBETA = TWO / BETA 
C-----------------------------
              IF (RNEW.LT.R_MAX) THEN
                BY(I,JCUT,K) =  (EXP(-RNEW**TWO)-C) / (ONE-C)
                BPHI =  BY(I,JCUT,K) * C_MT * A*RNEW**THREE /
     &               (A*RNEW**THREE + ONE) 
                IF (RNEW.NE.ZERO) THEN
                  BX(I,JCUT,K) =  BPHI*ZNEW/RNEW
                  BZ(I,JCUT,K) = -BPHI*XNEW/RNEW
                ENDIF
                CALL QROMB(FF,RNEW,R_MAX,DPR_TOT)
                DPR = OBETA/(ONE-C)**TWO*DPR_TOT -
     &               (BX(I,JCUT,K)**TWO + BY(I,JCUT,K)**TWO
     &               + BZ(I,JCUT,K)**TWO)*OBETA/TWO 
                PRE = RO(I,JCUT,K)*TT(I,JCUT,K)+DPR
                TT(I,JCUT,K) = PRE/RO(I,JCUT,K)
C-----note: RW is set to a negative!!! value because z is the depth!!!
                RW(I,JCUT,K) = -VZ0*RO(I,JCUT,K)
              ENDIF
            ENDDO
          ENDDO
          TT(I1:I2,J1:J2,K1:K2)=SPREAD(TT(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BX(I1:I2,J1:J2,K1:K2)=SPREAD(BX(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BY(I1:I2,J1:J2,K1:K2)=SPREAD(BY(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          BZ(I1:I2,J1:J2,K1:K2)=SPREAD(BZ(I1:I2,JCUT,K1:K2),2,J2-J1+1)
          RU = ZERO
          RV = ZERO
          RW(I1:I2,J1:J2,K1:K2)=SPREAD(RW(I1:I2,JCUT,K1:K2),2,J2-J1+1)
C**********************************************************************
C MODEL 4 : Layer with constant magnetic field
C**********************************************************************
        CASE (4)
C------------------
C Tube's Parameters
C------------------
          ZCENT = PAR(35) !center of the layer
          RMAX  = PAR(36) !half width of the layer
          CMT   = PAR(37) !may be used later for shearing the layer
          A     = PAR(38) !steepness of the border of layer
C------------------------------------------------------------
	JCUT=IY/2+1
        I1 = 2
        I2 = NX-IX+1
        J1 = 2
        J2 = NY-IY+1
        K1 = ILAP/2+1
        K2 = NZ-ILAP/2
        DO K=K1,K2
        DO I=I1,I2
        	ZNEW = ZEE(K)-ZCENT
              	RNEW = SQRT(ZNEW**TWO)
C-----------------------------
C Calculate the magnetic field
C N.B. OBETA = TWO / BETA !!!!
C-----------------------------
		BY(I,JCUT,K) = (1.00E00-TANH((RNEW-RMAX)/A))/2.0E00
		PRE = RO(I,JCUT,K)*TT(I,JCUT,K) - BY(I,JCUT,K)**TWO * OBETA/TWO
		RO(I,JCUT,K) = PRE/TT(I,JCUT,K)
	ENDDO
	ENDDO
	RO(I1:I2,J1:J2,K1:K2)=SPREAD(RO(I1:I2,JCUT,K1:K2),2,J2-J1+1)
        BY(I1:I2,J1:J2,K1:K2)=SPREAD(BY(I1:I2,JCUT,K1:K2),2,J2-J1+1)
        RU=ZERO
        RV=ZERO
        RW=ZERO
C**********************************************************************
C MODEL 5 : Layer of tubes with the same mag. flux as model=4
C           All the tubes identical, the same values of HH.
C           Maximum B=1.0, field strength determined by BETA.
C**********************************************************************
        CASE (5)
C----------------------------------------------------------------------
C  Tube parameters.
C----------------------------------------------------------------------
        ZCENT = PAR(35) !center of the layer
        RMAX  = PAR(36) !half width of the layer
        CMT   = PAR(37) != q the pitch = BT/(R*BY) = constant
        A     = PAR(38) !width of layer border taper
C
        HH    = PAR(55) !individual tube widths
C----------------------------------------------------------------------
C  Tube positions (ISTAG=0 tube rows aligned, ELSE tubes staggered).
C----------------------------------------------------------------------
        ISTAG=1
        N=1
        DCOL=XMAX/FLOAT(NCOL)
        DROW=TWO*RMAX/FLOAT(NROW)
        DO 10 K=1,NROW
        DO 10 I=1,NCOL
                IF (ISTAG.EQ.0) THEN
                        XN(N)=FLOAT(I)*DCOL-DCOL/TWO
                ELSE
                        XN(N)=FLOAT(I)*DCOL-FLOAT(K)*DCOL/TWO
                        IF (XN(N).LT.ZERO) THEN
                                XN(N)=XMAX+XN(N)
                        ENDIF
                ENDIF
                ZN(N)=ZCENT-RMAX-DROW/TWO+FLOAT(K)*DROW
                N=N+1
10      CONTINUE
C---------------------------------------------------------------------
C  Calculate magnetic field.
C---------------------------------------------------------------------
        DEN=-0.5E00/HH**2
C------------------------------------------------------------
C Calculate  2d tube in the plan (x,z), cut at JCUT=IY/2+1
C------------------------------------------------------------
	JCUT=IY/2+1
        I1 = 2
        I2 = NX-IX+1
	J1 = 2
        J2 = NY-IY+1
        K1 = ILAP/2+1
        K2 = NZ-ILAP/2
        DO 20 K=K1,K2
        DO 20 I=I1,I2
        DO 20 N=1,NTUBES
                RR1=ABS(EXX(I)-XN(N))
                RR2=ABS(XMAX-ABS(EXX(I)-XN(N)))
                IF (RR1.LE.XMAX/TWO) THEN
                        BY(I,JCUT,K)=BY(I,JCUT,K)+EXP(DEN*RR1**2)
     2                          *EXP(DEN*(ZEE(K)-ZN(N))**2)
                ENDIF
                IF (RR2.LT.XMAX/TWO) THEN
                        BY(I,JCUT,K)=BY(I,JCUT,K)+EXP(DEN*RR2**2)
     2                          *EXP(DEN*(ZEE(K)-ZN(N))**2)
                ENDIF
20      CONTINUE
        BMAX=MAXVAL(BY(I1:I2,JCUT,K1:K2))
        STOP 'TUBE: Communication update to MPI needed'
c	FIX WHEN ACTIVATING THIS SECTION
c        CALL BARRIER()
c        CALL SHMEM_REAL8_MAX_TO_ALL(R_BMAX,BMAX,1,0,0,NPE,PWORK1,ISYNC1)
c        CALL BARRIER()
        BMAX = R_BMAX
        BY=BY/BMAX
C----------------------------------------------------------------------
C  Calculate beta so that magnetic flux equals flux of uniform layer.
C  BETALAY = BETA of layer calculation (PAR(08)).
C----------------------------------------------------------------------
	BETALAY = 2.00E00/OBETA
        SURF=ZERO
        DO 30 N=1,NTUBES
                SURF=SURF+SQRT(TWO)*PI*HH**2*ERF(XMAX/TWO/SQRT(TWO)/HH)
     2                    *(ERF((ZCENT+RMAX-ZN(N))/SQRT(TWO)/HH)
     3                      -ERF((ZCENT-RMAX-ZN(N))/SQRT(TWO)/HH))
30      CONTINUE
	SURF=SURF/BMAX
        SURFLAY=TWO*RMAX*XMAX
        BETA=BETALAY*(SURF/SURFLAY)**2
        OBETA=2.00E00/BETA
C----------------------------------------------------------------------
C  Calculate density.
C----------------------------------------------------------------------
        DO 40 K=K1,K2
        DO 40 I=I1,I2
		PRE=RO(I,JCUT,K)*TT(I,JCUT,K)-BY(I,JCUT,K)**TWO*OBETA/TWO
              	RO(I,JCUT,K)=PRE/TT(I,JCUT,K)
40      CONTINUE
	RO(I1:I2,J1:J2,K1:K2)=SPREAD(RO(I1:I2,JCUT,K1:K2),2,J2-J1+1)
        BY(I1:I2,J1:J2,K1:K2)=SPREAD(BY(I1:I2,JCUT,K1:K2),2,J2-J1+1)
        RU=ZERO
        RV=ZERO
        RW=ZERO
C**********************************************************************
C End of the possible models
C**********************************************************************
        CASE DEFAULT
          STOP '3dmhdtube: Invalid NTUBE number'
        END SELECT
        END SUBROUTINE TUBE
C**********************************************************************
        FUNCTION FF(X)
C----------------------------------------------------------------------
        IMPLICIT REAL*8(A-H,O-Z)
c        IMPLICIT REAL(A-H,O-Z)
        COMMON /FFCOM/ A,C  
        COMMON /NUMBERS/ ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX
        FF=A**TWO * X**FIVE * (EXP(-X**TWO)-C)**TWO / 
     &       (A*X**THREE +ONE)**TWO
        END FUNCTION FF
C**********************************************************************
      subroutine getbackground(x,t,rho)
C----------------------------------------------------------------------
C  Written by M. Rempel (Fall 2001).
C----------------------------------------------------------------------
      implicit none

      integer kmax,kount,nvar,nok,nbad

      parameter (nvar=2)

      real*8 XP(100),YP(nvar,100),ystart(100),
     &     x,t,rho
      real*8 h1,hmin,x1,x2,dxsav,eps,krad

      parameter(eps=1d-8)

      external derivs1,rkqc

      COMMON /PATH1/ KMAX,KOUNT
      COMMON /PATH2/ DXSAV,XP,YP
c------------------------------------------------------------------------

      h1       = 1d-4
      hmin     = 1d-50
      dxsav    = 1d-4
      kmax     = 1

      x1       = 2.0
      x2       = x

      ystart(1)= 1.0
      ystart(2)= 1.0

      call ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS1,RKQC)

      t  =yp(1,kount)
      rho=yp(2,kount)/yp(1,kount)

      return
      end


c========================================================================
      subroutine derivs1(x,yt,dyt)
      implicit none

      real*8 x,yt(2),dyt(2),nrad,nad,os,krad,t,rho

c---------------------------------------------------------------------- 
      nad=0.41
      os=1.0

      t  =yt(1)
      rho=yt(2)/yt(1)
      call kappa(x,rho,t,krad)
      nrad=0.4/krad

      if(nrad.le.os*nad) then
         dyt(1)=nrad
         dyt(2)=yt(2)/yt(1)
      else
         dyt(1)=nad
         dyt(2)=yt(2)/yt(1)
      end if

      return
      end
c========================================================================
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RK
     *QC)
      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT REAL (A-H,O-Z)
      PARAMETER (MAXSTP=1000,NMAX=2,TWO=2.0,ZERO=0.0,TINY=1.E-30)
      COMMON /PATH1/ KMAX,KOUNT
      COMMON /PATH2/ DXSAV,XP(100),YP(2,100)

      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      EXTERNAL DERIVS, RKQC
C      PRINT*, 'odeint'

      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE
      XSAV=X-DXSAV*TWO
      DO 16 NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO 12 I=1,NVAR
          YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
12      CONTINUE
        IF(KMAX.GT.0)THEN
          IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1)THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO 13 I=1,NVAR
                YP(I,KOUNT)=Y(I)
13            CONTINUE
              XSAV=X
            ENDIF
          ENDIF
        ENDIF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO 14 I=1,NVAR
            YSTART(I)=Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XP(KOUNT)=X
            DO 15 I=1,NVAR
              YP(I,KOUNT)=Y(I)
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) PAUSE 'Stepsize smaller than minimum.'
        H=HNEXT
16    CONTINUE
      PAUSE 'Too many steps.'
      RETURN
      END

c======================================================================
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT REAL (A-H,O-Z)
      PARAMETER (NMAX=10,FCOR=1.0D0/15.0D0,
     .    ONE=1.,SAFETY=0.9,ERRCON=6.D-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20
      PSHRNK=-0.25
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV)PAUSE 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      RETURN
      END

c=======================================================================
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT REAL (A-H,O-Z)
      PARAMETER (NMAX=10)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      EXTERNAL DERIVS
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
C**********************************************************************
