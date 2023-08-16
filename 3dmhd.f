C******************************************************************************
	PROGRAM THREEDMHD
C------------------------------------------------------------------------------
C  This is a finite difference program for the Cray T3E to simulate three
C  dimensional magnetohydrodynamics in a rotating fully ionized
C  hydrogen ideal-gas of adiabatic stratification in the upper potion and
C  stable stratification in the lower.  Rotation and magnetic fields can be
C  turned on and off with logical parameters LROT and LMAG.  The code empolys 
C  irregular grids, second-order finite differencing, and fully explicit
C  time stepping.  Must be linked with 3dmhdset.f and 3dmhdsub.f.
C  Unless noted the code was written by Mark Peter Rast.
C  Original two-dimensional version finished 12/27/95.  MPP version 11/1/96.
C  Three-dimensional magnetohydrodynamic version 2/5/98.
C  Two dimensional decomposition of the domain by Matthias Rempel 03/21/02.
C------------------------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C
        CHARACTER*50 FIN0,FIN1,FF0,FF1,FLOG
	CHARACTER*80 STRNG,BLANKS
	CHARACTER*4  PESTRNG
C
        DIMENSION RU(NX,NY,NZ),RV(NX,NY,NZ),RW(NX,NY,NZ),RO(NX,NY,NZ)
     2		 ,TT(NX,NY,NZ)
	DIMENSION UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ)
	DIMENSION FU(NX,NY,NZ),FV(NX,NY,NZ),FW(NX,NY,NZ),FR(NX,NY,NZ)
     2		 ,FT(NX,NY,NZ)
	DIMENSION ZRU(NX,NY,NZ),ZRV(NX,NY,NZ),ZRW(NX,NY,NZ),ZRO(NX,NY,NZ)
     2		 ,ZTT(NX,NY,NZ)
	DIMENSION WW1(NX,NY,NZ),WW2(NX,NY,NZ),WW3(NX,NY,NZ)
	DIMENSION ILOC1(3),ILOC2(3),ILOC3(3)
C
	DIMENSION BX(NX,NY,NZ),BY(NX,NY,NZ),BZ(NX,NY,NZ)
	DIMENSION ZBX(NX,NY,NZ),ZBY(NX,NY,NZ),ZBZ(NX,NY,NZ)
	DIMENSION ILOC4(3)
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
        DIMENSION PAR1(96)
     	DIMENSION WMIN(4),WMOUT(4)
C
        COMMON/BIG/RU,SP1,RV,SP2,RW,SP3,RO,SP4,TT,SP5,UU,SP6,VV,SP7,WW
     2            ,SP8,FU,SP9,FV,SP10,FW,SP11,FR,SP12,FT,SP13
     3            ,ZRU,SP14,ZRV,SP15,ZRW,SP16,ZRO,SP17,ZTT
     4            ,SP18,WW1,SP19,WW2,SP20,WW3
C
     5            ,SP21,BX,SP22,BY,SP23,BZ,SP24,ZBX,SP25,ZBY,SP26,ZBZ
C
	COMMON/AJACOBI/EXX,DXXDX,D2XXDX2,DDX,WYY,DYYDY,D2YYDY2,DDY
     2		      ,ZEE,DZZDZ,D2ZZDZ2,DDZ
        COMMON/BCT/IXC,IYC,IZC,ITC,IBC
        COMMON/ITER/NTOTAL,NSTEP0,NIT
        COMMON/CPAR/CV,OCV,ORE,RE,REPR,THETA,GRAV,AMPT,SF,GAMMA
	COMMON/CROT/OMX,OMZ
	COMMON/CMAG/ORM,RM,OBETA,AMPB,BFH,BZP
	COMMON/CPEN/PZP,SIGMA,RKAPST,TB,RKAPA,DKAPA,RKAPM
	COMMON/CPER/TP,XP,YP,ZP,TC,QFH,HH
	COMMON/BOUNDS/XMAX,YMAX,ZMAX
        COMMON/GRID/DD,HX,H2X,HY,H2Y,HZ,H2Z,C13,C23,C43
        COMMON/CTIM/DT,TIMT,TIMC,TIMI
        COMMON/TRACE/UMACH
	COMMON/RUNGKU/GAM1,GAM2,GAM3,ZETA1,ZETA2
	COMMON/SPLINEX/KLOX,KHIX,HHX,SIGX,AAX,BBX,XPINV,XHH,ISEGX
	COMMON/SPLINEY/KLOY,KHIY,HHY,SIGY,AAY,BBY,YPINV,YHH,ISEGY
	COMMON/COMMUN/MYPE,MYPEY,MYPEZ,MPISIZE
	COMMON/RELAX/TSTART,TOFF,RLAX
C
	INCLUDE 'mpif.h'
C----------------------------------------------------------------------
C   Identify processor, and translate it into string.
C----------------------------------------------------------------------
	CALL MPI_INIT(IERR)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYPE,IERR)
	CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NN,IERR)
C----------------------------------------------------------------------
C  Define 2D processor geometry (M. Rempel).
C----------------------------------------------------------------------
	MYPEY=MOD(MYPE,NPEY)
	MYPEZ=MYPE/NPEY
	WRITE(*,'("MYPE ",I3," MYPEY ",I3," MYPEZ ",I3)') MYPE,MYPEY,
     2     MYPEZ
C
	IF (IWORD.EQ.8) THEN
c	IF (IWORD.EQ.2) THEN
		MPISIZE=MPI_DOUBLE_PRECISION
	ELSE
		MPISIZE=MPI_FLOAT
	ENDIF
C
c      	call system('hostname')
C
	CALL PEEXTN(PESTRNG)
C----------------------------------------------------------------------
C   Rename PARAMETERS for common blocks.
C----------------------------------------------------------------------
	IXC=IXCON
	IYC=IYCON
	IZC=IZCON
	ITC=ITCON
	IBC=IBCON
	SF=SFF
C
        C13=1.0E00/3.0E00
        C23=2.0E00*C13
        C43=4.0E00*C13
C
	GAM1=8.0E00/15.0E00
	GAM2=5.0E00/12.0E00
	GAM3=3.0E00/4.0E00
	ZETA1=-17.0E00/60.0E00
	ZETA2=-5.0E00/12.0E00
C
	HX=0.5E00*FLOAT(NPX-1)
        H2X=FLOAT(NPX-1)*FLOAT(NPX-1)
        HY=0.5E00*FLOAT(NPY-1)
        H2Y=FLOAT(NPY-1)*FLOAT(NPY-1)
        HZ=0.5E00*FLOAT(NPZ-1)
        H2Z=FLOAT(NPZ-1)*FLOAT(NPZ-1)
C----------------------------------------------------------------------
C  Call setup and initalise parameters.
C----------------------------------------------------------------------
        CALL SETUP(FINP,FOUT,IPAR,PAR)
C
	NCASE  = IPAR(01)	
	NCASEP = IPAR(02)
	NTOTAL = IPAR(03)
	NSTEP0 = IPAR(04)
	NSTART = IPAR(05)
C
	NIT    = 0
	NDUMP0 = 0
	NSW    = 0
C
	IPAR(06)  = IXC
	IPAR(07)  = IYC
	IPAR(08)  = IZC
	IPAR(09)  = ITC
	IPAR(10)  = IBC
	IPAR(11)  = NPX
	IPAR(12)  = NPY
	IPAR(13)  = NPZ
	IPAR(14)  = NPE
C	
	IPAR(17)  = NGRID
C
	IPAR(18)  = NCOL
	IPAR(19)  = NROW
C
	IPAR(20)  = ID
C
	IPAR(21)  = NPEY
C
	RE     = PAR(01)
	PR     = PAR(02)
	THETA  = PAR(03)
	GRAV   = PAR(04)
	RY     = PAR(05)
	ANG    = PAR(06)
	RM     = PAR(07)
	BETA   = PAR(08)
	PZP    = PAR(09)
	SIGMA  = PAR(10)
	POLYS  = PAR(11)
	TP     = PAR(12)
	TC     = PAR(13)
	XP     = PAR(14)
	YP     = PAR(15)
	ZP     = PAR(16)
	XMAX   = PAR(17)
	YMAX   = PAR(18)
	ZMAX   = PAR(19)
	AMPT   = PAR(20)
	AMPB   = PAR(21)
	BFH    = PAR(22)
	BZP    = PAR(23)
	QFH    = PAR(53)
	GAMMA  = PAR(54)
	HH     = PAR(55)
C
	TSTART = PAR(58)
        TOFF   = PAR(59)
        RLAX   = PAR(60)
C
        IF (LROT) THEN
                OMX=SIN(ANG)/RY
                OMZ=COS(ANG)/RY
        ELSE
                OMX=0.0E00
                OMZ=0.0E00
        ENDIF
C
	PAR(24) = XX1
	PAR(25) = XX2
	PAR(26) = YY1
	PAR(27) = YY2
	PAR(28) = ZZ1
	PAR(29) = ZZ2
	PAR(30) = SF
C
	PAR(39) = XA
	PAR(40) = XB
	PAR(41) = XC
	PAR(42) = XD
	PAR(43) = YA
	PAR(44) = YB
	PAR(45) = YC
	PAR(46) = YD
	PAR(47) = ZA
	PAR(48) = ZB
	PAR(49) = ZC
	PAR(50) = ZD
C
	PAR(56) = DH
	PAR(57) = DV
C	
	REPR=RE*PR
	ORE=1.0E00/RE
	IF (GRAV.EQ.0.0E00) THEN
		RKAPST=1.0E00
	ELSE
		RKAPST=(POLYS+1.0E00)*THETA/GRAV
	ENDIF
	ORM=1.0E00/RM
	OBETA=2.0E00/BETA
C----------------------------------------------------------------------
C  Zero clock.
C----------------------------------------------------------------------
	TIMC = 0.0E00
C----------------------------------------------------------------------
C  Calculate elements of grid transformation Jacobian.
C----------------------------------------------------------------------
	IF (NX.GT.IX+1) THEN
        	CALL XJACOBI
	ELSE
		EXX=0.0E00
		DXXDX=0.0E00
		D2XXDX=0.0E00
		DDX=0.0E00
	ENDIF
        CALL YJACOBI
        CALL ZJACOBI
C----------------------------------------------------------------------
C  Initialize spline common blocks for horizontal integrations.
C----------------------------------------------------------------------
	IF (NX.GT.IX+1) THEN
		CALL INIT_SPLINEX(EXX(IX:NX-IX+1))
	ENDIF
	CALL INIT_SPLINEY(WYY(IY:NY-IY+1))
C----------------------------------------------------------------------
C  Calculate thermal conductivity and its derivative with depth.
C  If LREM is true then RKAPA and DKAPA calculated in static (M. Rempel).
C----------------------------------------------------------------------
	IF (.NOT.LREM) THEN
	IF (PZP.EQ.0.0E00) THEN
                RKAPA=1.0E00
		DKAPA=0.0E00
        ELSE
                RKAPA=1.0E00
     2          +(RKAPST-1.0E00)/2.0E00*(1.0E00+TANH((ZEE-PZP)/SIGMA))
		DKAPA=(RKAPST-1.0E00)/2.0E00/SIGMA
     2				/COSH((ZEE-PZP)/SIGMA)**2.0D00
        ENDIF
	ENDIF
C----------------------------------------------------------------------
C  If NSTART = 0 call static.  If NSTART = else open files and copy 
C  arrays from old run.
C----------------------------------------------------------------------
        I1 = INDEX(FINP,'  ')-1
        I2 = INDEX(FOUT,'  ')-1
        FIN0 = FINP(1:I1)//'.dat0.'//PESTRNG
        FIN1 = FINP(1:I1)//'.par'
C
	IF (NSTART.EQ.0) THEN 
	    CALL STATIC
	    TIMI = 0.0E00
	ELSE
	    IF (LREM) THEN
C----------------------------------------------------------------------
C Call static to get RKAPA and DKAPA and overwrite the rest (M. Rempel).
C----------------------------------------------------------------------
		CALL STATIC
	    ENDIF
	    IF (MYPE.EQ.0) THEN
	    	OPEN(11,FILE=FIN1,STATUS='unknown',FORM='UNFORMATTED'
     2		,ACCESS='DIRECT',RECL=IWORD*96)
	    	READ(11,REC=1)PAR1
	    	CLOSE(11)
	    ENDIF
C
	    CALL MPI_BCAST(PAR1,96,MPISIZE,0,MPI_COMM_WORLD,IERR)
C
	    NX1=NINT(PAR1(11))
	    NY1=NINT(PAR1(12))/NINT(PAR1(21))
	    NZ1=NINT(PAR1(13))/NINT(PAR1(14)/PAR1(21))
	    IF (LMAG) THEN
	    	NSTRT0 = 8*(NSTART-1)
	    ELSE
		NSTRT0 = 5*(NSTART-1)
	    ENDIF
C
	    IF ((NX-IX.EQ.NX1).AND.(NY-IY.EQ.NY1)
     2					.AND.(NZ.EQ.NZ1+ILAP)) THEN
		IF (MYPE.EQ.0) THEN
		     OPEN(11,FILE=FIN1,STATUS='unknown'
     1                      ,FORM='UNFORMATTED'
     2                      ,ACCESS='DIRECT',RECL=IWORD*96)
		     READ(11,REC=NSTART)PAR1
            	     CLOSE(11)
            	ENDIF
C
		CALL MPI_BCAST(PAR1,96,MPISIZE,0,MPI_COMM_WORLD,IERR)
C
	    	TIMI = PAR1(64) + PAR1(65)
C
	    	OPEN(14,FILE=FIN0,STATUS='unknown'
     1                 ,FORM='UNFORMATTED'
     2		       ,ACCESS='DIRECT',RECL=IWORD*NX1*NY1*NZ1)
	    	READ(14,REC=NSTRT0+1)RU(2:NX-IX+1,2:NY-IY+1
     2						 ,ILAP/2+1:NZ-ILAP/2)
	    	READ(14,REC=NSTRT0+2)RV(2:NX-IX+1,2:NY-IY+1
     2						 ,ILAP/2+1:NZ-ILAP/2)
	    	READ(14,REC=NSTRT0+3)RW(2:NX-IX+1,2:NY-IY+1
     2						 ,ILAP/2+1:NZ-ILAP/2)
	    	READ(14,REC=NSTRT0+4)TT(2:NX-IX+1,2:NY-IY+1
     2						 ,ILAP/2+1:NZ-ILAP/2)
	    	READ(14,REC=NSTRT0+5)RO(2:NX-IX+1,2:NY-IY+1
     2						 ,ILAP/2+1:NZ-ILAP/2)
		IF (LMAG) THEN
			READ(14,REC=NSTRT0+6)BX(2:NX-IX+1,2:NY-IY+1
     2                                           ,ILAP/2+1:NZ-ILAP/2)
                	READ(14,REC=NSTRT0+7)BY(2:NX-IX+1,2:NY-IY+1
     2                                           ,ILAP/2+1:NZ-ILAP/2)
                	READ(14,REC=NSTRT0+8)BZ(2:NX-IX+1,2:NY-IY+1
     2                                           ,ILAP/2+1:NZ-ILAP/2)
		ENDIF
	    	CLOSE(14)
C
		IF (MYPEZ.EQ.NPEZ-1) THEN
			TB=TT(2,2,NZ-ILAP/2)
		ENDIF
C
	    ELSE
		WRITE(6,*)'MAIN:  Grid size mismatch, no interpolation'
		CALL MPI_FINALIZE(IERR)
		STOP
	    ENDIF
	ENDIF
C
        CALL COMMUNICATE
C----------------------------------------------------------------------
C  Calculate dt. Fully ionized H specific heats.
C  CV is nondimensional and is equal to Cv*mu/R = 1/(gamma-1)  
C  CP is nondimensional and is equal to Cp*mu/R 
C  GAMMA = CP/CV = Cp/Cv 
C  NOTE: CV is called Cv* in the description file
C  Time step for the thermal and magnetic diffusion terms reduced 
C  by a factor of 2 and those for the viscous diffusion by a factor 
C  of 2*4/3, based on instability arguments for pure diffusion equation,
C  allowing a safety factor increase from 0.4 to 0.8 (Rempel 03/01/02).
C----------------------------------------------------------------------
c	UU=SQRT((RU**2+RV**2+RW**2)/RO**2)
	CV=1.0E00/(GAMMA-1.0E00)   
	OCV=1.0E00/CV
	CP=CV*GAMMA                
c	VV=SQRT(CP*OCV*TT)
	UU(2:NX-IX+1,2:NY-IY+1,:)=(RU(2:NX-IX+1,2:NY-IY+1,:)**2
     2				  +RV(2:NX-IX+1,2:NY-IY+1,:)**2
     3				  +RW(2:NX-IX+1,2:NY-IY+1,:)**2)
     4					 /RO(2:NX-IX+1,2:NY-IY+1,:)**2
	VV(2:NX-IX+1,2:NY-IY+1,:)=GAMMA*TT(2:NX-IX+1,2:NY-IY+1,:)
	IF (LMAG) THEN
		VV(2:NX-IX+1,2:NY-IY+1,:)=VV(2:NX-IX+1,2:NY-IY+1,:)
     2					 +(BX(2:NX-IX+1,2:NY-IY+1,:)**2
     3                                    +BY(2:NX-IX+1,2:NY-IY+1,:)**2
     4                                    +BZ(2:NX-IX+1,2:NY-IY+1,:)**2)
     5                                  *OBETA/RO(2:NX-IX+1,2:NY-IY+1,:)
	ENDIF
	UU(2:NX-IX+1,2:NY-IY+1,:)=UU(2:NX-IX+1,2:NY-IY+1,:)
     2				 +VV(2:NX-IX+1,2:NY-IY+1,:)
     3				 +2.0E00*SQRT(UU(2:NX-IX+1,2:NY-IY+1,:)
     4					     *VV(2:NX-IX+1,2:NY-IY+1,:))
C
	ISW=1
	IF (ISW.EQ.0) THEN
C----------------------------------------------------------------------
C  Find global minimum timestep. Note that this approach yields a 
C  timestep possibly smaller than the minimum pointwise value.
C----------------------------------------------------------------------
	RMIN=MINVAL(RO(2:NX-IX+1,2:NY-IY+1,:))
	VMAX=SQRT(MAXVAL(UU(2:NX-IX+1,2:NY-IY+1,:)))
	IF (NX.GT.IX+1) THEN
        	DD=MIN(MINVAL(DDX(2:NX-IX+1)),MINVAL(DDY(2:NY-IY+1))
     2						          ,MINVAL(DDZ))
	ELSE
		DD=MIN(MINVAL(DDY(2:NY-IY+1)),MINVAL(DDZ))
	ENDIF
	RKAPM=MAXVAL(RKAPA)
C	
	WMIN(1)=RMIN
	WMIN(2)=DD
	CALL  MPI_ALLREDUCE(WMIN,WMOUT,2,MPISIZE,
     2			    MPI_MIN,MPI_COMM_WORLD,IERR)
	RMIN=WMOUT(1)
	DD=WMOUT(2)
C
	WMIN(1)=VMAX
	WMIN(2)=RKAPM
	CALL  MPI_ALLREDUCE(WMIN,WMOUT,2,MPISIZE,
     2			    MPI_MAX,MPI_COMM_WORLD,IERR)
	VMAX=WMOUT(1)
	RKAPM=WMOUT(2)
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
		DT4=0.0E00
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
2      	CONTINUE
	WMIN(1)=MINVAL(WW1(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	WMIN(2)=MINVAL(WW2(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	WMIN(3)=MINVAL(WW3(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	MINCNT=3
	IF (LMAG) THEN
		WMIN(4)=MINVAL(VV(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
		MINCNT=4
	ENDIF
C
	CALL  MPI_ALLREDUCE(WMIN,WMOUT,MINCNT,MPISIZE,
     2			    MPI_MIN,MPI_COMM_WORLD,IERR)
C
	DT1=WMOUT(1)
	DT2=WMOUT(2)
	DT3=WMOUT(3)
	IF (LMAG) THEN
		DT4=WMOUT(4)
		DT=SF*MIN(DT1,DT2,DT3,DT4)
	ELSE
		DT4=0.0E00
                DT=SF*MIN(DT1,DT2,DT3)
	ENDIF
	ENDIF
C
	NSOUND=NINT(1.0E00/DT)
C----------------------------------------------------------------------
C  Zero the fluxes.
C----------------------------------------------------------------------
	FU = 0.0E00
	FV = 0.0E00
	FW = 0.0E00
	FT = 0.0E00
	FR = 0.0E00
	WW1 = 0.0E00
	WW2 = 0.0E00
	WW3 = 0.0E00
C----------------------------------------------------------------------
C  Create the files for the results.
C----------------------------------------------------------------------
        FF0 = FOUT(1:I2)//'.dat0.'//PESTRNG
        FF1 = FOUT(1:I2)//'.par'
        FLOG = FOUT(1:I2)//'.lis'
C
	DO 10 I=1,80
		BLANKS(I:I) = ' '
10	CONTINUE
C
	DO 20 I=1,32
		PAR1(I) = FLOAT(IPAR(I))
20	CONTINUE
        DO 25 I=1,64
                PAR1(32+I) = PAR(I)
25      CONTINUE
C
	IF (NIT.EQ.0) THEN
	IF (MYPE.EQ.0) THEN
		IREC=0
		OPEN(UNIT=60,FILE=FLOG,STATUS='unknown',FORM='FORMATTED'
     2		   ,ACCESS='DIRECT',RECL=80)
C----------------------------------------------------------------------
C  Write words of greeting.
C----------------------------------------------------------------------
		IREC=IREC+1
		WRITE(STRNG,4000)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2		      		     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2		      		     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2		      		     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4002)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4003)NCASE
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		IF (NSTART.EQ.0) THEN
			WRITE(STRNG,4004)
		ELSEIF (NCASE.EQ.NCASEP) THEN
			WRITE(STRNG,4005)
		ELSE
			WRITE(STRNG,4006)NCASEP
		ENDIF
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                              				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				    				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				    				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4000)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4007)NPX,NPY,NPZ
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		IF (LROT) THEN
			WRITE(STRNG,4008)
		ELSE
			WRITE(STRNG,4009)
		ENDIF
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                               				CHAR(10)
		IREC=IREC+1
                IF (LMAG) THEN
                        WRITE(STRNG,4010)
                ELSE
                        WRITE(STRNG,4011)
                ENDIF
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                               				CHAR(10)
                IREC=IREC+1
                WRITE(STRNG,4012)
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
                WRITE(STRNG,4013)
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		IF (IZC.EQ.0) THEN
                	WRITE(STRNG,4014)
		ELSE IF (IZC.EQ.1) THEN
			WRITE(STRNG,4015)
		ENDIF
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		IF (ITC.EQ.0) THEN
                        WRITE(STRNG,4016)
		ELSE IF (ITC.EQ.1) THEN
                        WRITE(STRNG,4017)
		ELSE IF (ITC.EQ.2) THEN
                        WRITE(STRNG,4059)
		ELSE IF (ITC.EQ.3) THEN
                        WRITE(STRNG,4061)
		ELSE IF (ITC.EQ.4) THEN
                        WRITE(STRNG,4062)
                ELSE 
                        WRITE(STRNG,4060)
                ENDIF
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4018)
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4019)RE
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4020)PR
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4021)THETA
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4022)GRAV
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4023)RY
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4024)ANG
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
                WRITE(STRNG,4025)RM
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
                WRITE(STRNG,4026)BETA
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4027)PZP
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4028)SIGMA
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4029)POLYS
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
                IREC=IREC+1
                WRITE(STRNG,4030)TP
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
                IREC=IREC+1
                WRITE(STRNG,4031)XP
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
                WRITE(STRNG,4032)YP
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
                WRITE(STRNG,4089)ZP
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4033)XMAX
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
                IREC=IREC+1
                WRITE(STRNG,4034)YMAX
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4035)ZMAX
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4036)
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4037)SF
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
                IF (NSTART.EQ.0) THEN
                        IREC=IREC+1
                        WRITE(STRNG,4038)
                        CALL SLENGTH(STRNG,JLEN)
                        WRITE(60,3999,REC=IREC)STRNG(1:JLEN),
     2                                       BLANKS(1:79-JLEN),CHAR(10)
			IREC=IREC+1
			WRITE(STRNG,4039)AMPT
			CALL SLENGTH(STRNG,JLEN)
			WRITE(60,3999,REC=IREC)STRNG(1:JLEN),
     2				             BLANKS(1:79-JLEN),CHAR(10)
			IREC=IREC+1
			IF (LMAG) THEN
                                WRITE(STRNG,4040)AMPB
                        ELSE
                                WRITE(STRNG,4011)
                        ENDIF
                        CALL SLENGTH(STRNG,JLEN)
                        WRITE(60,3999,REC=IREC)STRNG(1:JLEN),
     2                                       BLANKS(1:79-JLEN),CHAR(10)
		ELSE                                                   
			IREC=IREC+1
			WRITE(STRNG,4041)NSTART,FINP
			CALL SLENGTH(STRNG,JLEN)
			WRITE(60,3999,REC=IREC)STRNG(1:JLEN),
     2				             BLANKS(1:79-JLEN),CHAR(10)
			IREC=IREC+1
			WRITE(STRNG,4042)TIMI
			CALL SLENGTH(STRNG,JLEN)
			WRITE(60,3999,REC=IREC)STRNG(1:JLEN),
     2				             BLANKS(1:79-JLEN),CHAR(10)
			IREC=IREC+1
			WRITE(STRNG,4043)NX1,NY1*NPEY,NZ1*NPEZ
			CALL SLENGTH(STRNG,JLEN)
			WRITE(60,3999,REC=IREC)STRNG(1:JLEN),
     2				             BLANKS(1:79-JLEN),CHAR(10)
		ENDIF
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4044)FOUT
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4045)NTOTAL
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4046)NSTEP0
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4047)DT1
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4048)DT2
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4049)DT3
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
                WRITE(STRNG,4050)DT4
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4051)DT
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4052)NSOUND
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
		IREC=IREC+1
		WRITE(STRNG,4000)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				     				CHAR(10)
                IREC=IREC+1
                WRITE(STRNG,4000)
                CALL SLENGTH(STRNG,JLEN)
                WRITE(60,3999,REC=IREC)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                                                          CHAR(10)
C
		CLOSE(60)
C                                  
	IREC=59
	ENDIF
	ENDIF
C----------------------------------------------------------------------
C  Ready to start.
C----------------------------------------------------------------------
	NBEG = NIT + 1
c----------------------------------------------------------------------
c  Stats.
c----------------------------------------------------------------------
c	call start()
c----------------------------------------------------------------------
c  Start timing (to time remove the cc's and link with lib_perf_mjp.a).
c----------------------------------------------------------------------
c  select counts to be made:
c  count machine cycles
cc        isel0 = 0
c  count flops
cc        isel1 = 10
c  count Scache misses, CBOX input 2
cc        isel2 = 15
c  clear and start the counters
cc        call  perfctr_start(isel0, isel1, isel2)
c
c!       stop 'loop'
	DO 1000 NK=NBEG,NTOTAL	
C
	CALL STEP
C
	PAR1(15) = FLOAT(NIT)
	PAR1(63) = DT
	IF (NK.EQ.NTOTAL) THEN
		PAR1(64) = TIMT
		PAR1(65)=0.0E00
	ELSE
		PAR1(64) = TIMI
		PAR1(65) = TIMC
	ENDIF
C
	IF (DT.LT.1.0E-08) GOTO 5010
C
	IF (MOD(NIT,NSTEP0).EQ.0) THEN
                KSTRT0 = NDUMP0
		IF (LMAG) THEN
			ISTRT0 = 8*NDUMP0
		ELSE
			ISTRT0 = 5*NDUMP0
		ENDIF
		NDUMP0 = NDUMP0+1     
C
		OPEN(10,FILE=FF0,STATUS='unknown'
     1                  ,FORM='UNFORMATTED'
     2			,ACCESS='DIRECT',RECL=IWORD*(NX-IX)*(NY-IY)
     3							 *(NZ-ILAP))
C
		WRITE(10,REC=ISTRT0+1)RU(2:NX-IX+1,2:NY-IY+1
     2						  ,ILAP/2+1:NZ-ILAP/2)
		WRITE(10,REC=ISTRT0+2)RV(2:NX-IX+1,2:NY-IY+1
     2						  ,ILAP/2+1:NZ-ILAP/2)
		WRITE(10,REC=ISTRT0+3)RW(2:NX-IX+1,2:NY-IY+1
     2						  ,ILAP/2+1:NZ-ILAP/2)
		WRITE(10,REC=ISTRT0+4)TT(2:NX-IX+1,2:NY-IY+1
     2						  ,ILAP/2+1:NZ-ILAP/2)
		WRITE(10,REC=ISTRT0+5)RO(2:NX-IX+1,2:NY-IY+1
     2						  ,ILAP/2+1:NZ-ILAP/2)
		IF (LMAG) THEN
			WRITE(10,REC=ISTRT0+6)BX(2:NX-IX+1,2:NY-IY+1
     2                                            ,ILAP/2+1:NZ-ILAP/2)
                	WRITE(10,REC=ISTRT0+7)BY(2:NX-IX+1,2:NY-IY+1
     2                                            ,ILAP/2+1:NZ-ILAP/2)
                	WRITE(10,REC=ISTRT0+8)BZ(2:NX-IX+1,2:NY-IY+1
     2                                            ,ILAP/2+1:NZ-ILAP/2)
		ENDIF
C
		CLOSE(10)
C
	IF (MYPE.EQ.0) THEN
                OPEN(11,FILE=FF1,STATUS='unknown',FORM='UNFORMATTED'
     2                 ,ACCESS='DIRECT',RECL=IWORD*96)
                WRITE(11,REC=KSTRT0+1)PAR1
                CLOSE(11)
C
                IREC=59+7*NIT/NSTEP0
C
		OPEN(60,FILE=FLOG,STATUS='unknown',FORM='FORMATTED'
     2		       ,ACCESS='DIRECT',RECL=80)
C
		WRITE(STRNG,4001)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC-6)STRNG(1:JLEN),
     2					 BLANKS(1:79-JLEN),CHAR(10)
C
		WRITE(STRNG,4053)NIT
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC-5)STRNG(1:JLEN),
     2					 BLANKS(1:79-JLEN),CHAR(10)
C
		WRITE(STRNG,4054)
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC-4)STRNG(1:JLEN),
     2					 BLANKS(1:79-JLEN),CHAR(10)
C
		WRITE(STRNG,4055)TIMT
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC-3)STRNG(1:JLEN),
     2					 BLANKS(1:79-JLEN),CHAR(10)
C
		WRITE(STRNG,4056)TIMC
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC-2)STRNG(1:JLEN),
     2					 BLANKS(1:79-JLEN),CHAR(10)
C
		WRITE(STRNG,4057)UMACH
		CALL SLENGTH(STRNG,JLEN)
		WRITE(60,3999,REC=IREC-1)STRNG(1:JLEN),
     2					 BLANKS(1:79-JLEN),CHAR(10)
C
		WRITE(STRNG,4058)DT
		CALL SLENGTH(STRNG,JLEN)
		WRITE (60,3999,REC=IREC)STRNG(1:JLEN),
     2					BLANKS(1:79-JLEN),CHAR(10)
C
		CLOSE(60)
	ENDIF
	ENDIF
C
1000	CONTINUE
c------------------------------------------------------------------------------
c  End stats.
c------------------------------------------------------------------------------
c	call summary()
c------------------------------------------------------------------------------
c  End timing stuff.
c------------------------------------------------------------------------------
c  stop, read counters (to time remove the cc's and link with lib_perf_mjp.a)
cc        call perfctr_stopread(ictr0, ictr1, ictr2)
cc        tsec=dble(ictr0)/3.0e08
cc        flop=dble(ictr1)
cc        speed=flop/tsec/1.0e06
cc        write(6,*)mype,' time: tsec=',tsec,' flop=',flop
cc        write(6,*)mype,' Mflops=',speed
C
	IF (MYPE.EQ.0) THEN
C
	OPEN(60,FILE=FLOG,STATUS='unknown',FORM='FORMATTED'
     2	      ,ACCESS='DIRECT',RECL=80)
	WRITE(STRNG,4000)
	CALL SLENGTH(STRNG,JLEN)
	WRITE(60,3999,REC=IREC+1)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				 CHAR(10)
	WRITE(STRNG,4090)
	CALL SLENGTH(STRNG,JLEN)
	WRITE(60,3999,REC=IREC+2)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				 CHAR(10)
	WRITE(STRNG,4000)
	CALL SLENGTH(STRNG,JLEN)
	WRITE(60,3999,REC=IREC+3)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				 CHAR(10)
	WRITE(STRNG,4000)
	CALL SLENGTH(STRNG,JLEN)
	WRITE(60,3999,REC=IREC+4)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				 CHAR(10)
	CLOSE(60)
C
	ENDIF
	CALL MPI_FINALIZE(IERR)
	STOP
C
5010	CONTINUE
	KSTRT0 = NDUMP0
        IF (LMAG) THEN
                ISTRT0 = 8*NDUMP0
        ELSE
                ISTRT0 = 5*NDUMP0
        ENDIF
        NDUMP0 = NDUMP0+1
C
        OPEN(10,FILE=FF0,STATUS='unknown'
     1          ,FORM='UNFORMATTED'
     2          ,ACCESS='DIRECT',RECL=IWORD*(NX-IX)*(NY-IY)
     3                                           *(NZ-ILAP))
C
        WRITE(10,REC=ISTRT0+1)RU(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
        WRITE(10,REC=ISTRT0+2)RV(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
        WRITE(10,REC=ISTRT0+3)RW(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
        WRITE(10,REC=ISTRT0+4)TT(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
        WRITE(10,REC=ISTRT0+5)RO(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
        IF (LMAG) THEN
                WRITE(10,REC=ISTRT0+6)BX(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
                WRITE(10,REC=ISTRT0+7)BY(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
                WRITE(10,REC=ISTRT0+8)BZ(2:NX-IX+1,2:NY-IY+1
     2                                    ,ILAP/2+1:NZ-ILAP/2)
        ENDIF
C
        CLOSE(10)
C----------------------------------------------------------------------
C  Find the pointwise minimum timestep.
C----------------------------------------------------------------------
	DO 5012 K=ILAP/2+1,NZ-ILAP/2
        DO 5012 J=2,NY-IY+1
        DO 5012 I=2,NX-IX+1
                DD=MIN(DDX(I),DDY(J),DDZ(K))
                WW1(I,J,K)=DD/SQRT(UU(I,J,K))
        IF (LREM) THEN
C----------------------------------------------------------------------
C  Thermal conductivity split into a constant and depth dependant part.
C----------------------------------------------------------------------
                WW2(I,J,K)=0.5*DD*DD*REPR*CV*RO(I,J,K)/(1.0+RKAPA(K))
        ELSE
                WW2(I,J,K)=0.5E00*DD*DD*REPR*CV*RO(I,J,K)/RKAPA(K)
        ENDIF
                WW3(I,J,K)=0.5E00*DD*DD*RE*RO(I,J,K)
                IF (LMAG) THEN
                        VV(I,J,K)=0.5E00*DD*DD*RM
                ENDIF
5012    CONTINUE
        WMIN1=MINVAL(WW1(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
        WMIN2=MINVAL(WW2(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
        WMIN3=MINVAL(WW3(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	ILOC1=MINLOC(WW1(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	ILOC2=MINLOC(WW2(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	ILOC3=MINLOC(WW3(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
	WRITE(6,*)MYPE,' Min DT1: ',WMIN1,ILOC1
	WRITE(6,*)MYPE,' Min DT2: ',WMIN2,ILOC2
	WRITE(6,*)MYPE,' Min DT3: ',WMIN3,ILOC3
	IF (LMAG) THEN
		WMIN4=MINVAL(VV(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
		ILOC4=MINLOC(VV(2:NX-IX+1,2:NY-IY+1,ILAP/2+1:NZ-ILAP/2))
		WRITE(6,*)MYPE,' Min DT4: ',WMIN4,ILOC4
	ENDIF
C
	IF (MYPE.EQ.0) THEN
	OPEN(11,FILE=FF1,STATUS='unknown',FORM='UNFORMATTED'
     2         ,ACCESS='DIRECT',RECL=IWORD*96)
        WRITE(11,REC=KSTRT0+1)PAR1
        CLOSE(11)
C
	IREC=58+7*NIT/NSTEP0
C
	OPEN(60,FILE=FLOG,STATUS='unknown',FORM='FORMATTED'
     2	      ,ACCESS='DIRECT',RECL=80)
	WRITE(STRNG,4000)
	CALL SLENGTH(STRNG,JLEN)
	WRITE(60,3999,REC=IREC+1)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				 CHAR(10)
	WRITE(STRNG,4100)DT,NIT-1
	CALL SLENGTH(STRNG,JLEN)
	WRITE(60,3999,REC=IREC+2)STRNG(1:JLEN),
     2			 	 BLANKS(1:79-JLEN),CHAR(10)
        WRITE(STRNG,4000)
        CALL SLENGTH(STRNG,JLEN)
        WRITE(60,3999,REC=IREC+3)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2                           CHAR(10)
	WRITE(STRNG,4000)
	CALL SLENGTH(STRNG,JLEN)
	WRITE(60,3999,REC=IREC+4)STRNG(1:JLEN),BLANKS(1:79-JLEN),
     2				 CHAR(10)
	CLOSE(60)
C
	ENDIF
	CALL MPI_FINALIZE(IERR)
	STOP
C
3999	FORMAT(A,A,A)
4000	FORMAT(78('-'),'@')
4001	FORMAT(78(' '),'@')
4002	FORMAT(5X,'THREE-D MAGNETOHYDRODYNAMICS','@')
4003	FORMAT(5X,'Case #',I3,'@')
4004	FORMAT(5X,'Static start','@')
4005	FORMAT(5X,'Continuation','@')
4006	FORMAT(5X,'Previous case #',I3,'@')
4007	FORMAT(5X,'Array size: NX=',I4,' NY=',I4,' NZ=',I4,'@')
4008    FORMAT(5X,'Rotating domain.','@')
4009    FORMAT(5X,'Nonrotating domain.','@')
4010    FORMAT(5X,'Magnetic fields included.','@')
4011    FORMAT(5X,'No magnetic fields.','@')
4012    FORMAT(5X,'Horizontally periodic.','@')
4013	FORMAT(5X,'Stress free and impenetrable lower boundary.','@')
4014	FORMAT(5X,'Constant temperature lower boundary.','@')
4015	FORMAT(5X,'Constant thermal flux lower boundary.','@')
4016	FORMAT(5X,'Constant temperature upper boundary.','@')
4017	FORMAT(5X,'Constant thermal flux upper boundary.','@')
4059	FORMAT(5X,'Embedded heat loss.','@')
4060	FORMAT(5X,'No temperature perturbation.','@')
4061	FORMAT(5X,'Embedded thermal.','@')
4062	FORMAT(5X,'Differential upper boundary temperature.','@')
4018	FORMAT(5X,'Parameters for this run:','@')
4019 	FORMAT(10X,'Re    = ',1P1D11.4,'@')
4020	FORMAT(10X,'Pr    = ',1P1D11.4,'@')
4021	FORMAT(10X,'Theta = ',1P1D11.4,'@')
4022	FORMAT(10X,'Grav  = ',1P1D11.4,'@')
4023	FORMAT(10X,'Ro    = ',1P1D11.4,'@')
4024	FORMAT(10X,'Angle = ',1P1D11.4,'@')
4025    FORMAT(10X,'Rm    = ',1P1D11.4,'@')
4026    FORMAT(10X,'Beta  = ',1P1D11.4,'@')
4027	FORMAT(10X,'Zpen  = ',1P1D11.4,'@')
4028	FORMAT(10X,'Sigma = ',1P1D11.4,'@')
4029	FORMAT(10X,'ms    = ',1P1D11.4,'@')
4030	FORMAT(10X,'Tp    = ',1P1D11.4,'@')
4031	FORMAT(10X,'Xp    = ',1P1D11.4,'@')
4032	FORMAT(10X,'Yp    = ',1P1D11.4,'@')
4089	FORMAT(10X,'Zp    = ',1P1D11.4,'@')
4033	FORMAT(10X,'Xmax  = ',1P1D11.4,'@')
4034	FORMAT(10X,'Ymax  = ',1P1D11.4,'@')
4035	FORMAT(10X,'Zmax  = ',1P1D11.4,'@')
4036	FORMAT(5X,'Time stepping safety factor:','@')
4037	FORMAT(10X,'SFF = ',1P1D10.4,'@')
4038	FORMAT(5X,'Calculation started from a static state.','@')
4039	FORMAT(5X,'Initial temperature perturbations: ',1P1D10.4,'@')
4040	FORMAT(5X,'Initial magnetic layer amplitude: ',1P1D10.4,'@')
4041	FORMAT(5X,'Started from dump ',I2,' of file ',A20,'@')
4042	FORMAT(5X,'Previous simulation time: ',1P1D10.4,'@')
4043    FORMAT(5X,'Previous size: NX=',I4,' NY=',I4,' NZ=',I4,'@')
4044	FORMAT(5X,'The results are stored in file ',A15,'@')
4045	FORMAT(5X,'Total number of iterations required = ',I8,'@')
4046	FORMAT(5X,'Variables dumped every ',I6,' iterations.','@')
4047	FORMAT(5X,'Advected fast-mode time: ',1P1D10.4,'@')
4048	FORMAT(5X,'Thermal diffusion time:  ',1P1D10.4,'@')
4049	FORMAT(5X,'Viscous diffusion time:  ',1P1D10.4,'@')
4050	FORMAT(5X,'Magnetic diffusion time: ',1P1D10.4,'@')
4051	FORMAT(5X,'Initial time-step value: ',1P1D10.4,'@')
4052	FORMAT(5X,'Iterations per crossing time: ',I8,'@')
4053  	FORMAT(5X,'Iteration ',I7,' completed','@')
4054	FORMAT(5X,'-----------------------------------','@')
4055	FORMAT(5X,'Total simulation time:   ',1P1D10.4,'@')
4056	FORMAT(5X,'Present simulation time: ',1P1D10.4,'@')
4057	FORMAT(5X,'Maximum Mach number:     ',1P1D10.4,'@')
4058 	FORMAT(5X,'Present time step:       ',1P1D10.4,'@')
4090	FORMAT(15X,'SIMULATION COMPLETE','@') 
4100	FORMAT(5X,'Time step ',1P1D11.4,'  Shut down at step ',I8,'@')
	END                        
C**********************************************************************
