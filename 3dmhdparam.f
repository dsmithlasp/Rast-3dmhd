        IMPLICIT REAL*8(A-H,O-Z)
c        IMPLICIT REAL(A-H,O-Z)
        LOGICAL LMAG,LROT,LPOT,LREM,LSHR
	CHARACTER*50 FINP,FOUT
C
	DIMENSION IPAR(32),PAR(64)
C
        PARAMETER(LROT=.FALSE.,LMAG=.FALSE.,LPOT=.FALSE.,LREM=.FALSE.,
     2		  LSHR=.FALSE.)
C	LROT: Rotation, LMAG: Magentic fields, LPOT: Ensure divB=0, LREM: Rempel relax, LSHR: Shear instability.
C  The minimum # of processors in the z direction is 2 at the moment!
	PARAMETER(NPEY=8,NPEZ=8,NPE=NPEY*NPEZ)
      	PARAMETER(NPX=512,NPY=512,NPZ=512)
C
        PARAMETER(IX=2,IY=2,IPAD=0,ILAP=2)
C For 2D problem set NX<IX
	PARAMETER(NRY=NPY/NPEY)
        PARAMETER(NRZ=NPZ/NPEZ)
        PARAMETER(NX=NPX+IX,NY=NRY+IY,NZ=NRZ+ILAP)
        PARAMETER(NGRID=0)
c	PARAMETER(XX1=-7.0E+00,XX2=7.0E+00)
c	PARAMETER(YY1=-7.0E+00,YY2=7.0E+00)
c	PARAMETER(ZZ1=-2.0E+00,ZZ2=1.0E-09)
c        PARAMETER(XA=0.5E+00,XB=0.48E+00,XC=28.0E+00,XD=0.025E+00)
c        PARAMETER(YA=0.5E+00,YB=0.48E+00,YC=28.0E+00,YD=0.025E+00)
c        PARAMETER(ZA=0.65E+00,ZB=1.0E+00,ZC=5.0E+00,ZD=0.2E+00)
	PARAMETER(XX1=0.0E+00,XX2=0.0E+00)
	PARAMETER(YY1=0.0E+00,YY2=0.0E+00)
	PARAMETER(ZZ1=0.0E+00,ZZ2=0.0E+00)
        PARAMETER(XA=0.0E+00,XB=0.0E+00,XC=0.0E+00,XD=0.0E+00)
        PARAMETER(YA=0.0E+00,YB=0.0E+00,YC=0.0E+00,YD=0.0E+00)
        PARAMETER(ZA=0.0E+00,ZB=0.0E+00,ZC=0.0E+00,ZD=0.0E+00)
        PARAMETER(IXCON=0,IYCON=0,IZCON=1,ITCON=4,IBCON=0)
c        PARAMETER(SFF=0.315E00)
        PARAMETER(SFF=0.8E00)
c        PARAMETER(SFF=0.4E00)
	PARAMETER(ID=0,DH=0.0E00,DV=0.0E00)
c        PARAMETER(IWORD=8)
        PARAMETER(IWORD=2)
C
	PARAMETER(NCOL=0,NROW=0)
	PARAMETER(NTUBES=NCOL*NROW)
C------------------------------------------------------------------------------
C NPX and NPZ are the number of grid points in each direction.  IX and
C IY must be at least two for second-order horizontal derivatives.  NRZ
C is the number of vertical grid points resident on each processor.
C ILAP=2 for second-order vertical derivatives. IPAD pads the memory
C location to increase perfomance.
C       
C The function used for the stretching of the grid is determined by the
C parameter NGRID. Putting NGRID=1 or NGRID=2 will make the code going
C through different pieces of the routine MKGRID which is called by
C XJACOBI, YJACOBI and ZJACOBI. To add a new grid definition just edit MKGRID. 
C
C NGRID=0: Uniform grid.
C        
C NGRID=1: Original arctan grid function.  XX1, XX2, YY1, YY2,
C and ZZ1, ZZ2 are the arctan parameters for the grid stretch
C (set to very small values, -1.0E-09, 1.0E-09 for no stretching.
C 
C NGRID=2: A function containing two ArcTan is used for the stretching.
C It produces a dip in grid spacing. 
C XA: determines the position of the dip. Takes a value between 0. and 1.
C XB: determines the half width of the dip. 0.2 means a width of 0.4.  
C XC: determines the sharpness of the walls. Set to 0 to have a regular grid
C XD: determines the depth of the dip in terms of DDX(center of dip)/DDX(0)
C     must be between 1.E-9 and 1-1.E-9
C
C IXCON and IYCON specify type of horizontal boundary condition (=0 periodic).
C NOTE: Horizontal periodicity hard-wired into border values in subroutine STEP.
C NOTE: Upper and lower boundaries hard-wired to be stress free and impenetrable
C in subroutine BCON.  IZCON specifies lower boundary temperature condition
C (=0 constant temperature, =1 constant thermal flux).  ITCON specifies plume
C source and/or upper boundary temperature condition (=0 fixed temperature,
C =1 fixed flux, =2 embedded heat loss, =3 thermal, =4 nonuniform boundary).  
C IBCON specifies the upper and lower boundary magnetic field conditions (=0 Bz=0).
C
C ID specifies the vertical variation in density diffusion (=0, no diffusion,
C =1, goes as 1/rhobar, =2 exponentially confined to upper and lower boundaries).
C------------------------------------------------------------------------------
