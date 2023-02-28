C**************************************************************
        SUBROUTINE SETUP(FINP,FOUT,IPAR,PAR)
C--------------------------------------------------------------
	INCLUDE '3dmhdparam.f'
C------------------------------------------------------------
        FINP = '1sicor3d1'
	FOUT = '1sicor3d1'
C--------------------------------------------------------------
	DO 10 I=1,32
10	IPAR(I) = 0
	DO 20 I=1,64
20	PAR(I) = 0.0E00
C--------------------------------------------------------------
C   01  NCASE  (case number)
C   02  NCASEP (previous case number if continuation run)
C   03  NTOTAL (total number of time steps required)
C   04  NSTEP0 (number of steps between field '.dat0','.par',and '.lis' dumps)
C   05  NSTART (0=new start, else= dump number from previous solution)
C   06  IXCON  (adjustable in PARAMETER statement, x boundary condition)
C   07  IYCON  (adjustable in PARAMETER statement, y boundary condition)
C   08  IZCON  (adjustable in PARAMETER statement, lower boundary condition)
C   09  ITCON  (adjustable in PARAMETER statement, upper temperature condition)
C   10  IBCON  (adjustable in PARAMETER statement, magnetic field condition)
C   11  NPX    (adjustable in PARAMETER statement, number x grid points)
C   12  NPY    (adjustable in PARAMETER statement, number x grid points)
C   13  NPZ    (adjustable in PARAMETER statement, number z grid points)
C   14  NPE    (adjustable in PARAMETER statement, total number of processors)
C   15  NIT    (nonadjustable, time step number)
C   16  NTUBE  (tube's model number (initial condition))
C   17  NGRID  (adjustable in PARAMETER statement, grid number (type of grid))
C   18  NCOL   (adjustable in PARAMETER statement, number of columns of tubes)
C   19  NROW   (adjustable in PARAMETER statement, number of rows of tubes)
C   20  ID     (adjustable in PARAMETER statement, density diffusion model)
C   21  NPEY   (adjustable in PARAMETER statement, processors in y direction)
C--------------------------------------------------------------
	IPAR(01)  = 1
	IPAR(02)  = 1
	IPAR(03)  = 1000
	IPAR(04)  = 1001
	IPAR(05)  = 0
C        
        IPAR(16)  = 0
C--------------------------------------------------------
C   01  RE    (soundspeed Reynolds number)
C   02  PR    (pseudo Prandtl number)
C   03  THETA (temperature gradient, THETA=ETA for LREM (M. Rempel))
C   04  GRAV  (nondimensional gravitational constant, = (M+1)*THETA; M=1.5)
C   05  RY    (soundspeed Rossby number)
C   06  ANG   (rotation axis angle from vertical)
C   07  RM    (magnetic Reynolds number)
C   08  BETA  (plasma beta)
C   09  PZP   (depth of transition to stable background)
C   10  SIGMA (width of transition region)
C   11  POLYS (polytropic index of lower stable layer)
C   12  TP    (plume temperature or flux, or boundary nonuniformity)
C   13  TC    (cooling time scale for perturbation onset; embedded onset time)
C   14  XP    (plume source x position)
C   15  YP    (plume source y position)
C   16  ZP    (plume source z position)
C   17  XMAX  (maximum x dimension)
C   18  YMAX  (maximum y dimension)
C   19  ZMAX  (maximum z dimension)
C   20  AMPT  (amplitude of new start random temperature fluctuations)
C   21  AMPB  (amplitude of new start magnetic field layer)
C   22  BFH   (full-width at half-maximum of new start magnetic field layer)
C   23  BZP   (depth of maximum of new start magnetic field layer)
C   24	XX1   (adjustable in PARAMETER statement, grid stretching parameter)
C   25	XX2   (adjustable in PARAMETER statement, grid stretching parameter)
C   26	YY1   (adjustable in PARAMETER statement, grid stretching parameter)
C   27	YY2   (adjustable in PARAMETER statement, grid stretching parameter)
C   28	ZZ1   (adjustable in PARAMETER statement, grid stretching parameter)
C   29	ZZ2   (adjustable in PARAMETER statement, grid stretching parameter)
C   30  SFF   (adjustable in PARAMETER statement, time stepping safety factor)
C   31  DT    (nonadjustable, time step)
C   32  TIMT  (nonadjustable, total simulation time)
C   33  TIMC  (nonadjustable, current simulation time)
C   34  XCENT (initial x position of the center of the magnetic tube)        
C   35  ZCENT (initial z position of the center of the magnetic tube)        
C   36  RMAX  (radius at which the magnetic field is set to zero)
C   37  CMT   (determine the pitch angle where B_phi has its maximum)
C   38  A     (determine the profile of B_phi)
C   39  XA    (adjustable in PARAMETER statement, grid stretching parameter)
C   40  XB    (adjustable in PARAMETER statement, grid stretching parameter)
C   41  XC    (adjustable in PARAMETER statement, grid stretching parameter)
C   42  XD    (adjustable in PARAMETER statement, grid stretching parameter)
C   43  YA    (adjustable in PARAMETER statement, grid stretching parameter)
C   44  YB    (adjustable in PARAMETER statement, grid stretching parameter)
C   45  YC    (adjustable in PARAMETER statement, grid stretching parameter)
C   46  YD    (adjustable in PARAMETER statement, grid stretching parameter)
C   47  ZA    (adjustable in PARAMETER statement, grid stretching parameter)
C   48  ZB    (adjustable in PARAMETER statement, grid stretching parameter)
C   49  ZC    (adjustable in PARAMETER statement, grid stretching parameter)
C   50  ZD    (adjustable in PARAMETER statement, grid stretching parameter)
C   51  LAMBDA(wavelength of the initial pertubation of the 3d tube)
C   52  VZ0   (amplitude  of the initial pertubation of the 3d tube)
C   53  QFH   (temporal full-width at half-maximum of embedded heat loss)
C   54  GAMMA (ratio of specific heats: Cp/Cv)
C   55  HH    (tube widths for layer of tubes, spatial width for plume)
C   56  DH    (adjustable in PARAMETER statement, horizontal density diffusion)
C   57  DV    (adjustable in PARAMETER statement, vertical density diffusion)
C   58  TSTART(time for heat flux relaxation start)
C   59  TOFF  (time scale for heat flux relaxation termination)
C   60  RLAX  (amplitude of heat flux relaxation)
C----------------------------------------------------------------
	PAR(01)  =  75.0E00
	PAR(02)  =  0.2E00
	PAR(03)  =  5.0E00
	PAR(04)  =  10.0E00
	PAR(05)  =  0.0E00
	PAR(06)  =  0.0E00
	PAR(07)  =  0.0E00
	PAR(08)  =  0.0E00
	PAR(09)  =  0.0E00
	PAR(10)  =  0.0E00
	PAR(11)  =  0.0E00
	PAR(12)  =  0.1E00
	PAR(13)  =  0.0E00
	PAR(14)  =  0.0E00
	PAR(15)  =  0.0E00
	PAR(16)  =  0.0E00
	PAR(17)  =  6.0E00
	PAR(18)  =  6.0E00
	PAR(19)  =  2.0E00
	PAR(20)  =  1.0E-03
	PAR(21)  =  0.0E00
	PAR(22)  =  0.0E00
        PAR(23)  =  0.0E00
C        
	PAR(34)  =  0.0E00
	PAR(35)  =  0.0E00
	PAR(36)  =  0.0E00
	PAR(37)  =  0.0E00
	PAR(38)  =  0.0E00
C        
        PAR(51)  =  0.0E00
        PAR(52)  =  0.0E00
C       
        PAR(53)  =  0.00E00
C       
        PAR(54)  =  5.00E00/3.00E00
C
	PAR(55)  =  0.0E00
C
	PAR(58)=0.0E00
	PAR(59)=0.0E00
	PAR(60)=0.0E00
C------------------------------------------------------------------  
	RETURN
	END
C**********************************************************************
c========================================================================
c  Specification of radiative conductivity (M. Rempel).
c========================================================================
      subroutine kappa(z,rho,t,krad)
      implicit none

      real*8 z,rho,t,alpha,beta,krad

c      beta=0.0
c      krad=t**6*sqrt(t)/rho**2
c      krad=(beta+krad)/(beta+1.0)
c      krad=krad+(1.25-krad)*exp(-(z/0.25)**2)

      krad=rho
c      krad=(0.9**10+krad**10)**0.1

      krad=krad+(1.25-krad)*exp(-5.0*z)

      return     
      end
c========================================================================
