        SUBROUTINE COMMUNICATE_CPU
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
        CALL COMM_MPI_CPU(RU)
        CALL COMM_MPI_CPU(RV)
        CALL COMM_MPI_CPU(RW)
        CALL COMM_MPI_CPU(TT)
        CALL COMM_MPI_CPU(RO)
        IF(LMAG) THEN
           CALL COMM_MPI_CPU(BX)
           CALL COMM_MPI_CPU(BY)
           CALL COMM_MPI_CPU(BZ)
        ENDIF
C
        RETURN
        END

C**********************************************************************
        SUBROUTINE COMM_MPI_CPU(VAR)
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