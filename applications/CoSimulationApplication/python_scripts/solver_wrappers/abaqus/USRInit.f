C==============================================================================
C======================== GLOBAL DATA =========================================
C==============================================================================
      BLOCK DATA

#ifdef MPI
      INTEGER ID,IDENTIFIED
      COMMON /IDENT/ ID,IDENTIFIED
      SAVE /IDENT/

      DATA IDENTIFIED/-1/
#endif

      END

C==============================================================================
C======================== IDENTIFY SUBROUTINE =================================
C==============================================================================
#ifdef MPI
      SUBROUTINE IDENTIFY
      
      IMPLICIT NONE
      
      INCLUDE 'mpif.h'

      INTEGER ID,IDENTIFIED,IERROR
      COMMON /IDENT/ ID,IDENTIFIED
      SAVE /IDENT/
      
      IF (|cpus| > 1) THEN
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,ID,IERROR)
         IF (IERROR < 0) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem while identifying')
         END IF
      ELSE
         ID = 0
      END IF
       
      IDENTIFIED = 1
      
      RETURN
      END
#endif

C==============================================================================
C======================== DLOAD SUBROUTINE ====================================
C==============================================================================

      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     &   COORDS,JLTYP,SNAME)

      IMPLICIT NONE

      INTEGER D,S
      PARAMETER (D = |dimension|)
      PARAMETER (S = |surfaces|)

#ifdef MPI
      INTEGER ID,IDENTIFIED
      COMMON /IDENT/ ID,IDENTIFIED
      SAVE /IDENT/
#else 
      INTEGER ID
#endif

      DOUBLE PRECISION F,TIME(2),COORDS(D)
      CHARACTER(LEN=20) :: FMT_FACES
      CHARACTER(LEN=200) :: FILENAME
      CHARACTER(LEN=80) :: SNAME
      INTEGER KSTEP,KINC,NOEL,NPT,LAYER,KSPT,JLTYP,R,UNIT_FACES(S)

      IF (KINC == 1) THEN
         FMT_FACES = '(2I6,|dimension|ES27.17E2)'
         UNIT_FACES = (/ (100+R,R=1,S) /)
     
#ifdef MPI
         IF (IDENTIFIED < 0) THEN
            CALL IDENTIFY
         END IF
#else
         ID = 0
#endif

         R = INDEX(SNAME,'SURFACE')
         READ(SNAME((R+7):LEN(TRIM(SNAME))),'(I)') R
         R = R+1

         WRITE(FILENAME,'(A,A,I0,A,I0,A,I0,A)')
     &      '|PWD|',
     &      '/|CSM_dir|/CSM_Time',
     &      (KSTEP-1),'Surface',(R-1),'Cpu',ID,'FacesBis.dat'

         OPEN(UNIT=UNIT_FACES(R),FILE=FILENAME,POSITION='APPEND')

         WRITE(UNIT_FACES(R),FMT_FACES) NOEL,NPT,COORDS

         CLOSE(UNIT_FACES(R))
      ELSE IF (KINC > 1) THEN
         CALL STDB_ABQERR(-3,'USR-abort: normal termination')
      END IF
          
      F = 0
      
      RETURN
      END

C==============================================================================
C======================== UTRACLOAD SUBROUTINE ================================
C==============================================================================

      SUBROUTINE UTRACLOAD(ALPHA,T_USER,KSTEP,KINC,TIME,
     &   NOEL,NPT,COORDS,DIRCOS,JLTYP,SNAME)

      IMPLICIT NONE
      
      INTEGER D,S
      PARAMETER (D = |dimension|)
      PARAMETER (S = |surfaces|)

#ifdef MPI
      INTEGER ID,IDENTIFIED
      COMMON /IDENT/ ID,IDENTIFIED
      SAVE /IDENT/
#else 
      INTEGER ID
#endif

      DOUBLE PRECISION ALPHA,T_USER(D),TIME(2),COORDS(D),DIRCOS(D,D)
      CHARACTER(LEN=20) :: FMT_FACES
      CHARACTER(LEN=200) :: FILENAME
      CHARACTER(LEN=80) :: SNAME
      INTEGER KSTEP,KINC,NOEL,NPT,JLTYP,R,UNIT_FACES(S)

      IF (KINC == 1) THEN
         FMT_FACES = '(2I6,|dimension|ES27.17E2)'
         UNIT_FACES = (/ (200+R,R=1,S) /)
     
#ifdef MPI
         IF (IDENTIFIED < 0) THEN
            CALL IDENTIFY
         END IF
#else
         ID = 0
#endif
     
         R = INDEX(SNAME,'SURFACE')
         READ(SNAME((R+7):LEN(TRIM(SNAME))),'(I)') R
         R = R+1

         WRITE(FILENAME,'(A,A,I0,A,I0,A,I0,A)')
     &      '|PWD|',
     &      '/|CSM_dir|/CSM_Time',
     &      (KSTEP-1),'Surface',(R-1),'Cpu',ID,'Faces.dat'

         OPEN(UNIT=UNIT_FACES(R),FILE=FILENAME,POSITION='APPEND')

         WRITE(UNIT_FACES(R),FMT_FACES) NOEL,NPT,COORDS

         CLOSE(UNIT_FACES(R))
      ELSE IF (KINC > 1) THEN
         CALL STDB_ABQERR(-3,'USR-abort: normal termination')
      END IF
          
      T_USER(1) = 1
      T_USER(2:) = 0
      ALPHA = 0
      
      RETURN
      END
