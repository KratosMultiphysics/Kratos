!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmgslib (basic use)

PROGRAM main

  IMPLICIT NONE

  !> Include the mmgs library hader file
  ! if the header file is in the "include" directory
  ! #include "libmmgsf.h"

  ! if the header file is in "include/mmg/mmgs"
#include "mmg/mmgs/libmmgsf.h"


  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol
  INTEGER            :: ier,argc
  CHARACTER(len=300) :: exec_name,filename,fileout

  WRITE(*,*) "  -- TEST MMGSLIB"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)

  IF ( argc /=2 ) THEN
     PRINT*," Usage: ",TRIM(exec_name)," input_file_name output_file_name"
     CALL EXIT(1);
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, filename)
  CALL get_command_argument(2, fileout)

  !> ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures
  !! args of InitMesh:
  !! MMG5_ARG_start: we start to give the args of a variadic func
  !! MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
  !! mmgMesh: your MMG5_pMesh (that store your mesh)
  !! MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
  !! mmgSol: your MMG5_pSol (that store your metric) */
  mmgMesh = 0
  mmgSol  = 0

  CALL MMGS_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMGS_loadMesh function that will read a .mesh(b)
  !!    file formatted or manually set your mesh using the MMGS_Set* functions

  !> with MMGS_loadMesh function
  CALL MMGS_loadMesh(mmgMesh,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMGS_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMGS_Set* functions

  !> With MMGS_loadSol function
  CALL MMGS_loadSol(mmgMesh,mmgSol,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMGS_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) CALL EXIT(105)

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
  !! Remark: %val(0) allow to pass the value 0 (i.e. NULL) instead of a pointer
  !! toward NULL.
  CALL MMGS_mmgslib(mmgMesh,mmgSol,ier)
  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMGSLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMGS_saveMesh/MMGS_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMGS_getMesh/MMGS_getSol functions

  !> 1) Automatically save the mesh
  CALL MMGS_saveMesh(mmgMesh,TRIM(ADJUSTL(fileout)),LEN(TRIM(ADJUSTL(fileout))),ier)
  IF ( ier /= 1 ) CALL EXIT(106)

  !> 2) Automatically save the solution
  CALL MMGS_saveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(fileout)),LEN(TRIM(ADJUSTL(fileout))),ier)
  IF ( ier /= 1 ) CALL EXIT(107)

  !> 3) Free the MMGS5 structures
  CALL MMGS_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main
