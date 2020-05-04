!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmg2dlib (basic use)

PROGRAM main

  IMPLICIT NONE

  !> Include here the mmg3d library hader file
! if the header file is in the "include" directory
! #include "libmmg2df.h"

! if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2df.h"

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol
  INTEGER            :: ier,argc
  CHARACTER(len=300) :: exec_name,filename,fileout

  PRINT*,"  -- TEST MMG2DLIB"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)

  IF ( argc /=2 ) THEN
     PRINT*," Usage: ",TRIM(ADJUSTL(exec_name))," input_file_name output_file_name"
     CALL EXIT(1);
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, filename)
  CALL get_command_argument(2, fileout)


  !> ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures
  !!   args of InitMesh:
  !! MMG5_ARG_start: we start to give the args of a variadic func
  !! MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
  !! mmgMesh: your MMG5_pMesh (that store your mesh)
  !! MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
  !! mmgSol: your MMG5_pSol (that store your metric) */

  mmgMesh = 0
  mmgSol  = 0

  CALL MMG2D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
  !! file formatted or manually set your mesh using the MMG2D_Set* functions

  !> with MMG2D_loadMesh function
  CALL MMG2D_loadMesh(mmgMesh,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG2D_loadMet function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG2D_Set* functions

  !> With MMG2D_loadSol function
  CALL MMG2D_loadSol(mmgMesh,mmgSol,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMG2D_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) CALL EXIT(107)

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
 ! NULLIFY(va)
  CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)
  IF ( ier == MMG5_STRONGFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH"
    STOP MMG5_STRONGFAILURE
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB"
  ELSE
     PRINT*,"MMG2DLIB SUCCEED"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMG2D_getMesh/MMG2D_getSol functions

  !> 1) Automatically save the mesh
  CALL MMG2D_saveMesh(mmgMesh,TRIM(ADJUSTL(fileout)),LEN(TRIM(ADJUSTL(fileout))),ier)
  IF ( ier /= 1 ) CALL EXIT(106)

  !> 2) Automatically save the solution
  CALL MMG2D_saveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(fileout)),LEN(TRIM(ADJUSTL(fileout))),ier)
  IF ( ier /= 1 ) CALL EXIT(107)

  !> 3) Free the MMG2D5 structures
  CALL MMG2D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main
