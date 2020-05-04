!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmglib (basic use)

PROGRAM main

  IMPLICIT NONE

!> Include here the mmg library hader file
! if the "include/mmg" dir is in your include path
!#include "libmmgf.h"

! if your include path do not contain the "mmg/mmg" subdirectories
#include "mmg/libmmgf.h"

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol
  INTEGER            :: ier,argc
  CHARACTER(len=350) :: tdfile,tdfileout, bdfile, bdfileout,exec_name,fileout

  WRITE(*,*) "  -- TEST MMGLIB"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)

  IF ( argc /=3 ) THEN
     PRINT*,argc," Usage: ",TRIM(exec_name),&
          " 2d_file_name 3d_file_name output_file_name"
     CALL EXIT(1);
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, bdfile)
  CALL get_command_argument(2, tdfile)
  CALL get_command_argument(3, fileout)

  !> ================== 2d remeshing using the mmg2d library
  !! ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures */
  !! args of InitMesh:
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
  !!    file formatted or manually set your mesh using the MMG2D_Set* functions

  !> with MMG2D_loadMesh function
  CALL MMG2D_loadMesh(mmgMesh,TRIM(ADJUSTL(bdfile)),&
       LEN(TRIM(ADJUSTL(bdfile))),ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG2D_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG2D_Set* functions
  CALL MMG2D_loadSol(mmgMesh,mmgSol,TRIM(ADJUSTL(bdfile)),LEN(TRIM(ADJUSTL(bdfile))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
  CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMG2D_getMesh/MMG2D_getSol functions

  !> 1) Automatically save the mesh
  CALL MMG2D_saveMesh(mmgMesh,TRIM(ADJUSTL(fileout)) // ".2d.mesh",&
       LEN(TRIM(ADJUSTL(fileout)) // ".2d.mesh"),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(106)
  ENDIF

  !> 2) Automatically save the solution

  !! b) save the solution in a file named bdfileout
  CALL MMG2D_saveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(fileout)) // ".2d.mesh",&
       LEN(TRIM(ADJUSTL(fileout)) // ".2d.mesh"),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(107)
  ENDIF

  !> 3) Free the MMG2D5 structures
  CALL MMG2D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  !> ================== surface remeshing using the mmgs library

  !! ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures */
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
  CALL MMGS_loadMesh(mmgMesh,TRIM(ADJUSTL(tdfile)),&
       LEN(TRIM(ADJUSTL(tdfile))),ier);
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMGS_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMGS_Set* functions

  !> With MMGS_loadSol function
  CALL MMGS_loadSol(mmgMesh,mmgSol,TRIM(ADJUSTL(tdfile)),&
       LEN(TRIM(ADJUSTL(tdfile))),ier);
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMGS_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /=1 ) CALL EXIT(105)

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
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
  CALL MMGS_saveMesh(mmgMesh,TRIM(ADJUSTL(fileout)) // ".s.mesh",&
       LEN(TRIM(ADJUSTL(fileout)) // ".s.mesh"),ier);
  IF ( ier /= 1 ) THEN
     CALL EXIT(108)
  ENDIF

  !> 2) Automatically save the solution
  CALL MMGS_saveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(fileout))//".s.sol",&
       LEN(TRIM(ADJUSTL(fileout))//".s.sol"),ier);
  IF ( ier /= 1 ) THEN
     CALL EXIT(107)
  ENDIF

  !> 3) Free the MMGS5 structures
  CALL MMGS_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  !> ================== 3d remeshing using the mmg3d library

  !! ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures */
  !! args of InitMesh:
  !! MMG5_ARG_start: we start to give the args of a variadic func
  !! MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
  !! mmgMesh: your MMG5_pMesh (that store your mesh)
  !! MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
  !! mmgSol: your MMG5_pSol (that store your metric) */

  mmgMesh = 0
  mmgSol  = 0

  CALL MMG3D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
  !!    file formatted or manually set your mesh using the MMG3D_Set* functions

  !> with MMG3D_loadMesh function
  CALL MMG3D_loadMesh(mmgMesh,TRIM(ADJUSTL(tdfile)),&
       LEN(TRIM(ADJUSTL(tdfile))),ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG3D_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG3D_Set* functions

  !> With MMG3D_loadSol function
  CALL MMG3D_loadSol(mmgMesh,mmgSol,TRIM(ADJUSTL(tdfile)),LEN(TRIM(ADJUSTL(tdfile))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMG3D_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) CALL EXIT(105)

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
  CALL MMG3D_mmg3dlib(mmgMesh,mmgSol,ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG3DLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMG3D_saveMesh/MMG3D_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMG3D_getMesh/MMG3D_getSol functions

  !> 1) Automatically save the mesh
  CALL MMG3D_saveMesh(mmgMesh,TRIM(ADJUSTL(fileout)) // ".3d.mesh",&
       LEN(TRIM(ADJUSTL(fileout)) // ".3d.mesh"),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(108)
  ENDIF

  !> 2) Automatically save the solution
  CALL MMG3D_saveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(fileout)) // ".3d.sol",&
       LEN(TRIM(ADJUSTL(fileout)) // ".3d.sol"),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(107)
  ENDIF


  !> 3) Free the MMG3D5 structures
  CALL MMG3D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main
