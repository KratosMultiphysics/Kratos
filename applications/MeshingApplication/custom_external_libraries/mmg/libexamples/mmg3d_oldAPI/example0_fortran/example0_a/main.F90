!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmg3dlib (basic use)

!> Include the mmg3d library hader file
! if the header file is in the "include" directory
! #include "libmmg3df.h"

! if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"

PROGRAM main

  IMPLICIT NONE

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol
  INTEGER            :: ier
  CHARACTER(len=255) :: pwd
  CHARACTER(len=300) :: filename

  WRITE(*,*) "  -- TEST MMG3DLIB"

  ! Name and path of the mesh file
  CALL getenv("PWD",pwd)
  WRITE(filename,*) TRIM(pwd),"/../libexamples/mmg3d/example0_fortran/example0_a/cube"

  !> ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures
  !! args of InitMesh: mesh=&mmgMesh, sol=&mmgSol, input mesh name, input sol name,
  !! output mesh name
  mmgMesh = 0
  mmgSol  = 0
  !! Remark: %val(0) allow to pass the value 0 (i.e. NULL) instead of a pointer
  !! toward NULL.
  CALL MMG5_Init_mesh(mmgMesh,mmgSol,%val(0))

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG5_loadMesh function that will read a .mesh(b)
  !!    file formatted or manually set your mesh using the MMG5_Set* functions

  !> with MMG5_loadMesh function
  !! a) (not mandatory): give the mesh name
  !!   (by default, the "mesh.mesh" file is oppened)
  CALL MMG5_Set_inputMeshName(mmgMesh,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier == 0 ) THEN
     CALL EXIT(101)
  ENDIF

  !> b) function calling
  CALL MMG5_loadMesh(mmgMesh,ier)
  IF ( ier == 0 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG5_loadMet function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG5_Set* functions

  !> With MMG5_loadMet function
  !! a) (not mandatory): give the sol name
  !!   (by default, the "mesh.sol" file is oppened)
  CALL MMG5_Set_inputSolName(mmgMesh,mmgSol,TRIM(ADJUSTL(filename)),&
       LEN(TRIM(ADJUSTL(filename))),ier)
  IF ( ier ==0 ) THEN
     CALL EXIT(103)
  ENDIF

  !> b) function calling
  CALL MMG5_loadMet(mmgMesh,mmgSol,ier)
  IF ( ier ==0 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMG5_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier ==0 ) CALL EXIT(105)

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
  CALL MMG5_mmg3dlib(mmgMesh,mmgSol,ier)
  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG3DLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMG5_saveMesh/MMG5_saveMet functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMG5_getMesh/MMG5_getSol functions

  !> 1) Automatically save the mesh
  !! a)  (not mandatory): give the ouptut mesh name using MMG5_Set_outputMeshName
  !!   (by default, the mesh is saved in the "mesh.o.mesh" file
  !!call MMG5_Set_outputMeshName(mmgMesh,"output.mesh",len("output.mesh"),ier)
  !! b) function calling
  CALL MMG5_saveMesh(mmgMesh,ier)

  !> 2) Automatically save the solution
  !! a)  (not mandatory): give the ouptut sol name using MMG5_Set_outputSolName
  !!   (by default, the mesh is saved in the "mesh.o.sol" file
  !!call MMG5_Set_outputSolName(mmgMesh,mmgSol,"output.sol",len("output.sol"),ier)
  !! b) function calling
  CALL MMG5_saveMet(mmgMesh,mmgSol,ier)

  !> 3) Free the MMG3D5 structures
  CALL MMG5_Free_all(mmgMesh,mmgSol,%val(0))

END PROGRAM main
