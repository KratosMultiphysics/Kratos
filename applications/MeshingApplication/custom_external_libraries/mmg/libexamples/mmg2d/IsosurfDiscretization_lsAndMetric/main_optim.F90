!! =============================================================================
!!  This file is part of the mmg software package for the tetrahedral
!!  mesh modification.
!!  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
!!
!!  mmg is free software: you can redistribute it and/or modify it
!!  under the terms of the GNU Lesser General Public License as published
!!  by the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  mmg is distributed in the hope that it will be useful, but WITHOUT
!!  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!!  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
!!  License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public
!!  License and of the GNU General Public License along with mmg (in
!!  files COPYING.LESSER and COPYING). If not, see
!!  <http://www.gnu.org/licenses/>. Please read their terms carefully and
!!  use this copy of the mmg distribution only if you accept them.
!! =============================================================================
!!

!>
!> Example of use of the mmg2dls function of the mmg2d library (basic use of
!> level-set discretization option): here the user only provide the level-set
!>
!> @author Charles Dapogny (LJLL, UPMC)
!> @author Pascal Frey (LJLL, UPMC)
!> @author Algiane Froehly (Inria / IMB, Université de Bordeaux)
!> @version 5
!> @copyright GNU Lesser General Public License.
!>

PROGRAM main

  IMPLICIT NONE

!> Include here the mmg2d library hader file
! if the header file is in the "include" directory
! #include "libmmg2df.h"

! if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2df.h"

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgLs,mmgMet
  INTEGER            :: ier,argc
  CHARACTER(len=300) :: exec_name,inname,outname,lsname
  !> to cast integers into MMG5F_INT integers
  INTEGER,PARAMETER :: immg = MMG5F_INT

  WRITE(*,*) "  -- TEST MMG2DLIB"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)

  IF ( argc /=3 ) THEN
     PRINT*," Usage: ",TRIM(exec_name)," meshfile lsfile meshout"
     CALL EXIT(1);
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, inname)
  CALL get_command_argument(2, lsname)
  CALL get_command_argument(3, outname)

  !> ------------------------------ STEP   I --------------------------
  !! 1) Initialisation of mesh and sol structures
  !! args of InitMesh:
  !! MMG5_ARG_start: we start to give the args of a variadic func
  !! MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
  !! mmgMesh: your MMG5_pMesh (that store your mesh)
  !! MMG5_ARG_ppLs: next arg will be a pointer over a MMG5_pSol storing a level-set
  !! &mmgLs: pointer toward your MMG5_pSol (that store your level-set)

  mmgMesh = 0
  mmgLs   = 0
  mmgMet  = 0

  CALL MMG2D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppLs,mmgLs, &
       MMG5_ARG_ppMet,mmgMet, &
       MMG5_ARG_end)

  !!------------------- Level set discretization option ---------------------
  ! Ask for level set discretization: note that it is important to do this step
  ! here because in iso mode, some filters are applied at mesh loading
  CALL MMG2D_Set_iparameter(mmgMesh,0_8,MMG2D_IPARAM_iso, 1_immg,ier)
  IF ( ier == 0 )  CALL EXIT(101)

  ! Ask for optim mode: compute the mean of input edge lengths.
  CALL MMG2D_Set_iparameter(mmgMesh,0_8,MMG2D_IPARAM_optim, 1_immg,ier)
  IF ( ier == 0 )  CALL EXIT(101)

  ! Ask to do this with anisotropic metric (unit tensor metric is computed).
  CALL MMG2D_Set_iparameter(mmgMesh,0_8,MMG2D_IPARAM_anisosize, 1_immg,ier)
  IF ( ier == 0 )  CALL EXIT(101)

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
  !!    file formatted or manually set your mesh using the MMG2D_Set* functions

  !> with MMG2D_loadMesh function
  CALL MMG2D_loadMesh(mmgMesh,TRIM(ADJUSTL(inname)),&
       LEN(TRIM(ADJUSTL(inname))),ier)
  IF ( ier == 0 )  CALL EXIT(102)

  !> 3) Build the level-set in MMG5 format
  !! Two solutions: just use the MMG2D_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your level-set using the MMG2D_Set* functions

  !> With MMG2D_loadSol function
  CALL MMG2D_loadSol(mmgMesh,mmgLs,TRIM(ADJUSTL(lsname)),&
       LEN(TRIM(ADJUSTL(lsname))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(104)
  ENDIF

  !> 4) (not mandatory): check the mesh and levelset sizes
  CALL MMG2D_Chk_meshData(mmgMesh,mmgLs,ier)
  IF ( ier /= 1 ) CALL EXIT(105)

  !> ------------------------------ STEP  II --------------------------
  !! isovalue discretization: as we don't want to impose an input metric we
  !! pass %val(0_8) instead of the metric structure as function argument
  CALL MMG2D_mmg2dls(mmgMesh,mmgLs,%val(0_8),ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH"
     STOP 2
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results
  !! Two solutions: just use the MMG2D_saveMesh function
  !!    that will write .mesh(b) formatted file or manually get your mesh
  !!    using the MMG2D_getMesh function

  !> 1) Automatically save the mesh
  CALL MMG2D_saveMesh(mmgMesh,TRIM(ADJUSTL(outname)),&
       LEN(TRIM(ADJUSTL(outname))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(106)
  ENDIF

  !> 2) (Not mandatory) Automatically save the output metric
  CALL MMG2D_saveSol(mmgMesh,mmgMet,TRIM(ADJUSTL(outname)),&
       LEN(TRIM(ADJUSTL(outname))),ier)
  IF ( ier /= 1 ) THEN
     CALL EXIT(107)
  ENDIF

   !> 3) Free the MMG2D5 structures
  CALL MMG2D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppLs,mmgLs, &
       MMG5_ARG_ppMet,mmgMet, &
       MMG5_ARG_end)

END PROGRAM
