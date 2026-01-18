!> @author
!> Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!> Test IO API of Mmg2d



PROGRAM main
  IMPLICIT NONE

!> Include the mmg2d library hader file */
! if the header file is in the "include" directory
! #include "libmmg2df.h"
! if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2df.h"

  MMG5_DATA_PTR_T  :: mmgMesh
  MMG5_DATA_PTR_T  :: mmgSol
  INTEGER          :: ier,argc,by_array
  CHARACTER(len=300) :: exec_name,filein,fileout,option
  !> to cast integers into MMG5F_INT integers
  INTEGER,PARAMETER :: immg = MMG5F_INT

  PRINT*,"  -- TEST MMG2DLIB"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)

  IF ( argc /=3 ) THEN
     PRINT*," Usage: ",TRIM(exec_name)," input_file_name output_file_name test_array_api"
     PRINT*,"     test_array_api = 0 to Get/Set the mesh field by field, 1 by array"
     CALL EXIT(1);
  ENDIF

 ! Name and path of the mesh file
  CALL get_command_argument(1, filein)
  CALL get_command_argument(2, fileout)
  CALL get_command_argument(3, option)

  READ(option, '(I2)') by_array

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

  call loadmesh(mmgMesh,trim(filein),by_array)

  !> ------------------------------ STEP  II --------------------------
  !! Disable corner detection to check the corner setting
  call mmg2d_set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_angle,0_immg,ier)

  !! remesh function
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
  !! save results */
  call writemesh(mmgMesh,trim(fileout),by_array)

  !> 3) Free the MMG2D structures
  CALL MMG2D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main

subroutine loadmesh(mesh,filename,by_array)

  implicit NONE

  MMG5_DATA_PTR_T, intent(inout) :: mesh
  character(len=*), intent(in)   :: filename
  integer, intent(in) :: by_array

  INTEGER :: ier
  INTEGER(MMG5F_INT) :: i,np,nc,ne,nt,nq
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vert
  INTEGER(MMG5F_INT), dimension(:,:), allocatable :: edg,tri,quad
  INTEGER(MMG5F_INT), dimension(:)  , allocatable :: vert_ref,edg_ref,tri_ref,quad_ref,corner


  !! Parse file
  open(12, file = TRIM(ADJUSTL(filename)), status = 'old')


  do i=1,3
     READ(12,*) ! MeshVersionFormatted 2
  ENDDO

  do i=1,3
     READ(12,*) ! Dimension
  ENDDO

  READ(12,*) !Vertices
  READ(12,*) np
  allocate(vert(2,np))
  allocate(vert_ref(np))
  do i=1,np
     READ(12,*) vert(:,i),vert_ref(i)
  ENDDO

  READ(12,*)
  READ(12,*)
  READ(12,*) !Corners
  READ(12,*) nc
  allocate(corner(nc))
  do i=1,nc
     READ(12,*) corner(i)
  ENDDO

  READ(12,*)
  READ(12,*)
  READ(12,*) !Edges
  READ(12,*) ne
  allocate(edg(2,ne))
  allocate(edg_ref(ne))
  do i=1,ne
     READ(12,*) edg(:,i), edg_ref(i)
  ENDDO

  READ(12,*)
  READ(12,*)
  READ(12,*) !Tria
  READ(12,*) nt
  allocate(tri(3,nt))
  allocate(tri_ref(nt))
  do i=1,nt
     READ(12,*) tri(:,i), tri_ref(i)
  ENDDO

  READ(12,*)
  READ(12,*)
  READ(12,*) !Quads
  READ(12,*) nq
  allocate(quad(4,nq))
  allocate(quad_ref(nq))
  do i=1,nq
     READ(12,*) quad(:,i), quad_ref(i)
  ENDDO

  close(12)

  !! Set mesh
  CALL MMG2D_Set_meshSize(mesh,np,nt,nq,ne,ier)
  IF ( ier /= 1 ) CALL EXIT(101)

  if ( 0 == by_array ) then
     do i= 1,np
        CALL MMG2D_Set_vertex(mesh, vert(1,i), vert(2,i), vert_ref(i),  i,ier)
        IF ( ier /= 1 ) THEN
           print*, "Fail to set vertex"
           CALL EXIT(102)
        ENDIF
     ENDDO

     do i= 1,nc
        CALL MMG2D_Set_corner(mesh, corner(i),ier)
        IF ( ier /= 1 ) THEN
           print*, "Fail to set corner"
           CALL EXIT(103)
        ENDIF
     ENDDO

     do i= 1,ne
        CALL MMG2D_Set_edge(mesh, edg(1,i), edg(2,i), edg_ref(i),  i,ier)
        IF ( ier /= 1 ) THEN
           print*, "Fail to set edge"
           CALL EXIT(102)
        ENDIF
     ENDDO

     do i= 1,nt
        CALL MMG2D_Set_triangle(mesh, tri(1,i), tri(2,i),tri(3,i),tri_ref(i),  i,ier)
        IF ( ier /= 1 ) THEN
           print*, "Fail to set tria"
           CALL EXIT(102)
        ENDIF
     ENDDO

     do i= 1,nq
        CALL MMG2D_Set_quadrilateral(mesh,quad(1,i),quad(2,i),quad(3,i),quad(4,i),quad_ref(i),  i,ier)
        IF ( ier /= 1 ) THEN
           print*, "Fail to set quad"
           CALL EXIT(102)
        ENDIF
     ENDDO
  ELSE
     CALL MMG2D_Set_vertices(mesh, vert, vert_ref,ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to set vertices"
        CALL EXIT(102)
     ENDIF

     do i= 1,nc
        CALL MMG2D_Set_corner(mesh, corner(i),ier)
        IF ( ier /= 1 ) THEN
           print*, "Fail to set corner"
           CALL EXIT(103)
        ENDIF
     ENDDO

     CALL MMG2D_Set_edges(mesh, edg, edg_ref,ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to set edges"
        CALL EXIT(102)
     ENDIF

     CALL MMG2D_Set_triangles(mesh, tri, tri_ref,ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to set trias"
        CALL EXIT(102)
     ENDIF

     CALL MMG2D_Set_quadrilaterals(mesh,quad,quad_ref,ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to set quads"
        CALL EXIT(102)
     ENDIF
  ENDIF

  DEALLOCATE(vert,vert_ref,corner,edg,edg_ref,tri,tri_ref,quad,quad_ref)

end subroutine loadmesh

subroutine writemesh(mesh,filename,by_array)

  MMG5_DATA_PTR_T, intent(inout) :: mesh
  character(len=*), intent(in)   :: filename
  integer, intent(in) :: by_array

  INTEGER(MMG5F_INT) :: k,np,nc,ne,nt,nq,nreq
  INTEGER :: inm=13,ier
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vert
  INTEGER(MMG5F_INT), dimension(:,:), allocatable :: edg,tri,quad
  INTEGER(MMG5F_INT), dimension(:)  , allocatable :: vert_ref,edg_ref,tri_ref,quad_ref
  INTEGER, DIMENSION(:), ALLOCATABLE :: corner,required

  !> a) get the size of the mesh: vertices, tetra, triangles, edges
  CALL MMG2D_Get_meshSize(mesh,np,nt,nq,ne,ier)
  IF ( ier /= 1 ) CALL EXIT(108)

  allocate(vert(2,np))
  allocate(vert_ref(np))
  ! Table to know if a vertex is corner
  ALLOCATE(corner(np))

  ! Table to know if a vertex/tetra/tria/edge is required
  ALLOCATE(required(np))

  allocate(edg(2,ne))
  allocate(edg_ref(ne))

  allocate(tri(3,nt))
  allocate(tri_ref(nt))

  allocate(quad(4,nq))
  allocate(quad_ref(nq))


  nreq = 0; nc = 0
  WRITE(inm,*) "Vertices"
  WRITE(inm,*) np

  if ( 0== by_array ) then
     DO k=1, np
        !> b) Vertex recovering
        !! Note that coordinates must be in double precision to match with the coordinate
        !! size in the C-library
        CALL MMG2D_Get_vertex(mesh,vert(1,k),vert(2,k),&
             vert_ref(k),corner(k),required(k),ier)
        IF ( ier /= 1 ) CALL EXIT(109)
        IF ( corner(k)/=0 )  nc=nc+1
        IF ( required(k)/=0 )  nreq=nreq+1
     ENDDO

     DO k=1, ne
        CALL MMG2D_Get_edge(mesh,edg(1,k),edg(2,k),&
             edg_ref(k),%val(0),%val(0),ier)
        IF ( ier /= 1 ) CALL EXIT(110)
     ENDDO

     DO k=1, nt
        CALL MMG2D_Get_triangle(mesh,tri(1,k),tri(2,k),tri(3,k),&
             tri_ref(k),%val(0),ier)
        IF ( ier /= 1 ) CALL EXIT(111)
     ENDDO

     DO k=1, nq
        CALL MMG2D_Get_quadrilateral(mesh,quad(1,k),quad(2,k),quad(3,k),quad(4,k),&
             quad_ref(k),%val(0),ier)
        IF ( ier /= 1 ) CALL EXIT(112)
     ENDDO
  ELSE

     CALL MMG2D_Get_vertices(mesh, vert, vert_ref,corner,required,ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to get vertices"
        CALL EXIT(102)
     ENDIF
     !! count corners and required points
     DO k=1, np
        IF ( corner(k)/=0 )  nc=nc+1
        IF ( required(k)/=0 )  nreq=nreq+1
     ENDDO

     CALL MMG2D_Get_edges(mesh, edg, edg_ref,%val(0),%val(0),ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to get edges"
        CALL EXIT(102)
     ENDIF

     CALL MMG2D_Get_triangles(mesh, tri, tri_ref,%val(0),ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to get trias"
        CALL EXIT(102)
     ENDIF

     CALL MMG2D_Get_quadrilaterals(mesh,quad,quad_ref,%val(0),ier)
     IF ( ier /= 1 ) THEN
        print*, "Fail to get quads"
        CALL EXIT(102)
     ENDIF

  ENDIf



  OPEN(unit=inm,file=TRIM(ADJUSTL(filename))//".mesh",form="formatted",status="replace")
  WRITE(inm,*) "MeshVersionFormatted 2"
  WRITE(inm,*)
  WRITE(inm,*)

  WRITE(inm,*) "Dimension 2"
  WRITE(inm,*)
  WRITE(inm,*)

  WRITE(inm,*) "Vertices"
  WRITE(inm,*) np

  DO k=1, np
     WRITE(inm,*) vert(:,k),vert_ref(k)
  ENDDO

  WRITE(inm,*)
  WRITE(inm,*)
  WRITE(inm,*) "Corners"
  WRITE(inm,*)  nc
  DO k=1, np
    IF ( corner(k)/=0 )  WRITE(inm,*) k
  ENDDO

  WRITE(inm,*)
  WRITE(inm,*)
  WRITE(inm,*) "RequiredVertices"
  WRITE(inm,*)  nreq
  DO k=1,np
    IF ( required(k)/=0 ) WRITE(inm,*) k
  ENDDO
  DEALLOCATE(corner,required)


  WRITE(inm,*)
  WRITE(inm,*)
  WRITE(inm,*) "Edges"
  WRITE(inm,*) ne
  DO k=1,ne
     WRITE(inm,*) Edg(:,k),Edg_ref(k)
  ENDDO

  WRITE(inm,*)
  WRITE(inm,*)
  WRITE(inm,*) "Triangles"
  WRITE(inm,*) nt
  DO k=1,nt
     WRITE(inm,*) Tri(:,k),tri_ref(k)
  ENDDO

  WRITE(inm,*)
  WRITE(inm,*)
  WRITE(inm,*) "Quadrilaterals"
  WRITE(inm,*) nq
  DO k=1,nq
     WRITE(inm,*) quad(:,k),quad_ref(k)
  ENDDO

  WRITE(inm,*) "End"
  CLOSE(inm)

  DEALLOCATE(vert,vert_ref,edg,edg_ref,tri,tri_ref,quad,quad_ref)


end subroutine writemesh
