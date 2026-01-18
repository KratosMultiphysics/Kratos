!!> @author
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

  MMG5_DATA_PTR_T  :: mmgMesh
  MMG5_DATA_PTR_T  :: mmgSol
  INTEGER          :: ier,argc
  CHARACTER(len=300) :: exec_name,fileout

  !> To save final mesh in a file
  INTEGER           :: inm=10
  !> To manually recover the mesh
  INTEGER(MMG5F_INT):: k, np, nt, na, nc, nr, nreq
  INTEGER           :: typEntity, typSol
  INTEGER(MMG5F_INT):: ref, Tria(3), Edge(2)
  DOUBLE PRECISION  :: Point(3),Sol
  INTEGER, DIMENSION(:), ALLOCATABLE :: corner, required, ridge
  CHARACTER(LEN=31) :: FMT="(E14.8,1X,E14.8,1X,E14.8,1X,I3)"
  !> to cast integers into MMG5F_INT integers
  INTEGER,PARAMETER :: immg = MMG5F_INT

  PRINT*,"  -- TEST MMGSLIB"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)

  IF ( argc /=1 ) THEN
     PRINT*," Usage: ",TRIM(exec_name)," output_file_name"
     CALL EXIT(1);
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, fileout)

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

  CALL MMGS_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMGS_loadMesh function that will read a .mesh(b)
  !! file formatted or manually set your mesh using the MMGS_Set* functions

  !> Manually set of the mesh
  !! a) give the size of the mesh: 12 vertices, 20 triangles, 0 edges
  np = 12
  nt = 20
  na = 0
  CALL MMGS_Set_meshSize(mmgMesh,np,nt,na,ier)
  IF ( ier /= 1 ) CALL EXIT(101)

  !> b) give the vertices: for each vertex, give the coordinates, the reference
  !!    and the position in mesh of the vertex
  !! Note that coordinates must be in double precision to match with the coordinate
  !! size in the C-library
  ref = 0
  CALL MMGS_Set_vertex(mmgMesh, 0.0D0, 0.0D0, 0.0D0, ref, 1_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 0.5D0, 0.0D0, 0.0D0, ref, 2_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 0.5D0, 0.0D0, 1.0D0, ref, 3_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 0.0D0, 0.0D0, 1.0D0, ref, 4_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 0.0D0, 1.0D0, 0.0D0, ref, 5_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 0.5D0, 1.0D0, 0.0D0, ref, 6_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 0.5D0, 1.0D0, 1.0D0, ref, 7_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 0.0D0, 1.0D0, 1.0D0, ref, 8_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 1.0D0, 0.0D0, 0.0D0, ref, 9_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 1.0D0, 1.0D0, 0.0D0, ref, 10_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 1.0D0, 0.0D0, 1.0D0, ref, 11_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMGS_Set_vertex(mmgMesh, 1.0D0, 1.0D0, 1.0D0, ref, 12_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(102)

  !> d) give the triangles (not mandatory): for each triangle,
  !!    give the vertices index, the reference and the position of the triangle
  CALL MMGS_Set_triangle(mmgMesh,  1_immg,  4_immg,  8_immg, 3_immg, 1_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  1_immg,  2_immg,  4_immg, 3_immg, 2_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  8_immg,  3_immg,  7_immg, 3_immg, 3_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  5_immg,  8_immg,  6_immg, 3_immg, 4_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  5_immg,  6_immg,  2_immg, 3_immg, 5_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  5_immg,  2_immg,  1_immg, 3_immg, 6_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  5_immg,  1_immg,  8_immg, 3_immg, 7_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  7_immg,  6_immg,  8_immg, 3_immg, 8_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  4_immg,  3_immg,  8_immg, 3_immg, 9_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  2_immg,  3_immg,  4_immg, 3_immg,10_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  9_immg,  3_immg,  2_immg, 4_immg,11_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh, 11_immg,  9_immg, 12_immg, 4_immg,12_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  7_immg, 11_immg, 12_immg, 4_immg,13_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  6_immg,  7_immg, 10_immg, 4_immg,14_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  6_immg, 10_immg,  9_immg, 4_immg,15_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  6_immg,  9_immg,  2_immg, 4_immg,16_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh, 12_immg, 10_immg,  7_immg, 4_immg,17_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh, 12_immg,  9_immg, 10_immg, 4_immg,18_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  3_immg, 11_immg,  7_immg, 4_immg,19_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMGS_Set_triangle(mmgMesh,  9_immg, 11_immg,  3_immg, 4_immg,20_immg,ier)
  IF ( ier /= 1 ) CALL EXIT(104)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMGS_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMGS_Set* functions

  !> Manually set of the sol
  !! a) give info for the sol structure: sol applied on vertex entities,
  !!    number of vertices=12, the sol is scalar
  CALL MMGS_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar,ier)
  IF ( ier /= 1 ) CALL EXIT(105)

  !> b) give solutions values and positions
  DO k=1,12
     CALL MMGS_Set_scalarSol(mmgSol,0.5D0,k,ier)
     IF ( ier /= 1 ) CALL EXIT(106)
  ENDDO

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMGS_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) CALL EXIT(107)

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
  CALL MMGS_mmgslib(mmgMesh,mmgSol,ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH"
    STOP MMG5_STRONGFAILURE
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMGSLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results */
  !! Two solutions: just use the MMGS_saveMesh/MMGS_saveSol functions
  !! that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !! using the MMGS_getMesh/MMGS_getSol functions

  !> 1) Manually get the mesh
  OPEN(unit=inm,file=TRIM(ADJUSTL(fileout))//".mesh",form="formatted",status="replace")
  WRITE(inm,*) "MeshVersionFormatted 2"
  WRITE(inm,*) "Dimension 3"

  !> a) get the new size of the mesh: vertices, triangles, edges (overwriting of old var)
  CALL MMGS_Get_meshSize(mmgMesh,np,nt,na,ier)
  IF ( ier /= 1 ) CALL EXIT(108)

  ! Table to know if a vertex is corner
  ALLOCATE(corner(np))

  ! Table to know if a vertex/tetra/tria/edge is required
  ALLOCATE(required(MAX(np,MAX(nt,na))))

  ! Table to know if a coponant is corner and/or required
  ALLOCATE(ridge(na))

  nreq = 0; nc = 0
  WRITE(inm,*)
  WRITE(inm,*) "Vertices"
  WRITE(inm,*) np

  DO k=1, np
     !> b) Vertex recovering
     !! Note that coordinates must be in double precision to match with the coordinate
     !! size in the C-library
     CALL MMGS_Get_vertex(mmgMesh,Point(1),Point(2),Point(3),&
          ref,corner(k),required(k),ier)
     IF ( ier /= 1 ) CALL EXIT(109)

     WRITE(inm,FMT) Point(1),Point(2),Point(3),ref
     IF ( corner(k)/=0 )  nc=nc+1
     IF ( required(k)/=0 )  nreq=nreq+1
  ENDDO

  WRITE(inm,*)
  WRITE(inm,*) "Corners"
  WRITE(inm,*)  nc

  DO k=1, np
    IF ( corner(k)/=0 )  WRITE(inm,*) k
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*) "RequiredVertices"
  WRITE(inm,*)  nreq

  DO k=1,np
    IF ( required(k)/=0 ) WRITE(inm,*) k
  ENDDO
  WRITE(inm,*)
  DEALLOCATE(corner)

  nreq = 0;
  WRITE(inm,*) "Triangles"
  WRITE(inm,*) nt

  DO k=1,nt
    !> d) Triangles recovering
     CALL MMGS_Get_triangle(mmgMesh,Tria(1),Tria(2),Tria(3),ref,required(k),ier)
     IF ( ier /= 1 ) CALL EXIT(110)
     WRITE(inm,*) Tria(1),Tria(2),Tria(3),ref
     IF ( required(k)/=0 )  nreq=nreq+1;
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*) "RequiredTriangles"
  WRITE(inm,*) nreq
  DO k=1,nt
    IF ( required(k)/=0 ) WRITE(inm,*) k
  ENDDO
  WRITE(inm,*)

  nreq = 0;nr = 0;
  WRITE(inm,*) "Edges"
  WRITE(inm,*) na
  DO k=1,na
     !> e) Edges recovering
     CALL MMGS_Get_edge(mmgMesh,Edge(1),Edge(2),ref,ridge(k),required(k),ier)
     IF ( ier /= 1 ) CALL EXIT(111)
     WRITE(inm,*) Edge(1),Edge(2),ref
     IF ( ridge(k)/=0 )     nr = nr+1
     IF ( required(k)/=0 )  nreq = nreq+1
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*) "RequiredEdges"
  WRITE(inm,*) nreq
  DO k=1,na
    IF ( required(k) /=0 ) WRITE(inm,*) k
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*) "Ridges"
  WRITE(inm,*) nr
  DO k=1,na
    IF ( ridge(k) /=0 ) WRITE(inm,*) k
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*) "End"
  CLOSE(inm)

  DEALLOCATE(required)
  DEALLOCATE(ridge)

  !> 2) Manually get the solution
  OPEN(unit=inm,file=TRIM(ADJUSTL(fileout))//".sol", &
       & form="formatted",status="replace")
  WRITE(inm,*) "MeshVersionFormatted 2"
  WRITE(inm,*) "Dimension 3"
  WRITE(inm,*)

  !> a) get the size of the sol: type of entity (SolAtVertices,...),
  !!    number of sol, type of solution (scalar, tensor...)
  CALL MMGS_Get_solSize(mmgMesh,mmgSol,typEntity,np,typSol,ier)
  IF ( ier /= 1 ) CALL EXIT(113)

  IF ( ( typEntity /= MMG5_Vertex ) .OR. ( typSol /= MMG5_Scalar ) ) THEN
     CALL EXIT(114);
  ENDIF

  WRITE(inm,*) "SolAtVertices"
  WRITE(inm,*) np
  WRITE(inm,*) "1 1"
  WRITE(inm,*)
  DO k=1,np
    !> b) Vertex recovering
     CALL MMGS_Get_scalarSol(mmgSol,Sol,ier)
     IF ( ier /= 1 ) CALL EXIT(115)
     WRITE(inm,*) Sol
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*) "End"
  CLOSE(inm)

  !> 3) Free the MMGS5 structures
  CALL MMGS_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)
END PROGRAM main
