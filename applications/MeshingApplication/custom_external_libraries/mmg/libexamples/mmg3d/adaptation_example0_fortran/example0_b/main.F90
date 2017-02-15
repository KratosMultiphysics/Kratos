!!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmg3dlib (basic use)

PROGRAM main

  IMPLICIT NONE

!> Include the mmg3d library hader file
! if the header file is in the "include" directory
! #include "libmmg3df.h"
! if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3df.h"

  MMG5_DATA_PTR_T  :: mmgMesh
  MMG5_DATA_PTR_T  :: mmgSol
  INTEGER          :: ier,k

  !> To save final mesh in a file
  INTEGER          :: inm=10
  !> To manually recover the mesh
  INTEGER          :: np, ne, nt, na, nc, nr, nreq, typEntity, typSol
  INTEGER          :: ref, Tetra(4), Tria(3), Edge(2)
  DOUBLE PRECISION :: Point(3),Sol
  INTEGER, DIMENSION(:), ALLOCATABLE :: corner, required, ridge
  CHARACTER(LEN=31) :: FMT="(E14.8,1X,E14.8,1X,E14.8,1X,I3)"

  PRINT*,"  -- TEST MMG3DLIB"

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

  CALL MMG3D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)
  CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_verbose,5,ier)

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
  !! file formatted or manually set your mesh using the MMG3D_Set* functions

  !> Manually set of the mesh
  !! a) give the size of the mesh: 12 vertices, 12 tetra,0 prisms, 20 triangles,
  !! 0 quads, 0 edges
  CALL MMG3D_Set_meshSize(mmgMesh,12,12,0,20,0,0,ier)
  IF ( ier /= 1 ) CALL EXIT(101)

  !> b) give the vertices: for each vertex, give the coordinates, the reference
  !!    and the position in mesh of the vertex
  !! Note that coordinates must be in double precision to match with the coordinate
  !! size in the C-library

  CALL MMG3D_Set_vertex(mmgMesh, 0.0D0, 0.0D0, 0.0D0, 0,  1,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 0.5D0, 0.0D0, 0.0D0, 0,  2,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 0.5D0, 0.0D0, 1.0D0, 0,  3,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 0.0D0, 0.0D0, 1.0D0, 0,  4,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 0.0D0, 1.0D0, 0.0D0, 0,  5,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 0.5D0, 1.0D0, 0.0D0, 0,  6,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 0.5D0, 1.0D0, 1.0D0, 0,  7,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 0.0D0, 1.0D0, 1.0D0, 0,  8,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 1.0D0, 0.0D0, 0.0D0, 0,  9,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 1.0D0, 1.0D0, 0.0D0, 0, 10,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 1.0D0, 0.0D0, 1.0D0, 0, 11,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG3D_Set_vertex(mmgMesh, 1.0D0, 1.0D0, 1.0D0, 0, 12,ier)
  IF ( ier /= 1 ) CALL EXIT(102)

  !> c) give the tetrahedras: for each tetrahedra,
  !!    give the vertices index, the reference and the position of the tetra
  CALL MMG3D_Set_tetrahedron(mmgMesh,  1,  4,  2,  8,1, 1,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  8,  3,  2,  7,1, 2,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  5,  2,  6,  8,1, 3,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  5,  8,  1,  2,1, 4,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  7,  2,  8,  6,1, 5,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  2,  4,  3,  8,1, 6,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  9,  2,  3,  7,2, 7,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  7, 11,  9, 12,2, 8,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  6,  9, 10,  7,2, 9,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  6,  7,  2,  9,2,10,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh, 12,  9,  7, 10,2,11,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG3D_Set_tetrahedron(mmgMesh,  9,  3, 11,  7,2,12,ier)
  IF ( ier /= 1 ) CALL EXIT(103)

  !> d) give the triangles (not mandatory): for each triangle,
  !!    give the vertices index, the reference and the position of the triangle
  CALL MMG3D_Set_triangle(mmgMesh,  1,  4,  8, 3, 1,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  1,  2,  4, 3, 2,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  8,  3,  7, 3, 3,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  5,  8,  6, 3, 4,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  5,  6,  2, 3, 5,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  5,  2,  1, 3, 6,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  5,  1,  8, 3, 7,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  7,  6,  8, 3, 8,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  4,  3,  8, 3, 9,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  2,  3,  4, 3,10,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  9,  3,  2, 4,11,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh, 11,  9, 12, 4,12,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  7, 11, 12, 4,13,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  6,  7, 10, 4,14,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  6, 10,  9, 4,15,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  6,  9,  2, 4,16,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh, 12, 10,  7, 4,17,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh, 12,  9, 10, 4,18,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  3, 11,  7, 4,19,ier)
  IF ( ier /= 1 ) CALL EXIT(104)
  CALL MMG3D_Set_triangle(mmgMesh,  9, 11,  3, 4,20,ier)
  IF ( ier /= 1 ) CALL EXIT(104)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG3D_loadSol function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG3D_Set* functions

  !> Manually set of the sol
  !! a) give info for the sol structure: sol applied on vertex entities,
  !!    number of vertices=12, the sol is scalar
  CALL MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,12,MMG5_Scalar,ier)
  IF ( ier /= 1 ) CALL EXIT(105)

  !> b) give solutions values and positions
  DO k=1,12
     CALL MMG3D_Set_scalarSol(mmgSol,0.5D0,k,ier)
     IF ( ier /= 1 ) CALL EXIT(106)
  ENDDO

  !> 4) (not mandatory): check if the number of given entities match with mesh size
  CALL MMG3D_Chk_meshData(mmgMesh,mmgSol,ier)
  IF ( ier /= 1 ) CALL EXIT(107)

  !> ------------------------------ STEP  II --------------------------
  !! remesh function
  CALL MMG3D_mmg3dlib(mmgMesh,mmgSol,ier)

  IF ( ier == MMG5_STRONGFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
    STOP MMG5_STRONGFAILURE
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG3DLIB"
  ENDIF

  !> ------------------------------ STEP III --------------------------
  !! get results */
  !! Two solutions: just use the MMG3D_saveMesh/MMG3D_saveSol functions
  !! that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !! using the MMG3D_getMesh/MMG3D_getSol functions

  !> 1) Manually get the mesh (in this example we show how to save the mesh
  !!    in the mesh.o.mesh file)
  OPEN(unit=inm,file="mesh.o.mesh",form="formatted",status="replace")
  WRITE(inm,*),"MeshVersionFormatted 2"
  WRITE(inm,*),"Dimension 3"

  !> a) get the size of the mesh: vertices, tetra,prisms, triangles, quads,edges
  CALL MMG3D_Get_meshSize(mmgMesh,np,ne,%val(0),nt,%val(0),na,ier)
  IF ( ier /= 1 ) CALL EXIT(108)

  ! Table to know if a vertex is corner
  ALLOCATE(corner(np))

  ! Table to know if a vertex/tetra/tria/edge is required
  ALLOCATE(required(MAX(MAX(np,ne),MAX(nt,na))))

  ! Table to know if a coponant is corner and/or required
  ALLOCATE(ridge(na))

  nreq = 0; nc = 0
  WRITE(inm,*)
  WRITE(inm,*),"Vertices"
  WRITE(inm,*),np

  DO k=1, np
     !> b) Vertex recovering
     !! Note that coordinates must be in double precision to match with the coordinate
     !! size in the C-library
     CALL MMG3D_Get_vertex(mmgMesh,Point(1),Point(2),Point(3),&
          ref,corner(k),required(k),ier)
     IF ( ier /= 1 ) CALL EXIT(109)

     WRITE(inm,FMT),Point(1),Point(2),Point(3),ref
     IF ( corner(k)/=0 )  nc=nc+1
     IF ( required(k)/=0 )  nreq=nreq+1
  ENDDO

  WRITE(inm,*)
  WRITE(inm,*),"Corners"
  WRITE(inm,*), nc

  DO k=1, np
    IF ( corner(k)/=0 )  WRITE(inm,*) ,k
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*),"RequiredVertices"
  WRITE(inm,*), nreq

  DO k=1,np
    IF ( required(k)/=0 ) WRITE(inm,*),k
  ENDDO
  WRITE(inm,*)
  DEALLOCATE(corner)

  nreq = 0;
  WRITE(inm,*),"Triangles"
  WRITE(inm,*),nt

  DO k=1,nt
    !> d) Triangles recovering
     CALL MMG3D_Get_triangle(mmgMesh,Tria(1),Tria(2),Tria(3),ref,required(k),ier)
     IF ( ier /= 1 ) CALL EXIT(110)
     WRITE(inm,*),Tria(1),Tria(2),Tria(3),ref
     IF ( required(k)/=0 )  nreq=nreq+1;
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*),"RequiredTriangles"
  WRITE(inm,*),nreq
  DO k=1,nt
    IF ( required(k)/=0 ) WRITE(inm,*),k
  ENDDO
  WRITE(inm,*)

  nreq = 0;nr = 0;
  WRITE(inm,*),"Edges"
  WRITE(inm,*),na
  DO k=1,na
     !> e) Edges recovering
     CALL MMG3D_Get_edge(mmgMesh,Edge(1),Edge(2),ref,ridge(k),required(k),ier)
     IF ( ier /= 1 ) CALL EXIT(111)
     WRITE(inm,*),Edge(1),Edge(2),ref
     IF ( ridge(k)/=0 )     nr = nr+1
     IF ( required(k)/=0 )  nreq = nreq+1
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*),"RequiredEdges"
  WRITE(inm,*),nreq
  DO k=1,na
    IF ( required(k) /=0 ) WRITE(inm,*),k
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*),"Ridges"
  WRITE(inm,*),nr
  DO k=1,na
    IF ( ridge(k) /=0 ) WRITE(inm,*),k
  ENDDO
  WRITE(inm,*)

  nreq = 0;
  WRITE(inm,*),"Tetrahedra"
  WRITE(inm,*),ne
  DO k=1,ne
    !> c) Tetra recovering
     CALL MMG3D_Get_tetrahedron(mmgMesh,Tetra(1),Tetra(2),Tetra(3),Tetra(4),&
          ref,required(k),ier)
    IF ( ier /= 1 ) CALL EXIT(112)
    WRITE(inm,*),Tetra(1),Tetra(2),Tetra(3),Tetra(4),ref
    IF ( required(k) /= 0 )  nreq = nreq+1
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*),"RequiredTetrahedra"
  WRITE(inm,*),nreq
  DO k=1,ne
    IF ( required(k) /= 0 ) WRITE(inm,*),k
  ENDDO

  WRITE(inm,*),"End"
  CLOSE(inm)

  DEALLOCATE(required)
  DEALLOCATE(ridge)

  !> 2) Manually get the solution (in this example we show how to save the
  !!    solution in the mesh.o.sol file)
  OPEN(unit=inm,file="mesh.o.sol",form="formatted",status="replace")
  WRITE(inm,*),"MeshVersionFormatted 2"
  WRITE(inm,*),"Dimension 3"
  WRITE(inm,*)

  !> a) get the size of the sol: type of entity (SolAtVertices,...),
  !!    number of sol, type of solution (scalar, tensor...)
  CALL MMG3D_Get_solSize(mmgMesh,mmgSol,typEntity,np,typSol,ier)
  IF ( ier /= 1 ) CALL EXIT(113)

  IF ( ( typEntity /= MMG5_Vertex ) .OR. ( typSol /= MMG5_Scalar ) ) THEN
     CALL EXIT(114);
  ENDIF

  WRITE(inm,*),"SolAtVertices"
  WRITE(inm,*),np
  WRITE(inm,*),"1 1"
  WRITE(inm,*)
  DO k=1,np
    !> b) Vertex recovering
     CALL MMG3D_Get_scalarSol(mmgSol,Sol,ier)
     IF ( ier /= 1 ) CALL EXIT(115)
     WRITE(inm,*),Sol
  ENDDO
  WRITE(inm,*)

  WRITE(inm,*),"End"
  CLOSE(inm)

  !> 3) Free the MMG3D5 structures
  CALL MMG3D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)
END PROGRAM main
