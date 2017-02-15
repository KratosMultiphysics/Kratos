!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmg2dlib (basic use)

PROGRAM main
  IMPLICIT NONE

  !> Include here the mmg2d library hader file */
  ! if the header file is in the "include" directory
  ! #include "libmmg2df.h"
  ! if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2df.h"


  MMG5_DATA_PTR_T  :: mmgMesh
  MMG5_DATA_PTR_T  :: mmgSol
  INTEGER          :: ier,k

  !> To save final mesh in a file
  INTEGER          :: inm=10
  !> To manually recover the mesh
  INTEGER          :: np, nt, na, nc, nr, nreq, typEntity, typSol
  INTEGER          :: ref, Tria(3), Edge(2)
  DOUBLE PRECISION :: Point(2),Sol
  INTEGER, DIMENSION(:), ALLOCATABLE :: corner, required, ridge
  CHARACTER(LEN=31) :: FMT="(E14.8,1X,E14.8,1X,I3)"
  CHARACTER(len=300) :: filename

  PRINT*,"  -- TEST MMG2DLIB"

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

  !> Manually set of the mesh
  !! a) give the size of the mesh: 4 vertices, 2 triangles, 0 edges
  CALL MMG2D_Set_meshSize(mmgMesh,4,2,0,ier)
  IF ( ier /= 1 ) CALL EXIT(101)

  !> b) give the vertices: for each vertex, give the coordinates, the reference
  !!    and the position in mesh of the vertex
  !! Note that coordinates must be in double precision to match with the coordinate
  !! size in the C-library

  CALL MMG2D_Set_vertex(mmgMesh, 0.0D0, 0.0D0, 0,  1,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG2D_Set_vertex(mmgMesh, 1.0D0, 0.0D0, 0,  2,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG2D_Set_vertex(mmgMesh, 0.0D0, 1.0D0, 0,  3,ier)
  IF ( ier /= 1 ) CALL EXIT(102)
  CALL MMG2D_Set_vertex(mmgMesh, 1.0D0, 1.0D0, 0,  4,ier)
  IF ( ier /= 1 ) CALL EXIT(102)


  !> c) give the triangles: for each triangle,
  !!    give the vertices index, the reference and the position of the triangle
  CALL MMG2D_Set_triangle(mmgMesh,  1,  2,  3,1, 1,ier)
  IF ( ier /= 1 ) CALL EXIT(103)
  CALL MMG2D_Set_triangle(mmgMesh,  2,  4,  3,1, 2,ier)
  IF ( ier /= 1 ) CALL EXIT(103)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG2D_loadMet function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG2D_Set* functions

  !> Manually set of the sol
  !! a) give info for the sol structure: sol applied on vertex entities,
  !!    number of vertices=12, the sol is scalar
  CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,4,MMG5_Scalar,ier)
  IF ( ier /= 1 ) CALL EXIT(105)

  !> b) give solutions values and positions
  DO k=1,4
     CALL MMG2D_Set_scalarSol(mmgSol,0.1D0,k,ier)
     IF ( ier /= 1 ) CALL EXIT(106)
  ENDDO

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
  !! get results */
  !! Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
  !! that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !! using the MMG2D_getMesh/MMG2D_getSol functions

  filename="sortie.mesh"
  !> First solution: Manually get the mesh (in this example we show how to save the mesh
  !!    in the sortie.mesh file)
  OPEN(unit=inm,file=filename,form="formatted",status="replace")
  WRITE(inm,*),"MeshVersionFormatted 2"
  WRITE(inm,*),"Dimension 2"

  !> a) get the size of the mesh: vertices triangles, edges

  CALL MMG2D_Get_meshSize(mmgMesh,np,nt,na,ier)
  IF ( ier /= 1 ) CALL EXIT(108)

  ! Table to know if a vertex is corner (comment if you don't want to get it)
  !!$ ALLOCATE(corner(np))

  ! Table to know if a vertex/tetra/tria/edge is required (comment if you don't
  ! want to get it)
  !!$ ALLOCATE(required(MAX(np,MAX(nt,na))))

  ! Table to know if a conponant is corner and/or required (comment if you don't
  ! want to get it)
  !!$ ALLOCATE(ridge(na))

  nreq = 0; nc = 0
  WRITE(inm,*)
  WRITE(inm,*),"Vertices"
  WRITE(inm,*),np

  DO k=1, np
     !> b) Vertex recovering
     !! Note that coordinates must be in double precision to match with the coordinate
     !! size in the C-library

     !!$ CALL MMG2D_Get_vertex(mmgMesh,Point(1),Point(2),&
     !!$      ref,corner(k),required(k),ier)

     ! Or, if you don't want to recover the corners and required points:
     CALL MMG2D_Get_vertex(mmgMesh,Point(1),Point(2),&
          ref,%val(0),%val(0),ier)


     IF ( ier /= 1 ) CALL EXIT(109)

     WRITE(inm,FMT),Point(1),Point(2),ref
     ! To comment if you don't want to get the corners/required
     !!$ IF ( corner(k)/=0 )  nc=nc+1
     !!$ IF ( required(k)/=0 )  nreq=nreq+1
  ENDDO

  ! Comment if you don't want to get the corners
  !!$ WRITE(inm,*)
  !!$ WRITE(inm,*),"Corners"
  !!$ WRITE(inm,*), nc
  !!$
  !!$ DO k=1, np
  !!$   IF ( corner(k)/=0 )  WRITE(inm,*) ,k
  !!$ ENDDO
  !!$ WRITE(inm,*)
  !!$
  ! Comment if you don't want to get the required vertices
  !!$ WRITE(inm,*),"RequiredVertices"
  !!$ WRITE(inm,*), nreq
  !!$
  !!$ DO k=1,np
  !!$   IF ( required(k)/=0 ) WRITE(inm,*),k
  !!$ ENDDO
  !!$ WRITE(inm,*)
  !!$ DEALLOCATE(corner)

  nreq = 0;
  WRITE(inm,*),"Triangles"
  WRITE(inm,*),nt

  DO k=1,nt
    !> d) Triangles recovering
     !!$ CALL MMG2D_Get_triangle(mmgMesh,Tria(1),Tria(2),Tria(3),ref,required(k),ier)

     ! Or if you don't want to get the required triangles
     CALL MMG2D_Get_triangle(mmgMesh,Tria(1),Tria(2),Tria(3),ref,%val(0),ier)

     IF ( ier /= 1 ) CALL EXIT(110)
     WRITE(inm,*),Tria(1),Tria(2),Tria(3),ref
     ! Comment if you don't want to get the required triangles
     !!$ IF ( required(k)/=0 )  nreq=nreq+1;
  ENDDO
  WRITE(inm,*)

  ! Comment if you don't want to get the required triangles
  !!$ WRITE(inm,*),"RequiredTriangles"
  !!$ WRITE(inm,*),nreq
  !!$ DO k=1,nt
  !!$   IF ( required(k)/=0 ) WRITE(inm,*),k
  !!$ ENDDO
  !!$ WRITE(inm,*)

  nreq = 0;nr = 0;
  WRITE(inm,*),"Edges"
  WRITE(inm,*),na
  DO k=1,na
     !> e) Edges recovering
     !!$ CALL MMG2D_Get_edge(mmgMesh,Edge(1),Edge(2),ref,ridge(k),required(k),ier)
     ! Or, if you don't want to get the required and sharp edges:
     CALL MMG2D_Get_edge(mmgMesh,Edge(1),Edge(2),ref,%val(0),%val(0),ier)

     IF ( ier /= 1 ) CALL EXIT(111)
     WRITE(inm,*),Edge(1),Edge(2),ref

     ! Comment if you don't want to get the ridges
     !!$ IF ( ridge(k)/=0 )     nr = nr+1
     ! Comment if you don't want to get the required edges
     !!$ IF ( required(k)/=0 )  nreq = nreq+1
  ENDDO
  WRITE(inm,*)

  ! Comment if you don't want to get the required edges
  !!$ WRITE(inm,*),"RequiredEdges"
  !!$ WRITE(inm,*),nreq
  !!$ DO k=1,na
  !!$   IF ( required(k) /=0 ) WRITE(inm,*),k
  !!$ ENDDO
  !!$ WRITE(inm,*)

  ! Comment if you don't want to get the ridges
  !!$ WRITE(inm,*),"Ridges"
  !!$ WRITE(inm,*),nr
  !!$ DO k=1,na
  !!$   IF ( ridge(k) /=0 ) WRITE(inm,*),k
  !!$ ENDDO
  !!$ WRITE(inm,*)

  WRITE(inm,*),"End"
  CLOSE(inm)

  ! Comment if you don't want to get the required edges
  !!$ DEALLOCATE(required)
  ! Comment if you don't want to get the ridges
  !!$ DEALLOCATE(ridge)

  !> Second solution: Use the API function to save the mesh
  !  CALL MMG2D_SaveMesh(mmgMesh,TRIM(ADJUSTL(filename)), &
  !       LEN(TRIM(ADJUSTL(filename))),ier)
  !  IF ( ier /= 1 ) CALL EXIT(108)
  !
  !> 3) Free the MMG2D5 structures
  CALL MMG2D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main
