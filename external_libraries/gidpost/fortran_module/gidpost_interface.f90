MODULE gidpost

! *************************************************************************

! *************************************************************************
!
!   AUTHORS: Carlos Labra (clabra@cimne.upc.edu)
!
! *************************************************************************

! **************************** Edit History *******************************
!
! Sep-01-2009  Initial creation based on the Fortran 2003 standard
!              Compiler accepted: Intel (>=11.0), SunStudio (>=12.1)
! Feb-18-2010  Minor changes, fix special cases
!              Add Compiler: gfortran-4.4
! Feb-24-2010  Minor changes.
!              Add Compiler: g95 (not allow derived type with BIND(C) and PRIVATE)
! Apr-22-2010  Fortran string to C string now uses explicit temporal variable.
!              Avoid creation of temporal array in run-time.
! Jun-11-2012  Add GID_REAL_KIND and GID_INT_KIND types, and routines GID_INT, GID_REAL.
!              Compatibility of diff real and integer kinds (to ensure no conflict with arguments).
!              GiD routines must receive integer(GID_INT_KIND) and real(GID_REAL_KIND) types,
!              or converted with GID_INT() and GID_REAL()
! Mar-13-2013  GiD_PostHDF5 added (gidpost must be compiled with this option)
!              Added routines GiD_MeshUnit and GiD_ResultUnit
! Jul-01-2013  Solved bug with VALUE attributes on interfaces
!
! *************************************************************************

  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  PRIVATE

!--------------------------------------------------------------------------
! PUBLIC PRECISION PARAMETERS
!--------------------------------------------------------------------------

  INTEGER, PARAMETER, PUBLIC :: GID_REAL_KIND = C_DOUBLE
  INTEGER, PARAMETER, PUBLIC :: GID_INT_KIND  = C_INT

!--------------------------------------------------------------------------
! PRIVATE PARAMETERS
!--------------------------------------------------------------------------

  INTEGER(C_INT),           PARAMETER :: GID_DUMMY_NULL     = -1  ! default (null) value on GiD types
  INTEGER,                  PARAMETER :: GID_MAX_CHAR_LEN   = 64  ! max length of strings
  INTEGER(C_INT),           PARAMETER :: GID_FILE_NULL_UNIT = -1  ! default (null) value on GiD file unit
  INTEGER(C_INT),           PARAMETER :: GID_NO_FILE = 0  ! default (null) value on GiD file unit
  CHARACTER(1,KIND=C_CHAR), PARAMETER :: C_CHAR_ARRAY(1) = (/ C_NULL_CHAR /)

!--------------------------------------------------------------------------
! GID NULL STRING PARAMETER
!--------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: GID_STRING_NULL = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)

!--------------------------------------------------------------------------
! GID OUTPUT FILE TYPES
!--------------------------------------------------------------------------

  TYPE, PUBLIC, BIND(C) :: GiD_PostMode
    !PRIVATE
    INTEGER(C_INT) :: dummy = GID_DUMMY_NULL
  END TYPE

  TYPE(GiD_PostMode), PARAMETER, PUBLIC :: GiD_PostAscii       = GiD_PostMode(0) !ascii
  TYPE(GiD_PostMode), PARAMETER, PUBLIC :: GiD_PostAsciiZipped = GiD_PostMode(1) !compresed ascii
  TYPE(GiD_PostMode), PARAMETER, PUBLIC :: GiD_PostBinary      = GiD_PostMode(2) !compresed binary
  TYPE(GiD_PostMode), PARAMETER, PUBLIC :: GiD_PostHDF5        = GiD_PostMode(3) !HDF5

!--------------------------------------------------------------------------
! GID MESH DIMENSIONS
!--------------------------------------------------------------------------

  TYPE, PUBLIC, BIND(C) :: GiD_Dimension
    !PRIVATE
    INTEGER(C_INT) :: dummy = GID_DUMMY_NULL
  END TYPE

  TYPE(GiD_Dimension), PARAMETER, PUBLIC :: GiD_2D = GiD_Dimension(2)
  TYPE(GiD_Dimension), PARAMETER, PUBLIC :: GiD_3D = GiD_Dimension(3)

!--------------------------------------------------------------------------
! GID ELEMENT TYPES
!--------------------------------------------------------------------------

  TYPE, PUBLIC, BIND(C) :: GiD_ElementType
    !PRIVATE
    INTEGER(C_INT) :: dummy = GID_DUMMY_NULL
  END TYPE

  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_NoElement     = GiD_ElementType(0)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Point         = GiD_ElementType(1)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Linear        = GiD_ElementType(2)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Triangle      = GiD_ElementType(3)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Quadrilateral = GiD_ElementType(4)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Tetrahedra    = GiD_ElementType(5)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Hexahedra     = GiD_ElementType(6)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Prism         = GiD_ElementType(7)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Piramid       = GiD_ElementType(8)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Sphere        = GiD_ElementType(9)
  TYPE(GiD_ElementType), PARAMETER, PUBLIC :: GiD_Circle        = GiD_ElementType(10)

!--------------------------------------------------------------------------
! GID RESULT TYPES
!--------------------------------------------------------------------------

  TYPE, PUBLIC, BIND(C) :: GiD_ResType
    !PRIVATE
    INTEGER(C_INT) :: dummy = GID_DUMMY_NULL
  END TYPE

  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_Scalar         = GiD_ResType(0)
  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_Vector         = GiD_ResType(1)
  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_Matrix         = GiD_ResType(2)
  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_PlainDefMatrix = GiD_ResType(3)
  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_MainMatrix     = GiD_ResType(4)
  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_LocalAxes      = GiD_ResType(5)
  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_ComplexScalar  = GiD_ResType(6)
  TYPE(GiD_ResType), PARAMETER, PUBLIC :: GiD_ComplexVector  = GiD_ResType(7)

!--------------------------------------------------------------------------
! GID RESULT LOCATION TYPES
!--------------------------------------------------------------------------

  TYPE, PUBLIC, BIND(C) :: GiD_ResLoc
    !PRIVATE
    INTEGER(C_INT) :: dummy = GID_DUMMY_NULL
  END TYPE

  TYPE(GiD_ResLoc), PARAMETER, PUBLIC :: GiD_onNodes = GiD_ResLoc(0)
  TYPE(GiD_ResLoc), PARAMETER, PUBLIC :: GiD_onGaussPoint = GiD_ResLoc(1)

!--------------------------------------------------------------------------
! GID FILE TYPE
!--------------------------------------------------------------------------

  TYPE, PUBLIC, BIND(C) :: GiD_File
    !PRIVATE
    INTEGER(C_INT) :: file_unit  = GID_FILE_NULL_UNIT
  END TYPE GiD_File
  
!--------------------------------------------------------------------------
! AUXILIAR VARIABLES
!--------------------------------------------------------------------------

  CHARACTER(LEN=1,KIND=C_CHAR), TARGET :: TempStringC(1024)
  INTEGER(C_INT), PARAMETER :: MAX_STRING_LEN = 1024

!--------------------------------------------------------------------------
! GENERIC INTERFACES
!--------------------------------------------------------------------------

  INTERFACE
  
!--------------------------------------------------------------------------
! GLOBAL ROUTINES
!--------------------------------------------------------------------------

    SUBROUTINE GiD_PostInit() BIND(C,NAME='GiD_PostInit')
    END SUBROUTINE GiD_PostInit
! ---
    SUBROUTINE GiD_PostDone() BIND(C,NAME='GiD_PostDone')
    END SUBROUTINE GiD_PostDone

!--------------------------------------------------------------------------
! MESH ROUTINES
!--------------------------------------------------------------------------

    !SUBROUTINE GiD_OpenPostMeshFile -> REDEFINED
    !SUBROUTINE GiD_fOpenPostMeshFile -> REDEFINED
! ---
    SUBROUTINE GiD_ClosePostMeshFile() BIND(C,NAME='GiD_ClosePostMeshFile')
    END SUBROUTINE GiD_ClosePostMeshFile
    SUBROUTINE GiD_fClosePostMeshFile( fd ) BIND(C,NAME='GiD_fClosePostMeshFile')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fClosePostMeshFile
! ---
    !SUBROUTINE GiD_BeginMeshGroup -> REDEFINED
    !SUBROUTINE GiD_fBeginMeshGroup -> REDEFINED             
! ---
    SUBROUTINE GiD_EndMeshGroup() BIND(C,NAME='GiD_EndMeshGroup')
    END SUBROUTINE GiD_EndMeshGroup
    SUBROUTINE GiD_fEndMeshGroup(fd) BIND(C,NAME='GiD_fEndMeshGroup')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndMeshGroup
! ---
    !SUBROUTINE GiD_MeshUnit -> REDEFINED
    !SUBROUTINE GiD_fMeshUnit -> REDEFINED
! ---
    !SUBROUTINE GiD_BeginMesh -> REDEFINED
    !SUBROUTINE GiD_fBeginMesh -> REDEFINED
! ---
    !SUBROUTINE GiD_BeginMeshColor -> REDEFINED
    !SUBROUTINE GiD_fBeginMeshColor -> REDEFINED
! ---
    SUBROUTINE GiD_EndMesh() BIND(C,NAME='GiD_EndMesh')
    END SUBROUTINE GiD_EndMesh
    SUBROUTINE GiD_fEndMesh(fd) BIND(C,NAME='GiD_fEndMesh')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndMesh
! ---
    SUBROUTINE GiD_BeginCoordinates() BIND(C,NAME='GiD_BeginCoordinates')
    END SUBROUTINE GiD_BeginCoordinates
    SUBROUTINE GiD_fBeginCoordinates(fd) BIND(C,NAME='GiD_fBeginCoordinates')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fBeginCoordinates
! ---
    SUBROUTINE GiD_EndCoordinates() BIND(C,NAME='GiD_EndCoordinates')
    END SUBROUTINE GiD_EndCoordinates
    SUBROUTINE GiD_fEndCoordinates(fd) BIND(C,NAME='GiD_fEndCoordinates')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndCoordinates
! ---
    SUBROUTINE GiD_WriteCoordinates2D(id,x,y) BIND(C,NAME='GiD_WriteCoordinates2D')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y
    END SUBROUTINE GiD_WriteCoordinates2D
    SUBROUTINE GiD_fWriteCoordinates2D(fd,id,x,y) BIND(C,NAME='GiD_fWriteCoordinates2D')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y
    END SUBROUTINE GiD_fWriteCoordinates2D
! ---
    ! SUBROUTINE GiD_WriteCoordinates  -> SPECIAL INTERFACE
    ! SUBROUTINE GiD_fWriteCoordinates  -> SPECIAL INTERFACE
! ---
    SUBROUTINE GiD_BeginElements() BIND(C,NAME='GiD_BeginElements')
    END SUBROUTINE GiD_BeginElements
    SUBROUTINE GiD_fBeginElements(fd) BIND(C,NAME='GiD_fBeginElements')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fBeginElements
! ---
    SUBROUTINE GiD_EndElements() BIND(C,NAME='GiD_EndElements')
    END SUBROUTINE GiD_EndElements
    SUBROUTINE GiD_fEndElements(fd) BIND(C,NAME='GiD_fEndElements')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndElements
! ---
    SUBROUTINE GiD_WriteElement(id,conec) BIND(C,NAME='GiD_WriteElement')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT),        INTENT(IN) :: conec(*)
    END SUBROUTINE GiD_WriteElement
    SUBROUTINE GiD_fWriteElement(fd,id,conec) BIND(C,NAME='GiD_fWriteElement')
      IMPORT :: C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT),        INTENT(IN) :: conec(*)
    END SUBROUTINE GiD_fWriteElement
! ---
    !SUBROUTINE GiD_WriteElementMat(id,conec) -> SPECIAL INTERFACE
    !SUBROUTINE GiD_fWriteElementMat(id,conec) -> SPECIAL INTERFACE
! ---
    SUBROUTINE GiD_WriteCircle(id,conec,radius,normal_x,normal_y,normal_z) BIND(C,NAME='GiD_WriteCircle')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: normal_x, normal_y, normal_z
    END SUBROUTINE GiD_WriteCircle
    SUBROUTINE GiD_fWriteCircle(fd,id,conec,radius,normal_x,normal_y,normal_z) BIND(C,NAME='GiD_fWriteCircle')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: normal_x, normal_y, normal_z
    END SUBROUTINE GiD_fWriteCircle
! ---
    SUBROUTINE GiD_WriteCircleMat(id,conec,radius,normal_x,normal_y,normal_z,material) BIND(C,NAME='GiD_WriteCircleMat')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: normal_x, normal_y, normal_z
      INTEGER(C_INT), VALUE, INTENT(IN) :: material
    END SUBROUTINE GiD_WriteCircleMat
    SUBROUTINE GiD_fWriteCircleMat(fd,id,conec,radius,normal_x,normal_y,normal_z,material) BIND(C,NAME='GiD_fWriteCircleMat')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: normal_x, normal_y, normal_z
      INTEGER(C_INT), VALUE, INTENT(IN) :: material
    END SUBROUTINE GiD_fWriteCircleMat
! ---
    SUBROUTINE GiD_WriteSphere(id,conec,radius) BIND(C,NAME='GiD_WriteSphere')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
    END SUBROUTINE GiD_WriteSphere
    SUBROUTINE GiD_fWriteSphere(fd,id,conec,radius) BIND(C,NAME='GiD_fWriteSphere')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
    END SUBROUTINE GiD_fWriteSphere
! ---
    SUBROUTINE GiD_WriteSphereMat(id,conec,radius,material) BIND(C,NAME='GiD_WriteSphereMat')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
      INTEGER(C_INT), VALUE, INTENT(IN) :: material
    END SUBROUTINE GiD_WriteSphereMat
    SUBROUTINE GiD_fWriteSphereMat(fd,id,conec,radius,material) BIND(C,NAME='GiD_fWriteSphereMat')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT), VALUE, INTENT(IN) :: conec
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: radius
      INTEGER(C_INT), VALUE, INTENT(IN) :: material
    END SUBROUTINE GiD_fWriteSphereMat

!--------------------------------------------------------------------------
! RESULTS ROUTINES
!--------------------------------------------------------------------------

    !SUBROUTINE GiD_OpenPostResultFile -> REDEFINED
    !SUBROUTINE GiD_fOpenPostResultFile -> REDEFINED

    SUBROUTINE GiD_ClosePostResultFile() BIND(C,NAME='GiD_ClosePostResultFile')
    END SUBROUTINE GiD_ClosePostResultFile
    SUBROUTINE GiD_fClosePostResultFile(fd) BIND(C,NAME='GiD_fClosePostResultFile')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fClosePostResultFile

    !SUBROUTINE GiD_BeginGaussPoint -> REDEFINED
    !SUBROUTINE GiD_fBeginGaussPoint -> REDEFINED

    SUBROUTINE GiD_EndGaussPoint() BIND(C,NAME='GiD_EndGaussPoint')
    END SUBROUTINE GiD_EndGaussPoint
    SUBROUTINE GiD_fEndGaussPoint(fd) BIND(C,NAME='GiD_fEndGaussPoint')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndGaussPoint

    SUBROUTINE GiD_WriteGaussPoint2D(x,y) BIND(C,NAME='GiD_WriteGaussPoint2D')
      IMPORT :: C_DOUBLE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y
    END SUBROUTINE GiD_WriteGaussPoint2D
    SUBROUTINE GiD_fWriteGaussPoint2D(fd,x,y) BIND(C,NAME='GiD_fWriteGaussPoint2D')
      IMPORT :: C_DOUBLE, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y
    END SUBROUTINE GiD_fWriteGaussPoint2D

    SUBROUTINE GiD_WriteGaussPoint3D(x,y,z) BIND(C,NAME='GiD_WriteGaussPoint3D')
      IMPORT :: C_DOUBLE
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z
    END SUBROUTINE GiD_WriteGaussPoint3D
    SUBROUTINE GiD_fWriteGaussPoint3D(fd,x,y,z) BIND(C,NAME='GiD_fWriteGaussPoint3D')
      IMPORT :: C_DOUBLE, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z
    END SUBROUTINE GiD_fWriteGaussPoint3D

    !SUBROUTINE GiD_BeginRangeTable -> REDEFINED
    !SUBROUTINE GiD_fBeginRangeTable -> REDEFINED

    SUBROUTINE GiD_EndRangeTable() BIND(C,NAME='GiD_EndRangeTable')
    END SUBROUTINE GiD_EndRangeTable
    SUBROUTINE GiD_fEndRangeTable(fd) BIND(C,NAME='GiD_fEndRangeTable')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndRangeTable

    !SUBROUTINE GiD_WriteMinRange -> REDEFINED
    !SUBROUTINE GiD_fWriteMinRange -> REDEFINED

    !SUBROUTINE GiD_WriteRange -> REDEFINED
    !SUBROUTINE GiD_fWriteRange -> REDEFINED

    !SUBROUTINE GiD_WriteMaxRange -> REDEFINED
    !SUBROUTINE GiD_fWriteMaxRange -> REDEFINED

    !SUBROUTINE GiD_BeginScalarResult -> REDEFINED
    !SUBROUTINE GiD_fBeginScalarResult -> REDEFINED

    !SUBROUTINE GiD_BeginVectorResult -> REDEFINED
    !SUBROUTINE GiD_fBeginVectorResult -> REDEFINED

    !SUBROUTINE GiD_Begin2DMatResult -> REDEFINED
    !SUBROUTINE GiD_fBegin2DMatResult -> REDEFINED

    !SUBROUTINE GiD_Begin3DMatResult -> REDEFINED
    !SUBROUTINE GiD_fBegin3DMatResult -> REDEFINED

    !SUBROUTINE GiD_BeginPDMMatResult -> REDEFINED
    !SUBROUTINE GiD_fBeginPDMMatResult -> REDEFINED

    !SUBROUTINE GiD_BeginMainMatResult -> REDEFINED
    !SUBROUTINE GiD_fBeginMainMatResult -> REDEFINED

    !SUBROUTINE GiD_BeginLAResult -> REDEFINED
    !SUBROUTINE GiD_fBeginLAResult -> REDEFINED
    
    !SUBROUTINE GiD_BeginComplexScalarResult -> REDEFINED
    !SUBROUTINE GiD_fBeginComplexScalarResult -> REDEFINED
    
    !SUBROUTINE GiD_BeginComplexVectorResult -> REDEFINED
    !SUBROUTINE GiD_fBeginComplexVectorResult -> REDEFINED

    !SUBROUTINE GiD_BeginResultHeader -> REDEFINED
    !SUBROUTINE GiD_fBeginResultHeader -> REDEFINED

    !SUBROUTINE GiD_ScalarComp -> REDEFINED
    !SUBROUTINE GiD_fScalarComp -> REDEFINED

    !SUBROUTINE GiD_VectorComp -> REDEFINED
    !SUBROUTINE GiD_fVectorComp -> REDEFINED

    !SUBROUTINE GiD_2DMatrixComp -> REDEFINED
    !SUBROUTINE GiD_f2DMatrixComp -> REDEFINED

    !SUBROUTINE GiD_3DMatrixComp -> REDEFINED
    !SUBROUTINE GiD_f3DMatrixComp -> REDEFINED

    !SUBROUTINE GiD_PDMComp -> REDEFINED
    !SUBROUTINE GiD_fPDMComp -> REDEFINED

    !SUBROUTINE GiD_MainMatrixComp -> REDEFINED
    !SUBROUTINE GiD_fMainMatrixComp -> REDEFINED

    !SUBROUTINE GiD_LAComponents -> REDEFINED
    !SUBROUTINE GiD_fLAComponents -> REDEFINED
    
    !SUBROUTINE GiD_ComplexScalarComp -> REDEFINED
    !SUBROUTINE GiD_fComplexScalarComp -> REDEFINED
    
    !SUBROUTINE GiD_ResultUnit -> REDEFINED
    !SUBROUTINE GiD_fResultUnit -> REDEFINED

    !SUBROUTINE GiD_BeginResultGroup -> REDEFINED
    !SUBROUTINE GiD_fBeginResultGroup -> REDEFINED

    !SUBROUTINE GiD_ResultDescriptor -> REDEFINED
    !SUBROUTINE GiD_fResultDescriptor -> REDEFINED

    SUBROUTINE GiD_ResultValues() BIND(C,NAME='GiD_ResultValues')
    END SUBROUTINE GiD_ResultValues
    SUBROUTINE GiD_fResultValues(fd) BIND(C,NAME='GiD_fResultValues')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fResultValues

    SUBROUTINE GiD_EndResult() BIND(C,NAME='GiD_EndResult')
    END SUBROUTINE GiD_EndResult
    SUBROUTINE GiD_fEndResult(fd) BIND(C,NAME='GiD_fEndResult')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndResult

    !SUBROUTINE GiD_BeginOnMeshGroup -> REDEFINED
    !SUBROUTINE GiD_fBeginOnMeshGroup -> REDEFINED

    SUBROUTINE GiD_EndOnMeshGroup() BIND(C,NAME='GiD_EndOnMeshGroup')
    END SUBROUTINE GiD_EndOnMeshGroup
    SUBROUTINE GiD_fEndOnMeshGroup(fd) BIND(C,NAME='GiD_fEndOnMeshGroup')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fEndOnMeshGroup

    SUBROUTINE GiD_FlushPostFile() BIND(C,NAME='GiD_FlushPostFile')
    END SUBROUTINE GiD_FlushPostFile
    SUBROUTINE GiD_fFlushPostFile(fd) BIND(C,NAME='GiD_fFlushPostFile')
      IMPORT :: GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    END SUBROUTINE GiD_fFlushPostFile

    SUBROUTINE GiD_WriteScalar(id,v) BIND(C,NAME='GiD_WriteScalar')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: v
    END SUBROUTINE GiD_WriteScalar
    SUBROUTINE GiD_fWriteScalar(fd,id,v) BIND(C,NAME='GiD_fWriteScalar')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: v
    END SUBROUTINE GiD_fWriteScalar

    SUBROUTINE GiD_Write2DVector(id,x,y) BIND(C,NAME='GiD_Write2DVector')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x, y
    END SUBROUTINE GiD_Write2DVector
    SUBROUTINE GiD_fWrite2DVector(fd,id,x,y) BIND(C,NAME='GiD_fWrite2DVector')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x, y
    END SUBROUTINE GiD_fWrite2DVector

    ! SUBROUTINE GiD_WriteVector -> SPECIAL INTERFACE
    ! SUBROUTINE GiD_fWriteVector -> SPECIAL INTERFACE

    ! SUBROUTINE GiD_WriteVectorModule -> SPECIAL INTERFACE
    ! SUBROUTINE GiD_fWriteVectorModule -> SPECIAL INTERFACE

    SUBROUTINE GiD_Write2DMatrix(id,Sxx,Syy,Sxy) BIND(C,NAME='GiD_Write2DMatrix')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Sxx,Syy,Sxy
    END SUBROUTINE GiD_Write2DMatrix
    SUBROUTINE GiD_fWrite2DMatrix(fd,id,Sxx,Syy,Sxy) BIND(C,NAME='GiD_fWrite2DMatrix')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Sxx,Syy,Sxy
    END SUBROUTINE GiD_fWrite2DMatrix

    SUBROUTINE GiD_Write3DMatrix(id,Sxx,Syy,Szz,Sxy,Syz,Sxz) BIND(C,NAME='GiD_Write3DMatrix')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Sxx,Syy,Szz,Sxy,Syz,Sxz
    END SUBROUTINE GiD_Write3DMatrix
    SUBROUTINE GiD_fWrite3DMatrix(fd,id,Sxx,Syy,Szz,Sxy,Syz,Sxz) BIND(C,NAME='GiD_fWrite3DMatrix')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Sxx,Syy,Szz,Sxy,Syz,Sxz
    END SUBROUTINE GiD_fWrite3DMatrix

    SUBROUTINE GiD_WritePlainDefMatrix(id,Sxx,Syy,Sxy,Szz) BIND(C,NAME='GiD_WritePlainDefMatrix')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Sxx,Syy,Sxy,Szz
    END SUBROUTINE GiD_WritePlainDefMatrix
    SUBROUTINE GiD_fWritePlainDefMatrix(fd,id,Sxx,Syy,Sxy,Szz) BIND(C,NAME='GiD_fWritePlainDefMatrix')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Sxx,Syy,Sxy,Szz
    END SUBROUTINE GiD_fWritePlainDefMatrix

    SUBROUTINE GiD_WriteMainMatrix(id,Si,Sii,Siii,Vix,Viy,Viz,Viix,Viiy,Viiz,Viiix,Viiiy,Viiiz) BIND(C,NAME='GiD_WriteMainMatrix')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Si,Sii,Siii,Vix,Viy,Viz,Viix,Viiy,Viiz,Viiix,Viiiy,Viiiz
    END SUBROUTINE GiD_WriteMainMatrix
    SUBROUTINE GiD_fWriteMainMatrix(fd,id,Si,Sii,Siii,Vix,Viy,Viz,Viix,Viiy,Viiz,Viiix,Viiiy,Viiiz) &
      BIND(C,NAME='GiD_fWriteMainMatrix')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Si,Sii,Siii,Vix,Viy,Viz,Viix,Viiy,Viiz,Viiix,Viiiy,Viiiz
    END SUBROUTINE GiD_fWriteMainMatrix

    SUBROUTINE GiD_WriteLocalAxes(id,euler_1,euler_2,euler_3) BIND(C,NAME='GiD_WriteLocalAxes')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: euler_1,euler_2,euler_3
    END SUBROUTINE GiD_WriteLocalAxes
    SUBROUTINE GiD_fWriteLocalAxes(fd,id,euler_1,euler_2,euler_3) BIND(C,NAME='GiD_fWriteLocalAxes')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: euler_1,euler_2,euler_3
    END SUBROUTINE GiD_fWriteLocalAxes
    
    SUBROUTINE GiD_WriteComplexScalar(id,Re,Im) BIND(C,NAME='GiD_WriteComplexScalar')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Re,Im
    END SUBROUTINE GiD_WriteComplexScalar
    SUBROUTINE GiD_fWriteComplexScalar(fd,id,Re,Im) BIND(C,NAME='GiD_fWriteComplexScalar')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Re, Im
    END SUBROUTINE GiD_fWriteComplexScalar
    
    SUBROUTINE GiD_WriteComplexVector(id,Rex,Imx,Rey,Imy,Rez,Imz) BIND(C,NAME='GiD_WriteComplexVector')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Rex,Imx,Rey,Imy,Rez,Imz
    END SUBROUTINE GiD_WriteComplexVector
    SUBROUTINE GiD_fWriteComplexVector(fd,id,Rex,Imx,Rey,Imy,Rez,Imz) BIND(C,NAME='GiD_fWriteComplexVector')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Rex,Imx,Rey,Imy,Rez,Imz
    END SUBROUTINE GiD_fWriteComplexVector

  END INTERFACE

!--------------------------------------------------------------------------
! SPECIAL INTERFACES
!--------------------------------------------------------------------------
  
  INTERFACE GiD_WriteCoordinates
    MODULE PROCEDURE GiD_WriteCoordinates_f
    SUBROUTINE GiD_WriteCoordinates_cpp(id,x,y,z) BIND(C,NAME='GiD_WriteCoordinates')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z
    END SUBROUTINE GiD_WriteCoordinates_cpp
  END INTERFACE
  INTERFACE GiD_fWriteCoordinates
    MODULE PROCEDURE GiD_fWriteCoordinates_f
    SUBROUTINE GiD_fWriteCoordinates_cpp(fd,id,x,y,z) BIND(C,NAME='GiD_fWriteCoordinates')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z
    END SUBROUTINE GiD_fWriteCoordinates_cpp
  END INTERFACE

  INTERFACE GiD_WriteElementMat
    MODULE PROCEDURE GiD_WriteElementMat_f
    SUBROUTINE GiD_WriteElementMat_cpp(id,conec) BIND(C,NAME='GiD_WriteElementMat')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT),        INTENT(IN) :: conec(*)
    END SUBROUTINE GiD_WriteElementMat_cpp
  END INTERFACE
  INTERFACE GiD_fWriteElementMat
    MODULE PROCEDURE GiD_fWriteElementMat_f
    SUBROUTINE GiD_fWriteElementMat_cpp(fd,id,conec) BIND(C,NAME='GiD_fWriteElementMat')
      IMPORT :: C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      INTEGER(C_INT),        INTENT(IN) :: conec(*)
    END SUBROUTINE GiD_fWriteElementMat_cpp
  END INTERFACE

  INTERFACE GiD_WriteVector
    MODULE PROCEDURE GiD_WriteVector_f
    SUBROUTINE GiD_WriteVector_cpp(id,x,y,z) BIND(C,NAME='GiD_WriteVector')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z
    END SUBROUTINE GiD_WriteVector_cpp
  END INTERFACE
  INTERFACE GiD_fWriteVector
    MODULE PROCEDURE GiD_fWriteVector_f
    SUBROUTINE GiD_fWriteVector_cpp(fd,id,x,y,z) BIND(C,NAME='GiD_fWriteVector')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z
    END SUBROUTINE GiD_fWriteVector_cpp
  END INTERFACE

  INTERFACE GiD_WriteVectorModule
    MODULE PROCEDURE GiD_WriteVectorModule_f
    SUBROUTINE GiD_WriteVectorModule_cpp(id,x,y,z,mod) BIND(C,NAME='GiD_WriteVectorModule')
      IMPORT :: C_DOUBLE, C_INT
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z,mod
    END SUBROUTINE GiD_WriteVectorModule_cpp
  END INTERFACE
  INTERFACE GiD_fWriteVectorModule
    MODULE PROCEDURE GiD_fWriteVectorModule_f
    SUBROUTINE GiD_fWriteVectorModule_cpp(fd,id,x,y,z,mod) BIND(C,NAME='GiD_fWriteVectorModule')
      IMPORT :: C_DOUBLE, C_INT, GiD_File
      TYPE(GiD_File), VALUE, INTENT(IN) :: fd
      INTEGER(C_INT), VALUE, INTENT(IN) :: id
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x,y,z,mod
    END SUBROUTINE GiD_fWriteVectorModule_cpp
  END INTERFACE

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE GiD_PostMode_Equal
  END INTERFACE

  !INTERFACE ASSIGNMENT(=)
  !  MODULE PROCEDURE Copy_GiD_File
  !END INTERFACE
  !PUBLIC :: ASSIGNMENT(=)

! AUXILIAR INTERFACES
  INTERFACE GiD_Int
    MODULE PROCEDURE GiD_Int_1
    MODULE PROCEDURE GiD_Int_2
    MODULE PROCEDURE GiD_Int_4
    MODULE PROCEDURE GiD_Int_8
  END INTERFACE
  
  INTERFACE GiD_Real
    MODULE PROCEDURE GiD_Real_4
    MODULE PROCEDURE GiD_Real_8
  END INTERFACE


!--------------------------------------------------------------------------
! PUBLIC ROUTINES
!--------------------------------------------------------------------------

  public :: GiD_File_is_null
  PUBLIC :: GiD_PostInit, GiD_PostDone
  PUBLIC :: GiD_OpenPostMeshFile, GiD_ClosePostMeshFile
  PUBLIC :: GiD_fOpenPostMeshFile, GiD_fOpenPostMeshFile_utf8, GiD_fClosePostMeshFile
  PUBLIC :: GiD_BeginMeshGroup, GiD_EndMeshGroup
  PUBLIC :: GiD_fBeginMeshGroup, GiD_fEndMeshGroup
  PUBLIC :: GiD_BeginMesh, GiD_BeginMeshColor, GiD_EndMesh
  PUBLIC :: GiD_fBeginMesh, GiD_fBeginMeshColor, GiD_fEndMesh
  PUBLIC :: GiD_MeshUnit, GiD_ResultUnit
  PUBLIC :: GiD_fMeshUnit, GiD_fResultUnit
  PUBLIC :: GiD_BeginCoordinates, GiD_EndCoordinates, GiD_WriteCoordinates, GiD_WriteCoordinates2D, GiD_WriteCoordinates3D
  PUBLIC :: GiD_fBeginCoordinates, GiD_fEndCoordinates, GiD_fWriteCoordinates, GiD_fWriteCoordinates2D, GiD_fWriteCoordinates3D
  PUBLIC :: GiD_BeginElements, GiD_EndElements, GiD_WriteElement, GiD_WriteElementMat
  PUBLIC :: GiD_fBeginElements, GiD_fEndElements, GiD_fWriteElement, GiD_fWriteElementMat
  PUBLIC :: GiD_WriteCircle, GiD_WriteCircleMat, GiD_WriteSphere, GiD_WriteSphereMat
  PUBLIC :: GiD_fWriteCircle, GiD_fWriteCircleMat, GiD_fWriteSphere, GiD_fWriteSphereMat
  PUBLIC :: GiD_OpenPostResultFile, GiD_ClosePostResultFile
  PUBLIC :: GiD_fOpenPostResultFile, GiD_fOpenPostResultFile_utf8, GiD_fClosePostResultFile
  PUBLIC :: GiD_BeginGaussPoint, GiD_EndGaussPoint, GiD_WriteGaussPoint2D, GiD_WriteGaussPoint3D
  PUBLIC :: GiD_fBeginGaussPoint, GiD_fEndGaussPoint, GiD_fWriteGaussPoint2D, GiD_fWriteGaussPoint3D
  PUBLIC :: GiD_BeginRangeTable, GiD_EndRangeTable, GiD_WriteMinRange, GiD_WriteRange, GiD_WriteMaxRange
  PUBLIC :: GiD_fBeginRangeTable, GiD_fEndRangeTable, GiD_fWriteMinRange, GiD_fWriteRange, GiD_fWriteMaxRange
  PUBLIC :: GiD_BeginResult, GiD_BeginScalarResult, GiD_BeginVectorResult, GiD_Begin2DMatResult, GiD_Begin3DMatResult
  PUBLIC :: GiD_fBeginResult, GiD_fBeginScalarResult, GiD_fBeginVectorResult, GiD_fBegin2DMatResult, GiD_fBegin3DMatResult
  PUBLIC :: GiD_BeginPDMMatResult, GiD_BeginMainMatResult, GiD_BeginLAResult, GiD_BeginResultHeader
  PUBLIC :: GiD_BeginComplexScalarResult, GiD_BeginComplexVectorResult
  PUBLIC :: GiD_fBeginPDMMatResult, GiD_fBeginMainMatResult, GiD_fBeginLAResult, GiD_fBeginResultHeader
  PUBLIC :: GiD_fBeginComplexScalarResult, GiD_fBeginComplexVectorResult
  PUBLIC :: GiD_ResultComponents, GiD_ScalarComp, GiD_VectorComp, GiD_2DMatrixComp, GiD_3DMatrixComp, GiD_PDMComp
  PUBLIC :: GiD_fResultComponents, GiD_fScalarComp, GiD_fVectorComp, GiD_f2DMatrixComp, GiD_f3DMatrixComp, GiD_fPDMComp
  PUBLIC :: GiD_MainMatrixComp, GiD_LAComponents, GiD_ComplexScalarComp
  PUBLIC :: GiD_fMainMatrixComp, GiD_fLAComponents, GiD_fComplexScalarComp
  PUBLIC :: GiD_BeginResultGroup, GiD_ResultDescription, GiD_ResultValues, GiD_EndResult
  PUBLIC :: GiD_fBeginResultGroup, GiD_fResultDescription, GiD_fResultValues, GiD_fEndResult
  PUBLIC :: GiD_BeginOnMeshGroup, GiD_EndOnMeshGroup
  PUBLIC :: GiD_fBeginOnMeshGroup, GiD_fEndOnMeshGroup
  PUBLIC :: GiD_FlushPostFile
  PUBLIC :: GiD_fFlushPostFile
  PUBLIC :: GiD_WriteScalar, GiD_WriteVector, GiD_Write2DVector, GiD_Write3DVector, GiD_WriteVectorModule, GiD_Write2DMatrix
  PUBLIC :: GiD_fWriteScalar, GiD_fWriteVector, GiD_fWrite2DVector, GiD_fWrite3DVector, GiD_fWriteVectorModule
  PUBLIC :: GiD_fWrite2DMatrix
  PUBLIC :: GiD_Write3DMatrix, GiD_WritePlainDefMatrix, GiD_WriteMainMatrix, GiD_WriteLocalAxes
  PUBLIC :: GiD_WriteComplexScalar, GiD_WriteComplexVector
  PUBLIC :: GiD_fWrite3DMatrix, GiD_fWritePlainDefMatrix, GiD_fWriteMainMatrix, GiD_fWriteLocalAxes
  PUBLIC :: GiD_fWriteComplexScalar, GiD_fWriteComplexVector
  PUBLIC :: OPERATOR(==)
  PUBLIC :: GiD_Int, GiD_Real

!--------------------------------------------------------------------------

CONTAINS

  
    FUNCTION GiD_File_is_null( fd ) RESULT(res)
  
      TYPE(GiD_File), INTENT(IN) :: fd
      LOGICAL :: res
      res = .FALSE.
      if ( fd%file_unit == GID_NO_FILE) then
          res = .TRUE.
      end if
      
    END FUNCTION GiD_File_is_null
    
!--------------------------------------------------------------------------
! REDEFINED SUBROUTINES WITH STRING ARGUMENTS
!--------------------------------------------------------------------------

  SUBROUTINE GiD_OpenPostMeshFile( filename, postmode )

    CHARACTER(LEN=*),   INTENT(IN) :: filename
    TYPE(GiD_PostMode), INTENT(IN) :: postmode

    CHARACTER(1,KIND=C_CHAR) :: filename_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_OpenPostMeshFile_cpp( filename_cpp, postmode ) BIND(C,NAME='GiD_OpenPostMeshFile')
        IMPORT:: GiD_PostMode, C_CHAR
        CHARACTER(1,KIND=C_CHAR),  INTENT(IN) :: filename_cpp(1)
        TYPE(GiD_PostMode), VALUE, INTENT(IN) :: postmode
      END SUBROUTINE GiD_OpenPostMeshFile_cpp
    END INTERFACE

    filename_cpp = C_STRING(filename)
    CALL GiD_OpenPostMeshFile_cpp( filename_cpp, postmode )

  END SUBROUTINE GiD_OpenPostMeshFile

  
  
  FUNCTION GiD_fOpenPostMeshFile( filename, postmode ) RESULT(fd)

    CHARACTER(LEN=*),   INTENT(IN) :: filename
    TYPE(GiD_PostMode), INTENT(IN) :: postmode
    TYPE(GiD_File) :: fd

    CHARACTER(1,KIND=C_CHAR) :: filename_cpp(MAX_STRING_LEN)

    INTERFACE
      INTEGER(C_INT) FUNCTION GiD_fOpenPostMeshFile_cpp( filename_cpp, postmode ) BIND(C,NAME='GiD_fOpenPostMeshFile')
        IMPORT:: GiD_PostMode, C_CHAR, GiD_File, C_INT
        CHARACTER(1,KIND=C_CHAR),  INTENT(IN) :: filename_cpp(1)
        TYPE(GiD_PostMode), VALUE, INTENT(IN) :: postmode
      END FUNCTION GiD_fOpenPostMeshFile_cpp
    END INTERFACE

    filename_cpp = C_STRING(filename)
    fd%file_unit = GiD_fOpenPostMeshFile_cpp( filename_cpp, postmode )

  END FUNCTION GiD_fOpenPostMeshFile

  FUNCTION GiD_fOpenPostMeshFile_utf8( filename, postmode ) RESULT(fd)

    CHARACTER(LEN=*),   INTENT(IN) :: filename
    TYPE(GiD_PostMode), INTENT(IN) :: postmode
    TYPE(GiD_File) :: fd

    CHARACTER(1,KIND=C_CHAR) :: filename_cpp(MAX_STRING_LEN)

    INTERFACE
      INTEGER(C_INT) FUNCTION GiD_fOpenPostMeshFile_utf8_cpp( filename_cpp, postmode ) BIND(C,NAME='GiD_fOpenPostMeshFile_utf8')
        IMPORT:: GiD_PostMode, C_CHAR, GiD_File, C_INT
        CHARACTER(1,KIND=C_CHAR),  INTENT(IN) :: filename_cpp(1)
        TYPE(GiD_PostMode), VALUE, INTENT(IN) :: postmode
      END FUNCTION GiD_fOpenPostMeshFile_utf8_cpp
    END INTERFACE

    filename_cpp = C_STRING(filename)
    fd%file_unit = GiD_fOpenPostMeshFile_utf8_cpp( filename_cpp, postmode )

  END FUNCTION GiD_fOpenPostMeshFile_utf8

!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_BeginMeshGroup( label )

    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_BeginMeshGroup_cpp( label_cpp ) BIND(C,NAME='GiD_BeginMeshGroup')
        IMPORT :: C_CHAR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_BeginMeshGroup_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_BeginMeshGroup_cpp( label_cpp )

  END SUBROUTINE GiD_BeginMeshGroup
  
  SUBROUTINE GiD_fBeginMeshGroup( fd, label )

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fBeginMeshGroup_cpp( fd, label_cpp ) BIND(C,NAME='GiD_fBeginMeshGroup')
        IMPORT :: C_CHAR, GiD_File
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_fBeginMeshGroup_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_fBeginMeshGroup_cpp( fd, label_cpp )

  END SUBROUTINE GiD_fBeginMeshGroup
  
!--------------------------------------------------------------------------
   
  SUBROUTINE GiD_MeshUnit( unit_name )

    CHARACTER(LEN=*), INTENT(IN) :: unit_name
    
    CHARACTER(1,KIND=C_CHAR) :: unit_name_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_MeshUnit_cpp( unit_name_cpp ) BIND(C,NAME='GiD_MeshUnit')      
        IMPORT :: C_CHAR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: unit_name_cpp(1)
      END SUBROUTINE GiD_MeshUnit_cpp
    END INTERFACE

    unit_name_cpp = C_STRING(unit_name)
    call GiD_MeshUnit_cpp(unit_name_cpp)

  END SUBROUTINE GiD_MeshUnit

  SUBROUTINE GiD_fMeshUnit( fd, unit_name )

    TYPE(GiD_File), VALUE, INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: unit_name     
    
    CHARACTER(1,KIND=C_CHAR) :: unit_name_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fMeshUnit_cpp( fd, unit_name_cpp ) BIND(C,NAME='GiD_fMeshUnit')      
        IMPORT :: GiD_File, C_CHAR
        TYPE(GiD_File), VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: unit_name_cpp(1)
      END SUBROUTINE GiD_fMeshUnit_cpp
    END INTERFACE

    unit_name_cpp = C_STRING(unit_name)
    call GiD_fMeshUnit_cpp(fd,unit_name_cpp)

  END SUBROUTINE GiD_fMeshUnit

!--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginMesh(label,dim,elemtype,nnode)

    CHARACTER(LEN=*),      INTENT(IN) :: label
    TYPE(GiD_Dimension),   INTENT(IN) :: dim
    TYPE(GiD_ElementType), INTENT(IN) :: elemtype
    INTEGER(C_INT),        INTENT(IN) :: nnode

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_BeginMesh_cpp(label_cpp,dim,elemtype,nnode) BIND(C,NAME='GiD_BeginMesh')
        IMPORT:: GiD_Dimension, GiD_ElementType, C_INT, C_CHAR
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: label_cpp(1)
        TYPE(GiD_Dimension),   VALUE, INTENT(IN) :: dim
        TYPE(GiD_ElementType), VALUE, INTENT(IN) :: elemtype
        INTEGER(C_INT),        VALUE, INTENT(IN) :: nnode
      END SUBROUTINE GiD_BeginMesh_cpp
    END INTERFACE
  
    label_cpp = C_STRING(label)
    CALL GiD_BeginMesh_cpp(label_cpp,dim,elemtype,nnode)

  END SUBROUTINE GiD_BeginMesh
    
  SUBROUTINE GiD_fBeginMesh(fd,label,dim,elemtype,nnode)

    TYPE(GiD_File),        INTENT(IN) :: fd
    CHARACTER(LEN=*),      INTENT(IN) :: label
    TYPE(GiD_Dimension),   INTENT(IN) :: dim
    TYPE(GiD_ElementType), INTENT(IN) :: elemtype
    INTEGER(C_INT),        INTENT(IN) :: nnode

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fBeginMesh_cpp(fd,label_cpp,dim,elemtype,nnode) BIND(C,NAME='GiD_fBeginMesh')
        IMPORT:: GiD_Dimension, GiD_ElementType, C_INT, C_CHAR, GiD_File
        TYPE(GiD_File),        VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: label_cpp(1)
        TYPE(GiD_Dimension),   VALUE, INTENT(IN) :: dim
        TYPE(GiD_ElementType), VALUE, INTENT(IN) :: elemtype
        INTEGER(C_INT),        VALUE, INTENT(IN) :: nnode
      END SUBROUTINE GiD_fBeginMesh_cpp
    END INTERFACE
  
    label_cpp = C_STRING(label)
    CALL GiD_fBeginMesh_cpp(fd,label_cpp,dim,elemtype,nnode)

  END SUBROUTINE GiD_fBeginMesh
  
!--------------------------------------------------------------------------
    
  SUBROUTINE GiD_BeginMeshColor(label,dim,elemtype,nnode,red,green,blue)

    CHARACTER(LEN=*),      INTENT(IN) :: label
    TYPE(GiD_Dimension),   INTENT(IN) :: dim
    TYPE(GiD_ElementType), INTENT(IN) :: elemtype
    INTEGER(C_INT),        INTENT(IN) :: nnode
    REAL(C_DOUBLE),        INTENT(IN) :: red, green, blue

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_BeginMeshColor_cpp(label_cpp,dim,elemtype,nnode,red,green,blue) BIND(C,NAME='GiD_BeginMeshColor')
        IMPORT:: C_DOUBLE, GiD_Dimension, GiD_ElementType, C_INT, C_CHAR
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: label_cpp(1)
        TYPE(GiD_Dimension),   VALUE, INTENT(IN) :: dim
        TYPE(GiD_ElementType), VALUE, INTENT(IN) :: elemtype
        INTEGER(C_INT),        VALUE, INTENT(IN) :: nnode
        REAL(C_DOUBLE),        VALUE, INTENT(IN) :: red, green, blue
      END SUBROUTINE GiD_BeginMeshColor_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_BeginMeshColor_cpp(label_cpp,dim,elemtype,nnode,red,green,blue)

  END SUBROUTINE GiD_BeginMeshColor
    
  SUBROUTINE GiD_fBeginMeshColor(fd,label,dim,elemtype,nnode,red,green,blue)

    TYPE(GiD_File),        INTENT(IN) :: fd
    CHARACTER(LEN=*),      INTENT(IN) :: label
    TYPE(GiD_Dimension),   INTENT(IN) :: dim
    TYPE(GiD_ElementType), INTENT(IN) :: elemtype
    INTEGER(C_INT),        INTENT(IN) :: nnode
    REAL(C_DOUBLE),        INTENT(IN) :: red, green, blue

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fBeginMeshColor_cpp(fd,label_cpp,dim,elemtype,nnode,red,green,blue) BIND(C,NAME='GiD_fBeginMeshColor')
        IMPORT:: C_DOUBLE, GiD_Dimension, GiD_ElementType, C_INT, C_CHAR, GiD_File
        TYPE(GiD_File),        VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: label_cpp(1)
        TYPE(GiD_Dimension),   VALUE, INTENT(IN) :: dim
        TYPE(GiD_ElementType), VALUE, INTENT(IN) :: elemtype
        INTEGER(C_INT),        VALUE, INTENT(IN) :: nnode
        REAL(C_DOUBLE),        VALUE, INTENT(IN) :: red, green, blue
      END SUBROUTINE GiD_fBeginMeshColor_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_fBeginMeshColor_cpp(fd,label_cpp,dim,elemtype,nnode,red,green,blue)

  END SUBROUTINE GiD_fBeginMeshColor
  
!--------------------------------------------------------------------------
    
  SUBROUTINE GiD_OpenPostResultFile(filename,postmode)

    CHARACTER(LEN=*),   INTENT(IN) :: filename
    TYPE(GiD_PostMode), INTENT(IN) :: postmode

    CHARACTER(1,KIND=C_CHAR) :: filename_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_OpenPostResultFile_cpp(filename_cpp,postmode) BIND(C,NAME='GiD_OpenPostResultFile')
        IMPORT:: GiD_PostMode, C_CHAR
        CHARACTER(1,KIND=C_CHAR),  INTENT(IN) :: filename_cpp(1)
        TYPE(GiD_PostMode), VALUE, INTENT(IN) :: postmode
      END SUBROUTINE GiD_OpenPostResultFile_cpp
    END INTERFACE

    filename_cpp = C_STRING(filename)
    CALL GiD_OpenPostResultFile_cpp(filename_cpp,postmode)

  END SUBROUTINE GiD_OpenPostResultFile

  FUNCTION GiD_fOpenPostResultFile(filename,postmode) RESULT(fd)

    CHARACTER(LEN=*),   INTENT(IN) :: filename
    TYPE(GiD_PostMode), INTENT(IN) :: postmode
    TYPE(GiD_File) :: fd

    CHARACTER(1,KIND=C_CHAR) :: filename_cpp(MAX_STRING_LEN)

    INTERFACE
      INTEGER(C_INT) FUNCTION GiD_fOpenPostResultFile_cpp(filename_cpp,postmode) BIND(C,NAME='GiD_fOpenPostResultFile')
        IMPORT:: GiD_PostMode, C_CHAR, GiD_File, C_INT
        CHARACTER(1,KIND=C_CHAR),  INTENT(IN) :: filename_cpp(1)
        TYPE(GiD_PostMode), VALUE, INTENT(IN) :: postmode
      END FUNCTION GiD_fOpenPostResultFile_cpp
    END INTERFACE

    filename_cpp = C_STRING(filename)
    fd%file_unit = GiD_fOpenPostResultFile_cpp(filename_cpp,postmode)

  END FUNCTION GiD_fOpenPostResultFile
  
  FUNCTION GiD_fOpenPostResultFile_utf8(filename,postmode) RESULT(fd)

    CHARACTER(LEN=*),   INTENT(IN) :: filename
    TYPE(GiD_PostMode), INTENT(IN) :: postmode
    TYPE(GiD_File) :: fd

    CHARACTER(1,KIND=C_CHAR) :: filename_cpp(MAX_STRING_LEN)

    INTERFACE
      INTEGER(C_INT) FUNCTION GiD_fOpenPostResultFile_utf8_cpp(filename_cpp,postmode) BIND(C,NAME='GiD_fOpenPostResultFile_utf8')
        IMPORT:: GiD_PostMode, C_CHAR, GiD_File, C_INT
        CHARACTER(1,KIND=C_CHAR),  INTENT(IN) :: filename_cpp(1)
        TYPE(GiD_PostMode), VALUE, INTENT(IN) :: postmode
      END FUNCTION GiD_fOpenPostResultFile_utf8_cpp
    END INTERFACE

    filename_cpp = C_STRING(filename)
    fd%file_unit = GiD_fOpenPostResultFile_utf8_cpp(filename_cpp,postmode)

  END FUNCTION GiD_fOpenPostResultFile_utf8
  
!--------------------------------------------------------------------------
    
  SUBROUTINE GiD_BeginGaussPoint(label,elemtype,meshname,ngauss,nodeincluded,internalcoord)
    
    CHARACTER(LEN=*),      INTENT(IN) :: label
    TYPE(GiD_ElementType), INTENT(IN) :: elemtype
    CHARACTER(LEN=*),      INTENT(IN) :: meshname
    INTEGER(C_INT),        INTENT(IN) :: ngauss
    INTEGER(C_INT),        INTENT(IN) :: nodeincluded
    INTEGER(C_INT),        INTENT(IN) :: internalcoord

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN), meshname_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_BeginGaussPoint_cpp(label_cpp,elemtype,meshname_cpp,ngauss,nodeincluded,internalcoord) &
        BIND(C,NAME='GiD_BeginGaussPoint')
        IMPORT:: GiD_ElementType, C_INT, C_CHAR
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: label_cpp(1)
        TYPE(GiD_ElementType), VALUE, INTENT(IN) :: elemtype
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: meshname_cpp(1)
        INTEGER(C_INT),        VALUE, INTENT(IN) :: ngauss
        INTEGER(C_INT),        VALUE, INTENT(IN) :: nodeincluded
        INTEGER(C_INT),        VALUE, INTENT(IN) :: internalcoord
      END SUBROUTINE GiD_BeginGaussPoint_cpp
    END INTERFACE

    label_cpp    = C_STRING(label)
    meshname_cpp = C_STRING(meshname)
    CALL GiD_BeginGaussPoint_cpp(label_cpp,elemtype,meshname_cpp,ngauss,nodeincluded,internalcoord)

  END SUBROUTINE GiD_BeginGaussPoint
    
  SUBROUTINE GiD_fBeginGaussPoint(fd,label,elemtype,meshname,ngauss,nodeincluded,internalcoord)
    
    TYPE(GiD_File),        INTENT(IN) :: fd
    CHARACTER(LEN=*),      INTENT(IN) :: label
    TYPE(GiD_ElementType), INTENT(IN) :: elemtype
    CHARACTER(LEN=*),      INTENT(IN) :: meshname
    INTEGER(C_INT),        INTENT(IN) :: ngauss
    INTEGER(C_INT),        INTENT(IN) :: nodeincluded
    INTEGER(C_INT),        INTENT(IN) :: internalcoord

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN), meshname_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fBeginGaussPoint_cpp(fd,label_cpp,elemtype,meshname_cpp,ngauss,nodeincluded,internalcoord) &
        BIND(C,NAME='GiD_fBeginGaussPoint')
        IMPORT:: GiD_ElementType, C_INT, C_CHAR, GiD_File
        TYPE(GiD_File),        VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: label_cpp(1)
        TYPE(GiD_ElementType), VALUE, INTENT(IN) :: elemtype
        CHARACTER(1,KIND=C_CHAR),     INTENT(IN) :: meshname_cpp(1)
        INTEGER(C_INT),        VALUE, INTENT(IN) :: ngauss
        INTEGER(C_INT),        VALUE, INTENT(IN) :: nodeincluded
        INTEGER(C_INT),        VALUE, INTENT(IN) :: internalcoord
      END SUBROUTINE GiD_fBeginGaussPoint_cpp
    END INTERFACE

    label_cpp    = C_STRING(label)
    meshname_cpp = C_STRING(meshname)
    CALL GiD_fBeginGaussPoint_cpp(fd,label_cpp,elemtype,meshname_cpp,ngauss,nodeincluded,internalcoord)

  END SUBROUTINE GiD_fBeginGaussPoint
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_BeginRangeTable(label)
  
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_BeginRangeTable_cpp(label_cpp) BIND(C,NAME='GiD_BeginRangeTable')
        IMPORT :: C_CHAR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_BeginRangeTable_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_BeginRangeTable_cpp(label_cpp)

  END SUBROUTINE GiD_BeginRangeTable
  
  SUBROUTINE GiD_fBeginRangeTable(fd,label)
  
    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fBeginRangeTable_cpp(fd,label_cpp) BIND(C,NAME='GiD_fBeginRangeTable')
        IMPORT :: C_CHAR, GiD_File
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_fBeginRangeTable_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_fBeginRangeTable_cpp(fd,label_cpp)

  END SUBROUTINE GiD_fBeginRangeTable
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_WriteMinRange(minvalue,label)

    REAL(C_DOUBLE),   INTENT(IN) :: minvalue
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_WriteMinRange_cpp(minvalue,label_cpp) BIND(C,NAME='GiD_WriteMinRange')
        IMPORT :: C_DOUBLE, C_CHAR
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: minvalue
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_WriteMinRange_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_WriteMinRange_cpp(minvalue,label_cpp)

  END SUBROUTINE GiD_WriteMinRange

  SUBROUTINE GiD_fWriteMinRange(fd,minvalue,label)

    TYPE(GiD_File),   INTENT(IN) :: fd
    REAL(C_DOUBLE),   INTENT(IN) :: minvalue
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fWriteMinRange_cpp(fd,minvalue,label_cpp) BIND(C,NAME='GiD_fWriteMinRange')
        IMPORT :: C_DOUBLE, C_CHAR, GiD_File
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: minvalue
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_fWriteMinRange_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_fWriteMinRange_cpp(fd,minvalue,label_cpp)

  END SUBROUTINE GiD_fWriteMinRange

!--------------------------------------------------------------------------

  SUBROUTINE GiD_WriteRange(minvalue,maxvalue,label)

    REAL(C_DOUBLE),   INTENT(IN) :: minvalue, maxvalue
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_WriteRange_cpp(minvalue,maxvalue,label_cpp) BIND(C,NAME='GiD_WriteRange')
        IMPORT :: C_DOUBLE, C_CHAR
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: minvalue, maxvalue
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_WriteRange_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_WriteRange_cpp(minvalue,maxvalue,label_cpp)

  END SUBROUTINE GiD_WriteRange

  SUBROUTINE GiD_fWriteRange(fd,minvalue,maxvalue,label)

    TYPE(GiD_File),   INTENT(IN) :: fd
    REAL(C_DOUBLE),   INTENT(IN) :: minvalue, maxvalue
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fWriteRange_cpp(fd,minvalue,maxvalue,label_cpp) BIND(C,NAME='GiD_fWriteRange')
        IMPORT :: C_DOUBLE, C_CHAR, GiD_File
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: minvalue, maxvalue
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_fWriteRange_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_fWriteRange_cpp(fd,minvalue,maxvalue,label_cpp)

  END SUBROUTINE GiD_fWriteRange

!--------------------------------------------------------------------------

  SUBROUTINE GiD_WriteMaxRange(maxvalue,label)

    REAL(C_DOUBLE),   INTENT(IN) :: maxvalue
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_WriteMaxRange_cpp(maxvalue,label_cpp) BIND(C,NAME='GiD_WriteMaxRange')
        IMPORT :: C_DOUBLE, C_CHAR
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: maxvalue
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_WriteMaxRange_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_WriteMaxRange_cpp(maxvalue,label_cpp)

  END SUBROUTINE GiD_WriteMaxRange

  SUBROUTINE GiD_fWriteMaxRange(fd,maxvalue,label)

    TYPE(GiD_File),   INTENT(IN) :: fd
    REAL(C_DOUBLE),   INTENT(IN) :: maxvalue
    CHARACTER(LEN=*), INTENT(IN) :: label

    CHARACTER(1,KIND=C_CHAR) :: label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fWriteMaxRange_cpp(fd,maxvalue,label_cpp) BIND(C,NAME='GiD_fWriteMaxRange')
        IMPORT :: C_DOUBLE, C_CHAR, GiD_File
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: maxvalue
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: label_cpp(1)
      END SUBROUTINE GiD_fWriteMaxRange_cpp
    END INTERFACE

    label_cpp = C_STRING(label)
    CALL GiD_fWriteMaxRange_cpp(fd,maxvalue,label_cpp)

  END SUBROUTINE GiD_fWriteMaxRange

!--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginScalarResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp

    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr
    TYPE(C_PTR) :: RangeTable_ptr
    TYPE(C_PTR) :: Comp_ptr

    INTERFACE
      SUBROUTINE GiD_BeginScalarResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_cpp,RangeTable_cpp,Comp_cpp) &
        BIND(C,NAME='GiD_BeginScalarResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result_cpp(1)       ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis_cpp(1)     ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: Step                !
        TYPE(GiD_ResLoc),  VALUE, INTENT(IN) :: Where               ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: GaussPointsName_cpp ! gauss point name (just ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: RangeTable_cpp
        TYPE(C_PTR),       VALUE, INTENT(IN) :: Comp_cpp
      END SUBROUTINE GiD_BeginScalarResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr = C_PTR_STRING(Comp,Comp_cpp)
    Result_cpp = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginScalarResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr,Comp_ptr)

  END SUBROUTINE GiD_BeginScalarResult

  SUBROUTINE GiD_fBeginScalarResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(1)

    Comp_f(1) = Comp

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_Scalar,Where,GaussPointsName,RangeTable,1,Comp_f)

  END SUBROUTINE GiD_fBeginScalarResult

!--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginVectorResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4

    CHARACTER(1,KIND=C_CHAR), TARGET  :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: Comp_cpp(GID_MAX_CHAR_LEN,4)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr
    TYPE(C_PTR) :: RangeTable_ptr
    TYPE(C_PTR) :: Comp_ptr(4)

    INTERFACE
      SUBROUTINE GiD_BeginVectorResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4) &
        BIND(C,NAME='GiD_BeginVectorResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Comp1, Comp2, Comp3, Comp4
      END SUBROUTINE GiD_BeginVectorResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1) = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3) = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4) = C_PTR_STRING(Comp4,Comp_cpp(:,4))
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginVectorResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr, &
                                   Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4))

  END SUBROUTINE GiD_BeginVectorResult

  SUBROUTINE GiD_fBeginVectorResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(4)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3;  Comp_f(4) = Comp4

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_Vector,Where,GaussPointsName,RangeTable,4,Comp_f)

  END SUBROUTINE GiD_fBeginVectorResult
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_Begin2DMatResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3

    CHARACTER(1,KIND=C_CHAR), TARGET  :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: Comp_cpp(GID_MAX_CHAR_LEN,3)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR)  :: GaussPointsName_ptr
    TYPE(C_PTR)  :: RangeTable_ptr
    TYPE(C_PTR), TARGET :: Comp_ptr(3)

    INTERFACE
      SUBROUTINE GiD_Begin2DMatResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3) &
        BIND(C,NAME='GiD_Begin2DMatResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Comp1, Comp2, Comp3
      END SUBROUTINE GiD_Begin2DMatResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1) = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3) = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_Begin2DMatResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr, &
                                  Comp_ptr(1),Comp_ptr(2),Comp_ptr(3))

  END SUBROUTINE GiD_Begin2DMatResult
  
  SUBROUTINE GiD_fBegin2DMatResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(3)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_Matrix,Where,GaussPointsName,RangeTable,3,Comp_f)

  END SUBROUTINE GiD_fBegin2DMatResult
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_Begin3DMatResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4, Comp5, Comp6

    CHARACTER(1,KIND=C_CHAR), TARGET  :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: Comp_cpp(GID_MAX_CHAR_LEN,6)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR)  :: RangeTable_ptr
    TYPE(C_PTR)  :: GaussPointsName_ptr
    TYPE(C_PTR), TARGET :: Comp_ptr(6)

    INTERFACE
      SUBROUTINE GiD_Begin3DMatResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,&
                                          Comp4,Comp5,Comp6) BIND(C,NAME='GiD_Begin3DMatResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Comp1, Comp2, Comp3, Comp4, Comp5, Comp6
      END SUBROUTINE GiD_Begin3DMatResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1) = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3) = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4) = C_PTR_STRING(Comp4,Comp_cpp(:,4))
    Comp_ptr(5) = C_PTR_STRING(Comp5,Comp_cpp(:,5))
    Comp_ptr(6) = C_PTR_STRING(Comp6,Comp_cpp(:,6))
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_Begin3DMatResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr, &
                                  Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4),Comp_ptr(5),Comp_ptr(6))

  END SUBROUTINE GiD_Begin3DMatResult
  
  SUBROUTINE GiD_fBegin3DMatResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4, Comp5, Comp6

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(6)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3
    Comp_f(4) = Comp4;  Comp_f(5) = Comp5;  Comp_f(6) = Comp6

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_Matrix,Where,GaussPointsName,RangeTable,6,Comp_f)

  END SUBROUTINE GiD_fBegin3DMatResult
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginPDMMatResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4

    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,4)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr
    TYPE(C_PTR) :: RangeTable_ptr
    TYPE(C_PTR) :: Comp_ptr(4)

    INTERFACE
      SUBROUTINE GiD_BeginPDMMatResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4) &
        BIND(C,NAME='GiD_BeginPDMMatResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Comp1, Comp2, Comp3, Comp4
      END SUBROUTINE GiD_BeginPDMMatResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1) = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3) = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4) = C_PTR_STRING(Comp4,Comp_cpp(:,4))
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginPDMMatResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr, &
                                   Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4))

  END SUBROUTINE GiD_BeginPDMMatResult

  SUBROUTINE GiD_fBeginPDMMatResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(4)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3;  Comp_f(4) = Comp4

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_PlainDefMatrix,Where,GaussPointsName,RangeTable,4,Comp_f)

  END SUBROUTINE GiD_fBeginPDMMatResult
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_BeginMainMatResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6, &
                                    Comp7,Comp8,Comp9,Comp10,Comp11,Comp12)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12

    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,12)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr
    TYPE(C_PTR) :: RangeTable_ptr
    TYPE(C_PTR) :: Comp_ptr(12)

    INTERFACE
      SUBROUTINE GiD_BeginMainMatResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4,Comp5, &
                                            Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12) BIND(C,NAME='GiD_BeginMainMatResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Comp1, Comp2, Comp3, Comp4,  Comp5,  Comp6
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Comp7, Comp8, Comp9, Comp10, Comp11, Comp12
      END SUBROUTINE GiD_BeginMainMatResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1)  = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2)  = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3)  = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4)  = C_PTR_STRING(Comp4,Comp_cpp(:,4))
    Comp_ptr(5)  = C_PTR_STRING(Comp5,Comp_cpp(:,5))
    Comp_ptr(6)  = C_PTR_STRING(Comp6,Comp_cpp(:,6))
    Comp_ptr(7)  = C_PTR_STRING(Comp7,Comp_cpp(:,7))
    Comp_ptr(8)  = C_PTR_STRING(Comp8,Comp_cpp(:,8))
    Comp_ptr(9)  = C_PTR_STRING(Comp9,Comp_cpp(:,9))
    Comp_ptr(10) = C_PTR_STRING(Comp10,Comp_cpp(:,10))
    Comp_ptr(11) = C_PTR_STRING(Comp11,Comp_cpp(:,11))
    Comp_ptr(12) = C_PTR_STRING(Comp12,Comp_cpp(:,12))
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginMainMatResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr,  &
                                    Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4),Comp_ptr(5),Comp_ptr(6),   &
                                    Comp_ptr(7),Comp_ptr(8),Comp_ptr(9),Comp_ptr(10),Comp_ptr(11),Comp_ptr(12) )

  END SUBROUTINE GiD_BeginMainMatResult
  
  SUBROUTINE GiD_fBeginMainMatResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6, &
                                    Comp7,Comp8,Comp9,Comp10,Comp11,Comp12)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(12)

    Comp_f(1) = Comp1;  Comp_f(2)  = Comp2;  Comp_f(3)  = Comp3;  Comp_f(4)  = Comp4
    Comp_f(5) = Comp5;  Comp_f(6)  = Comp6;  Comp_f(7)  = Comp7;  Comp_f(8)  = Comp8
    Comp_f(9) = Comp9;  Comp_f(10) = Comp10; Comp_f(11) = Comp11; Comp_f(12) = Comp12

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_MainMatrix,Where,GaussPointsName,RangeTable,12,Comp_f)

  END SUBROUTINE GiD_fBeginMainMatResult

!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_BeginLAResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3

    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,3)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR)  :: GaussPointsName_ptr
    TYPE(C_PTR)  :: RangeTable_ptr
    TYPE(C_PTR), TARGET :: Comp_ptr(3)

    INTERFACE
      SUBROUTINE GiD_BeginLAResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3) &
        BIND(C,NAME='GiD_BeginLAResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Comp1, Comp2, Comp3
      END SUBROUTINE GiD_BeginLAResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1) = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3) = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginLAResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr, &
                               Comp_ptr(1),Comp_ptr(2),Comp_ptr(3))

  END SUBROUTINE GiD_BeginLAResult
  
  SUBROUTINE GiD_fBeginLAResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Comp1,Comp2,Comp3)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(3)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_LocalAxes,Where,GaussPointsName,RangeTable,3,Comp_f)

  END SUBROUTINE GiD_fBeginLAResult
  

  !--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginComplexScalarResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Re,Im)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Re, Im

    CHARACTER(1,KIND=C_CHAR), TARGET  :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: Comp_cpp(GID_MAX_CHAR_LEN,2)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr
    TYPE(C_PTR) :: RangeTable_ptr
    TYPE(C_PTR) :: Comp_ptr(2)

    INTERFACE
      SUBROUTINE GiD_BeginComplexScalarResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Re,Im) &
        BIND(C,NAME='GiD_BeginComplexScalarResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Re, Im
      END SUBROUTINE GiD_BeginComplexScalarResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1) = C_PTR_STRING(Re,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Im,Comp_cpp(:,2))    
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginComplexScalarResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr, &
                                   Comp_ptr(1),Comp_ptr(2))

  END SUBROUTINE GiD_BeginComplexScalarResult

  SUBROUTINE GiD_fBeginComplexScalarResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Re,Im)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Re, Im

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(2)

    Comp_f(1) = Re;  Comp_f(2) = Im;

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_ComplexScalar,Where,GaussPointsName,RangeTable,2,Comp_f)

  END SUBROUTINE GiD_fBeginComplexScalarResult
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginComplexVectorResult(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Rex,Imx,Rey,Imy,Rez,Imz)

    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Rex, Imx, Rey, Imy, Rez, Imz

    CHARACTER(1,KIND=C_CHAR), TARGET  :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: Comp_cpp(GID_MAX_CHAR_LEN,6)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr
    TYPE(C_PTR) :: RangeTable_ptr
    TYPE(C_PTR) :: Comp_ptr(6)

    INTERFACE
      SUBROUTINE GiD_BeginComplexVectorResult_cpp(Result,Analysis,Step,Where,GaussPointsName,RangeTable,Rex,Imx,Rey,Imy,Rez,Imz) &
        BIND(C,NAME='GiD_BeginComplexVectorResult')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),   VALUE,  INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc), VALUE,  INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: RangeTable
        TYPE(C_PTR),      VALUE,  INTENT(IN) :: Rex, Imx, Rey, Imy, Rez, Imz
      END SUBROUTINE GiD_BeginComplexVectorResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Comp_ptr(1) = C_PTR_STRING(Rex,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Imx,Comp_cpp(:,2))    
    Comp_ptr(3) = C_PTR_STRING(Rey,Comp_cpp(:,3))
    Comp_ptr(4) = C_PTR_STRING(Imy,Comp_cpp(:,4))    
    Comp_ptr(5) = C_PTR_STRING(Rez,Comp_cpp(:,5))
    Comp_ptr(6) = C_PTR_STRING(Imz,Comp_cpp(:,6))    

    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginComplexVectorResult_cpp(Result_cpp,Analysis_cpp,Step,Where,GaussPointsName_ptr,RangeTable_ptr, &
                                   Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4),Comp_ptr(5),Comp_ptr(6))

  END SUBROUTINE GiD_BeginComplexVectorResult

  SUBROUTINE GiD_fBeginComplexVectorResult(fd,Result,Analysis,Step,Where,GaussPointsName,RangeTable,Rex,Imx,Rey,Imy,Rez,Imz)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: RangeTable
    CHARACTER(LEN=*), INTENT(IN) :: Rex, Imx, Rey, Imy, Rez, Imz

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(6)

    Comp_f(1) = Rex;  Comp_f(2) = Imx;
    Comp_f(2) = Rey;  Comp_f(3) = Imy;
    Comp_f(4) = Rez;  Comp_f(5) = Imz;

    CALL GiD_fBeginResult(fd,Result,Analysis,Step,GiD_ComplexScalar,Where,GaussPointsName,RangeTable,6,Comp_f)

  END SUBROUTINE GiD_fBeginComplexVectorResult
  
!--------------------------------------------------------------------------
  
    SUBROUTINE GiD_BeginResultHeader(Result,Analysis,Step,Type,Where,GaussPointsName)

    CHARACTER(LEN=*),  INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*),  INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),    INTENT(IN) :: Step           !
    TYPE(GiD_ResType), INTENT(IN) :: Type           ! result tipe (scalar,vector,matrix)
    TYPE(GiD_ResLoc),  INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*),  INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)

    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr

    INTERFACE
      SUBROUTINE GiD_BeginResultHeader_cpp(Result,Analysis,Step,Type,Where,GaussPointsName) &
        BIND(C,NAME='GiD_BeginResultHeader')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, GiD_ResType, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: Step               !
        TYPE(GiD_ResType), VALUE, INTENT(IN) :: Type               ! result tipe (scalar,vector,matrix)
        TYPE(GiD_ResLoc),  VALUE, INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
      END SUBROUTINE GiD_BeginResultHeader_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginResultHeader_cpp(Result_cpp,Analysis_cpp,Step,Type,Where,GaussPointsName_ptr)

  END SUBROUTINE GiD_BeginResultHeader
  
  SUBROUTINE GiD_fBeginResultHeader(fd,Result,Analysis,Step,Type,Where,GaussPointsName)

    TYPE(GiD_File),    INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*),  INTENT(IN) :: Result         ! result name
    CHARACTER(LEN=*),  INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),    INTENT(IN) :: Step           !
    TYPE(GiD_ResType), INTENT(IN) :: Type           ! result tipe (scalar,vector,matrix)
    TYPE(GiD_ResLoc),  INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*),  INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)

    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr

    INTERFACE
      SUBROUTINE GiD_fBeginResultHeader_cpp(fd,Result,Analysis,Step,Type,Where,GaussPointsName) &
        BIND(C,NAME='GiD_fBeginResultHeader')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, GiD_ResType, GiD_File, C_PTR
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd                 ! file
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result(1)          ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: Step               !
        TYPE(GiD_ResType), VALUE, INTENT(IN) :: Type               ! result tipe (scalar,vector,matrix)
        TYPE(GiD_ResLoc),  VALUE, INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
      END SUBROUTINE GiD_fBeginResultHeader_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_fBeginResultHeader_cpp(fd,Result_cpp,Analysis_cpp,Step,Type,Where,GaussPointsName_ptr)

  END SUBROUTINE GiD_fBeginResultHeader
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_ScalarComp(Comp)

    CHARACTER(LEN=*), INTENT(IN) :: Comp
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN)
    TYPE(C_PTR) :: Comp_ptr

    INTERFACE
      SUBROUTINE GiD_ScalarComp_cpp(Comp) BIND(C,NAME='GiD_ScalarComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp
      END SUBROUTINE GiD_ScalarComp_cpp
    END INTERFACE

    ! Optional values:
    Comp_ptr = C_PTR_STRING(Comp,Comp_cpp)

    CALL GiD_ScalarComp_cpp(Comp_ptr)

  END SUBROUTINE GiD_ScalarComp  

  SUBROUTINE GiD_fScalarComp(fd,Comp)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Comp

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(1)

    Comp_f(1)  = Comp

    CALL GiD_fResultComponents(fd,1,Comp_f)

  END SUBROUTINE GiD_fScalarComp  

!--------------------------------------------------------------------------

  SUBROUTINE GiD_VectorComp(Comp1,Comp2,Comp3,Comp4)

    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,4)
    TYPE(C_PTR) :: Comp_ptr(4)

    INTERFACE
      SUBROUTINE GiD_VectorComp_cpp(Comp1,Comp2,Comp3,Comp4) BIND(C,NAME='GiD_VectorComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp1, Comp2, Comp3, Comp4
      END SUBROUTINE GiD_VectorComp_cpp
    END INTERFACE
    
    Comp_ptr(1)  = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2)  = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3)  = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4)  = C_PTR_STRING(Comp4,Comp_cpp(:,4))

    CALL GiD_VectorComp_cpp(Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4))

  END SUBROUTINE GiD_VectorComp  

  SUBROUTINE GiD_fVectorComp(fd,Comp1,Comp2,Comp3,Comp4)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(4)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3;  Comp_f(4) = Comp4

    CALL GiD_fResultComponents(fd,4,Comp_f)

  END SUBROUTINE GiD_fVectorComp
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_2DMatrixComp(Comp1,Comp2,Comp3)

    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,3)
    TYPE(C_PTR) :: Comp_ptr(3)

    INTERFACE
      SUBROUTINE GiD_2DMatrixComp_cpp(Comp1,Comp2,Comp3) BIND(C,NAME='GiD_2DMatrixComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp1, Comp2, Comp3
      END SUBROUTINE GiD_2DMatrixComp_cpp
    END INTERFACE
    
    Comp_ptr(1)  = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2)  = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3)  = C_PTR_STRING(Comp3,Comp_cpp(:,3))

    CALL GiD_2DMatrixComp_cpp(Comp_ptr(1),Comp_ptr(2),Comp_ptr(3))

  END SUBROUTINE GiD_2DMatrixComp
  
  SUBROUTINE GiD_f2DMatrixComp(fd,Comp1,Comp2,Comp3)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3
    
    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(3)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3

    CALL GiD_fResultComponents(fd,3,Comp_f)

  END SUBROUTINE GiD_f2DMatrixComp
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_3DMatrixComp(Comp1,Comp2,Comp3,Comp4,Comp5,Comp6)

    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4, Comp5, Comp6
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,6)
    TYPE(C_PTR) :: Comp_ptr(6)

    INTERFACE
      SUBROUTINE GiD_3DMatrixComp_cpp(Comp1,Comp2,Comp3,Comp4,Comp5,Comp6) BIND(C,NAME='GiD_3DMatrixComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp1, Comp2, Comp3, Comp4, Comp5, Comp6
      END SUBROUTINE GiD_3DMatrixComp_cpp
    END INTERFACE
    
    Comp_ptr(1) = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3) = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4) = C_PTR_STRING(Comp4,Comp_cpp(:,4))
    Comp_ptr(5) = C_PTR_STRING(Comp5,Comp_cpp(:,5))
    Comp_ptr(6) = C_PTR_STRING(Comp6,Comp_cpp(:,6))

    CALL GiD_3DMatrixComp_cpp(Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4),Comp_ptr(5),Comp_ptr(6))

  END SUBROUTINE GiD_3DMatrixComp
  
  SUBROUTINE GiD_f3DMatrixComp(fd,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4, Comp5, Comp6
    
    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(6)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3;
    Comp_f(4) = Comp4;  Comp_f(5) = Comp5;  Comp_f(6) = Comp6

    CALL GiD_fResultComponents(fd,6,Comp_f)

  END SUBROUTINE GiD_f3DMatrixComp
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_PDMComp(Comp1,Comp2,Comp3,Comp4)

    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,4)
    TYPE(C_PTR) :: Comp_ptr(4)

    INTERFACE
      SUBROUTINE GiD_PDMComp_cpp(Comp1,Comp2,Comp3,Comp4) BIND(C,NAME='GiD_PDMComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp1, Comp2, Comp3, Comp4
      END SUBROUTINE GiD_PDMComp_cpp
    END INTERFACE
    
    Comp_ptr(1) = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2) = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3) = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4) = C_PTR_STRING(Comp4,Comp_cpp(:,4))

    CALL GiD_PDMComp_cpp(Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4))

  END SUBROUTINE GiD_PDMComp

  SUBROUTINE GiD_fPDMComp(fd,Comp1,Comp2,Comp3,Comp4)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4
    
    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(4)

    Comp_f(1) = Comp1;  Comp_f(2) = Comp2;  Comp_f(3) = Comp3;  Comp_f(4) = Comp4

    CALL GiD_fResultComponents(fd,4,Comp_f)

  END SUBROUTINE GiD_fPDMComp
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_MainMatrixComp(Comp1,Comp2,Comp3,Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12)

    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,12)
    TYPE(C_PTR) :: Comp_ptr(12)

    INTERFACE
      SUBROUTINE GiD_MainMatrixComp_cpp(Comp1,Comp2,Comp3,Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12) &
        BIND(C,NAME='GiD_MainMatrixComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp1, Comp2, Comp3, Comp4,  Comp5,  Comp6
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp7, Comp8, Comp9, Comp10, Comp11, Comp12
      END SUBROUTINE GiD_MainMatrixComp_cpp
    END INTERFACE
    
    Comp_ptr(1)  = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2)  = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3)  = C_PTR_STRING(Comp3,Comp_cpp(:,3))
    Comp_ptr(4)  = C_PTR_STRING(Comp4,Comp_cpp(:,4))
    Comp_ptr(5)  = C_PTR_STRING(Comp5,Comp_cpp(:,5))
    Comp_ptr(6)  = C_PTR_STRING(Comp6,Comp_cpp(:,6))
    Comp_ptr(7)  = C_PTR_STRING(Comp7,Comp_cpp(:,7))
    Comp_ptr(8)  = C_PTR_STRING(Comp8,Comp_cpp(:,8))
    Comp_ptr(9)  = C_PTR_STRING(Comp9,Comp_cpp(:,9))
    Comp_ptr(10) = C_PTR_STRING(Comp10,Comp_cpp(:,10))
    Comp_ptr(11) = C_PTR_STRING(Comp11,Comp_cpp(:,11))
    Comp_ptr(12) = C_PTR_STRING(Comp12,Comp_cpp(:,12))

    CALL GiD_MainMatrixComp_cpp(Comp_ptr(1),Comp_ptr(2),Comp_ptr(3),Comp_ptr(4),Comp_ptr(5),Comp_ptr(6),   &
                                Comp_ptr(7),Comp_ptr(8),Comp_ptr(9),Comp_ptr(10),Comp_ptr(11),Comp_ptr(12) )

  END SUBROUTINE GiD_MainMatrixComp  
  
  SUBROUTINE GiD_fMainMatrixComp(fd,Comp1,Comp2,Comp3,Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3, Comp4,Comp5,Comp6,Comp7,Comp8,Comp9,Comp10,Comp11,Comp12
    
    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(12)

    Comp_f(1) = Comp1;  Comp_f(2)  = Comp2;  Comp_f(3)  = Comp3;  Comp_f(4)  = Comp4
    Comp_f(5) = Comp5;  Comp_f(6)  = Comp6;  Comp_f(7)  = Comp7;  Comp_f(8)  = Comp8
    Comp_f(9) = Comp9;  Comp_f(10) = Comp10; Comp_f(11) = Comp11; Comp_f(12) = Comp12

    CALL GiD_fResultComponents(fd,12,Comp_f)

  END SUBROUTINE GiD_fMainMatrixComp  

!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_LAComponents(Comp1,Comp2,Comp3)

    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,3)
    TYPE(C_PTR) :: Comp_ptr(3)

    INTERFACE
      SUBROUTINE GiD_LAComponents_cpp(Comp1,Comp2,Comp3) BIND(C,NAME='GiD_LAComponents')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Comp1, Comp2, Comp3
      END SUBROUTINE GiD_LAComponents_cpp
    END INTERFACE
    
    Comp_ptr(1)  = C_PTR_STRING(Comp1,Comp_cpp(:,1))
    Comp_ptr(2)  = C_PTR_STRING(Comp2,Comp_cpp(:,2))
    Comp_ptr(3)  = C_PTR_STRING(Comp3,Comp_cpp(:,3))

    CALL GiD_LAComponents_cpp(Comp_ptr(1),Comp_ptr(3),Comp_ptr(3))

  END SUBROUTINE GiD_LAComponents
  
  SUBROUTINE GiD_fLAComponents(fd,Comp1,Comp2,Comp3)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Comp1, Comp2, Comp3
    
    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(3)

    Comp_f(1) = Comp1; Comp_f(2) = Comp2; Comp_f(3) = Comp3

    CALL GiD_fResultComponents(fd,3,Comp_f)

  END SUBROUTINE GiD_fLAComponents
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_ComplexScalarComp(Re,Im)

    CHARACTER(LEN=*), INTENT(IN) :: Re,Im
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Re_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Im_cpp(GID_MAX_CHAR_LEN)
    TYPE(C_PTR) :: Re_ptr, Im_ptr

    INTERFACE
      SUBROUTINE GiD_ComplexScalarComp_cpp(Re,Im) BIND(C,NAME='GiD_ComplexScalarComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Re,Im
      END SUBROUTINE GiD_ComplexScalarComp_cpp
    END INTERFACE

    ! Optional values:
    Re_ptr = C_PTR_STRING(Re,Re_cpp)
    Im_ptr = C_PTR_STRING(Im,Im_cpp)

    CALL GiD_ComplexScalarComp_cpp(Re_ptr,Im_ptr)

  END SUBROUTINE GiD_ComplexScalarComp  

  SUBROUTINE GiD_fComplexScalarComp(fd,Re,Im)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Re,Im

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(2)

    Comp_f(1)  = Re
    Comp_f(2)  = Im

    CALL GiD_fResultComponents(fd,2,Comp_f)

  END SUBROUTINE GiD_fComplexScalarComp  
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_ComplexVectorComp(Rex,Imx,Rey,Imy,Rez,Imz)

    CHARACTER(LEN=*), INTENT(IN) :: Rex,Imx,Rey,Imy,Rez,Imz
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: Rex_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Imx_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Rey_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Imy_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Rez_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Imz_cpp(GID_MAX_CHAR_LEN)
    TYPE(C_PTR) :: Rex_ptr, Imx_ptr, Rey_ptr, Imy_ptr, Rez_ptr, Imz_ptr

    INTERFACE
      SUBROUTINE GiD_ComplexVectorComp_cpp(Rex,Imx,Rey,Imy,Rez,Imz) BIND(C,NAME='GiD_ComplexVectorComp')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE, INTENT(IN) :: Rex,Imx,Rey,Imy,Rez,Imz
      END SUBROUTINE GiD_ComplexVectorComp_cpp
    END INTERFACE

    ! Optional values:
    Rex_ptr = C_PTR_STRING(Rex,Rex_cpp)
    Imx_ptr = C_PTR_STRING(Imx,Imx_cpp)
    Rey_ptr = C_PTR_STRING(Rey,Rey_cpp)
    Imy_ptr = C_PTR_STRING(Imy,Imy_cpp)
    Rez_ptr = C_PTR_STRING(Rez,Rez_cpp)
    Imz_ptr = C_PTR_STRING(Imz,Imz_cpp)

    CALL GiD_ComplexVectorComp_cpp(Rex_ptr, Imx_ptr, Rey_ptr, Imy_ptr, Rez_ptr, Imz_ptr)

  END SUBROUTINE GiD_ComplexVectorComp  

  SUBROUTINE GiD_fComplexVectorComp(fd,Rex,Imx,Rey,Imy,Rez,Imz)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Rex,Imx,Rey,Imy,Rez,Imz

    CHARACTER(GID_MAX_CHAR_LEN) :: Comp_f(6)

    Comp_f(1)  = Rex
    Comp_f(2)  = Imx
    Comp_f(3)  = Rey
    Comp_f(4)  = Imy
    Comp_f(5)  = Rez
    Comp_f(6)  = Imz

    CALL GiD_fResultComponents(fd,6,Comp_f)

  END SUBROUTINE GiD_fComplexVectorComp    
    
!--------------------------------------------------------------------------
    
  SUBROUTINE GiD_ResultUnit( unit_name )
    
    CHARACTER(LEN=*), INTENT(IN) :: unit_name
    
    CHARACTER(1,KIND=C_CHAR) :: unit_name_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_ResultUnit_cpp( unit_name_cpp ) BIND(C,NAME='GiD_ResultUnit')      
        IMPORT :: C_CHAR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: unit_name_cpp(1)
      END SUBROUTINE GiD_ResultUnit_cpp
    END INTERFACE

    unit_name_cpp = C_STRING(unit_name)
    call GiD_ResultUnit_cpp(unit_name_cpp)

  END SUBROUTINE GiD_ResultUnit

  SUBROUTINE GiD_fResultUnit( fd, unit_name )

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: unit_name
    
    CHARACTER(1,KIND=C_CHAR) :: unit_name_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fResultUnit_cpp( fd, unit_name_cpp ) BIND(C,NAME='GiD_fResultUnit')      
        IMPORT :: GiD_File, C_CHAR
        TYPE(GiD_File), VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: unit_name_cpp(1)
      END SUBROUTINE GiD_fResultUnit_cpp
    END INTERFACE

    unit_name_cpp = C_STRING(unit_name)
    call GiD_fResultUnit_cpp(fd,unit_name_cpp)

  END SUBROUTINE GiD_fResultUnit

!--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginResultGroup(Analysis,Step,Where,GaussPointsName)

    CHARACTER(LEN=*), INTENT(IN) :: Analysis        ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step            !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where           ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR) :: Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr

    INTERFACE
      SUBROUTINE GiD_BeginResultGroup_cpp(Analysis,Step,Where,GaussPointsName) BIND(C,NAME='GiD_BeginResultGroup')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc),  VALUE, INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
      END SUBROUTINE GiD_BeginResultGroup_cpp
    END INTERFACE
    
    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_BeginResultGroup_cpp(Analysis_cpp,Step,Where,GaussPointsName_ptr)

  END SUBROUTINE GiD_BeginResultGroup
  
  SUBROUTINE GiD_fBeginResultGroup(fd,Analysis,Step,Where,GaussPointsName)

    TYPE(GiD_File),   INTENT(IN) :: fd             ! file
    CHARACTER(LEN=*), INTENT(IN) :: Analysis       ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),   INTENT(IN) :: Step           !
    TYPE(GiD_ResLoc), INTENT(IN) :: Where          ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*), INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    
    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR) :: Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR) :: GaussPointsName_ptr

    INTERFACE
      SUBROUTINE GiD_fBeginResultGroup_cpp(fd,Analysis,Step,Where,GaussPointsName) BIND(C,NAME='GiD_fBeginResultGroup')
        IMPORT :: C_DOUBLE, GiD_ResLoc, C_CHAR, GiD_File, C_PTR
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd                 ! file
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis(1)        ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: Step               !
        TYPE(GiD_ResLoc),  VALUE, INTENT(IN) :: Where              ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: GaussPointsName    ! gauss point name (just ongausspoint)
      END SUBROUTINE GiD_fBeginResultGroup_cpp
    END INTERFACE
    
    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    Analysis_cpp = C_STRING(Analysis)

    CALL GiD_fBeginResultGroup_cpp(fd,Analysis_cpp,Step,Where,GaussPointsName_ptr)

  END SUBROUTINE GiD_fBeginResultGroup
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_ResultDescription(Result,Type)

    CHARACTER(LEN=*),  INTENT(IN) :: Result         ! result name
    TYPE(GiD_ResType), INTENT(IN) :: Type           ! result tipe (scalar,vector,matrix)
    
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_ResultDescription_cpp(Result_cpp,Type) BIND(C,NAME='GiD_ResultDescription')
        IMPORT :: C_CHAR, GiD_ResType
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result_cpp(1)      ! result name
        TYPE(GiD_ResType), VALUE, INTENT(IN) :: Type               ! result tipe (scalar,vector,matrix)
      END SUBROUTINE GiD_ResultDescription_cpp
    END INTERFACE

    Result_cpp   = C_STRING(Result)
    CALL GiD_ResultDescription_cpp(Result_cpp,Type)

  END SUBROUTINE GiD_ResultDescription
  
  SUBROUTINE GiD_fResultDescription(fd,Result,Type)

    TYPE(GiD_File),    INTENT(IN) :: fd             ! GiD file
    CHARACTER(LEN=*),  INTENT(IN) :: Result         ! result name
    TYPE(GiD_ResType), INTENT(IN) :: Type           ! result tipe (scalar,vector,matrix)

    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fResultDescription_cpp(fd,Result_cpp,Type) BIND(C,NAME='GiD_fResultDescription')
        IMPORT :: C_CHAR, GiD_ResType, GiD_File
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result_cpp(1)      ! result name
        TYPE(GiD_ResType), VALUE, INTENT(IN) :: Type               ! result tipe (scalar,vector,matrix)
      END SUBROUTINE GiD_fResultDescription_cpp
    END INTERFACE

    Result_cpp   = C_STRING(Result)
    CALL GiD_fResultDescription_cpp(fd,Result_cpp,Type)

  END SUBROUTINE GiD_fResultDescription
  
!--------------------------------------------------------------------------

    SUBROUTINE GiD_BeginOnMeshGroup(Label)

    CHARACTER(LEN=*), INTENT(IN) :: Label

    CHARACTER(1,KIND=C_CHAR) :: Label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_BeginOnMeshGroup_cpp(Label_cpp) BIND(C,NAME='GiD_BeginOnMeshGroup')
        IMPORT :: C_CHAR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Label_cpp(1)
      END SUBROUTINE GiD_BeginOnMeshGroup_cpp
    END INTERFACE

    Label_cpp = C_STRING(Label)
    CALL GiD_BeginOnMeshGroup_cpp(Label_cpp)

  END SUBROUTINE GiD_BeginOnMeshGroup
  
  SUBROUTINE GiD_fBeginOnMeshGroup(fd,Label)

    TYPE(GiD_File),   INTENT(IN) :: fd
    CHARACTER(LEN=*), INTENT(IN) :: Label

    CHARACTER(1,KIND=C_CHAR) :: Label_cpp(MAX_STRING_LEN)

    INTERFACE
      SUBROUTINE GiD_fBeginOnMeshGroup_cpp(fd,Label_cpp) BIND(C,NAME='GiD_fBeginOnMeshGroup')
        IMPORT :: C_CHAR, GiD_File
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Label_cpp(1)
      END SUBROUTINE GiD_fBeginOnMeshGroup_cpp
    END INTERFACE

    Label_cpp = C_STRING(Label)
    CALL GiD_fBeginOnMeshGroup_cpp(fd,label_cpp)

  END SUBROUTINE GiD_fBeginOnMeshGroup

!--------------------------------------------------------------------------
!  AUXILIAR ROUTINES
!--------------------------------------------------------------------------
  
  PURE FUNCTION C_STRING(f_string)

    CHARACTER(LEN=*), INTENT(IN) :: f_string
    CHARACTER(1,KIND=C_CHAR) :: C_STRING(LEN_TRIM(f_string)+1)
    
    C_STRING = TRANSFER(TRIM(f_string)//C_NULL_CHAR,C_CHAR_ARRAY)

  END FUNCTION C_STRING

  FUNCTION C_PTR_STRING(string_for,string_cpp)

    CHARACTER(LEN=*),                     INTENT(IN)    :: string_for                   ! fortran string
    CHARACTER(LEN=1,KIND=C_CHAR), TARGET, INTENT(INOUT) :: string_cpp(GID_MAX_CHAR_LEN) ! cpp format string
    TYPE(C_PTR) :: C_PTR_STRING                                                         ! pointer to cpp format string

    C_PTR_STRING = C_NULL_PTR
    IF( string_for .NE. GID_STRING_NULL ) THEN
      string_cpp   = C_STRING(string_for)
      C_PTR_STRING = C_LOC(string_cpp)
    END IF

  END FUNCTION C_PTR_STRING

!--------------------------------------------------------------------------
! SPECIAL INTERFACES
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_WriteCoordinates3D(id,x,y,z)

    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: x, y, z

    CALL GiD_WriteCoordinates_cpp(id,x,y,z)

  END SUBROUTINE GiD_WriteCoordinates3D
  
  SUBROUTINE GiD_fWriteCoordinates3D(fd,id,x,y,z)

    TYPE(GiD_File), INTENT(IN) :: fd
    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: x, y, z

    CALL GiD_fWriteCoordinates_cpp(fd,id,x,y,z)

  END SUBROUTINE GiD_fWriteCoordinates3D
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_WriteCoordinates_f(id,xyz)

    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: xyz(:)

    SELECT CASE (SIZE(xyz))
    CASE (2)
      CALL GiD_WriteCoordinates2D(id,xyz(1),xyz(2))
    CASE (3)
      CALL GiD_WriteCoordinates_cpp(id,xyz(1),xyz(2),xyz(3))
    END SELECT

  END SUBROUTINE GiD_WriteCoordinates_f

  SUBROUTINE GiD_fWriteCoordinates_f(fd,id,xyz)

    TYPE(GiD_File), INTENT(IN) :: fd
    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: xyz(:)

    SELECT CASE (SIZE(xyz))
    CASE (2)
      CALL GiD_fWriteCoordinates2D(fd,id,xyz(1),xyz(2))
    CASE (3)
      CALL GiD_fWriteCoordinates_cpp(fd,id,xyz(1),xyz(2),xyz(3))
    END SELECT

  END SUBROUTINE GiD_fWriteCoordinates_f

!--------------------------------------------------------------------------

  SUBROUTINE GiD_WriteElementMat_f(id,conec,material)
    
    INTEGER(C_INT), INTENT(IN) :: id
    INTEGER(C_INT), INTENT(IN) :: conec(:)
    INTEGER(C_INT), INTENT(IN) :: material
   
    INTEGER(C_INT) :: conec_cpp(SIZE(conec)+1)

    conec_cpp(:SIZE(conec))  = conec
    conec_cpp(SIZE(conec)+1) = material
    CALL GiD_WriteElementMat_cpp(id,conec_cpp)

  END SUBROUTINE GiD_WriteElementMat_f

  SUBROUTINE GiD_fWriteElementMat_f(fd,id,conec,material)
    
    TYPE(GiD_File), INTENT(IN) :: fd
    INTEGER(C_INT), INTENT(IN) :: id
    INTEGER(C_INT), INTENT(IN) :: conec(:)
    INTEGER(C_INT), INTENT(IN) :: material
   
    INTEGER(C_INT) :: conec_cpp(SIZE(conec)+1)

    conec_cpp(:SIZE(conec))  = conec
    conec_cpp(SIZE(conec)+1) = material
    CALL GiD_fWriteElementMat_cpp(fd,id,conec_cpp)

  END SUBROUTINE GiD_fWriteElementMat_f
  
!--------------------------------------------------------------------------
  
  SUBROUTINE GiD_Write3DVector(id,x,y,z)

    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: x,y,z

    CALL GiD_WriteVector_cpp(id,x,y,z)

  END SUBROUTINE GiD_Write3DVector
  
  SUBROUTINE GiD_fWrite3DVector(fd,id,x,y,z)

    TYPE(GiD_File), INTENT(IN) :: fd
    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: x,y,z

    CALL GiD_fWriteVector_cpp(fd,id,x,y,z)

  END SUBROUTINE GiD_fWrite3DVector
  
!--------------------------------------------------------------------------

  SUBROUTINE GiD_WriteVector_f(id,vector)

    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: vector(:)

    SELECT CASE (SIZE(vector))
    CASE (2)
      CALL GiD_Write2DVector(id,vector(1),vector(2))
    CASE (3)
      CALL GiD_WriteVector_cpp(id,vector(1),vector(2),vector(3))
    END SELECT

  END SUBROUTINE GiD_WriteVector_f

  SUBROUTINE GiD_fWriteVector_f(fd,id,vector)

    TYPE(GiD_File), INTENT(IN) :: fd
    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: vector(:)

    SELECT CASE (SIZE(vector))
    CASE (2)
      CALL GiD_fWrite2DVector(fd,id,vector(1),vector(2))
    CASE (3)
      CALL GiD_fWriteVector_cpp(fd,id,vector(1),vector(2),vector(3))
    END SELECT

  END SUBROUTINE GiD_fWriteVector_f

!--------------------------------------------------------------------------

  SUBROUTINE GiD_WriteVectorModule_f(id,vector)

    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: vector(4)

    CALL GiD_WriteVectorModule_cpp(id,vector(1),vector(2),vector(3),vector(4))

  END SUBROUTINE GiD_WriteVectorModule_f

  SUBROUTINE GiD_fWriteVectorModule_f(fd,id,vector)

    TYPE(GiD_File), INTENT(IN) :: fd
    INTEGER(C_INT), INTENT(IN) :: id
    REAL(C_DOUBLE), INTENT(IN) :: vector(4)

    CALL GiD_fWriteVectorModule_cpp(fd,id,vector(1),vector(2),vector(3),vector(4))

  END SUBROUTINE GiD_fWriteVectorModule_f

!--------------------------------------------------------------------------

  SUBROUTINE GiD_BeginResult(Result,Analysis,Step,Type,Where,GaussPointsName,RangeTable,nComp,Comp)

    CHARACTER(LEN=*),  INTENT(IN) :: Result     ! result name
    CHARACTER(LEN=*),  INTENT(IN) :: Analysis   ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),    INTENT(IN) :: Step       !
    TYPE(GiD_ResType), INTENT(IN) :: Type       ! result tipe (scalar,vector,matrix)
    TYPE(GiD_ResLoc),  INTENT(IN) :: Where      ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*),  INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*),  INTENT(IN) :: RangeTable ! 
    INTEGER(C_INT),    INTENT(IN) :: nComp 
    CHARACTER(LEN=*),  INTENT(IN) :: Comp(nComp)

    CHARACTER(1,KIND=C_CHAR), TARGET :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,nComp)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR)  :: GaussPointsName_ptr
    TYPE(C_PTR)  :: RangeTable_ptr
    TYPE(C_PTR), TARGET :: Comp_ptr(nComp)
    INTEGER :: i

    INTERFACE
      SUBROUTINE GiD_BeginResult_cpp(Result_cpp,Analysis_cpp,Step,Type,Where,GaussPointsName,RangeTable,nComp,Comp) &
        BIND(C,NAME='GiD_BeginResult')
        IMPORT :: GiD_File, GiD_ResType, GiD_ResLoc, C_CHAR, C_INT, C_DOUBLE, C_PTR
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result_cpp(1)   ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis_cpp(1) ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: Step       !
        TYPE(GiD_ResType), VALUE, INTENT(IN) :: Type       ! result tipe (scalar,vector,matrix)
        TYPE(GiD_ResLoc),  VALUE, INTENT(IN) :: Where      ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: RangeTable ! 
        INTEGER(C_INT),    VALUE, INTENT(IN) :: nComp
        TYPE(C_PTR),       VALUE, INTENT(IN) :: Comp
      END SUBROUTINE GiD_BeginResult_cpp
    END INTERFACE
    
    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    DO i = 1, nComp
      IF (Comp(i).NE.GID_STRING_NULL) THEN
        Comp_cpp(:,i) = C_STRING(Comp(i)); Comp_ptr(i)  = C_LOC(Comp_cpp(:,i))
      ELSE
        Comp_cpp(:,i) = C_NULL_CHAR; Comp_ptr(i) = C_LOC(Comp_cpp(:,i))
      END IF
    END DO
   
    Analysis_cpp = C_STRING(Analysis)
    Result_cpp   = C_STRING(Result)
    call GiD_BeginResult_cpp(Result_cpp,Analysis_cpp,Step,Type,Where,RangeTable_ptr,GaussPointsName_ptr,&
                             nComp,C_LOC(Comp_ptr))

  END SUBROUTINE GiD_BeginResult

  SUBROUTINE GiD_fBeginResult(fd,Result,Analysis,Step,Type,Where,GaussPointsName,RangeTable,nComp,Comp)

    TYPE(GiD_File),    INTENT(IN) :: fd         ! file
    CHARACTER(LEN=*),  INTENT(IN) :: Result     ! result name
    CHARACTER(LEN=*),  INTENT(IN) :: Analysis   ! analysis type (timestep,load,...)
    REAL(C_DOUBLE),    INTENT(IN) :: Step       !
    TYPE(GiD_ResType), INTENT(IN) :: Type       ! result tipe (scalar,vector,matrix)
    TYPE(GiD_ResLoc),  INTENT(IN) :: Where      ! result location (onnodes,ongausspoint)
    CHARACTER(LEN=*),  INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
    CHARACTER(LEN=*),  INTENT(IN) :: RangeTable ! 
    INTEGER(C_INT),    INTENT(IN) :: nComp 
    CHARACTER(LEN=*),  INTENT(IN) :: Comp(nComp)

    CHARACTER(1,KIND=C_CHAR), TARGET  :: GaussPointsName_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: RangeTable_cpp(GID_MAX_CHAR_LEN)
    CHARACTER(1,KIND=C_CHAR), TARGET  :: Comp_cpp(GID_MAX_CHAR_LEN,nComp)
    CHARACTER(1,KIND=C_CHAR) :: Result_cpp(MAX_STRING_LEN), Analysis_cpp(MAX_STRING_LEN)
    TYPE(C_PTR)  :: GaussPointsName_ptr
    TYPE(C_PTR)  :: RangeTable_ptr
    TYPE(C_PTR), TARGET :: Comp_ptr(nComp)
    INTEGER :: i

    INTERFACE
      SUBROUTINE GiD_fBeginResult_cpp(fd,Result_cpp,Analysis_cpp,Step,Type,Where,GaussPointsName,RangeTable,nComp,Comp) &
        BIND(C,NAME='GiD_fBeginResult')
        IMPORT :: GiD_File, GiD_ResType, GiD_ResLoc, C_CHAR, C_INT, C_DOUBLE, C_PTR
        TYPE(GiD_File),    VALUE, INTENT(IN) :: fd
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Result_cpp(1)   ! result name
        CHARACTER(1,KIND=C_CHAR), INTENT(IN) :: Analysis_cpp(1) ! analysis type (timestep,load,...)
        REAL(C_DOUBLE),    VALUE, INTENT(IN) :: Step            !
        TYPE(GiD_ResType), VALUE, INTENT(IN) :: Type            ! result tipe (scalar,vector,matrix)
        TYPE(GiD_ResLoc),  VALUE, INTENT(IN) :: Where           ! result location (onnodes,ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: GaussPointsName ! gauss point name (just ongausspoint)
        TYPE(C_PTR),       VALUE, INTENT(IN) :: RangeTable
        INTEGER(C_INT),    VALUE, INTENT(IN) :: nComp
        TYPE(C_PTR),       VALUE, INTENT(IN) :: Comp
      END SUBROUTINE GiD_fBeginResult_cpp
    END INTERFACE

    ! Optional values:
    GaussPointsName_ptr = C_PTR_STRING(GaussPointsName,GaussPointsName_cpp)
    RangeTable_ptr      = C_PTR_STRING(RangeTable,RangeTable_cpp)
    Result_cpp   = C_STRING(Result)
    Analysis_cpp = C_STRING(Analysis)
    DO i = 1, nComp
      IF (Comp(i).NE.GID_STRING_NULL) THEN
        Comp_cpp(:,i) = C_STRING(Comp(i)); Comp_ptr(i)  = C_LOC(Comp_cpp(:,i))
      ELSE
        Comp_cpp(:,i) = C_NULL_CHAR; Comp_ptr(i) = C_LOC(Comp_cpp(:,i))
      END IF
    END DO
    
    call GiD_fBeginResult_cpp(fd,Result_cpp,Analysis_cpp,Step,Type,Where,RangeTable_ptr,GaussPointsName_ptr,&
                              nComp,C_LOC(Comp_ptr))

  END SUBROUTINE GiD_fBeginResult

!--------------------------------------------------------------------------

  SUBROUTINE GiD_ResultComponents(nComp,Comp)

    INTEGER(C_INT),   INTENT(IN) :: nComp
    CHARACTER(LEN=*), INTENT(IN) :: Comp(nComp)

    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,nComp)
    TYPE(C_PTR), TARGET :: Comp_ptr(nComp)
    INTEGER :: i

    INTERFACE
      SUBROUTINE GiD_ResultComponents_cpp(nComp,Comp) BIND(C,NAME='GiD_ResultComponents')
        IMPORT :: C_INT, C_PTR
        INTEGER(C_INT), VALUE, INTENT(IN) :: nComp
        TYPE(C_PTR),    VALUE, INTENT(IN) :: Comp
      END SUBROUTINE GiD_ResultComponents_cpp
    END INTERFACE

    DO i = 1, nComp
      Comp_ptr(i) = C_PTR_STRING(Comp(i),Comp_cpp(:,i))
    END DO
    
    call GiD_ResultComponents_cpp(nComp,C_LOC(Comp_ptr))

  END SUBROUTINE GiD_ResultComponents

  SUBROUTINE GiD_fResultComponents(fd,nComp,Comp)

    TYPE(GiD_File),   INTENT(IN) :: fd         ! file
    INTEGER(C_INT),   INTENT(IN) :: nComp
    CHARACTER(LEN=*), INTENT(IN) :: Comp(nComp)

    CHARACTER(1,KIND=C_CHAR), TARGET :: Comp_cpp(GID_MAX_CHAR_LEN,nComp)
    TYPE(C_PTR), TARGET :: Comp_ptr(nComp)
    INTEGER :: i

    INTERFACE
      SUBROUTINE GiD_fResultComponents_cpp(fd,nComp,Comp) BIND(C,NAME='GiD_fResultComponents')
        IMPORT :: GiD_File, C_INT, C_PTR
        TYPE(GiD_File), VALUE, INTENT(IN) :: fd
        INTEGER(C_INT), VALUE, INTENT(IN) :: nComp
        TYPE(C_PTR),    VALUE, INTENT(IN) :: Comp
      END SUBROUTINE GiD_fResultComponents_cpp
    END INTERFACE

    DO i = 1, nComp
      Comp_ptr(i) = C_PTR_STRING(Comp(i),Comp_cpp(:,i))
    END DO
    
    call GiD_fResultComponents_cpp(fd,nComp,C_LOC(Comp_ptr))

  END SUBROUTINE GiD_fResultComponents

!--------------------------------------------------------------------------
! TYPE COMPARISONS
!--------------------------------------------------------------------------

PURE LOGICAL FUNCTION GiD_PostMode_Equal(PostMode_1,PostMode_2)
    TYPE(GiD_PostMode), INTENT(IN) :: PostMode_1, PostMode_2
    GiD_PostMode_Equal = ( PostMode_1%dummy == PostMode_2%dummy )
END FUNCTION GiD_PostMode_Equal

!--------------------------------------------------------------------------
! AUXILIAR ROUTINES
!--------------------------------------------------------------------------

ELEMENTAL INTEGER(GID_INT_KIND) FUNCTION GiD_Int_1(number)
  INTEGER(1), INTENT(IN) :: number
  GiD_Int_1 = int(number,GID_INT_KIND)
END FUNCTION GiD_Int_1
ELEMENTAL INTEGER(GID_INT_KIND) FUNCTION GiD_Int_2(number)
  INTEGER(2), INTENT(IN) :: number
  GiD_Int_2 = int(number,GID_INT_KIND)
END FUNCTION GiD_Int_2
ELEMENTAL INTEGER(GID_INT_KIND) FUNCTION GiD_Int_4(number)
  INTEGER(4), INTENT(IN) :: number
  GiD_Int_4 = int(number,GID_INT_KIND)
END FUNCTION GiD_Int_4
ELEMENTAL INTEGER(GID_INT_KIND) FUNCTION GiD_Int_8(number)
  INTEGER(8), INTENT(IN) :: number
  GiD_Int_8 = int(number,GID_INT_KIND)
END FUNCTION GiD_Int_8

!--

ELEMENTAL REAL(GID_REAL_KIND) FUNCTION GiD_Real_4(number)
  REAL(4), INTENT(IN) :: number
  GiD_Real_4 = int(number,GID_REAL_KIND)
END FUNCTION GiD_Real_4
ELEMENTAL REAL(GID_REAL_KIND) FUNCTION GiD_Real_8(number)
  REAL(8), INTENT(IN) :: number
  GiD_Real_8 = int(number,GID_REAL_KIND)
END FUNCTION GiD_Real_8

END MODULE gidpost
