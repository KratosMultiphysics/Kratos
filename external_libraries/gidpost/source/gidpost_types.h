#ifndef _GIDPOST_TYPES_
#define _GIDPOST_TYPES_

/* opening mode */

typedef enum {
  GiD_PostUndefined = -1,
  GiD_PostAscii=0, 
  GiD_PostAsciiZipped, 
  GiD_PostBinary,
  GiD_PostHDF5
} GiD_PostMode;

#define GiD_PostBinary2 GiD_PostHDF5 /* back compatibility */

/* domain dimension */

typedef enum { GiD_2D = 2, GiD_3D = 3 } GiD_Dimension;

/* post element types */

typedef enum {
  GiD_NoElement = 0,
  GiD_Point,
  GiD_Linear,
  GiD_Triangle,
  GiD_Quadrilateral,
  GiD_Tetrahedra,
  GiD_Hexahedra,
  GiD_Prism,
  GiD_Pyramid,
  GiD_Sphere,
  GiD_Circle
} GiD_ElementType;

typedef enum {
  GiD_Scalar = 0,
  GiD_Vector,
  GiD_Matrix,
  GiD_PlainDeformationMatrix,
  GiD_MainMatrix,
  GiD_LocalAxes,
  GiD_ComplexScalar,
  GiD_ComplexVector,
  GiD_ComplexMatrix
} GiD_ResultType;

typedef enum { 
  GiD_OnNodes=0, 
  GiD_OnGaussPoints, 
  GiD_OnNurbsLine, 
  GiD_OnNurbsSurface, 
  GiD_OnNurbsVolume 
} GiD_ResultLocation;

typedef unsigned int GiD_FILE;

#define GP_OK                   0
#define GP_ERROR_INVALID_STATE -1
#define GP_ERROR_NOMEM         -2
#define GP_ERROR_FILEOPENED    -3
#define GP_ERROR_OPENFAILED    -4
#define GP_ERROR_HANDLEFAIL    -5
#define GP_ERROR_WRITESTRING   -6
#define GP_ERROR_WRITEPOINT    -7
#define GP_ERROR_SCOPE         -8
#define GP_ERROR_NULLSTRING    -9
#define GP_ERROR_ZEROCOMPONENTS -10
#define GP_ERROR_NULLFILE       -11
#define GP_ERROR_NOTINGROUP     -12

#endif // #ifndef _GIDPOST_TYPES_
