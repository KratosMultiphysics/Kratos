/* gidpost 1.7 */
/*
 *  gidpost.cc--
 *
 *    This file implement a C interface for generating postprocessing
 *    results in the 'New postprocess format' of GiD. See declaration
 *    in gidpost.h
 *
 */

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "gidpostInt.h"
#include "gidpost.h"

#ifdef WIN32
#define snprintf _snprintf
#endif

/* ---------------------------------------------------------------------------
 *
 *  Global files
 *
 * ---------------------------------------------------------------------------
 */

static CPostFile * MeshFile = NULL;
static CPostFile * ResultFile = NULL;

static  CPostFile * outputMesh = NULL;

typedef enum {
  POST_UNDEFINED,
  POST_S0,           /* TOP level */
  POST_MESH_S0,      /* MESH header */
  POST_MESH_COORD0,  /* inside a Coordinate block */
  POST_MESH_COORD1,  /* after a Coordinate block but inside a MESH */
  POST_MESH_ELEM,    /* inside an Element block */
  POST_GAUSS_S0,     /* GAUSS point block: implicit */
  POST_GAUSS_GIVEN,  /* GAUSS point block: explicit */
  POST_RANGE_S0,     /* RANGE table block */
  POST_RESULT_SINGLE,    /* Result block */
  POST_RESULT_GROUP, /* Result group block */
  POST_RESULT_DESC,  /* Result description block */
  POST_RESULT_VALUES /* writing values */
} post_state;

static post_state level_mesh = POST_UNDEFINED;
static post_state level_res  = POST_UNDEFINED;

static GP_CONST char * level_desc[] = {
  "UNDEFINED level",
  "TOP level",
  "MESH header",
  "inside a Coordinate block",
  "after a Coordinate block but inside a MESH",
  "inside an Element block",
  "GAUSS point block: implicit",
  "GAUSS point block: explicit",
  "RANGE table block",
  "Result block",
  "Result group block",
  "Result description block",
  "writing values",
  "unknown"
};

static GP_CONST char * GetStateDesc(int i)
{
  static int last = sizeof(level_desc)/sizeof(level_desc[0]) - 1;
  return (i<0 || i >= last) ? level_desc[last] : level_desc[i]; 
} 

static int GP_number_check = 0;   /* number of gauss points to be written */ 
static int gauss_written = 0;      /* number of gauss points written */
static int flag_isgroup = 0;      /* is this a result group */
static int flag_begin_values = 0; /* Is values written ? */
/* not used --
static int active_type;
static int active_dim;
*/

static CBufferValues buffer_values;

inline CPostFile * GetMeshFile()
{
  return (outputMesh = MeshFile ? MeshFile : ResultFile);
}

#define NMAX_DIMS 4
typedef struct {
  GP_CONST char * str;
  size_t dims[NMAX_DIMS];
} SResultTypeInfo;

static SResultTypeInfo _ResultTypeInfo[] = {
  {"Scalar",                  {1,  0, 0, 0}},
  {"Vector",                  {2,  3, 4, 0}},
  {"Matrix",                  {3,  6, 0, 0}},
  {"PlainDeformationMatrix",  {4,  0, 0, 0}},
  {"MainMatrix",              {12, 0, 0, 0}},
  {"LocalAxes",               {3,  0, 0, 0}}
};

GP_CONST char * GetResultTypeName(GiD_ResultType type, size_t s)
{
  static char buffer[255];
  char * ptr;
  int i;
  
  strcpy(buffer, _ResultTypeInfo[int(type)].str);
  ptr = &(buffer[0]) + strlen(buffer);
  if (s) {
    for (i = 0; i <  NMAX_DIMS && _ResultTypeInfo[i].dims[i]; i++) {
      if (s == _ResultTypeInfo[i].dims[i])
        break;
    }
    if (i == NMAX_DIMS)
      printf("Invalid dimension %u for type %s\n", s, buffer);
    else
      sprintf(ptr, ":%u", s);
  }
  return buffer;
}

void GetResultTypeMinMaxValues(GiD_ResultType type, size_t &min, size_t &max)
{
  int i = 0;
  
  min = _ResultTypeInfo[i].dims[0];
  for (i = 1; i < NMAX_DIMS && _ResultTypeInfo[i].dims[i]; i++)
    ;
  max = _ResultTypeInfo[i].dims[i-1];
}

/*int GetResultTypeMinValues(GiD_ResultType type)
{
  return _ResultTypeInfo[int(type)].min_values;
}

int GetResultTypeMaxValues(GiD_ResultType type)
{
  return _ResultTypeInfo[int(type)].max_values;
}*/

CPostFile * NewFile( GiD_PostMode Mode )
{
  switch ( Mode ) {
    case GiD_PostAscii:
      return new CPostAscii;
    case GiD_PostAsciiZipped:
      return new CPostAsciiZ;      
    case GiD_PostBinary:
      return new CPostBinary;      
  }
  return NULL;
};

static GP_CONST char * strElementType[]= {
  "",
  "Point",
  "Linear",
  "Triangle",
  "Quadrilateral",
  "Tetrahedra",
  "Hexahedra",
  "Prism",
  "Pyramid",
  "Sphere",
  "Circle"
};

GP_CONST char * GetElementTypeName( GiD_ElementType type )
{
  return strElementType[int(type)];
}

static int CheckState(post_state s_req, post_state s_cur)
{
  if (s_req!=s_cur) {
    printf("invalid state '%s' should be '%s'\n", GetStateDesc(s_cur), GetStateDesc(s_req));
    return 0;
  }
  return 1;
}

static int ValidateConnectivity(GiD_ElementType etype , int NNode)
{
  int error;

  switch (etype) {
  case GiD_Point:
  case GiD_Sphere:
  case GiD_Circle:
    error = (NNode != 1);
    break;
  case GiD_Linear:
    error = (NNode != 2 && NNode != 3);
    break;
  case GiD_Triangle:
    error = (NNode != 3 && NNode != 6);
    break;
  case GiD_Quadrilateral:
    error = (NNode != 4 && NNode != 8 && NNode != 9);
    break;
  case GiD_Tetrahedra:
    error = (NNode != 4 && NNode != 10);
    break;
  case GiD_Hexahedra:
    error = (NNode != 8 && NNode != 20 && NNode != 27);
    break;
  case GiD_Prism:
    error = (NNode != 6 && NNode != 15);
    break;
  case GiD_Pyramid:
    error = (NNode != 5 && NNode != 13);
    break;
   default:
    printf("invalid type of element code %d", etype);
    return 0;
  }
  if (error) {
    printf("invalid number of nodes '%d' for element type %s\n",
	   NNode, GetElementTypeName(etype));
    return 0;
  }
  return 1;
}

// let's change the double quotes to simple ones, to avoid potential problems
static char *change_quotes(char *str) {
  if ( str && *str) {
    for ( unsigned int in = 0; in < strlen( str); in++) {
      if ( str[ in] == '"')
	str[ in] = '\'';
    }
  }
  return str;
}

/* ---------------------------------------------------------------------------
 *
 *  Post Mesh Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post mesh file
 */

int GiD_OpenPostMeshFile(GP_CONST char * FileName, GiD_PostMode Mode )
{
  /* Binary mode not allowed in mesh file!!! */
  assert(Mode!=GiD_PostBinary);
  assert(!MeshFile);
  level_mesh = POST_UNDEFINED;
  if (MeshFile) {
    /* must be closed */
    return 1;
  }
  if (!(MeshFile = NewFile(Mode)))
    /* not enough memory */
    return 2;
  GetMeshFile();
  /* now outputMesh points to MeshFile */
  if (MeshFile->Open(FileName)) {
    /* Open failed */
    return 4;
  }
  /*
  if (MeshFile->WritePostHeader()) {
    // WritePostHeader failed
    GiD_ClosePostMeshFile();
    return 5;
  }
  */
  level_mesh = POST_S0;
  return 0;
}

/*
 *  Close the current post mesh file
 */

int GiD_ClosePostMeshFile()
{
  int fail = 1;

  assert(MeshFile);
  assert(CheckState(POST_S0, level_mesh));    
  
  if (MeshFile) {
    delete MeshFile;
    fail = CPostFile::fail;
    CPostFile::fail = 0;
    MeshFile = NULL;
    /* reset outpuMesh */
    GetMeshFile();
    level_mesh = POST_UNDEFINED;
  }
  return fail;
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile
 */

int GiD_BeginMesh(GP_CONST char * MeshName, GiD_Dimension Dim,
		  GiD_ElementType EType, int NNode)
{
  CPostFile * mesh;
  int fail = 1;

  assert(CheckState(POST_S0,level_mesh));
  /* here we sould validate EType & NNode */
  assert(ValidateConnectivity(EType,NNode));

  mesh = GetMeshFile();
  if (mesh) {
    char line[LINE_SIZE];
    char *mesh_name = change_quotes( strdup( MeshName));

    snprintf(line, LINE_SIZE-1,
	     "MESH \"%s\" dimension %d ElemType %s Nnode %d",
	     mesh_name, int(Dim), GetElementTypeName(EType), NNode);
    free( mesh_name);
    if ( !(fail = mesh->WriteString(line)) ) 
      mesh->SetConnectivity(NNode);
    level_mesh = POST_MESH_S0;
  }
  return fail;
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile. With this function you can
 *    specify a color for the mesh by its RGB components where each
 *    component take values on the interval [0,1].
 */

int GiD_BeginMeshColor(GP_CONST char * MeshName, GiD_Dimension Dim,
		       GiD_ElementType EType, int NNode,
		       double Red, double Green, double Blue)
{
  CPostFile * mesh;
  int fail ;

  mesh = GetMeshFile();
  if (!(fail=GiD_BeginMesh(MeshName, Dim, EType, NNode))) {
    char line[LINE_SIZE];
    snprintf(line, LINE_SIZE-1,
	     "# color %g %g %g", Red, Green, Blue);
    if ( (fail = mesh->WriteString(line)) ) 
      level_mesh = POST_UNDEFINED;
  }
  return fail;
}

/*
 *  End current mesh
 */

int GiD_EndMesh()
{
  /* check & change state */

  /*
  assert(CheckState(POST_MESH_ELEM,level_mesh));
  level_mesh = POST_S0;
  */

  return 0;
}

/*
 *  Start a coordinate block in the current mesh
 */


int GiD_BeginCoordinates()
{
  /* state checking */
  assert(CheckState(POST_MESH_S0, level_mesh));
  level_mesh = POST_MESH_COORD0;
  return outputMesh->BeginCoordinates();
}

/*
 *  Close the current coordinate block
 */

int GiD_EndCoordinates()
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, level_mesh));
  level_mesh = POST_MESH_COORD1;
  return outputMesh->WriteString("End Coordinates");
}

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes withing the group.
 */

int GiD_BeginMeshGroup(GP_CONST char* Name)
{
  CPostFile * mesh;
  int fail = 1;
  
  assert(CheckState(POST_S0,level_mesh));
  
  mesh = GetMeshFile();
  if (mesh) {
    char line[LINE_SIZE];
    char *name = change_quotes( strdup( Name));
    
    snprintf(line, LINE_SIZE-1,
             "Group \"%s\"", name);
    free( name);
    fail = mesh->WriteString(line);
    level_mesh = POST_S0;
  }
  return fail;
}

/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */

int GiD_EndMeshGroup()
{
  return outputMesh->WriteString("End Group");
}

/*
 *  Write a coordinate member at the current Coordinates Block 
 */

int GiD_WriteCoordinates( int id, double x, double y, double z )
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, level_mesh));
  /* keep in the same level */
  return outputMesh->WriteValues(id, 3, x, y, z);
}

int GiD_WriteCoordinates2D(int id, double x, double y)
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, level_mesh));
  /* keep in the same level */
  return outputMesh->IsBinary() ? outputMesh->WriteValues(id, 3, x, y, 0.0) :
    outputMesh->WriteValues(id, 2, x, y);
}

/*
 *  Start a elements block in the current mesh
 */

int GiD_BeginElements()
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD1, level_mesh));  
  level_mesh = POST_MESH_ELEM;
  return outputMesh->BeginElements();
}

/*
 *  Close the current elements block
 */

int GiD_EndElements()
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));    
  level_mesh = POST_S0;
  return outputMesh->WriteString("End Elements");
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh.
 *  
 */

int GiD_WriteElement( int id, int nid[] )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));
  /* keep on the same state */
  return outputMesh->WriteElement(id, outputMesh->GetConnectivity(), nid);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */

int GiD_WriteElementMat( int id, int nid[] )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));    
  /* keep on the same state */
  return outputMesh->WriteElement(id, outputMesh->GetConnectivity()+1, nid);
}

/*
 *  Write an sphere element member at the current Elements Block.
 *  An sphere element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the sphere element.
 *  
 */
int GiD_WriteSphere( int id, int nid, double r )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));    
  /* keep on the same state */
  outputMesh->WriteInteger(id, 0);
  outputMesh->WriteInteger(nid,1);
  outputMesh->WriteDouble (r,  2);
  if (outputMesh->IsBinary()) {
    outputMesh->WriteInteger(1,1);    
  }
  return 0;
}

/*
 *  Write an sphere element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  An sphere element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the sphere element.
 *
 *     mat: material identification.
 *  
 */
int GiD_WriteSphereMat( int id, int nid, double r, int mat )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));    
  /* keep on the same state */
  outputMesh->WriteInteger(id,  0);
  outputMesh->WriteInteger(nid, 1);
  outputMesh->WriteDouble (r,   1);
  outputMesh->WriteInteger(mat, 2);
  return 0;
}

/*
 *  Write a circle element member at the current Elements Block.
 *  A circle element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the circle element.      
 *
 *     nx, ny, nz : normal to the plane containing the circle.
 *  
 */
int GiD_WriteCircle( int id, int nid, double r,
                     double nx, double ny, double nz )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));    
  /* keep on the same state */
  outputMesh->WriteInteger(id,  0);
  outputMesh->WriteInteger(nid, 1);
  outputMesh->WriteDouble (r,   1);
  outputMesh->WriteDouble (nx,  1);
  outputMesh->WriteDouble (ny,  1);
  outputMesh->WriteDouble (nz,  2);
  if (outputMesh->IsBinary()) {
    outputMesh->WriteInteger(1,1);    
  }
  return 0;
}

/*
 *  Write a circle element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  A circle element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the circle element.      
 *
 *     nx, ny, nz : normal to the plane containing the circle.
 *  
 */
int GiD_WriteCircleMat( int id, int nid, double r,
                        double nx, double ny, double nz, int mat )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));    
  /* keep on the same state */
  outputMesh->WriteInteger(id,  0);
  outputMesh->WriteInteger(nid, 1);
  outputMesh->WriteDouble (r,   1);
  outputMesh->WriteDouble (nx,  1);
  outputMesh->WriteDouble (ny,  1);
  outputMesh->WriteDouble (nz,  1);
  outputMesh->WriteInteger(mat, 2);
  return 0;
}

/* ---------------------------------------------------------------------------
 *
 *  Post Result Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post result file
 */

int GiD_OpenPostResultFile(GP_CONST char * FileName, GiD_PostMode Mode )
{
  assert(ResultFile==NULL);
  level_res = POST_UNDEFINED;
  if (!MeshFile)
    level_mesh = POST_UNDEFINED;
  if (ResultFile) {
    /* must be closed */
    return 1;
  }
  if ( !(ResultFile = NewFile(Mode)) )
    /* not enough memory */
    return 2;
  if (ResultFile->Open(FileName)) {
    /* could not open file */
    return 4;
  }
  if (ResultFile->WritePostHeader()) {
    /* WritePostHeader failed */
    GiD_ClosePostResultFile();
    return 5;
  }
  level_res = POST_S0;
  if (!MeshFile)
    level_mesh = POST_S0;
  buffer_values.SetFile(ResultFile);
  return 0;
}

/*
 *  Close the current post result file
 */

int GiD_ClosePostResultFile()
{
  int fail = 1;

  assert(ResultFile!=NULL);
  assert(CheckState(POST_S0, level_res));    

  if (ResultFile) {
    delete ResultFile;
    fail = CPostFile::fail;
    CPostFile::fail = 0;
    ResultFile = NULL;
    /* reset outputMesh pointer */
    GetMeshFile();
  }
  level_res = POST_UNDEFINED;
  if (!MeshFile)
    level_mesh = POST_UNDEFINED;
  return fail;
}

/*
 *  Begin Gauss Points definition
 */

int GiD_BeginGaussPoint(GP_CONST char * name, GiD_ElementType EType, GP_CONST char * MeshName,
			 int GP_number, int NodesIncluded, int InternalCoord )
{
  /* check state & validation */
  assert(CheckState(POST_S0, level_res));    

  char line[LINE_SIZE];

  char *gp_name = change_quotes( strdup( name));

  snprintf( line, LINE_SIZE-1,
	    "GaussPoints \"%s\" ElemType %s", gp_name, GetElementTypeName(EType));
  if ( MeshName && *MeshName ) {
    char *mesh_name = change_quotes( strdup( MeshName));
    strcat(line, " \"");
    strcat(line, mesh_name);
    strcat(line, "\"");
    free( mesh_name);
  }
  free( gp_name); 
  gp_name = NULL;
  if ( ResultFile->WriteString(line) )
    return 1;
  snprintf(line, LINE_SIZE, "Number Of Gauss Points: %d", GP_number);
  if ( ResultFile->WriteString(line) )
    return 1;
  /* here we could save the number of GP in order to check at
     EndGaussPoint */
  GP_number_check = GP_number;
  if ( EType == GiD_Linear ) {
    if ( NodesIncluded ) {
      if ( ResultFile->WriteString("  Nodes included") )
	return 1;
    } else
      if ( ResultFile->WriteString("  Nodes not included") )
	return 1;
  }
  if (InternalCoord) {
    if (ResultFile->WriteString("Natural Coordinates: Internal"))
      return 1;
    level_res = POST_GAUSS_S0;
  } else {
    if (ResultFile->WriteString("Natural Coordinates: Given"))
      return 1;
    level_res = POST_GAUSS_GIVEN;    
    /* here we can save the size of the coordinates to check later
       in WriteGaussPointXX*/
  }
  return 0;
}

/*
 *  End current Gauss Points definition
 */

static int CheckGaussPointEnd()
{
  if (level_res != POST_GAUSS_S0 && level_res != POST_GAUSS_GIVEN) {
    printf("Invalid call of GiD_EndGaussPoint. Current state is '%s' and should be '%s' or '%s'\n",
	   GetStateDesc(level_res), GetStateDesc(POST_GAUSS_S0), GetStateDesc(POST_GAUSS_GIVEN));
    return 0;
  }
  return 1;
}

static int CheckGaussPointGiven()
{
  if (gauss_written !=  GP_number_check) {
    printf("missmatch in gauss point given, written %d and %d were requiered",
	   gauss_written, GP_number_check);
    return 0;
  }
  return 1;
}

int GiD_EndGaussPoint()
{
  /* check state */
  assert(CheckGaussPointEnd());
#ifndef NDEBUG
  if (level_res == POST_GAUSS_GIVEN)
    assert(CheckGaussPointGiven());
#endif
  level_res = POST_S0;
  GP_number_check = gauss_written = 0;
  return ResultFile->WriteString("End GaussPoints");
}

/*
 *  Write internal gauss point coordinate.
 */

int GiD_WriteGaussPoint2D( double x, double y )
{
  /* check state */
  assert(CheckState(POST_GAUSS_GIVEN, level_res));    
  if (!ResultFile->Write2D(x, y)) {
    ++gauss_written;
    return 0;
  }
  return 1;
}

int GiD_WriteGaussPoint3D( double x, double y, double z )
{
  /* check state */
  assert(CheckState(POST_GAUSS_GIVEN, level_res));
  if (!ResultFile->Write3D(x, y, z)) {
    ++gauss_written;
    return 0;    
  }
  return 1;
}

/*
 *  Begin a Range Table definition
 */

int GiD_BeginRangeTable(GP_CONST char * name )
{
  char line[LINE_SIZE];
  /* check & update state */
  assert(CheckState(POST_S0, level_res));
  

  char *rt_name = change_quotes( strdup( name));
  snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", rt_name);
  free( rt_name);
  if (!ResultFile->WriteString(line)) {
    level_res = POST_RANGE_S0;
    return 0;
  }
  return 1;
}

/*
 *  End a Range Table definition
 */

int GiD_EndRangeTable()
{
  /* check & update state */
  assert(CheckState(POST_RANGE_S0, level_res));

  if (!ResultFile->WriteString("End ResultRangesTable")) {
    level_res = POST_S0;
    return 0;
  }
  return 1;
}

/*
 *  Write Range functions --
 *
 *   WriteMinRange : write a range with an implicit minimum value, the
 *   minimum absolute in the result set.
 *
 *   WriteRange : write an explicit range.
 *
 *   WritemaxRange: write a range with an implicit maximum value, the
 *   maximum absolute in the result set.
 */

int GiD_WriteMinRange(double max, GP_CONST char * name)
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckState(POST_RANGE_S0, level_res));

  char *tmp_name = change_quotes( strdup( name));
  snprintf(line, LINE_SIZE-1, " - %g : \"%s\"", max, tmp_name);
  free( tmp_name);
  return ResultFile->WriteString(line);
}

int GiD_WriteRange(double min, double max, GP_CONST char * name)
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckState(POST_RANGE_S0, level_res));
  
  char *tmp_name = change_quotes( strdup( name));
  snprintf(line, LINE_SIZE-1, " %g - %g : \"%s\"", min, max, tmp_name);
  free( tmp_name);
  return ResultFile->WriteString(line);
}

int GiD_WriteMaxRange( double min, GP_CONST char * name )
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckState(POST_RANGE_S0, level_res));
  
  char *tmp_name = change_quotes( strdup( name));
  snprintf(line, LINE_SIZE-1, "%g - : \"%s\"", min, tmp_name);
  free( tmp_name);
  return ResultFile->WriteString(line);
}

/*
 *  Begin Result Block
 */

int GiD_BeginResult( GP_CONST char * Result, GP_CONST char * Analysis, double step,
		     GiD_ResultType Type, GiD_ResultLocation Where,
		     GP_CONST char * GaussPointsName, GP_CONST char * RangeTable, 
		     int compc, GP_CONST char * compv[])
{
  char line[LINE_SIZE];
  const char * loc;
  int i;

  /* check & change state */
  assert(CheckState(POST_S0, level_res));

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  char *res_name = change_quotes( strdup( Result));
  char *analysis_name = change_quotes( strdup( Analysis));
  snprintf(line, LINE_SIZE-1, "Result \"%s\" \"%s\" %g %s %s",
	   res_name, analysis_name, step, GetResultTypeName(Type), loc);
  free( res_name);
  free( analysis_name);
  if (Where == GiD_OnGaussPoints) {
    assert(GaussPointsName);
    char *tmp_name = change_quotes( strdup( GaussPointsName));
    strcat(line, " \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
  }
  if (ResultFile->WriteString(line))
    return 1;
  if (RangeTable) {
    assert(RangeTable);
    char *tmp_name = change_quotes( strdup( RangeTable));
    snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", tmp_name);
    free( tmp_name);
    if  (ResultFile->WriteString(line))
      return 1;
  }
  if ( compc > 0 ) {
    snprintf(line, LINE_SIZE-1, "ComponentNames");
    for ( i = 0; i < compc; i++ ) {  
      char *tmp_name = change_quotes( strdup( compv[ i]));
      strcat(line, " ");
      strcat(line, "\"");
      strcat(line, tmp_name);
      strcat(line, "\"");
      free( tmp_name);
    }
    if (ResultFile->WriteString(line))
      return 1;
  }
  flag_isgroup = 0;
  flag_begin_values = 1;
  level_res = POST_RESULT_VALUES;
  return ResultFile->BeginValues();
}

int GiD_BeginResultHeader( GP_CONST char * Result, GP_CONST char * Analysis, double step,
			   GiD_ResultType Type, GiD_ResultLocation Where,
			   GP_CONST char * GaussPointsName)
{
  char line[LINE_SIZE];
  const char * loc;

  /* check & change state */
  assert(CheckState(POST_S0, level_res));
  assert(Result);
  assert(Analysis);

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  char *res_name = change_quotes( strdup( Result));
  char *analysis_name = change_quotes( strdup( Analysis));
  snprintf(line, LINE_SIZE-1, "Result \"%s\" \"%s\" %g %s %s",
	   res_name, analysis_name, step, GetResultTypeName(Type), loc);
  free( res_name);
  free( analysis_name);
  if (Where == GiD_OnGaussPoints) {
    assert(GaussPointsName);
    char *tmp_name = change_quotes( strdup( GaussPointsName));
    strcat(line, " \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
  }
  if (ResultFile->WriteString(line)) {
    /* could not write result header */
    return 1;
  }
  level_res = POST_RESULT_SINGLE;
  flag_isgroup = 0;
  flag_begin_values = 0;
  return 0;
}

static int CheckResultHeaderState()
{
  if (level_res == POST_RESULT_SINGLE || level_res == POST_RESULT_DESC)
    return 1;
  printf("Invalid result state '%s'. Should be: '%s' or '%s'\n",
	 GetStateDesc(level_res),
	 GetStateDesc(POST_RESULT_SINGLE),
	 GetStateDesc(POST_RESULT_DESC));
  return 0;
}

int GiD_ResultRange(GP_CONST char * RangeTable)
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckResultHeaderState());
  assert(RangeTable);
  
  if (RangeTable) {
    char *tmp_name = change_quotes( strdup( RangeTable));
    snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", tmp_name);
    free( tmp_name);
    return ResultFile->WriteString(line);
  }
  return 1;
}

int GiD_ResultComponents(int compc, GP_CONST char * compv[])
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckResultHeaderState());
  assert(compc>0);
  
  if (compc > 0) {
    snprintf(line, LINE_SIZE-1, "ComponentNames");
    for (int i = 0; i < compc; i++) {
      char *tmp_name = change_quotes( strdup( compv[ i]));
      strcat(line, " ");
      strcat(line, "\"");
      strcat(line, tmp_name);
      strcat(line, "\"");
      free( tmp_name);
    }
    return ResultFile->WriteString(line);
  }
  return 1;
}

int GiD_ResultUnit(GP_CONST char * UnitName)
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckResultHeaderState());
  assert(UnitName);

  if (UnitName) {
    char *tmp_name = change_quotes( strdup( UnitName));
    snprintf(line, LINE_SIZE-1, "UnitName \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
    return ResultFile->WriteString(line);
  }
  return 1;
}

int GiD_BeginResultGroup(GP_CONST char * Analysis, double step, GiD_ResultLocation Where,
			 GP_CONST char * GaussPointsName)
{
  char line[LINE_SIZE];
  const char * loc;

  /* check & change state */
  assert(CheckState(POST_S0, level_res));
  assert(Analysis);

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  char *analysis_name = change_quotes( strdup( Analysis));
  snprintf(line, LINE_SIZE-1, "ResultGroup \"%s\" %g %s", analysis_name, step, loc);
  free( analysis_name);
  if (Where == GiD_OnGaussPoints) {
    assert(GaussPointsName);
    char *tmp_name = change_quotes( strdup( GaussPointsName));
    strcat(line, " \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
  }
  if (ResultFile->WriteString(line)) {
    /* could not write result header */
    return 1;
  }
  level_res = POST_RESULT_GROUP;
  /* initialize values counter */
  flag_isgroup = 1;
  flag_begin_values = 0;
  buffer_values.OnBeginResultGroup();
  return 0;
}

static int CheckStateDesc()
{
  if (level_res == POST_RESULT_GROUP || level_res == POST_RESULT_DESC)
    return 1;
  printf("Invalid result state '%s'. Should be: '%s' or '%s'\n",
	 GetStateDesc(level_res),
	 GetStateDesc(POST_RESULT_GROUP),
	 GetStateDesc(POST_RESULT_DESC));
  return 0;
}

static
int _GiD_ResultDescription_(GP_CONST char * Result, GiD_ResultType Type, size_t s = 0)
{
  char line[LINE_SIZE];

  /* check & change state */
  assert(CheckStateDesc());
  assert(Result);

  line[0] = '\0';
  char *tmp_name = change_quotes( strdup( Result));
  snprintf(line, LINE_SIZE-1, "ResultDescription \"%s\" %s", tmp_name, GetResultTypeName(Type,s));
  free( tmp_name);
  if (ResultFile->WriteString(line)) {
    /* could not write result description */
    return 1;
  }
  level_res = POST_RESULT_DESC;
  /* update number of values to check per location */
  buffer_values.OnNewType(Type);
  return 0;
}

int GiD_ResultDescription( GP_CONST char * Result, GiD_ResultType Type)
{
  return _GiD_ResultDescription_(Result, Type);
}

int GiD_ResultDescriptionDim( GP_CONST char * Result, GiD_ResultType Type, size_t s)
{
  return _GiD_ResultDescription_(Result, Type, s);
}

int GiD_ResultValues()
{
  return 0;
}

/*
 *  End Result Block
 */

int GiD_EndResult()
{
  int _fail = 0;

  /* check & change state */
  assert(CheckState(POST_RESULT_VALUES, level_res));
  if (flag_isgroup) {
    //buffer_values.FlushValues();
    assert(buffer_values.IsEmpty());
  }

  _fail = ResultFile->EndValues();
  ResultFile->ResetLastID();
  level_res = POST_S0;
  return _fail;
}

/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */

int GiD_BeginOnMeshGroup(char * Name)
{
  char line[LINE_SIZE];

  snprintf(line, LINE_SIZE-1, "OnGroup \"%s\"", Name);
  return ResultFile->WriteString(line);
}

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */

int GiD_EndOnMeshGroup()
{
  return ResultFile->WriteString("End OnGroup");
}

/*
 * Flushes all pending output into the compressed file.
 */

int GiD_FlushPostFile()
{
  return ResultFile->Flush();
}

/*
 *  Write result functions
 */

inline int GiD_EnsureBeginValues()
{
  if (!flag_begin_values) {
    if (!ResultFile->BeginValues()) {
      level_res = POST_RESULT_VALUES;
      if (flag_isgroup) {
        buffer_values.OnBeginValues();
      }
      flag_begin_values = 1;
      return 0;
    }
  }
  return 1;
}

int GiD_WriteScalar( int id, double v )
{
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  return flag_isgroup ?  buffer_values.WriteValues(GiD_Scalar, id, 1, v) :
    ResultFile->WriteValues(id, 1, v);
}

int GiD_Write2DVector( int id, double x, double y)
{
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  if (flag_isgroup || !ResultFile->IsBinary()) {
    return flag_isgroup ? buffer_values.WriteValues(GiD_Vector, id, 2, x, y) :
      ResultFile->WriteValues(id, 2, x, y);
  } else {
    /* single result & binary */
    double mod = sqrt(x*x + y*y);
  
    return ResultFile->WriteValues(id, 4, x, y, 0.0, mod);
  }  
}

int GiD_WriteVector( int id, double x, double y, double z )
{
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  if (flag_isgroup || !ResultFile->IsBinary()) {
    return flag_isgroup ? buffer_values.WriteValues(GiD_Vector, id, 3, x, y, z) :
      ResultFile->WriteValues(id, 3, x, y, z);
  } else {
    /* single result & binary */
    double mod = sqrt(x*x + y*y + z*z);
  
    return ResultFile->WriteValues(id, 4, x, y, z, mod);
  }  
}

int GiD_WriteVectorModule( int id, double x, double y, double z, double mod )
{
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  /* 4-vectors can not be written on RG-ASCII */
  return flag_isgroup ? buffer_values.WriteValues(GiD_Vector, id, 4, x, y, z, mod) :
    ResultFile->WriteValues(id, 4, x, y, z, mod );
}

int GiD_Write2DMatrix( int id, double Sxx, double Syy, double Sxy )
{
  if (ResultFile->IsBinary()) {
    return GiD_Write3DMatrix(id,
                             Sxx, Syy,         0.0 /*Szz*/,
                             Sxy, 0.0 /*Syz*/, 0.0 /*Sxz*/);
  }
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));
  return flag_isgroup ?
    buffer_values.WriteValues(GiD_Matrix, id, 3, Sxx, Syy, Sxy) :
    ResultFile->WriteValues(              id, 3, Sxx, Syy, Sxy);
}

int GiD_Write3DMatrix( int id, double Sxx, double Syy, double Szz,
		   double Sxy, double Syz, double Sxz )
{
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));
  
  return flag_isgroup ?
    buffer_values.WriteValues(GiD_Matrix, id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz) :
    ResultFile->WriteValues(id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz);
}

int GiD_WritePlainDefMatrix( int id, double Sxx, double Syy, double Sxy, double Szz )
{
  if (ResultFile->IsBinary()) {
    return GiD_Write3DMatrix(id,
                             Sxx, Syy,         Szz,
                             Sxy, 0.0 /*Syz*/, 0.0 /*Sxz*/);
  }
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));
   
  return flag_isgroup ? 
    buffer_values.WriteValues(GiD_PlainDeformationMatrix, id, 4, Sxx, Syy, Sxy, Szz) :
    ResultFile->WriteValues(id, 4, Sxx, Syy, Sxy, Szz);
}

int GiD_WriteMainMatrix( int id,
		     double Si, double Sii, double Siii,
		     double Vix, double Viy, double Viz,
		     double Viix, double Viiy, double Viiz,
		     double Viiix, double Viiiy, double Viiiz )
{
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));
  
  return flag_isgroup ? 
    buffer_values.WriteValues(GiD_MainMatrix, id, 12, Si, Sii, Siii,
			      Vix, Viy, Viz,
			      Viix, Viiy, Viiz,
			      Viiix, Viiiy, Viiiz) :
    ResultFile->WriteValues(id, 12, Si, Sii, Siii,
			    Vix, Viy, Viz,
			    Viix, Viiy, Viiz,
			    Viiix, Viiiy, Viiiz);
}

int GiD_WriteLocalAxes( int id, double euler_1, double euler_2, double euler_3 )
{
  GiD_EnsureBeginValues();
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));
  
  return flag_isgroup ?
    buffer_values.WriteValues(GiD_LocalAxes, id, 3, euler_1, euler_2, euler_3) :
    ResultFile->WriteValues(id, 3, euler_1, euler_2, euler_3);
}
