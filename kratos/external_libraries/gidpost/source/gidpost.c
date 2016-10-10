/* gidpost 2.0 */
/*
 *  gidpost.c--
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
#include <ctype.h>

#include "gidpostInt.h"
#include "gidpost.h"
#include "gidpostHash.h"
  
#ifdef HDF5
  #include "gidpostHDF5.h"
#endif


#ifdef WIN32
#define snprintf _snprintf
#endif
  
/* let's change the double quotes to simple ones, to avoid potential problems*/
static char *change_quotes(char *str) {
  unsigned int in;
  
  if ( str && *str) {
    for ( in = 0; in < strlen( str); in++) {
      if ( str[ in] == '"')
	str[ in] = '\'';
    }
  }
  return str;
}

/* to check if a string must be enclosed in quotes*/
int string_hasspace(const char *str){
  unsigned char *tstr;
  /*
    assert(str);
  */
  if ( !str || !*str) return 0;
  tstr=(unsigned char*)str;
  while(*tstr){
    if(isspace(*tstr))
      return 1;
    tstr++;
  }  
  return 0;
}


/* ---------------------------------------------------------------------------
 *
 *  Global files
 *
 * ---------------------------------------------------------------------------
 */

static CPostFile *MeshFile = NULL;
static CPostFile *ResultFile = NULL;

static  CPostFile *outputMesh = NULL;

/*
  hay que encapsular HDF5 en las estructuras internas de GiDPost
*/
GiD_PostMode PostMode;

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

#define FD2FILE(_fd_, _File_)                     \
  do {                                            \
    assert(_fd_);                                 \
    _File_ = (CPostFile*)GiD_HashFind(_fd_);      \
    if (!_File_) {                                \
      /* file handler does not exists */          \
      return -8;                                  \
    }                                             \
  } while (0);

CPostFile *GetMeshFile()
{
  return (outputMesh = MeshFile ? MeshFile : ResultFile);
}

CPostFile * NewFile(GiD_PostMode Mode)
{
  GiD_PostInit();
  switch ( Mode ) {
    case GiD_PostAscii:
      return CPostAscii_Create();
    case GiD_PostAsciiZipped:
      return CPostAsciiZ_Create();      
    case GiD_PostBinary:
      return CPostBinary_Create();
    case GiD_PostHDF5:
      return NULL;
    default:
      return NULL;
  }
}

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
  "Circle",
  "Point"
};

GP_CONST char * GetElementTypeName( GiD_ElementType type )
{
  return strElementType[(int)(type)];
}

static
int CheckState(post_state s_req, post_state s_cur)
{
  if (s_req!=s_cur) {
    printf("invalid state '%s' should be '%s'\n", GetStateDesc(s_cur), GetStateDesc(s_req));
    return 0;
  }
  return 1;
}

static
int ValidateConnectivity(GiD_ElementType etype , int NNode)
{
  int error;

  switch (etype) {
  case GiD_Point:
  case GiD_Sphere:
  case GiD_Circle:
  case GiD_Cluster:
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


int GiD_PostInit()
{
  return GiD_HashInit();
}

int GiD_PostDone()
{
  return GiD_HashDone();
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
  PostMode=Mode;

  if(PostMode==GiD_PostHDF5){
#ifdef HDF5
      return GiD_OpenPostMeshFile_HDF5(FileName);
#else
      return 5;
#endif
    }
  /* Binary mode not allowed in mesh file!!! */
  assert(Mode!=GiD_PostBinary);
  assert(!MeshFile);
  if (MeshFile) {
    /* must be closed */
    return 1;
  }
  if (!(MeshFile = NewFile(Mode)))
    /* not enough memory */
    return 2;
  GetMeshFile();
  /* now outputMesh points to MeshFile */
  if (CPostFile_Open(MeshFile, FileName)) {
    /* Open failed */
    CPostFile_Release(MeshFile);
    MeshFile = NULL;
    return 4;
  }
  MeshFile->level_mesh = POST_S0;
  return 0;
}

GiD_FILE GiD_fOpenPostMeshFile(GP_CONST char * FileName,
		               GiD_PostMode Mode)
{
  CPostFile *File = NULL;
  GiD_FILE fd;
  
  /* Binary mode not allowed in mesh file!!! */
  assert(Mode!=GiD_PostBinary);

  if (!(File=NewFile(Mode)))
    /* not enough memory = -2 */
    return 0;
  /* now open the File */
  if (CPostFile_Open(File, FileName)) {
    /* Open failed = -4 */
    return 0;
  }
  fd = GiD_HashAdd(File);
  if (!fd) {
    /* could not create a file handler = -6 */
    CPostFile_Release(File);
    return 0;
  }
  return fd;
}

/*
 *  Close the current post mesh file
 */
int GiD_ClosePostMeshFile()
{
  int fail = 1;

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ClosePostMeshFile_HDF5();
  }
#endif
  assert(MeshFile);
  assert(CheckState(POST_S0, MeshFile->level_mesh));    
  
  if (MeshFile) {
    fail = CPostFile_Release(MeshFile);
    MeshFile = NULL;
    /* reset outpuMesh */
    GetMeshFile();
    /* level_mesh = POST_UNDEFINED; */
  }
  return fail;
}

int GiD_fClosePostMeshFile(GiD_FILE fd)
{
  int fail = 1;
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  fail = CPostFile_Release(File);
  GiD_HashRemove(fd);
  
  return fail;
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile
 */

int _GiD_BeginMesh(CPostFile *File,
		   GP_CONST char * MeshName, GiD_Dimension Dim,
		   GiD_ElementType EType, int NNode)
{
  int fail = 1;
  char line[LINE_SIZE];
  char *mesh_name;
  
  /* here we sould validate EType & NNode */
  assert(ValidateConnectivity(EType,NNode));
  assert(File);
  /*
    assert(CheckState(POST_S0,File->level_mesh));
  */

  mesh_name = change_quotes(strdup(MeshName));
  snprintf(line, LINE_SIZE-1,
	   "MESH \"%s\" dimension %d ElemType %s Nnode %d",
	   mesh_name, (int)(Dim), GetElementTypeName(EType), NNode);
  free(mesh_name);
  if (!(fail = CPostFile_WriteString(File,line))) 
    CPostFile_SetConnectivity(File, NNode);
  File->level_mesh = POST_MESH_S0;
  return fail;
}

int GiD_BeginMesh(GP_CONST char * MeshName, GiD_Dimension Dim,
		  GiD_ElementType EType, int NNode)
{
  CPostFile *mesh;
  
  #ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginMesh_HDF5(MeshName, Dim,EType,NNode);
  }
#endif

  mesh = GetMeshFile();
  return _GiD_BeginMesh(mesh, MeshName, Dim, EType, NNode);
  
}

int GiD_fBeginMesh(GiD_FILE fd, GP_CONST char * MeshName,
		   GiD_Dimension Dim, GiD_ElementType EType, int NNode)
{
  CPostFile *mesh;

  FD2FILE(fd,mesh);
  return _GiD_BeginMesh(mesh, MeshName, Dim, EType, NNode);
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile. With this function you can
 *    specify a color for the mesh by its RGB components where each
 *    component take values on the interval [0,1].
 */

int _GiD_BeginMeshColor(CPostFile *File,
		        GP_CONST char * MeshName, GiD_Dimension Dim,
		        GiD_ElementType EType, int NNode,
		        double Red, double Green, double Blue)
{
  int fail ;
  char line[LINE_SIZE];
  
  if (!(fail=_GiD_BeginMesh(File, MeshName, Dim, EType, NNode))) {
    snprintf(line, LINE_SIZE-1, "# color %g %g %g", Red, Green, Blue);
    if ( (fail = CPostFile_WriteString(File, line)) ) 
      File->level_mesh = POST_UNDEFINED;
  }
  return fail;
}

int GiD_BeginMeshColor(GP_CONST char * MeshName, GiD_Dimension Dim,
		       GiD_ElementType EType, int NNode,
		       double Red, double Green, double Blue)
{
  CPostFile * mesh;

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginMeshColor_HDF5(MeshName, Dim,EType,NNode,Red,Green, Blue);
  }
#endif

  mesh = GetMeshFile();
  return _GiD_BeginMeshColor(mesh, MeshName, Dim, EType, NNode,
		             Red, Green, Blue);
}

int GiD_fBeginMeshColor(GiD_FILE fd, GP_CONST char * MeshName,
		        GiD_Dimension Dim, GiD_ElementType EType,
		        int NNode,
		        double Red, double Green, double Blue)
{
  CPostFile * mesh;

  FD2FILE(fd,mesh);
    
  return _GiD_BeginMeshColor(mesh, MeshName, Dim, EType, NNode,
		             Red, Green, Blue);
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
  
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndMesh_HDF5();
  }
#endif

  return 0;
}

int GiD_fEndMesh(GiD_FILE fd)
{
  return 0;
}

int _GiD_MeshUnit(CPostFile *File,GP_CONST char * UnitName)
{
  if (UnitName) {
    char line[LINE_SIZE];
    char *tmp_name;
    tmp_name = change_quotes( strdup( UnitName));
    snprintf(line, LINE_SIZE-1, "Unit \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
    return CPostFile_WriteString(File, line);
  }  
  return 1;
}

int GiD_MeshUnit(GP_CONST char * UnitName)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return  GiD_MeshUnit_HDF5(UnitName);
  }
#endif
  return _GiD_MeshUnit(outputMesh,UnitName);
  return 1;
}


int GiD_fMeshUnit(GiD_FILE fd,GP_CONST char * UnitName)
{
  CPostFile *File = NULL;
  FD2FILE(fd,File);
  return _GiD_MeshUnit(File,UnitName);
}

int GiD_MeshLocalAxes(GP_CONST char * Result, GP_CONST char * Analysis,double step)
{
#ifdef HDF5
    if(PostMode==GiD_PostHDF5){
    return GiD_MeshLocalAxes_HDF5(Result,Analysis,step);
  }
#endif
  return 1;
}

/*
 *  Start a coordinate block in the current mesh
 */


int _GiD_BeginCoordinates(CPostFile *File)
{
  /* state checking */
  assert(File);
  assert(CheckState(POST_MESH_S0, File->level_mesh));
  File->level_mesh = POST_MESH_COORD0;
  return CPostFile_BeginCoordinates(File);
}

int GiD_BeginCoordinates()
{

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginCoordinates_HDF5();
  }
#endif

  return _GiD_BeginCoordinates(outputMesh);
}

int GiD_fBeginCoordinates(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  return _GiD_BeginCoordinates(File);
}

/*
 *  Close the current coordinate block
 */

int _GiD_EndCoordinates(CPostFile *File)
{

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndCoordinates_HDF5();
  }
#endif

  /* state checking */
  assert(CheckState(POST_MESH_COORD0, File->level_mesh));
  File->level_mesh = POST_MESH_COORD1;
  CPostFile_ResetLastID(File); /* m_lastID is also checked when writting coordinates */
  return CPostFile_WriteString(File, "End Coordinates");
}

int GiD_EndCoordinates()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndCoordinates_HDF5();
  }
#endif

  return _GiD_EndCoordinates(outputMesh);
}

int GiD_fEndCoordinates(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_EndCoordinates(File);
}

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes withing the group.
 */

int _GiD_BeginMeshGroup( CPostFile *File, GP_CONST char* Name )
{
  int fail = 1;
  char line[LINE_SIZE];
  char *name;
    
  assert(File);
  assert(CheckState(POST_S0,File->level_mesh));
  
  name = change_quotes(strdup(Name));

  snprintf(line, LINE_SIZE-1, "Group \"%s\"", name);
  free(name);
  fail = CPostFile_WriteString(File, line);
  File->level_mesh = POST_S0;
  return fail;
}

int GiD_BeginMeshGroup( GP_CONST char* Name )
{
  CPostFile * Mesh;
  
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return -1;
  }
#endif

  Mesh = GetMeshFile();
  return _GiD_BeginMeshGroup(Mesh, Name);
}

int GiD_fBeginMeshGroup(GiD_FILE fd, GP_CONST char* Name)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_BeginMeshGroup(File, Name);
}

/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */

int GiD_EndMeshGroup()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return -1;
  }
#endif

  return CPostFile_WriteString(outputMesh, "End Group");
}

int GiD_fEndMeshGroup(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return CPostFile_WriteString(File, "End Group");
}

/*
 *  Write a coordinate member at the current Coordinates Block 
 */

int GiD_WriteCoordinates( int id, double x, double y, double z )
{
  int res = 0;
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    res = GiD_WriteCoordinates_HDF5(id,x,y,z);
    return res;
  }
#endif

  /* state checking */
  assert(CheckState(POST_MESH_COORD0, outputMesh->level_mesh));
  /* keep in the same level */
  res = CPostFile_WriteValuesVA(outputMesh, id, 3, x, y, z);
  CPostFile_ResetLastID( outputMesh); /* inside _WriteValuesVA there is some trick to write gaussian  results */
  return res;
}

int GiD_fWriteCoordinates(GiD_FILE fd, int id, double x, double y, double z)
{
  CPostFile *File = NULL;
  int res = 0;

  FD2FILE(fd,File);
  
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, File->level_mesh));
  /* keep in the same level */
  res = CPostFile_WriteValuesVA(File, id, 3, x, y, z);
  CPostFile_ResetLastID(File); /* inside _WriteValuesVA there is some trick to write gaussian  results */
  return res;
}

int _GiD_WriteCoordinates2D(CPostFile *File, int id, double x, double y)
{
  int res = 0;
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, File->level_mesh));
  /* keep in the same level */
  res = CPostFile_IsBinary(File)
    ? CPostFile_WriteValuesVA(File,id, 3, x, y, 0.0)
    : CPostFile_WriteValuesVA(File, id, 2, x, y);
  CPostFile_ResetLastID(File); /* inside _WriteValuesVA there is some trick to write gaussian  results */
  return res;
}

int GiD_WriteCoordinates2D(int id, double x, double y)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    int res = GiD_WriteCoordinates2D_HDF5(id,x,y);
    return res;
  }
#endif

  return _GiD_WriteCoordinates2D(outputMesh, id, x, y);
}

int GiD_fWriteCoordinates2D(GiD_FILE fd, int id, double x, double y)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WriteCoordinates2D(File, id, x, y);
}

/*
 *  Start a elements block in the current mesh
 */

int _GiD_BeginElements(CPostFile *File)
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD1, File->level_mesh));  
  File->level_mesh = POST_MESH_ELEM;
  return CPostFile_BeginElements(File);
}
int GiD_BeginElements()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginElements_HDF5();
  }
#endif

  return _GiD_BeginElements(outputMesh);
}

int GiD_fBeginElements(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  return _GiD_BeginElements(File);
}

/*
 *  Close the current elements block
 */

int _GiD_EndElements(CPostFile *File)
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  File->level_mesh = POST_S0;
  return CPostFile_WriteString(File, "End Elements");
}

int GiD_EndElements()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndElements_HDF5();
  }
#endif

  return _GiD_EndElements(outputMesh);
}

int GiD_fEndElements(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
  FD2FILE(fd,File);
  
  return _GiD_EndElements(File);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh.
 *  
 */

int _GiD_WriteElement(CPostFile *File, int id, int nid[])
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));
  /* keep on the same state */
  return CPostFile_WriteElement(File,
		                id, CPostFile_GetConnectivity(File), nid);
}
int GiD_WriteElement( int id, int nid[] )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteElement_HDF5(id,nid);
  }
#endif

  return _GiD_WriteElement(outputMesh, id, nid);
}

int GiD_fWriteElement(GiD_FILE fd, int id, int nid[])
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteElement(File, id, nid);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */

int _GiD_WriteElementMat(CPostFile *File, int id, int nid[])
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  /* keep on the same state */
  return CPostFile_WriteElement(File,
		                id, CPostFile_GetConnectivity(File)+1,
		                nid);
}

int GiD_WriteElementMat(int id, int nid[])
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteElementMat_HDF5(id,nid);
  }
#endif

  return _GiD_WriteElementMat(outputMesh, id, nid);
}

int GiD_fWriteElementMat(GiD_FILE fd, int id, int nid[])
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteElementMat(File, id, nid);
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

int _GiD_WriteSphere(CPostFile *File, int id, int nid, double r)
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id, 0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,  2);
  if (CPostFile_IsBinary(File)) {
    CPostFile_WriteInteger(File, 1, 1);    
  }
  return 0;
}

int GiD_WriteSphere(int id, int nid, double r)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteSphere_HDF5(id,nid,r);
  }
#endif

  return _GiD_WriteSphere(outputMesh, id, nid, r);
}

int GiD_fWriteSphere(GiD_FILE fd, int id, int nid, double r)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteSphere(File, id, nid, r);
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

int _GiD_WriteSphereMat(CPostFile * File, int id, int nid, double r, int mat)
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,   1);
  CPostFile_WriteInteger(File, mat, 2);
  return 0;
}

int GiD_WriteSphereMat(int id, int nid, double r, int mat)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteSphereMat_HDF5(id,nid,r,mat);
  }
#endif

  return _GiD_WriteSphereMat(outputMesh, id, nid, r, mat);
}

int GiD_fWriteSphereMat(GiD_FILE fd, int id, int nid, double r, int mat)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteSphereMat(File, id, nid, r, mat);
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

int _GiD_WriteCircle(CPostFile *File,
		     int id, int nid, double r,
		     double nx, double ny, double nz)
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,   1);
  CPostFile_WriteDouble (File, nx,  1);
  CPostFile_WriteDouble (File, ny,  1);
  CPostFile_WriteDouble (File, nz,  2);
  if (CPostFile_IsBinary(File)) {
    CPostFile_WriteInteger(File, 1, 1);    
  }
  return 0;
}

int GiD_WriteCircle(int id, int nid, double r,
		    double nx, double ny, double nz)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteCircle_HDF5(id,nid,r,nx,ny,nz);
  }
#endif

  return _GiD_WriteCircle(outputMesh, id, nid, r, nx, ny, nz);
}

int GiD_fWriteCircle(GiD_FILE fd, int id, int nid, double r,
		     double nx, double ny, double nz)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteCircle(File, id, nid, r, nx, ny, nz);
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

int _GiD_WriteCircleMat(CPostFile *File, int id, int nid, double r,
		        double nx, double ny, double nz, int mat)
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,   1);
  CPostFile_WriteDouble (File, nx,  1);
  CPostFile_WriteDouble (File, ny,  1);
  CPostFile_WriteDouble (File, nz,  1);
  CPostFile_WriteInteger(File, mat, 2);
  return 0;
}

int GiD_WriteCircleMat(int id, int nid, double r,
		       double nx, double ny, double nz, int mat)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteCircleMat_HDF5(id,nid,r,nx,ny,nz,mat);
  }
#endif
  
  return _GiD_WriteCircleMat(outputMesh, id, nid, r, nx, ny, nz, mat);
}

int GiD_fWriteCircleMat(GiD_FILE fd, int id, int nid, double r,
		        double nx, double ny, double nz, int mat)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteCircleMat(File, id, nid, r, nx, ny, nz, mat);
}

/*
 *  Write a cluster element member at the current Elements Block.
 *  A cluster element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *  
 */

int _GiD_WriteCluster(CPostFile *File, int id, int nid)
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id, 0);
  CPostFile_WriteInteger(File, nid, 1);
  if (CPostFile_IsBinary(File)) {
    CPostFile_WriteInteger(File, 1, 1);    
  }
  return 0;
}

int GiD_WriteCluster(int id, int nid)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteCluster_HDF5(id,nid);
  }
#endif

  return _GiD_WriteCluster(outputMesh, id, nid);
}

int GiD_fWriteCluster(GiD_FILE fd, int id, int nid)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteCluster(File, id, nid);
}

/*
 *  Write a cluster element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  A cluster element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     mat: material identification.
 *  
 */

int _GiD_WriteClusterMat(CPostFile * File, int id, int nid, int mat)
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, File->level_mesh));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteInteger(File, mat, 2);
  return 0;
}

int GiD_WriteClusterMat(int id, int nid, int mat)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteClusterMat_HDF5(id,nid,mat);
  }
#endif

  return _GiD_WriteClusterMat(outputMesh, id, nid, mat);
}

int GiD_fWriteClusterMat(GiD_FILE fd, int id, int nid, int mat)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteClusterMat(File, id, nid, mat);
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
  PostMode=Mode;

  if(PostMode==GiD_PostHDF5){
#ifdef HDF5
    return GiD_OpenPostResultFile_HDF5(FileName);
#else
    return 5;
#endif
  }

  assert(ResultFile==NULL);
  if (ResultFile) {
    /* must be closed */
    return 1;
  }
  if ( !(ResultFile = NewFile(Mode)) )
    /* not enough memory */
    return 2;
  if (CPostFile_Open(ResultFile, FileName)) {
    /* could not open file */
    CPostFile_Release(ResultFile);
    return 4;
  }
  if (!MeshFile)
    ResultFile->level_mesh = POST_UNDEFINED;
  if (CPostFile_WritePostHeader(ResultFile)) {
    /* WritePostHeader failed */
    GiD_ClosePostResultFile();
    return 5;
  }
  ResultFile->level_res = POST_S0;
  if (!MeshFile)
    ResultFile->level_mesh = POST_S0;
  return 0;
}

GiD_FILE GiD_fOpenPostResultFile(GP_CONST char * FileName, GiD_PostMode Mode)
{
  CPostFile *File = NULL;
  GiD_FILE fd;

  /* level_res = POST_UNDEFINED; */
  if (!(File=NewFile(Mode)))
    /* not enough memory = -2 */
    return 0;

  /* now open the File */
  if (CPostFile_Open(File, FileName)) {
    /* Open failed = -4 */
    return 0;
  }
  fd = GiD_HashAdd(File);
  if (!fd) {
    /* could not create a file handler = -6 */
    CPostFile_Release(File);
    return 0;
  }

  if (CPostFile_WritePostHeader(File)) {
    /* WritePostHeader failed = -5 */
    GiD_ClosePostResultFile();
    return 0;
  } else {
    File->level_res = POST_S0;
    File->level_mesh = POST_S0;
  }
  return fd;
}

/*
 *  Close the current post result file
 */

int GiD_ClosePostResultFile()
{
  int fail = 1;

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ClosePostResultFile_HDF5();
  }
#endif

  assert(ResultFile!=NULL);
  assert(CheckState(POST_S0, ResultFile->level_res));    

  if (ResultFile) {
    fail = CPostFile_Release(ResultFile);
    ResultFile = NULL;
    /* reset outputMesh pointer */
    GetMeshFile();
  }
  /*
  level_res = POST_UNDEFINED;
  if (!MeshFile)
    level_mesh = POST_UNDEFINED;
  */
  return fail;
}

int GiD_fClosePostResultFile(GiD_FILE fd)
{
  int fail = 1;
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  fail = CPostFile_Release(File);
  GiD_HashRemove(fd);
  
  return fail;
}

/*
 *  Begin Gauss Points definition
 */

int _GiD_BeginGaussPoint(CPostFile *File,
		         GP_CONST char * name, GiD_ElementType EType,
		         GP_CONST char * MeshName,
		         int GP_number, int NodesIncluded, int InternalCoord)
{
  char line[LINE_SIZE];
  char *gp_name;
  char *mesh_name;
  
  /* check state & validation */
  assert(CheckState(POST_S0, File->level_res));    

  gp_name = change_quotes(strdup( name));

  snprintf(line, LINE_SIZE-1,
	   "GaussPoints \"%s\" ElemType %s",
	   gp_name,
	   GetElementTypeName(EType));
  if (MeshName && *MeshName) {
    mesh_name = change_quotes(strdup(MeshName));
    strcat(line, " \"");
    strcat(line, mesh_name);
    strcat(line, "\"");
    free(mesh_name);
  }
  free(gp_name); 
  gp_name = NULL;
  if (CPostFile_WriteString(File, line))
    return 1;
  snprintf(line, LINE_SIZE, "Number Of Gauss Points: %d", GP_number);
  if (CPostFile_WriteString(File, line))
    return 1;
  /* here we could save the number of GP in order to check at
     EndGaussPoint */
  File->GP_number_check = GP_number;
  if (EType == GiD_Linear) {
    if (NodesIncluded) {
      if (CPostFile_WriteString(File, "  Nodes included"))
	return 1;
    } else
      if (CPostFile_WriteString(File, "  Nodes not included"))
	return 1;
  }
  if (InternalCoord) {
    if (CPostFile_WriteString(File, "Natural Coordinates: Internal"))
      return 1;
    File->level_res = POST_GAUSS_S0;
  } else {
    if (CPostFile_WriteString(File, "Natural Coordinates: Given"))
      return 1;
    File->level_res = POST_GAUSS_GIVEN;    
    /* here we can save the size of the coordinates to check later
       in WriteGaussPointXX*/
  }
  return 0;
}

int GiD_BeginGaussPoint(GP_CONST char * name, GiD_ElementType EType,
		        GP_CONST char * MeshName,
		        int GP_number, int NodesIncluded, int InternalCoord)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginGaussPoint_HDF5(name,EType,MeshName,GP_number,NodesIncluded,InternalCoord);
  }
#endif
  
  return _GiD_BeginGaussPoint(ResultFile,
		              name, EType, MeshName, GP_number,
		              NodesIncluded, InternalCoord);
}

int GiD_fBeginGaussPoint(GiD_FILE fd, GP_CONST char * name,
		         GiD_ElementType EType,
		         GP_CONST char * MeshName,
		         int GP_number, int NodesIncluded, int InternalCoord)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_BeginGaussPoint(File,
		              name, EType, MeshName, GP_number,
		              NodesIncluded, InternalCoord);
}


/*
 *  End current Gauss Points definition
 */

static int CheckGaussPointEnd(CPostFile* File)
{
  if (File->level_res != POST_GAUSS_S0 &&
      File->level_res != POST_GAUSS_GIVEN) {
    printf("Invalid call of GiD_EndGaussPoint. Current state is '%s' and should be '%s' or '%s'\n",
	   GetStateDesc(File->level_res),
	   GetStateDesc(POST_GAUSS_S0),
	   GetStateDesc(POST_GAUSS_GIVEN));
    return 0;
  }
  return 1;
}

static int CheckGaussPointGiven(int written, int check)
{
  if (written !=  check) {
    printf("missmatch in gauss point given, written %d and %d were requiered",
	   written, check);
    return 0;
  }
  return 1;
}

int _GiD_EndGaussPoint(CPostFile *File)
{
  /* check state */
  assert(CheckGaussPointEnd(File));
#ifndef NDEBUG
  if (File->level_res == POST_GAUSS_GIVEN)
    assert(CheckGaussPointGiven(File->gauss_written, File->GP_number_check));
#endif
  File->level_res = POST_S0;
  File->GP_number_check = File->gauss_written = 0;
  return CPostFile_WriteString(File, "End GaussPoints");
}

int GiD_EndGaussPoint()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndGaussPoint_HDF5();
  }
#endif

  return _GiD_EndGaussPoint(ResultFile);
}

int GiD_fEndGaussPoint(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_EndGaussPoint(File);
}

/*
 *  Write internal gauss point coordinate.
 */

int _GiD_WriteGaussPoint2D(CPostFile *File, double x, double y)
{
  /* check state */
  assert(CheckState(POST_GAUSS_GIVEN, File->level_res));    
  if (!CPostFile_Write2D(File, x, y)) {
    ++File->gauss_written;
    return 0;
  }
  return 1;
}

int GiD_WriteGaussPoint2D(double x, double y)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteGaussPoint2D_HDF5(x,y);
  }
#endif

  return _GiD_WriteGaussPoint2D(ResultFile, x, y);
}

int GiD_fWriteGaussPoint2D(GiD_FILE fd, double x, double y)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteGaussPoint2D(File, x, y);
}

int _GiD_WriteGaussPoint3D(CPostFile *File, double x, double y, double z)
{
  /* check state */
  assert(CheckState(POST_GAUSS_GIVEN, File->level_res));
  if (!CPostFile_Write3D(File, x, y, z)) {
    ++File->gauss_written;
    return 0;    
  }
  return 1;
}

int GiD_WriteGaussPoint3D(double x, double y, double z)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteGaussPoint3D_HDF5(x,y,z);
  }
#endif
  
  return _GiD_WriteGaussPoint3D(ResultFile, x, y, z);
}

int GiD_fWriteGaussPoint3D(GiD_FILE fd, double x, double y, double z)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteGaussPoint3D(File, x, y, z);
}

/*
 *  Begin a Range Table definition
 */

int _GiD_BeginRangeTable(CPostFile *File, GP_CONST char * name)
{
  char line[LINE_SIZE];
  char *rt_name;
  
  /* check & update state */
  assert(CheckState(POST_S0, File->level_res));
  
  rt_name = change_quotes(strdup( name));
  snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", rt_name);
  free(rt_name);
  if (!CPostFile_WriteString(File, line)) {
    File->level_res = POST_RANGE_S0;
    return 0;
  }
  return 1;
}

int GiD_BeginRangeTable(GP_CONST char * name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginRangeTable_HDF5(name);
  }
#endif
  
  return _GiD_BeginRangeTable(ResultFile, name); 
}

int GiD_fBeginRangeTable(GiD_FILE fd, GP_CONST char * name)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_BeginRangeTable(File, name);
}

/*
 *  End a Range Table definition
 */

int _GiD_EndRangeTable(CPostFile *File)
{
  /* check & update state */
  assert(CheckState(POST_RANGE_S0, File->level_res));

  if (!CPostFile_WriteString(File, "End ResultRangesTable")) { 
    File->level_res = POST_S0;
    return 0;
  }
  return 1;
}

int GiD_EndRangeTable()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndRangeTable_HDF5();
  }
#endif

  return _GiD_EndRangeTable(ResultFile);
}

int GiD_fEndRangeTable(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_EndRangeTable(File);
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

int _GiD_WriteMinRange(CPostFile* File, double max, GP_CONST char * name)
{
  char line[LINE_SIZE];
  char *tmp_name;
  
  /* check state */
  assert(CheckState(POST_RANGE_S0, File->level_res));

  tmp_name = change_quotes(strdup( name));
  snprintf(line, LINE_SIZE-1, " - %g : \"%s\"", max, tmp_name);
  free(tmp_name);
  return CPostFile_WriteString(File, line);
}

int GiD_WriteMinRange(double max, GP_CONST char * name)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteMinRange_HDF5(max,name);
  }
#endif

  return _GiD_WriteMinRange(ResultFile, max, name);
}

int GiD_fWriteMinRange(GiD_FILE fd, double max, GP_CONST char * name)
{
  CPostFile *File = NULL;
  
  FD2FILE(fd,File);

  return _GiD_WriteMinRange(File, max, name);
}

int _GiD_WriteRange(CPostFile* File,
		    double min, double max, GP_CONST char * name)
{
  char line[LINE_SIZE];
  char *tmp_name;
    
/* check state */
  assert(CheckState(POST_RANGE_S0, File->level_res));
  
  tmp_name = change_quotes(strdup( name));
  snprintf(line, LINE_SIZE-1, " %g - %g : \"%s\"", min, max, tmp_name);
  free( tmp_name);
  return CPostFile_WriteString(File, line);
}

int GiD_WriteRange(double min, double max, GP_CONST char * name)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteRange_HDF5(min,max,name);
  }
#endif
  
  return _GiD_WriteRange(ResultFile, min, max, name);
}

int GiD_fWriteRange(GiD_FILE fd, double min, double max, GP_CONST char * name)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_WriteRange(File, min, max, name);
}

int _GiD_WriteMaxRange(CPostFile *File, double min, GP_CONST char * name)
{
  char line[LINE_SIZE];
  char *tmp_name;
  
  /* check state */
  assert(CheckState(POST_RANGE_S0, File->level_res));
  
  tmp_name = change_quotes( strdup( name));
  snprintf(line, LINE_SIZE-1, "%g - : \"%s\"", min, tmp_name);
  free( tmp_name);
  return CPostFile_WriteString(File, line);
}

int GiD_WriteMaxRange( double min, GP_CONST char * name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteMaxRange_HDF5(min,name);
  }
#endif

  return _GiD_WriteMaxRange(ResultFile, min, name);
}

int GiD_fWriteMaxRange(GiD_FILE fd, double min, GP_CONST char * name)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WriteMaxRange(File, min, name);
}

/*
 *  Begin Result Block
 */

int _GiD_BeginResult(CPostFile *File,
		     GP_CONST char     *Result,
		     GP_CONST char     *Analysis,
		     double             step,
		     GiD_ResultType     Type,
		     GiD_ResultLocation Where,
		     GP_CONST char     *GaussPointsName,
		     GP_CONST char     *RangeTable, 
		     int                compc,
		     GP_CONST char      *compv[])
{
  char line[LINE_SIZE];
  const char * loc;
  char *res_name, *analysis_name, *tmp_name;
  int i;

  /* check & change state */
  assert(CheckState(POST_S0, File->level_res));

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  res_name = change_quotes(strdup( Result));
  analysis_name = change_quotes(strdup( Analysis));
  snprintf(line, LINE_SIZE-1, "Result \"%s\" \"%s\" %.16g %s %s",
	   res_name, analysis_name, step, GetResultTypeName(Type, 0), loc);
  free(res_name);
  free(analysis_name);
  if (Where == GiD_OnGaussPoints) {
    assert(GaussPointsName);
    tmp_name = change_quotes(strdup(GaussPointsName));
    strcat(line, " \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free(tmp_name);
  }
  if (CPostFile_WriteString(File, line))
    return 1;
  if (RangeTable) {
    assert(RangeTable);
    tmp_name = change_quotes(strdup(RangeTable));
    snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", tmp_name);
    free(tmp_name);
    if (CPostFile_WriteString(File, line))
      return 1;
  }
  if (compc > 0) {
    snprintf(line, LINE_SIZE-1, "ComponentNames");
    for (i = 0; i < compc; i++) {  
      tmp_name = change_quotes(strdup(compv[ i]));
      strcat(line, " ");
      strcat(line, "\"");
      strcat(line, tmp_name);
      strcat(line, "\"");
      free(tmp_name);
    }
    if (CPostFile_WriteString(File, line))
      return 1;
  }
  File->flag_isgroup = 0;
  File->flag_begin_values = 1;
  File->level_res = POST_RESULT_VALUES;
  return CPostFile_BeginValues(File);
}

int GiD_BeginResult(GP_CONST char     *Result,
		    GP_CONST char     *Analysis,
		    double             step,
		    GiD_ResultType     Type,
		    GiD_ResultLocation Where,
		    GP_CONST char     *GaussPointsName,
		    GP_CONST char     *RangeTable, 
		    int                compc,
		    GP_CONST char     *compv[])
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginResult_HDF5(Result,Analysis,step,Type,Where,GaussPointsName,RangeTable,compc,compv);
  }
#endif
  
  return _GiD_BeginResult(ResultFile, Result, Analysis, step, Type, Where,
  GaussPointsName, RangeTable, compc, compv);
}

int GiD_fBeginResult(GiD_FILE fd, GP_CONST char * Result,
		     GP_CONST char * Analysis,
		     double step,
		     GiD_ResultType Type, GiD_ResultLocation Where,
		     GP_CONST char * GaussPointsName,
		     GP_CONST char * RangeTable, 
		     int compc, GP_CONST char * compv[])
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_BeginResult(File, Result, Analysis, step, Type, Where,
		          GaussPointsName, RangeTable, compc, compv);
}

int _GiD_BeginResultHeader(CPostFile         *File,
		           GP_CONST char     *Result,
		           GP_CONST char     *Analysis,
		           double             step,
		           GiD_ResultType     Type,
		           GiD_ResultLocation Where,
		           GP_CONST char     *GaussPointsName)
{
  char line[LINE_SIZE];
  const char * loc;
  char *res_name, *analysis_name, *tmp_name;
  
  /* check & change state */
  assert(File);
  assert(CheckState(POST_S0, File->level_res));
  assert(Result);
  assert(Analysis);

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  res_name = change_quotes(strdup(Result));
  analysis_name = change_quotes(strdup(Analysis));
  snprintf(line, LINE_SIZE-1, "Result \"%s\" \"%s\" %.16g %s %s",
	   res_name, analysis_name, step, GetResultTypeName(Type,0), loc);
  free(res_name);
  free(analysis_name);
  if (Where == GiD_OnGaussPoints) {
    assert(GaussPointsName);
    tmp_name = change_quotes(strdup(GaussPointsName));
    strcat(line, " \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
  }
  if (CPostFile_WriteString(File, line)) {
    /* could not write result header */
    return 1;
  }
  File->level_res = POST_RESULT_SINGLE;
  File->flag_isgroup = 0;
  File->flag_begin_values = 0;
  return 0;
}

int GiD_BeginResultHeader(GP_CONST char     *Result,
		          GP_CONST char     *Analysis,
		          double             step,
		          GiD_ResultType     Type,
		          GiD_ResultLocation Where,
		          GP_CONST char     *GaussPointsName)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginResultHeader_HDF5(Result,Analysis,step,Type,Where,GaussPointsName);
  }
#endif
  
  return _GiD_BeginResultHeader(ResultFile, Result, Analysis, step, Type,
  Where, GaussPointsName);
}

int GiD_fBeginResultHeader(GiD_FILE fd, GP_CONST char * Result,
		           GP_CONST char * Analysis, 
		           double step,
		           GiD_ResultType Type, GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_BeginResultHeader(File, Result, Analysis, step, Type,
		                Where, GaussPointsName);
}

static int CheckResultHeaderState(CPostFile* File)
{
  if (File->level_res == POST_RESULT_SINGLE ||
      File->level_res == POST_RESULT_DESC)
    return 1;
  printf("Invalid result state '%s'. Should be: '%s' or '%s'\n",
	 GetStateDesc(File->level_res),
	 GetStateDesc(POST_RESULT_SINGLE),
	 GetStateDesc(POST_RESULT_DESC));
  return 0;
}

int _GiD_ResultRange(CPostFile *File, GP_CONST char * RangeTable)
{
  char line[LINE_SIZE];
  char *tmp_name;

  /* check state */
  assert(File);
  assert(CheckResultHeaderState(File));
  assert(RangeTable);
  
  if (RangeTable) {
    tmp_name = change_quotes(strdup(RangeTable));
    snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", tmp_name);
    free(tmp_name);
    return CPostFile_WriteString(File, line);
  }
  return 1;
}

int GiD_ResultRange(GP_CONST char * RangeTable)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultRange_HDF5(RangeTable);
  }
#endif
  
  return _GiD_ResultRange(ResultFile, RangeTable);
}

int GiD_fResultRange(GiD_FILE fd, GP_CONST char * RangeTable)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_ResultRange(File, RangeTable);
}

int _GiD_ResultComponents(CPostFile *File, int compc, GP_CONST char * compv[])
{
  char line[LINE_SIZE];
  char *tmp_name;
  int i;
  
  /* check state */
  assert(File);
  assert(CheckResultHeaderState(File));
  assert(compc>0);
  
  if (compc > 0) {
    snprintf(line, LINE_SIZE-1, "ComponentNames");
    for (i = 0; i < compc; i++) {
      tmp_name = change_quotes(strdup(compv[ i]));
      strcat(line, " ");
      strcat(line, "\"");
      strcat(line, tmp_name);
      strcat(line, "\"");
      free(tmp_name);
    }
    return CPostFile_WriteString(File, line);
  }
  return 1;
}

int GiD_ResultComponents(int compc, GP_CONST char * compv[])
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultComponents_HDF5(compc,compv);
  }
#endif

  return _GiD_ResultComponents(ResultFile, compc, compv);
}

int GiD_fResultComponents(GiD_FILE fd, int compc, GP_CONST char * compv[])
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_ResultComponents(File, compc, compv);
}

int _GiD_ResultUnit(CPostFile *File, GP_CONST char * UnitName)
{
  /* check state */
  assert(CheckResultHeaderState(ResultFile));
  assert(UnitName);

  if (UnitName) {
    char line[LINE_SIZE];
    char *tmp_name;
    tmp_name = change_quotes( strdup( UnitName));
    snprintf(line, LINE_SIZE-1, "Unit \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
    return CPostFile_WriteString(File, line);
  }
  return 1;
}

int GiD_ResultUnit(GP_CONST char * UnitName)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultUnit_HDF5(UnitName);
  }
#endif
  return _GiD_ResultUnit(ResultFile, UnitName);  
}

int GiD_fResultUnit(GiD_FILE fd, GP_CONST char * UnitName)
{
  CPostFile *File = NULL;
  FD2FILE(fd,File);
  return _GiD_ResultUnit(File, UnitName);  
}

int _GiD_ResultUserDefined(CPostFile *File, GP_CONST char * Name,GP_CONST char * Value)
  {  
  char line[LINE_SIZE];
  char* tmp_name=change_quotes(strdup(Name));
  char* tmp_value=change_quotes(strdup(Value));
  if(string_hasspace(Name)){
    if(string_hasspace(Value)){
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined \"%s\" %s",tmp_name,tmp_value);  
    } else {
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined \"%s\" \"%s\"",tmp_name,tmp_value);  
    }
  } else {
    if(string_hasspace(Value)){
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined %s %s",tmp_name,tmp_value);  
    } else {
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined %s \"%s\"",tmp_name,tmp_value);  
    }
  }
  free(tmp_name);
  free(tmp_value);
  return CPostFile_WriteString(File, line);
}

int GiD_ResultUserDefined(GP_CONST char * Name,GP_CONST char * Value)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultUserDefined_HDF5(Name,Value);
  }
#endif
  return _GiD_ResultUserDefined(ResultFile,Name,Value);
}

int GiD_fResultUserDefined(GiD_FILE fd, GP_CONST char * Name,GP_CONST char * Value)
{
  CPostFile *File = NULL;
  FD2FILE(fd,File);
  return _GiD_ResultUserDefined(File,Name,Value);
}

int _GiD_BeginResultGroup(CPostFile         *File,
		          GP_CONST char     *Analysis,
		          double             step,
		          GiD_ResultLocation Where,
		          GP_CONST char     *GaussPointsName)
{
  char line[LINE_SIZE];
  const char * loc;
  char *analysis_name, *tmp_name;
  
  /* check & change state */
  assert(File);
  assert(CheckState(POST_S0, File->level_res));
  assert(Analysis);

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  analysis_name = change_quotes(strdup(Analysis));
  snprintf(line, LINE_SIZE-1, "ResultGroup \"%s\" %.16g %s",
	   analysis_name, step, loc);
  free(analysis_name);
  if (Where == GiD_OnGaussPoints) {
    assert(GaussPointsName);
    tmp_name = change_quotes(strdup(GaussPointsName));
    strcat(line, " \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
  }
  if (CPostFile_WriteString(File, line)) {
    /* could not write result header */
    return 1;
  }
  File->level_res = POST_RESULT_GROUP;
  /* initialize values counter */
  File->flag_isgroup = 1;
  File->flag_begin_values = 0;
  CPostFile_ResultGroupOnBegin(File);
  return 0;
}

int GiD_BeginResultGroup(GP_CONST char     *Analysis,
		         double             step,
		         GiD_ResultLocation Where,
		         GP_CONST char     *GaussPointsName)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginResultGroup_HDF5(Analysis,step,Where,GaussPointsName);
  }
#endif

  return _GiD_BeginResultGroup(ResultFile, Analysis, step, Where,
		               GaussPointsName);
}

int GiD_fBeginResultGroup(GiD_FILE fd, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_BeginResultGroup(File, Analysis, step, Where,
		               GaussPointsName);
}

static int CheckStateDesc(CPostFile* File)
{
  if (File->level_res == POST_RESULT_GROUP ||
      File->level_res == POST_RESULT_DESC)
    return 1;
  printf("Invalid result state '%s'. Should be: '%s' or '%s'\n",
	 GetStateDesc(File->level_res),
	 GetStateDesc(POST_RESULT_GROUP),
	 GetStateDesc(POST_RESULT_DESC));
  return 0;
}

static
int _GiD_ResultDescription_(CPostFile     *File,
		            GP_CONST char *Result,
		            GiD_ResultType Type,
		            size_t         s)
{
  char line[LINE_SIZE];
  char *tmp_name;

  /* check & change state */
  assert(File);
  assert(CheckStateDesc(File));
  assert(Result);

  line[0] = '\0';
  tmp_name = change_quotes(strdup(Result));
  snprintf(line, LINE_SIZE-1, "ResultDescription \"%s\" %s",
	   tmp_name, GetResultTypeName(Type,s));
  free(tmp_name);
  if (CPostFile_WriteString(File, line)) {
    /* could not write result description */
    return 1;
  }
  File->level_res = POST_RESULT_DESC;
  /* update number of values to check per location */
  CPostFile_ResultGroupOnNewType(File, Type);
  return 0;
}

int GiD_ResultDescription(GP_CONST char * Result, GiD_ResultType Type)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultDescription_HDF5(Result,Type);
  }
#endif

  return _GiD_ResultDescription_(ResultFile, Result, Type, 0);
}

int GiD_fResultDescription(GiD_FILE fd,
		           GP_CONST char * Result, GiD_ResultType Type)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_ResultDescription_(File, Result, Type, 0);
}

int GiD_ResultDescriptionDim(GP_CONST char * Result, GiD_ResultType Type,
		             size_t s)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    /* disregard s in hdf5, because ResultGroup is really written in HDF5 as independent Results*/
    return GiD_ResultDescription_HDF5(Result,Type);
  }
#endif

  return _GiD_ResultDescription_(ResultFile, Result, Type, s);
}

int GiD_fResultDescriptionDim(GiD_FILE fd,
		              GP_CONST char * Result, GiD_ResultType Type,
		              size_t s)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_ResultDescription_(File, Result, Type, s);
}

int GiD_ResultLocalAxes(GP_CONST char * Result, GP_CONST char * Analysis,
  double step,double vx,double vy,double vz)
{
#ifdef HDF5
    if(PostMode==GiD_PostHDF5){
    return GiD_ResultLocalAxes_HDF5(Result,Analysis,step,vx,vy,vz);
  }
#endif
  return 1;
}

int GiD_ResultIsDeformationVector(int boolean)
{
#ifdef HDF5
    if(PostMode==GiD_PostHDF5){
    return GiD_ResultIsDeformationVector_HDF5(boolean);
  }
#endif
  return 1;
}

int GiD_ResultValues()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultValues_HDF5();
  }
#endif

  return 0;
}

int GiD_fResultValues(GiD_FILE fd)
{
  return 0;
}

/*
 *  End Result Block
 */

int _GiD_EndResult(CPostFile *File)
{
  int _fail = 0;
  int status;

  /* check & change state */
  assert(File);
  assert(CheckState(POST_RESULT_VALUES, File->level_res));
  if (File->flag_isgroup) {
    status = CPostFile_ResultGroupIsEmpty(File);
    assert(status);
  }

  _fail = CPostFile_EndValues(File);
  CPostFile_ResetLastID(File);
  File->level_res = POST_S0;
  return _fail;
}

int GiD_EndResult()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndResult_HDF5();
  }
#endif

  return _GiD_EndResult(ResultFile);
}

int GiD_fEndResult(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_EndResult(File);
}

/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */

int _GiD_BeginOnMeshGroup(CPostFile *File, char * Name)
{
  char line[LINE_SIZE];
  
  snprintf(line, LINE_SIZE-1, "OnGroup \"%s\"", Name);
  return CPostFile_WriteString(File, line);
}

int GiD_BeginOnMeshGroup(char * Name)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return -1;
  }
#endif

  return _GiD_BeginOnMeshGroup(ResultFile, Name);
}

int GiD_fBeginOnMeshGroup(GiD_FILE fd, char * Name)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return _GiD_BeginOnMeshGroup(File, Name);
}

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */

int GiD_EndOnMeshGroup()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return -1;
  }
#endif

  return CPostFile_WriteString(ResultFile, "End OnGroup");
}

int GiD_fEndOnMeshGroup(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return CPostFile_WriteString(File, "End OnGroup");
}

/*
 * Flushes all pending output into the compressed file.
 */

int GiD_FlushPostFile()
{
  int res1 = 0, res2 = 0;
  
#ifdef HDF5
  if( PostMode == GiD_PostHDF5 ) {
    return GiD_FlushPostFile_HDF5();
  }
#endif
  if( MeshFile ) res1 = CPostFile_Flush( MeshFile );
  if( !res1 && ResultFile ) res2 = CPostFile_Flush( ResultFile );
  return res1 || res2;
}

int GiD_fFlushPostFile(GiD_FILE fd)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);

  return CPostFile_Flush(File);
}

/*
 *  Write result functions
 */

static
int GiD_EnsureBeginValues(CPostFile *File)
{
  assert(File);
  
  if (!File->flag_begin_values) {
    if (!CPostFile_BeginValues(File)) {
      File->level_res = POST_RESULT_VALUES;
      if (File->flag_isgroup) {
	CPostFile_ResultGroupOnBeginValues(File);
      }
      File->flag_begin_values = 1;
      return 0;
    }
  }
  return 1;
}

int _GiD_WriteScalar(CPostFile *File, int id, double v )
{
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));

  return File->flag_isgroup ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_Scalar, id, 1, v) :
    CPostFile_WriteValues(File, id, 1, &v);
}

int GiD_WriteScalar(int id, double v)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteScalar_HDF5(id,v);
  }
#endif

  return _GiD_WriteScalar(ResultFile, id, v);
}

int GiD_fWriteScalar(GiD_FILE fd, int id, double v)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WriteScalar(File, id, v);
}

int _GiD_Write2DVector(CPostFile *File, int id, double x, double y)
{
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));

  if (File->flag_isgroup || !CPostFile_IsBinary(File)) {
    return File->flag_isgroup ?
      CPostFile_ResultGroupWriteValues(File,
		                       GiD_Vector, id, 2, x, y) :
      CPostFile_WriteValuesVA(File, id, 2, x, y);
  } else {
    /* single result & binary */
    double mod = sqrt(x*x + y*y);
  
    return CPostFile_WriteValuesVA(File, id, 4, x, y, 0.0, mod);
  }  
}

int GiD_Write2DVector(int id, double x, double y)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_Write2DVector_HDF5(id,x,y);
  }
#endif

  return _GiD_Write2DVector(ResultFile, id, x, y);
}

int GiD_fWrite2DVector(GiD_FILE fd, int id, double x, double y)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_Write2DVector(File, id, x, y);
}

int _GiD_WriteVector(CPostFile *File, int id, double x, double y, double z)
{
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));

  if (File->flag_isgroup || !CPostFile_IsBinary(File)) {
    return File->flag_isgroup
      ?
      CPostFile_ResultGroupWriteValues(File,
		                       GiD_Vector, id, 3, x, y, z)
      :
      CPostFile_WriteValuesVA(File, id, 3, x, y, z);
  } else {
    /* single result & binary */
    double mod = sqrt(x*x + y*y + z*z);
  
    return CPostFile_WriteValuesVA(File, id, 4, x, y, z, mod);
  }  
}

int GiD_WriteVector(int id, double x, double y, double z)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteVector_HDF5(id,x,y,z);
  }
#endif
  
  return _GiD_WriteVector(ResultFile, id, x, y, z);
}

int GiD_fWriteVector(GiD_FILE fd, int id, double x, double y, double z)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WriteVector(File, id, x, y, z);
}

int _GiD_WriteVectorModule(CPostFile *File,
		           int id, double x, double y, double z, double mod)
{
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));

  /* 4-vectors can not be written on RG-ASCII */
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_Vector, id, 4, x, y, z, mod)
    :
    CPostFile_WriteValuesVA(File, id, 4, x, y, z, mod );
}

int GiD_WriteVectorModule(int id, double x, double y, double z, double mod)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteVectorModule_HDF5(id,x,y,z,mod);
  }
#endif

  return _GiD_WriteVectorModule(ResultFile, id, x, y, z, mod);
}

int GiD_fWriteVectorModule(GiD_FILE fd,
		           int id, double x, double y, double z, double mod)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WriteVectorModule(File, id, x, y, z, mod);
}

int _GiD_Write3DMatrix(CPostFile *File,
		       int id, double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz);

int _GiD_Write2DMatrix(CPostFile *File,
		       int id, double Sxx, double Syy, double Sxy)
{
  assert(File);
  if (CPostFile_IsBinary(File)) {
    return _GiD_Write3DMatrix(File,
		              id,
		              Sxx, Syy,         0.0 /*Szz*/,
		              Sxy, 0.0 /*Syz*/, 0.0 /*Sxz*/);
  }
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_Matrix, id, 3, Sxx, Syy, Sxy)
    :
    CPostFile_WriteValuesVA(File, id, 3, Sxx, Syy, Sxy);
}

int GiD_Write2DMatrix(int id, double Sxx, double Syy, double Sxy)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_Write2DMatrix_HDF5(id,Sxx,Syy,Sxy);
  }
#endif

  return _GiD_Write2DMatrix(ResultFile, id, Sxx, Syy, Sxy);
}

int GiD_fWrite2DMatrix(GiD_FILE fd, int id, double Sxx, double Syy, double Sxy)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_Write2DMatrix(File, id, Sxx, Syy, Sxy);
}

int _GiD_Write3DMatrix(CPostFile *File,
		       int id, double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz)
{
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_Matrix,
		                     id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz)
    :
    CPostFile_WriteValuesVA(File, id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz);
}

int GiD_Write3DMatrix(int id,
		      double Sxx, double Syy, double Szz,
		      double Sxy, double Syz, double Sxz)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_Write3DMatrix_HDF5(id,Sxx,Syy,Szz,Sxy,Syz,Sxz);
  }
#endif

  return _GiD_Write3DMatrix(ResultFile, id,
		            Sxx, Syy, Szz,
		            Sxy, Syz, Sxz);
}

int GiD_fWrite3DMatrix(GiD_FILE fd, int id,
		       double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_Write3DMatrix(File, id,
		            Sxx, Syy, Szz,
		            Sxy, Syz, Sxz);
}

int _GiD_WritePlainDefMatrix(CPostFile *File, int id,
		             double Sxx, double Syy, double Sxy, double Szz )
{
#if 0
    if (CPostFile_IsBinary(File)) {
      return _GiD_Write3DMatrix(File, id,
                                Sxx, Syy,         Szz,
                                Sxy, 0.0 /*Syz*/, 0.0 /*Sxz*/);
    }
#endif
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));
   
  return File->flag_isgroup
    ? 
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_PlainDeformationMatrix,
		                     id, 4, Sxx, Syy, Sxy, Szz)
    :
    CPostFile_WriteValuesVA(File, id, 4, Sxx, Syy, Sxy, Szz);
}

int GiD_WritePlainDefMatrix(int id,
		            double Sxx, double Syy, double Sxy, double Szz )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WritePlainDefMatrix_HDF5(id,Sxx,Syy,Sxy,Szz);
  }
#endif
  
  return _GiD_WritePlainDefMatrix(ResultFile,
		                  id, Sxx, Syy, Sxy, Szz);
}

int GiD_fWritePlainDefMatrix(GiD_FILE fd, int id,
		             double Sxx, double Syy, double Sxy, double Szz )
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WritePlainDefMatrix(File,
		                  id, Sxx, Syy, Sxy, Szz);
}

int _GiD_WriteMainMatrix(CPostFile *File, int id,
		         double Si, double Sii, double Siii,
		         double Vix, double Viy, double Viz,
		         double Viix, double Viiy, double Viiz,
		         double Viiix, double Viiiy, double Viiiz)
{
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));
  
  return File->flag_isgroup
    ? 
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_MainMatrix, id, 12,
		                     Si,    Sii,   Siii,
		                     Vix,   Viy,   Viz,
		                     Viix,  Viiy,  Viiz,
		                     Viiix, Viiiy, Viiiz)
    :
    CPostFile_WriteValuesVA(File, id, 12,
		            Si,    Sii,   Siii,
		            Vix,   Viy,   Viz,
		            Viix,  Viiy,  Viiz,
		            Viiix, Viiiy, Viiiz);
}

int GiD_WriteMainMatrix(int id,
		        double Si, double Sii, double Siii,
		        double Vix, double Viy, double Viz,
		        double Viix, double Viiy, double Viiz,
		        double Viiix, double Viiiy, double Viiiz)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteMainMatrix_HDF5(id,Si,Sii,Siii,Vix,Viy,Viz,Viix,Viiy,Viiz,Viiix,Viiiy,Viiiz);
  }
#endif

  return _GiD_WriteMainMatrix(ResultFile, id, 
		              Si,  Sii,  Siii,
		              Vix,  Viy,  Viz,
		              Viix,  Viiy,  Viiz,
		              Viiix,  Viiiy,  Viiiz);
}

int GiD_fWriteMainMatrix(GiD_FILE fd, int id,
		         double Si, double Sii, double Siii,
		         double Vix, double Viy, double Viz,
		         double Viix, double Viiy, double Viiz,
		         double Viiix, double Viiiy, double Viiiz)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WriteMainMatrix(File, id, 
		              Si,  Sii,  Siii,
		              Vix,  Viy,  Viz,
		              Viix,  Viiy,  Viiz,
		              Viiix,  Viiiy,  Viiiz);
}

int _GiD_WriteLocalAxes(CPostFile *File,
		        int id, double euler_1, double euler_2, double euler_3)
{
  GiD_EnsureBeginValues(File);
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, File->level_res));
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_LocalAxes, id, 3,
		                     euler_1, euler_2, euler_3)
    :
    CPostFile_WriteValuesVA(File, id, 3, euler_1, euler_2, euler_3);
}

int GiD_WriteLocalAxes(int id, double euler_1, double euler_2, double euler_3)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteLocalAxes_HDF5(id,euler_1,euler_2,euler_3);
  }
#endif

  return _GiD_WriteLocalAxes(ResultFile, id, euler_1, euler_2, euler_3);
}

int GiD_fWriteLocalAxes(GiD_FILE fd,
		        int id, double euler_1, double euler_2, double euler_3)
{
  CPostFile *File = NULL;

  FD2FILE(fd,File);
  
  return _GiD_WriteLocalAxes(File, id, euler_1, euler_2, euler_3);
}

/*
 * Complex numbers
 */

int _GiD_WriteComplexScalar( CPostFile *File,
                             int id, double complex_real, double complex_imag) {
  GiD_EnsureBeginValues( File);
  /* check state */
  assert( CheckState( POST_RESULT_VALUES, File->level_res));
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues( File,
				      GiD_ComplexScalar, id, 2,
				      complex_real, complex_imag)
    :
    CPostFile_WriteValuesVA(File, id, 2, complex_real, complex_imag);
}

int GiD_WriteComplexScalar( int id, double complex_real, double complex_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_WriteComplexScalar_HDF5(id, complex_real, complex_imag);
  }
#endif
  return _GiD_WriteComplexScalar(ResultFile, id, complex_real, complex_imag);
}

int GiD_fWriteComplexScalar(GiD_FILE fd, int id, double complex_real, double complex_imag) {
  CPostFile *File = NULL;
  FD2FILE( fd, File);
  return _GiD_WriteComplexScalar(File, id, complex_real, complex_imag);
}

int _GiD_Write2DComplexVector( CPostFile *File, int id,
			       double x_real, double x_imag,
			       double y_real, double y_imag) {
  double mod2_r = x_real * x_real + y_real * y_real;
  double mod2_i = x_imag * x_imag + y_imag * y_imag;
  double mod_r = sqrt( mod2_r);
  double mod_i = sqrt( mod2_i);
  double mod = sqrt( mod2_r + mod2_i);

  GiD_EnsureBeginValues( File);
  /* check state */
  assert( CheckState( POST_RESULT_VALUES, File->level_res));
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues( File,
				      GiD_ComplexVector, id,  9,
				      x_real, x_imag, y_real, y_imag, 0.0, 0.0,
				      mod_r, mod_i, mod)
    :
    CPostFile_WriteValuesVA(File, id, 9, 
			    x_real, x_imag, y_real, y_imag, 0.0, 0.0,
			    mod_r, mod_i, mod);
}

int GiD_Write2DComplexVector( int id, 
			      double x_real, double x_imag,
			      double y_real, double y_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_Write2DComplexVector_HDF5(id, x_real, x_imag, y_real, y_imag);
  }
#endif
  return _GiD_Write2DComplexVector(ResultFile, id, x_real, x_imag, y_real, y_imag);
}

int GiD_fWrite2DComplexVector( GiD_FILE fd, int id,
			       double x_real, double x_imag,
			       double y_real, double y_imag) {
  CPostFile *File = NULL;
  FD2FILE( fd, File);
  return _GiD_Write2DComplexVector(File, id, x_real, x_imag, y_real, y_imag);
}

int _GiD_WriteComplexVector( CPostFile *File, int id,
			     double x_real, double x_imag,
			     double y_real, double y_imag,
			     double z_real, double z_imag) {
  double mod2_r = x_real * x_real + y_real * y_real + z_real * z_real;
  double mod2_i = x_imag * x_imag + y_imag * y_imag + z_imag * z_imag;
  double mod_r = sqrt( mod2_r);
  double mod_i = sqrt( mod2_i);
  double mod = sqrt( mod2_r + mod2_i);
  
  GiD_EnsureBeginValues( File);
  /* check state */
  assert( CheckState( POST_RESULT_VALUES, File->level_res));

  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues( File,
				      GiD_ComplexVector, id, 9,
				      x_real, x_imag, y_real, y_imag, z_real, z_imag,
				      mod_r, mod_i, mod)
    :
    CPostFile_WriteValuesVA(File, id, 9, 
			    x_real, x_imag, y_real, y_imag, z_real, z_imag,
			    mod_r, mod_i, mod);
}

int GiD_WriteComplexVector( int id, 
			    double x_real, double x_imag,
			    double y_real, double y_imag,
			    double z_real, double z_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_WriteComplexVector_HDF5(id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
  }
#endif
  return _GiD_WriteComplexVector(ResultFile, id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
}

int GiD_fWriteComplexVector( GiD_FILE fd, int id,
			     double x_real, double x_imag,
			     double y_real, double y_imag,
			     double z_real, double z_imag) {
  CPostFile *File = NULL;
  FD2FILE( fd, File);
  return _GiD_WriteComplexVector(File, id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
}
