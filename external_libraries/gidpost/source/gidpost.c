/* gidpost */
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
#include <assert.h>

#include "gidpostInt.h"
#include "gidpost.h"
#include "gidpostHash.h"
#include "gidpostFILES.h"

#ifdef ENABLE_HDF5
  #include "gidpostHDF5.h"
#endif

/* defined in gidpostFILES.c */
extern CPostFile *G_MeshFile;
extern CPostFile *G_ResultFile;
extern CPostFile *G_outputMesh;
// extern GiD_PostMode G_PostMode;

#if 0
#define FD2FILE(_fd_, _File_)                     \
  do {                                            \
    assert(_fd_);                                 \
    _File_ = (CPostFile*)GiD_HashFind(_fd_);      \
    if (!_File_) {                                \
      /* file handler does not exists */          \
      return -8;                                  \
    }                                             \
  } while (0);
#endif

CPostFile *FD2FILE( GiD_FILE fd) {
  CPostFile *File = NULL;
  assert( fd);
  File = ( CPostFile *)GiD_HashFind( fd);
  return File;
}

int GiD_PostInit()
{ // ensure GiD_PostInit is called only once if we are doing more things here
  // already doing it inside GiD_HashInit()
  /* const char *lib_version = */
  GiD_PostGetVersion(); // so that the 'static char *' are created (hdf5 may not be thread safe)
  const int ret = GiD_HashInit();
#ifndef NDEBUG
  printf( "Debug version: %s\n", GiD_PostGetVersion() );
#endif // NDEBUG
  return ret;
}

int GiD_PostDone()
{
  return GiD_HashDone();
}


GIDPOST_API GP_CONST char *GiD_PostGetVersion( void ) {
  static char *version_str = NULL;
  if ( !version_str ) {
    char buf[ 1024 ];
    snprintf( buf, 1024, "GiDpost %d.%d", GP_VERSION_MAJOR, GP_VERSION_MINOR );
#ifdef ENABLE_HDF5
    const char *hdf5_version = GiD_GetHDF5Version();
    if ( hdf5_version ) {
      size_t len_buf = strlen( buf );
      snprintf( &buf[ len_buf ], 1024 - len_buf, " - %s", hdf5_version );
    }
#endif // ENABLE_HDF5
    version_str = strdup( buf );
  }
  return version_str;
}

// to check if hdf5 is ThreadSafe to be used when threads/OpenMP are writing in hdf5
GIDPOST_API int GiD_PostIsThreadSafe( GiD_PostMode Mode ) {
  int is_thread_safe = 1;
  if ( Mode == GiD_PostHDF5 ) {
    is_thread_safe = -1;
#ifdef ENABLE_HDF5
    is_thread_safe = GiD_IsThreadSafe_HDF5();
#endif // ENABLE_HDF5
  }
  return is_thread_safe;
}

/* ---------------------------------------------------------------------------
 *
 *  Post Mesh Interface
 *
 * ---------------------------------------------------------------------------
 */

// Do not complain about implementing deprecated api
#ifdef WIN32
// #if defined( _MSC_VER )
// disable deprecated declarations
#pragma warning(disable:4996)
// #endif
#else // WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif // WIN32

/*
 *  Open a new post mesh file
 */

int GiD_OpenPostMeshFile(GP_CONST char * FileName, GiD_PostMode Mode )
{
  // G_PostMode=Mode;

  assert(!G_MeshFile);
  if (G_MeshFile) {
    /* must be closed */
    return GP_ERROR_FILEOPENED;
  }
  G_MeshFile = _GiDfiles_NewFile( Mode ); // inside does also this G_MeshFile->m_post_mode = Mode;
  if (!G_MeshFile )
    /* not enough memory */
    return GP_ERROR_NOMEM;
  _GiDfiles_GetMeshFile();

  if(Mode==GiD_PostHDF5){
#ifdef ENABLE_HDF5
    if ( G_MeshFile->m_hdf5_file != NULL ) {
      return GP_ERROR_FILEOPENED;
    }
    G_MeshFile->m_hdf5_file = new_CurrentHdf5WriteData();
    assert(G_MeshFile->m_hdf5_file );
    if ( !G_MeshFile->m_hdf5_file ) {
      return GP_ERROR_OPENFAILED;
    }
    CPostFile_PushState( G_MeshFile, POST_S0 );
    return GiD_OpenPostMeshFile_HDF5( G_MeshFile->m_hdf5_file, FileName);
#else
    return GP_ERROR_INVALID_STATE;
#endif
    }
  /* Binary mode not allowed in mesh file!!! */
  assert(Mode!=GiD_PostBinary);
  /* now G_outputMesh points to G_MeshFile */
  if (CPostFile_Open(G_MeshFile, FileName)) {
    /* Open failed */
    CPostFile_Release(G_MeshFile);
    G_MeshFile = NULL;
    return GP_ERROR_OPENFAILED;
  }
  
  //G_MeshFile->level_mesh = POST_S0;
  CPostFile_PushState( G_MeshFile, POST_S0 );
  return GP_OK;
}

GiD_FILE GiD_fOpenPostMeshFile( GP_CONST char * FileName,
		                GiD_PostMode Mode )
{
  CPostFile *File = NULL;
  File = _GiDfiles_NewFile( Mode ); // inside does also this File->m_post_mode = Mode;
  if ( !File ) {
    /* not enough memory = -2 */
    return 0;
  }

#ifdef ENABLE_HDF5
  if ( Mode == GiD_PostHDF5 ) {

    CurrentHdf5WriteData *my_hdf5_file = new_CurrentHdf5WriteData();
    assert( new_CurrentHdf5WriteData );
    if ( !my_hdf5_file ) {
      return 0;
    }
    int fail = GiD_OpenPostMeshFile_HDF5( my_hdf5_file, FileName );
    if ( fail ) {
      delete_CurrentHdf5WriteData( my_hdf5_file );
      CPostFile_Release( File );
      return 0;
    }
    File->m_hdf5_file = my_hdf5_file;

    GiD_FILE fd;
    fd = GiD_HashAdd( File );
    if ( !fd ) {
      /* could not create a file handler = -6 */
      GiD_ClosePostMeshFile_HDF5( my_hdf5_file );
      delete_CurrentHdf5WriteData( my_hdf5_file );
      CPostFile_Release( File );
      return 0;
    }
    CPostFile_PushState( File, POST_S0 );
    return fd;
  }
#endif
  
  /* Binary mode not allowed in mesh file!!! */
  assert( Mode != GiD_PostBinary );
  /* now open the File */
  if ( CPostFile_Open( File, FileName ) ) {
    /* Open failed = -4 */
    return 0;
  }

  GiD_FILE fd;
  fd = GiD_HashAdd( File );
  if (!fd) {
    /* could not create a file handler = -6 */
    CPostFile_Release(File);
    return 0;
  }
  CPostFile_PushState( File, POST_S0 );
  return fd;
}

/*
 *  Close the current post mesh file
 */
int GiD_ClosePostMeshFile()
{
  int fail = 1;
  assert( G_MeshFile );
  assert( _GiDfiles_CheckState( POST_S0, G_MeshFile ) );  

#ifdef ENABLE_HDF5
  if ( G_MeshFile->m_post_mode == GiD_PostHDF5 ) {
    if ( G_MeshFile->m_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    fail = GiD_ClosePostMeshFile_HDF5( G_MeshFile->m_hdf5_file );
    delete_CurrentHdf5WriteData( G_MeshFile->m_hdf5_file );
    G_MeshFile->m_hdf5_file = NULL;
    fail = CPostFile_Release( G_MeshFile );
    G_MeshFile = NULL;
    /* reset outpuMesh */
    _GiDfiles_GetMeshFile();
    return fail;
  }
#endif  
  
  if ( G_MeshFile ) {
    fail = CPostFile_Release( G_MeshFile );
    G_MeshFile = NULL;
    /* reset outpuMesh */
    _GiDfiles_GetMeshFile( );
  }
  return fail;
}

int GiD_fClosePostMeshFile( GiD_FILE fd )
{
  int fail = 1;

  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;

    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    fail = GiD_ClosePostMeshFile_HDF5( my_hdf5_file );
    delete_CurrentHdf5WriteData( my_hdf5_file );
    my_hdf5_file = NULL;
    File->m_hdf5_file = NULL;
    fail = CPostFile_Release( File );
    GiD_HashRemove( fd );
    return fail;
  }
#endif
  
  fail = CPostFile_Release( File );
  GiD_HashRemove( fd );
  
  return fail;
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if G_MeshFile is opened write using this file
 *    in other case write using G_ResultFile
 */

int GiD_BeginMesh(GP_CONST char * MeshName, GiD_Dimension Dim,
		  GiD_ElementType EType, int NNode)
{

  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginMesh_HDF5( mesh->m_hdf5_file, MeshName, Dim,EType,NNode);
  }
#endif
  return _GiDfiles_BeginMesh(mesh, MeshName, Dim, EType, NNode);
  
}

int GiD_fBeginMesh(GiD_FILE fd, GP_CONST char * MeshName,
		   GiD_Dimension Dim, GiD_ElementType EType, int NNode)
{  
  CPostFile *mesh = NULL;
  if ( ( mesh = FD2FILE( fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( mesh->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = mesh->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginMesh_HDF5( my_hdf5_file, MeshName, Dim, EType, NNode );
  }
#endif

  return _GiDfiles_BeginMesh(mesh, MeshName, Dim, EType, NNode);
}
/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if G_MeshFile is opened write using this file
 *    in other case write using G_ResultFile. With this function you can
 *    specify a color for the mesh by its RGB components where each
 *    component take values on the interval [0,1].
 */

int GiD_BeginMeshColor(GP_CONST char * MeshName, GiD_Dimension Dim,
		       GiD_ElementType EType, int NNode,
		       double Red, double Green, double Blue)
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();

#ifdef ENABLE_HDF5
  if ( mesh->m_post_mode == GiD_PostHDF5 ){
    return GiD_BeginMeshColor_HDF5( mesh->m_hdf5_file,  MeshName, Dim, EType, NNode, Red, Green, Blue );
  }
#endif

  return _GiDfiles_BeginMeshColor( mesh, MeshName, Dim, EType, NNode,
		             Red, Green, Blue );
}

int GiD_fBeginMeshColor(GiD_FILE fd, GP_CONST char * MeshName,
		        GiD_Dimension Dim, GiD_ElementType EType,
		        int NNode,
		        double Red, double Green, double Blue)
{  
  CPostFile *mesh = NULL;
  if ( ( mesh = FD2FILE( fd)) == NULL) { return -8;}
    
#ifdef ENABLE_HDF5
  if ( mesh->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = mesh->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginMeshColor_HDF5( my_hdf5_file, MeshName, Dim, EType, NNode, Red, Green, Blue );
  }
#endif

  return _GiDfiles_BeginMeshColor(mesh, MeshName, Dim, EType, NNode,
		             Red, Green, Blue);
}

/*
 *  End current mesh
 */

int GiD_EndMesh()
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();  
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_EndMesh_HDF5( mesh->m_hdf5_file );
  }
#endif

  return _GiDfiles_EndMesh( mesh );
}

int GiD_fEndMesh(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_EndMesh_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_EndMesh( File );
}

int GiD_MeshUnit(GP_CONST char * UnitName)
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return  GiD_MeshUnit_HDF5( mesh->m_hdf5_file, UnitName);
  }
#endif
  return _GiDfiles_MeshUnit( mesh, UnitName );
  // return 1;
}


int GiD_fMeshUnit(GiD_FILE fd,GP_CONST char * UnitName)
{ 
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_MeshUnit_HDF5( my_hdf5_file, UnitName );
  }
#endif

  return _GiDfiles_MeshUnit(File,UnitName);
}

int GiD_MeshLocalAxes( GP_CONST char *result_name, GP_CONST char *analysis_name, double step_value ) {
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
    if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_MeshLocalAxes_HDF5( mesh->m_hdf5_file, result_name, analysis_name, step_value );
  }
#endif
  return _GiDFiles_MeshLocalAxes( mesh, result_name, analysis_name, step_value );
}

int GiD_fMeshLocalAxes( GiD_FILE fd, GP_CONST char *result_name, GP_CONST char *analysis_name, double step_value ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_MeshLocalAxes_HDF5( my_hdf5_file, result_name, analysis_name, step_value );
  }
#endif
  return _GiDFiles_MeshLocalAxes( File, result_name, analysis_name, step_value );
}

/*
 *  Start a coordinate block in the current mesh
 */

int GiD_BeginCoordinates()
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginCoordinates_HDF5( mesh->m_hdf5_file );
  }
#endif

  return _GiDfiles_BeginCoordinates( mesh );
}

int GiD_fBeginCoordinates(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginCoordinates_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_BeginCoordinates( File );
}

/*
 *  Close the current coordinate block
 */

int GiD_EndCoordinates()
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_EndCoordinates_HDF5( mesh->m_hdf5_file );
  }
#endif

  return _GiDfiles_EndCoordinates( mesh );
}

int GiD_fEndCoordinates(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_EndCoordinates_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_EndCoordinates(File);
}

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes within the group.
 */

int GiD_BeginMeshGroup( GP_CONST char* Name )
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginMeshGroup_HDF5( mesh->m_hdf5_file,  Name);
  }
#endif

  return _GiDfiles_BeginMeshGroup( mesh, Name );
}

int GiD_fBeginMeshGroup(GiD_FILE fd, GP_CONST char* Name)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginMeshGroup_HDF5( my_hdf5_file, Name );
  }
#endif

  return _GiDfiles_BeginMeshGroup(File, Name);
}

/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */

int GiD_EndMeshGroup()
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode == GiD_PostHDF5 ) {
    return GiD_EndMeshGroup_HDF5( mesh->m_hdf5_file );
    }
#endif
  return _GiDfiles_EndMeshGroup( mesh );
}

int GiD_fEndMeshGroup(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_EndMeshGroup_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_EndMeshGroup( File );
}


/*
 *  Write a coordinate member at the current Coordinates Block 
 */

int GiD_WriteCoordinates( int id, double x, double y, double z )
{
  int res = 0;
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    res = GiD_WriteCoordinates_HDF5( mesh->m_hdf5_file, id,x,y,z);
    return res;
  }
#endif

  return _GiDfiles_WriteCoordinates( mesh, id, x, y, z );
}

int GiD_fWriteCoordinates(GiD_FILE fd, int id, double x, double y, double z)
{
  CPostFile *File = NULL;
  // int res = 0;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteCoordinates_HDF5( my_hdf5_file, id, x, y, z );
  }
#endif

  return _GiDfiles_WriteCoordinates( File, id, x, y, z );
}

int GiD_WriteCoordinates2D(int id, double x, double y)
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    int res = GiD_WriteCoordinates2D_HDF5( mesh->m_hdf5_file, id,x,y);
    return res;
  }
#endif

  return _GiDfiles_WriteCoordinates2D( mesh, id, x, y );
}

int GiD_fWriteCoordinates2D(GiD_FILE fd, int id, double x, double y)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteCoordinates2D_HDF5( my_hdf5_file, id, x, y );
  }
#endif
    
  return _GiDfiles_WriteCoordinates2D( File, id, x, y );
}

// GiD_fWriteCoordinatesBlock includes BeginCoordinates() and EndCoordinates()
int GiD_fWriteCoordinatesBlock( GiD_FILE fd, int num_points, GP_CONST double *xyz_array ) {
  CPostFile *File = NULL;
  // int res = 0;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteCoordinatesBlock_HDF5( my_hdf5_file, num_points, xyz_array );
  }
#endif

  return _GiDfiles_WriteCoordinatesBlock( File, num_points, xyz_array);
}

// GiD_fWriteCoordinatesBlock includes BeginCoordinates() and EndCoordinates()
int GiD_fWriteCoordinatesIdBlock( GiD_FILE fd, int num_points,GP_CONST int *list_ids, GP_CONST double *xyz_array ) {
  CPostFile *File = NULL;
  // int res = 0;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteCoordinatesIdBlock_HDF5( my_hdf5_file, num_points, list_ids, xyz_array );
  }
#endif

  return _GiDfiles_WriteCoordinatesIdBlock( File, num_points, list_ids, xyz_array );
}

/*
 *  Start a elements block in the current mesh
 */

int GiD_BeginElements()
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginElements_HDF5( mesh->m_hdf5_file );
  }
#endif

  return _GiDfiles_BeginElements( mesh );
}

int GiD_fBeginElements(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginElements_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_BeginElements( File );
}

/*
 *  Close the current elements block
 */

int GiD_EndElements()
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_EndElements_HDF5( mesh->m_hdf5_file );
  }
#endif

  return _GiDfiles_EndElements( mesh );
}

int GiD_fEndElements(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
    
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_EndElements_HDF5( my_hdf5_file );
  }
#endif
    
  return _GiDfiles_EndElements(File);
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
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteElement_HDF5( mesh->m_hdf5_file, id,nid);
  }
#endif

  return _GiDfiles_WriteElement( mesh, id, nid);
}

int GiD_fWriteElement(GiD_FILE fd, int id, int nid[])
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteElement_HDF5( my_hdf5_file, id, nid );
  }
#endif

  return _GiDfiles_WriteElement(File, id, nid);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */

int GiD_WriteElementMat(int id, int nid[])
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteElementMat_HDF5( mesh->m_hdf5_file, id,nid);
  }
#endif

  return _GiDfiles_WriteElementMat( mesh, id, nid);
}

int GiD_fWriteElementMat(GiD_FILE fd, int id, int nid[])
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteElementMat_HDF5( my_hdf5_file, id, nid );
  }
#endif

  return _GiDfiles_WriteElementMat(File, id, nid);
}


// GiD_fWriteElementsBlock includes BeginElements() and EndElements()
int GiD_fWriteElementsBlock( GiD_FILE fd, int num_elements, GP_CONST int *connectivities ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd ) ) == NULL ) {    return -8;  }

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteElementsBlock_HDF5( my_hdf5_file, num_elements, connectivities );
  }
#endif

  return _GiDfiles_WriteElementsBlock( File, num_elements, connectivities );
}

int GiD_fWriteElementsIdBlock( GiD_FILE fd, int num_elements, GP_CONST int *list_ids, GP_CONST int *connectivities ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd ) ) == NULL ) {
    return -8;
  }

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteElementsIdBlock_HDF5( my_hdf5_file, num_elements, list_ids, connectivities );
  }
#endif

  return _GiDfiles_WriteElementsIdBlock( File, num_elements, list_ids, connectivities );
}

int GiD_fWriteElementsMatBlock( GiD_FILE fd, int num_elements, GP_CONST int *connectivities,
                                GP_CONST int *lst_material_id ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd ) ) == NULL ) {
    return -8;
  }

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteElementsMatBlock_HDF5( my_hdf5_file, num_elements, connectivities, lst_material_id );
  }
#endif

  return _GiDfiles_WriteElementsMatBlock( File, num_elements, connectivities, lst_material_id );
}

int GiD_fWriteElementsIdMatBlock( GiD_FILE fd, int num_elements, GP_CONST int *list_ids, GP_CONST int *connectivities,
                                  GP_CONST int *lst_material_id ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd ) ) == NULL ) {
    return -8;
  }

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteElementsIdMatBlock_HDF5( my_hdf5_file, num_elements, list_ids, connectivities, lst_material_id );
  }
#endif

  return _GiDfiles_WriteElementsIdMatBlock( File, num_elements, list_ids, connectivities, lst_material_id );
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

int GiD_WriteSphere(int id, int nid, double r)
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteSphere_HDF5( mesh->m_hdf5_file, id,nid,r);
  }
#endif

  return _GiDfiles_WriteSphere( mesh, id, nid, r);
}

int GiD_fWriteSphere(GiD_FILE fd, int id, int nid, double r)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteSphere_HDF5( my_hdf5_file, id, nid, r );
  }
#endif

  return _GiDfiles_WriteSphere(File, id, nid, r);
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

int GiD_WriteSphereMat(int id, int nid, double r, int mat)
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteSphereMat_HDF5( mesh->m_hdf5_file, id,nid,r,mat);
  }
#endif

  return _GiDfiles_WriteSphereMat( mesh, id, nid, r, mat);
}

int GiD_fWriteSphereMat(GiD_FILE fd, int id, int nid, double r, int mat)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteSphereMat_HDF5( my_hdf5_file, id, nid, r, mat );
  }
#endif

  return _GiDfiles_WriteSphereMat(File, id, nid, r, mat);
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

int GiD_WriteCircle(int id, int nid, double r,
		    double nx, double ny, double nz)
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteCircle_HDF5( mesh->m_hdf5_file, id,nid,r,nx,ny,nz);
  }
#endif

  return _GiDfiles_WriteCircle( mesh , id, nid, r, nx, ny, nz);
}

int GiD_fWriteCircle(GiD_FILE fd, int id, int nid, double r,
		     double nx, double ny, double nz)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteCircle_HDF5( my_hdf5_file, id, nid, r, nx, ny, nz );
  }
#endif

  return _GiDfiles_WriteCircle(File, id, nid, r, nx, ny, nz);
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

int GiD_WriteCircleMat(int id, int nid, double r,
		       double nx, double ny, double nz, int mat)
{
  CPostFile *mesh = _GiDfiles_GetMeshFile();
#ifdef ENABLE_HDF5
  if( mesh->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteCircleMat_HDF5( mesh->m_hdf5_file, id,nid,r,nx,ny,nz,mat);
  }
#endif
  
  return _GiDfiles_WriteCircleMat( mesh, id, nid, r, nx, ny, nz, mat);
}

int GiD_fWriteCircleMat(GiD_FILE fd, int id, int nid, double r,
		        double nx, double ny, double nz, int mat)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteCircleMat_HDF5( my_hdf5_file, id, nid, r, nx, ny, nz, mat );
  }
#endif

  return _GiDfiles_WriteCircleMat(File, id, nid, r, nx, ny, nz, mat);
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

int GiD_OpenPostResultFile( GP_CONST char * FileName, GiD_PostMode Mode )
{
  // G_PostMode=Mode;

  if ( G_ResultFile ) {
    /* must be closed */
    return  GP_ERROR_FILEOPENED;
  }

  G_ResultFile = _GiDfiles_NewFile( Mode ); // inside does also this G_ResultFile->m_post_mode = Mode;
  if ( !G_ResultFile ) {
    /* not enough memory */
    return GP_ERROR_NOMEM;
  }

  if( Mode ==GiD_PostHDF5){
#ifdef ENABLE_HDF5
    if ( G_ResultFile->m_hdf5_file != NULL ) {
      return GP_ERROR_FILEOPENED;
    }
    G_ResultFile->m_hdf5_file = new_CurrentHdf5WriteData();
    assert( G_ResultFile->m_hdf5_file );
    if ( !G_ResultFile->m_hdf5_file ) {
      return GP_ERROR_OPENFAILED;
    }
    CPostFile_PushState( G_ResultFile, POST_S0 );
    return GiD_OpenPostResultFile_HDF5( G_ResultFile->m_hdf5_file, FileName);
#else
    return GP_ERROR_INVALID_STATE;
#endif
  }

  if ( CPostFile_Open( G_ResultFile, FileName ) )
    {
    /* could not open file */
    CPostFile_Release( G_ResultFile );
    return GP_ERROR_OPENFAILED;
    }

  if ( CPostFile_WritePostHeader( G_ResultFile ) ) {
    /* WritePostHeader failed */
    GiD_ClosePostResultFile();
    return GP_ERROR_WRITESTRING;
  }
  CPostFile_PushState( G_ResultFile, POST_S0 );
  return 0;
}

GiD_FILE GiD_fOpenPostResultFile(GP_CONST char * FileName, GiD_PostMode Mode)
{
  GiD_FILE fd;
  CPostFile *File = NULL;
  File = _GiDfiles_NewFile( Mode ); // inside does also this File->m_post_mode = Mode;
  if ( !File ) {
    /* not enough memory = -2 GP_ERROR_NOMEM */
    return 0;
  }

#ifdef ENABLE_HDF5
  if ( Mode == GiD_PostHDF5 ) {

    CurrentHdf5WriteData *my_hdf5_file = new_CurrentHdf5WriteData();
    assert( new_CurrentHdf5WriteData );
    if ( !my_hdf5_file ) {
      return 0;
    }
    int fail = GiD_OpenPostResultFile_HDF5( my_hdf5_file, FileName );
    if ( fail ) {
      delete_CurrentHdf5WriteData( my_hdf5_file );
      CPostFile_Release( File );
      return 0;
    }
    File->m_hdf5_file = my_hdf5_file;

    GiD_FILE fd;
    fd = GiD_HashAdd( File );
    if ( !fd ) {
      /* could not create a file handler = -6 */
      GiD_ClosePostResultFile_HDF5( my_hdf5_file );
      delete_CurrentHdf5WriteData( my_hdf5_file );
      CPostFile_Release( File );
      return 0;
    }
    CPostFile_PushState( File, POST_S0 );
    return fd;
  }
#endif

  /* now open the File */
  if ( CPostFile_Open( File, FileName ) ) 
    {
    /* Open failed = -4 GP_ERROR_OPENFAILED */
    return 0;
    }
  fd = GiD_HashAdd( File );
  if ( !fd ) 
    {
    /* could not create a file handler = -5 GP_ERROR_HANDLEFAIL */
    CPostFile_Release( File );
    return 0;
    }

  if ( CPostFile_WritePostHeader( File ) )
    {
    /* WritePostHeader failed = -6 GP_ERROR_WRITESTRING */
    GiD_ClosePostResultFile( );
    return 0;
    } 
  else 
    {
    CPostFile_PushState( File, POST_S0 );
    }
  return fd;
}

/*
 *  Close the current post result file
 */

int GiD_ClosePostResultFile( )
{
  int fail;

  assert( G_ResultFile != NULL );
  assert( _GiDfiles_CheckState( POST_S0, G_ResultFile ) );   

#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    if ( G_ResultFile->m_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    fail = GiD_ClosePostResultFile_HDF5( G_ResultFile->m_hdf5_file );
    delete_CurrentHdf5WriteData( G_ResultFile->m_hdf5_file );
    G_ResultFile->m_hdf5_file = NULL;
    fail = CPostFile_Release( G_ResultFile );
    G_ResultFile = NULL;
    return fail;
  }
#endif 

  if ( G_ResultFile ) 
    {
    fail = CPostFile_Release(G_ResultFile);
    G_ResultFile = NULL;
    /* reset G_outputMesh pointer */
    _GiDfiles_GetMeshFile();
    return fail;
    }
  return GP_ERROR_NULLFILE;
}

int GiD_fClosePostResultFile(GiD_FILE fd)
{
  int fail = 1;

  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;

    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    fail = GiD_ClosePostResultFile_HDF5( my_hdf5_file );
    delete_CurrentHdf5WriteData( my_hdf5_file );
    my_hdf5_file = NULL;
    File->m_hdf5_file = NULL;
    assert( _GiDfiles_CheckState( POST_S0, File ) );
    fail = CPostFile_Release( File );
    GiD_HashRemove( fd );
    return fail;
  }
#endif
  
  assert( _GiDfiles_CheckState( POST_S0, File ) );    
  fail = CPostFile_Release( File );
  GiD_HashRemove( fd );
  
  return fail;
}

/*
 *  Begin Gauss Points definition
 */

int GiD_BeginGaussPoint( GP_CONST char * name, GiD_ElementType EType,
                         GP_CONST char * MeshName,
                         int GP_number, int NodesIncluded, int InternalCoord )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginGaussPoint_HDF5( G_ResultFile->m_hdf5_file, name,EType,MeshName,GP_number,NodesIncluded,InternalCoord);
  }
#endif
  
  return _GiDfiles_BeginGaussPoint( G_ResultFile,
                               name, EType, MeshName, GP_number,
                               NodesIncluded, InternalCoord );
}

int GiD_fBeginGaussPoint( GiD_FILE fd, GP_CONST char * name,
                          GiD_ElementType EType,
                          GP_CONST char * MeshName,
                          int GP_number, int NodesIncluded, int InternalCoord )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_BeginGaussPoint_HDF5( my_hdf5_file, name, EType, MeshName, GP_number, NodesIncluded, InternalCoord );
  }
#endif

  return _GiDfiles_BeginGaussPoint( File,
                               name, EType, MeshName, GP_number,
                               NodesIncluded, InternalCoord );
}


/*
 *  End current Gauss Points definition
 */

int GiD_EndGaussPoint( )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_EndGaussPoint_HDF5( G_ResultFile->m_hdf5_file );
  }
#endif

  return _GiDfiles_EndGaussPoint( G_ResultFile );
}

int GiD_fEndGaussPoint( GiD_FILE fd )
{
  CPostFile *File = NULL;
  if ( ( File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_EndGaussPoint_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_EndGaussPoint( File );
}

/*
 *  Write internal gauss point coordinate.
 */

int GiD_WriteGaussPoint2D(double x, double y)
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteGaussPoint2D_HDF5( G_ResultFile->m_hdf5_file, x,y);
  }
#endif

  return _GiDfiles_WriteGaussPoint2D( G_ResultFile, x, y );
}

int GiD_fWriteGaussPoint2D(GiD_FILE fd, double x, double y)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteGaussPoint2D_HDF5( my_hdf5_file, x, y );
  }
#endif

  return _GiDfiles_WriteGaussPoint2D( File, x, y );
}

int GiD_WriteGaussPoint3D( double x, double y, double z )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteGaussPoint3D_HDF5( G_ResultFile->m_hdf5_file, x,y,z);
  }
#endif
  
  return _GiDfiles_WriteGaussPoint3D( G_ResultFile, x, y, z );
}

int GiD_fWriteGaussPoint3D(GiD_FILE fd, double x, double y, double z)
{
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteGaussPoint3D_HDF5( my_hdf5_file, x, y, z );
  }
#endif

  return _GiDfiles_WriteGaussPoint3D( File, x, y, z );
}

/*
 *  Begin a Range Table definition
 */

int GiD_BeginRangeTable( GP_CONST char * name )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginRangeTable_HDF5( G_ResultFile->m_hdf5_file, name);
  }
#endif
  
  return _GiDfiles_BeginRangeTable( G_ResultFile, name ); 
}

int GiD_fBeginRangeTable( GiD_FILE fd, GP_CONST char * name )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_BeginRangeTable_HDF5( my_hdf5_file, name );
  }
#endif

  return _GiDfiles_BeginRangeTable( File, name );
}

/*
 *  End a Range Table definition
 */

int GiD_EndRangeTable()
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_EndRangeTable_HDF5( G_ResultFile->m_hdf5_file );
  }
#endif

  return _GiDfiles_EndRangeTable( G_ResultFile );
}

int GiD_fEndRangeTable(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_EndRangeTable_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_EndRangeTable( File );
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

int GiD_WriteMinRange( double max, GP_CONST char * name )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteMinRange_HDF5( G_ResultFile->m_hdf5_file, max,name);
  }
#endif

  return _GiDfiles_WriteMinRange( G_ResultFile, max, name );
}

int GiD_fWriteMinRange( GiD_FILE fd, double max, GP_CONST char * name )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
    
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteMinRange_HDF5( my_hdf5_file, max, name );
  }
#endif

  return _GiDfiles_WriteMinRange( File, max, name );
}

int GiD_WriteRange( double min, double max, GP_CONST char * name )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteRange_HDF5( G_ResultFile->m_hdf5_file, min,max,name);
  }
#endif
  
  return _GiDfiles_WriteRange( G_ResultFile, min, max, name );
}

int GiD_fWriteRange( GiD_FILE fd, double min, double max, GP_CONST char * name )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteRange_HDF5( my_hdf5_file, min, max, name );
  }
#endif


  return _GiDfiles_WriteRange(File, min, max, name);
}

int GiD_WriteMaxRange( double min, GP_CONST char * name )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteMaxRange_HDF5( G_ResultFile->m_hdf5_file, min,name);
  }
#endif

  return _GiDfiles_WriteMaxRange(G_ResultFile, min, name);
}

int GiD_fWriteMaxRange(GiD_FILE fd, double min, GP_CONST char * name)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteMaxRange_HDF5( my_hdf5_file, min, name );
  }
#endif

  
  return _GiDfiles_WriteMaxRange(File, min, name);
}

/*
 *  Begin Result Block
 */

int GiD_BeginResult( GP_CONST char     *Result,
                     GP_CONST char     *Analysis,
                     double             step,
                     GiD_ResultType     Type,
                     GiD_ResultLocation Where,
                     GP_CONST char     *GaussPointsName,
                     GP_CONST char     *RangeTable, 
                     int                compc,
                     GP_CONST char     *compv[] )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginResult_HDF5( G_ResultFile->m_hdf5_file, Result,Analysis,step,Type,Where,GaussPointsName,RangeTable,compc,compv);
  }
#endif
  
  return _GiDfiles_BeginResult( G_ResultFile, Result, Analysis, step, Type, Where,
                           GaussPointsName, RangeTable, compc, compv );
}

int GiD_fBeginResult( GiD_FILE fd, GP_CONST char * Result,
                      GP_CONST char * Analysis,
                      double step,
                      GiD_ResultType Type, GiD_ResultLocation Where,
                      GP_CONST char * GaussPointsName,
                      GP_CONST char * RangeTable, 
		     int compc, GP_CONST char * compv[] )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginResult_HDF5( my_hdf5_file, Result, Analysis, step, Type, Where, GaussPointsName, RangeTable, compc, compv );
  }
#endif

  return _GiDfiles_BeginResult( File, Result, Analysis, step, Type, Where,
                           GaussPointsName, RangeTable, compc, compv );
}

int GiD_BeginResultHeader( GP_CONST char     *Result,
                           GP_CONST char     *Analysis,
                           double             step,
                           GiD_ResultType     Type,
                           GiD_ResultLocation Where,
                           GP_CONST char     *GaussPointsName )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginResultHeader_HDF5( G_ResultFile->m_hdf5_file, Result,Analysis,step,Type,Where,GaussPointsName);
  }
#endif
  
  return _GiDfiles_BeginResultHeader( G_ResultFile, Result, Analysis, step, Type,
                                 Where, GaussPointsName );
}

int GiD_fBeginResultHeader( GiD_FILE fd, GP_CONST char * Result,
                            GP_CONST char * Analysis, 
                            double step,
                            GiD_ResultType Type, GiD_ResultLocation Where,
                            GP_CONST char * GaussPointsName )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginResultHeader_HDF5( my_hdf5_file, Result, Analysis, step, Type, Where, GaussPointsName );
  }
#endif

  return _GiDfiles_BeginResultHeader( File, Result, Analysis, step, Type,
                                 Where, GaussPointsName );
}

int GiD_ResultRange( GP_CONST char * RangeTable )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_ResultRange_HDF5( G_ResultFile->m_hdf5_file, RangeTable);
  }
#endif

  return _GiDfiles_ResultRange( G_ResultFile, RangeTable );
}

int GiD_fResultRange( GiD_FILE fd, GP_CONST char * RangeTable )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_ResultRange_HDF5( my_hdf5_file, RangeTable );
  }
#endif

  return _GiDfiles_ResultRange( File, RangeTable );
}

int GiD_ResultComponents( int compc, GP_CONST char * compv[] )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_ResultComponents_HDF5( G_ResultFile->m_hdf5_file, compc,compv);
  }
#endif

  return _GiDfiles_ResultComponents( G_ResultFile, compc, compv );
}

int GiD_fResultComponents( GiD_FILE fd, int compc, GP_CONST char * compv[] )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_ResultComponents_HDF5( my_hdf5_file, compc, compv );
  }
#endif

  return _GiDfiles_ResultComponents( File, compc, compv );
}

int GiD_ResultUnit( GP_CONST char * UnitName )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_ResultUnit_HDF5( G_ResultFile->m_hdf5_file, UnitName);
  }
#endif
  return _GiDfiles_ResultUnit( G_ResultFile, UnitName );  
}

int GiD_fResultUnit(GiD_FILE fd, GP_CONST char * UnitName)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_ResultUnit_HDF5( my_hdf5_file, UnitName );
  }
#endif

  return _GiDfiles_ResultUnit( File, UnitName );  
}

int GiD_ResultUserDefined( GP_CONST char * Name,GP_CONST char * Value )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_ResultUserDefined_HDF5( G_ResultFile->m_hdf5_file, Name,Value);
  }
#endif
  return _GiDfiles_ResultUserDefined( G_ResultFile, Name, Value );
}

int GiD_fResultUserDefined( GiD_FILE fd, GP_CONST char * Name,GP_CONST char * Value )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_ResultUserDefined_HDF5( my_hdf5_file, Name, Value );
  }
#endif

  return _GiDfiles_ResultUserDefined( File, Name, Value );
}

int GiD_BeginResultGroup( GP_CONST char     *Analysis,
                          double             step,
                          GiD_ResultLocation Where,
                          GP_CONST char     *GaussPointsName )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginResultGroup_HDF5( G_ResultFile->m_hdf5_file, Analysis,step,Where,GaussPointsName);
  }
#endif

  return _GiDfiles_BeginResultGroup( G_ResultFile, Analysis, step, Where,
                                GaussPointsName );
}

int GiD_fBeginResultGroup( GiD_FILE fd, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
                           GP_CONST char * GaussPointsName )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginResultGroup_HDF5( my_hdf5_file, Analysis, step, Where, GaussPointsName );
  }
#endif

  return _GiDfiles_BeginResultGroup( File, Analysis, step, Where,
                                GaussPointsName );
}

int GiD_ResultDescription( GP_CONST char * Result, GiD_ResultType Type )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_ResultDescription_HDF5( G_ResultFile->m_hdf5_file, Result,Type);
  }
#endif

  return _GiDfiles_ResultDescription_( G_ResultFile, Result, Type, 0 );
}

int GiD_fResultDescription( GiD_FILE fd,
		            GP_CONST char * Result, GiD_ResultType Type )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_ResultDescription_HDF5( my_hdf5_file, Result, Type );
  }
#endif

  return _GiDfiles_ResultDescription_( File, Result, Type, 0 );
}

int GiD_ResultDescriptionDim( GP_CONST char * Result, GiD_ResultType Type,
                              size_t s )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    /* disregard s in hdf5, because ResultGroup is really written in hdf5 as independent Results*/
    return GiD_ResultDescription_HDF5( G_ResultFile->m_hdf5_file, Result,Type);
  }
#endif

  return _GiDfiles_ResultDescription_( G_ResultFile, Result, Type, s );
}

int GiD_fResultDescriptionDim( GiD_FILE fd,
                               GP_CONST char * Result, GiD_ResultType Type,
                               size_t s )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    /* disregard s in hdf5, because ResultGroup is really written in hdf5 as independent Results*/
    return GiD_ResultDescription_HDF5( my_hdf5_file, Result, Type );
  }
#endif

  return _GiDfiles_ResultDescription_( File, Result, Type, s );
}

int GiD_ResultLocalAxes( GP_CONST char *Result, GP_CONST char *Analysis, double step, double vx, double vy,
                         double vz ) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5 ) {
    return GiD_ResultLocalAxes_HDF5( G_ResultFile->m_hdf5_file, Result, Analysis, step, vx, vy, vz );
  }
#endif
  return _GiDFiles_ResultLocalAxes( G_ResultFile, Result, Analysis, step, vx, vy, vz );
}

int GiD_fResultLocalAxes( GiD_FILE fd, GP_CONST char *Result, GP_CONST char *Analysis, double step, double vx, double vy,
                         double vz ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd ) ) == NULL ) {
    return -8;
  }
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    return GiD_ResultLocalAxes_HDF5( File->m_hdf5_file, Result, Analysis, step, vx, vy, vz );
  }
#endif
  return _GiDFiles_ResultLocalAxes( File, Result, Analysis, step, vx, vy, vz );
}

int GiD_ResultIsDeformationVector( int is_deformation_vector ) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5 ) {
    return GiD_ResultIsDeformationVector_HDF5( G_ResultFile->m_hdf5_file, is_deformation_vector );
  }
#endif
  return _GiDFiles_ResultIsDeformationVector( G_ResultFile, is_deformation_vector );
}
int GiD_fResultIsDeformationVector( GiD_FILE fd, int is_deformation_vector ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd ) ) == NULL ) {
    return -8;
  }
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    return GiD_ResultIsDeformationVector_HDF5( File->m_hdf5_file, is_deformation_vector );
  }
#endif
  return _GiDFiles_ResultIsDeformationVector( File, is_deformation_vector );
}

int GiD_ResultValues()
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_ResultValues_HDF5( G_ResultFile->m_hdf5_file );
  }
#endif

  return GP_OK;
}

int GiD_fResultValues( GiD_FILE fd )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_ResultValues_HDF5( my_hdf5_file );
  }
#endif

  return GP_OK;
}

/*
 *  End Result Block
 */

int GiD_EndResult( )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_EndResult_HDF5( G_ResultFile->m_hdf5_file );
  }
#endif

  return _GiDfiles_EndResult( G_ResultFile );
}

int GiD_fEndResult( GiD_FILE fd )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_EndResult_HDF5( my_hdf5_file );
  }
#endif

  return _GiDfiles_EndResult( File );
}


// GiD_fWriteResultBlock includes BeginResult()/BeginValues() and EndValues()/EndResult()
GIDPOST_API
int GiD_fWriteResultBlock( GiD_FILE fd, GP_CONST char *result_name, GP_CONST char *analysis_name, double step_value,
                           GiD_ResultType result_type, GiD_ResultLocation result_location, GP_CONST char *gauss_point_name,
                           GP_CONST char *range_table_name, int num_component_names, GP_CONST char **list_component_names,
                           GP_CONST char *unit_name,
			               int num_result_values, GP_CONST int *list_result_ids, 
                           int num_component_values, GP_CONST double *list_component_values ) {
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd ) ) == NULL ) {
    return -8;
  }

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteResultBlock_HDF5( my_hdf5_file, result_name, analysis_name, step_value, result_type, result_location,
                                      gauss_point_name, range_table_name, num_component_names, list_component_names,
                                      unit_name,
                                      num_result_values, list_result_ids, num_component_values, list_component_values );
  }
#endif

  return _GiD_WriteResultBlock( File, result_name, analysis_name, step_value, result_type, result_location, gauss_point_name,
                                range_table_name, num_component_names, list_component_names, unit_name,
                                num_result_values, list_result_ids, num_component_values, list_component_values );
}

/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */

int GiD_BeginOnMeshGroup( char * Name )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_BeginOnMeshGroup_HDF5( G_ResultFile->m_hdf5_file,  Name);
  }
#endif

  return _GiDfiles_BeginOnMeshGroup( G_ResultFile, Name );
}

int GiD_fBeginOnMeshGroup(GiD_FILE fd, char * Name)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_BeginOnMeshGroup_HDF5( my_hdf5_file, Name );
  }
#endif

  return _GiDfiles_BeginOnMeshGroup( File, Name );
}

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */

int GiD_EndOnMeshGroup()
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_EndOnMeshGroup_HDF5( G_ResultFile->m_hdf5_file );
  }
#endif
  return _GiDfiles_EndOnMeshGroup( G_ResultFile );
}

int GiD_fEndOnMeshGroup(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_EndOnMeshGroup_HDF5( my_hdf5_file );
  }
#endif
  
  return _GiDfiles_EndOnMeshGroup( File );
}

/*
 * Flushes all pending output into the compressed file.
 */

int GiD_FlushPostFile() {  
  int res1=(G_MeshFile)?CPostFile_Flush(G_MeshFile):0;
  int res2=(G_ResultFile)?CPostFile_Flush(G_ResultFile):0;
  int res=(res1||res2)?1:0;
  return res;
}

int GiD_fFlushPostFile(GiD_FILE fd)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    if ( my_hdf5_file == NULL ) {
      return GP_ERROR_NULLFILE; // or may be GP_ERROR_HANDLEFAIL
    }
    return GiD_FlushPostFile_HDF5( my_hdf5_file );
  }
#endif

  return CPostFile_Flush( File );
}

/*
 *  Write result functions
 */

// static
// int GiD_EnsureBeginValues( CPostFile *File )
// {
// 
//   assert( File );
//   
//   if ( !File->flag_begin_values ) 
//     {
//     if ( !CPostFile_BeginValues( File ) ) 
//       {
//       // post_state st = CPostFile_TopState( File );
//       
//       CPostFile_PushState( File, POST_RESULT_VALUES );
//       if ( File->flag_isgroup )
//         {
// 	CPostFile_ResultGroupOnBeginValues( File );
//         }
//       File->flag_begin_values = 1;
//       return 0;
//       }
//     }
//   else
//     {
//     assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
//     }
//   return 1;
// }

int GiD_WriteScalar(int id, double v)
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteScalar_HDF5( G_ResultFile->m_hdf5_file, id,v);
  }
#endif

  return _GiDfiles_WriteScalar( G_ResultFile, id, v );
}

int GiD_fWriteScalar(GiD_FILE fd, int id, double v)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteScalar_HDF5( my_hdf5_file, id, v );
  }
#endif

  
  return _GiDfiles_WriteScalar( File, id, v );
}

int GiD_Write2DVector( int id, double x, double y )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_Write2DVector_HDF5( G_ResultFile->m_hdf5_file,  id, x, y );
  }
#endif

  return _GiDfiles_Write2DVector( G_ResultFile, id, x, y );
}

int GiD_fWrite2DVector( GiD_FILE fd, int id, double x, double y )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_Write2DVector_HDF5( my_hdf5_file, id, x, y );
  }
#endif

  return _GiDfiles_Write2DVector( File, id, x, y );
}

int GiD_WriteVector( int id, double x, double y, double z )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteVector_HDF5( G_ResultFile->m_hdf5_file, id,x,y,z);
  }
#endif
  
  return _GiDfiles_WriteVector( G_ResultFile, id, x, y, z );
}

int GiD_fWriteVector( GiD_FILE fd, int id, double x, double y, double z )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteVector_HDF5( my_hdf5_file, id, x, y, z );
  }
#endif

  return _GiDfiles_WriteVector( File, id, x, y, z );
}

int GiD_WriteVectorModule(int id, double x, double y, double z, double mod)
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteVectorModule_HDF5( G_ResultFile->m_hdf5_file, id,x,y,z,mod);
  }
#endif

  return _GiDfiles_WriteVectorModule(G_ResultFile, id, x, y, z, mod);
}

int GiD_fWriteVectorModule(GiD_FILE fd,
		           int id, double x, double y, double z, double mod)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteVectorModule_HDF5( my_hdf5_file, id, x, y, z, mod );
  }
#endif

  return _GiDfiles_WriteVectorModule(File, id, x, y, z, mod);
}

int GiD_Write2DMatrix(int id, double Sxx, double Syy, double Sxy)
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_Write2DMatrix_HDF5( G_ResultFile->m_hdf5_file, id,Sxx,Syy,Sxy);
  }
#endif

  return _GiDfiles_Write2DMatrix(G_ResultFile, id, Sxx, Syy, Sxy);
}

int GiD_fWrite2DMatrix(GiD_FILE fd, int id, double Sxx, double Syy, double Sxy)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_Write2DMatrix_HDF5( my_hdf5_file, id, Sxx, Syy, Sxy );
  }
#endif

  return _GiDfiles_Write2DMatrix(File, id, Sxx, Syy, Sxy);
}

int GiD_Write3DMatrix(int id,
		      double Sxx, double Syy, double Szz,
		      double Sxy, double Syz, double Sxz)
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_Write3DMatrix_HDF5( G_ResultFile->m_hdf5_file, id,Sxx,Syy,Szz,Sxy,Syz,Sxz);
  }
#endif

  return _GiDfiles_Write3DMatrix(G_ResultFile, id,
		            Sxx, Syy, Szz,
		            Sxy, Syz, Sxz);
}

int GiD_fWrite3DMatrix(GiD_FILE fd, int id,
		       double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_Write3DMatrix_HDF5( my_hdf5_file, id, Sxx, Syy, Szz, Sxy, Syz, Sxz );
  }
#endif

  return _GiDfiles_Write3DMatrix(File, id,
		            Sxx, Syy, Szz,
		            Sxy, Syz, Sxz);
}

int GiD_WritePlainDefMatrix(int id,
		            double Sxx, double Syy, double Sxy, double Szz )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WritePlainDefMatrix_HDF5( G_ResultFile->m_hdf5_file, id,Sxx,Syy,Sxy,Szz);
  }
#endif
  
  return _GiDfiles_WritePlainDefMatrix(G_ResultFile,
		                  id, Sxx, Syy, Sxy, Szz);
}

int GiD_fWritePlainDefMatrix(GiD_FILE fd, int id,
		             double Sxx, double Syy, double Sxy, double Szz )
{
  CPostFile *File = NULL;
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WritePlainDefMatrix_HDF5( my_hdf5_file, id, Sxx, Syy, Sxy, Szz );
  }
#endif

  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
  return _GiDfiles_WritePlainDefMatrix(File,
		                  id, Sxx, Syy, Sxy, Szz);
}

int GiD_WriteMainMatrix(int id,
		        double Si, double Sii, double Siii,
		        double Vix, double Viy, double Viz,
		        double Viix, double Viiy, double Viiz,
		        double Viiix, double Viiiy, double Viiiz)
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteMainMatrix_HDF5( G_ResultFile->m_hdf5_file, id,Si,Sii,Siii,Vix,Viy,Viz,Viix,Viiy,Viiz,Viiix,Viiiy,Viiiz);
  }
#endif

  return _GiDfiles_WriteMainMatrix(G_ResultFile, id, 
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
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteMainMatrix_HDF5( my_hdf5_file, id, Si, Sii, Siii, Vix, Viy, Viz, Viix, Viiy, Viiz, Viiix, Viiiy, Viiiz );
  }
#endif

  return _GiDfiles_WriteMainMatrix(File, id, 
		              Si,  Sii,  Siii,
		              Vix,  Viy,  Viz,
		              Viix,  Viiy,  Viiz,
		              Viiix,  Viiiy,  Viiiz);
}

int GiD_WriteLocalAxes(int id, double euler_1, double euler_2, double euler_3)
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteLocalAxes_HDF5( G_ResultFile->m_hdf5_file, id,euler_1,euler_2,euler_3);
  }
#endif

  return _GiDfiles_WriteLocalAxes(G_ResultFile, id, euler_1, euler_2, euler_3);
}

int GiD_fWriteLocalAxes(GiD_FILE fd,
		        int id, double euler_1, double euler_2, double euler_3)
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteLocalAxes_HDF5( my_hdf5_file, id, euler_1, euler_2, euler_3 );
  }
#endif

  
  return _GiDfiles_WriteLocalAxes(File, id, euler_1, euler_2, euler_3);
}

/*
 * Complex numbers
 */

int GiD_WriteComplexScalar( int id, double complex_real, double complex_imag) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5) {
    return GiD_WriteComplexScalar_HDF5( G_ResultFile->m_hdf5_file, id, complex_real, complex_imag);
  }
#endif
  return _GiDfiles_WriteComplexScalar(G_ResultFile, id, complex_real, complex_imag);
}

int GiD_fWriteComplexScalar(GiD_FILE fd, int id, double complex_real, double complex_imag) {
  CPostFile *File = NULL;
  if ( (  File = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteComplexScalar_HDF5( my_hdf5_file, id, complex_real, complex_imag );
  }
#endif

  return _GiDfiles_WriteComplexScalar(File, id, complex_real, complex_imag);
}

int GiD_Write2DComplexVector( int id, 
			      double x_real, double x_imag,
			      double y_real, double y_imag) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5) {
    return GiD_Write2DComplexVector_HDF5( G_ResultFile->m_hdf5_file, id, x_real, x_imag, y_real, y_imag);
  }
#endif
  return _GiDfiles_Write2DComplexVector(G_ResultFile, id, x_real, x_imag, y_real, y_imag);
}

int GiD_fWrite2DComplexVector( GiD_FILE fd, int id,
			       double x_real, double x_imag,
			       double y_real, double y_imag) {
  CPostFile *File = NULL;  
  if ( (  File = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_Write2DComplexVector_HDF5( my_hdf5_file, id, x_real, x_imag, y_real, y_imag );
  }
#endif

  return _GiDfiles_Write2DComplexVector(File, id, x_real, x_imag, y_real, y_imag);
}

int GiD_WriteComplexVector( int id, 
			    double x_real, double x_imag,
			    double y_real, double y_imag,
			    double z_real, double z_imag) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5) {
    return GiD_WriteComplexVector_HDF5( G_ResultFile->m_hdf5_file, id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
  }
#endif
  return _GiDfiles_WriteComplexVector(G_ResultFile, id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
}

int GiD_fWriteComplexVector( GiD_FILE fd, int id,
			     double x_real, double x_imag,
			     double y_real, double y_imag,
			     double z_real, double z_imag) {
  CPostFile *File = NULL;  
  if ( (  File = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteComplexVector_HDF5( my_hdf5_file, id, x_real, x_imag, y_real, y_imag, z_real, z_imag );
  }
#endif

  return _GiDfiles_WriteComplexVector(File, id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
}


int GiD_Write2DComplexMatrix(int id, double Sxx_real, double Syy_real, double Sxy_real,
                             double Sxx_imag, double Syy_imag, double Sxy_imag) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5) {
    return GiD_Write2DComplexMatrix_HDF5( G_ResultFile->m_hdf5_file, id, Sxx_real, Syy_real, Sxy_real,
                                     Sxx_imag, Syy_imag, Sxy_imag);
  }
#endif
  return _GiDfiles_Write2DComplexMatrix(G_ResultFile, id, Sxx_real, Syy_real, Sxy_real,
                                        Sxx_imag, Syy_imag, Sxy_imag);
}

int GiD_fWrite2DComplexMatrix(GiD_FILE fd, int id,
                              double Sxx_real, double Syy_real, double Sxy_real,
                              double Sxx_imag, double Syy_imag, double Sxy_imag) {
  CPostFile *File = NULL;  
  if ( (  File = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_Write2DComplexMatrix_HDF5( my_hdf5_file, id, Sxx_real, Syy_real, Sxy_real,
                                          Sxx_imag, Syy_imag, Sxy_imag );
  }
#endif

  return _GiDfiles_Write2DComplexMatrix(File, id, Sxx_real, Syy_real, Sxy_real,
                                        Sxx_imag, Syy_imag, Sxy_imag);
}

int GiD_Write3DComplexMatrix(int id,
                             double Sxx_real, double Syy_real, double Szz_real,
                             double Sxy_real, double Syz_real, double Sxz_real,
                             double Sxx_imag, double Syy_imag, double Szz_imag,
                             double Sxy_imag, double Syz_imag, double Sxz_imag) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5) {
    return GiD_WriteComplexMatrix_HDF5( G_ResultFile->m_hdf5_file, id,
                                       Sxx_real, Syy_real, Szz_real,
                                       Sxy_real, Syz_real, Sxz_real,
                                       Sxx_imag, Syy_imag, Szz_imag,
                                       Sxy_imag, Syz_imag, Sxz_imag);
  }
#endif
  return _GiDfiles_WriteComplexMatrix(G_ResultFile, id, 
                                      Sxx_real, Syy_real, Szz_real,
                                      Sxy_real, Syz_real, Sxz_real,
                                      Sxx_imag, Syy_imag, Szz_imag,
                                      Sxy_imag, Syz_imag, Sxz_imag);
}

int GiD_fWrite3DComplexMatrix(GiD_FILE fd, int id,
                              double Sxx_real, double Syy_real, double Szz_real,
                              double Sxy_real, double Syz_real, double Sxz_real,
                              double Sxx_imag, double Syy_imag, double Szz_imag,
                              double Sxy_imag, double Syz_imag, double Sxz_imag) {
  CPostFile *File = NULL;  
  if ( (  File = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteComplexMatrix_HDF5( my_hdf5_file, id,
                                        Sxx_real, Syy_real, Szz_real,
                                        Sxy_real, Syz_real, Sxz_real,
                                        Sxx_imag, Syy_imag, Szz_imag,
                                        Sxy_imag, Syz_imag, Sxz_imag );
  }
#endif

  return _GiDfiles_WriteComplexMatrix(File, id, 
                                      Sxx_real, Syy_real, Szz_real,
                                      Sxy_real, Syz_real, Sxz_real,
                                      Sxx_imag, Syy_imag, Szz_imag,
                                      Sxy_imag, Syz_imag, Sxz_imag);
}


/* 
* Nurbs 
*/

int GiD_WriteNurbsSurface( int id, int n, double* v )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteNurbsSurface_HDF5( G_ResultFile->m_hdf5_file, id,n,v);
  }
#endif
  
  return _GiDfiles_WriteNurbsSurface( G_ResultFile, id, n, v );
}

int GiD_fWriteNurbsSurface( GiD_FILE fd, int id, int n, double* v )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteNurbsSurface_HDF5( my_hdf5_file, id, n, v );
  }
#endif

  return _GiDfiles_WriteNurbsSurface( File, id, n, v );
}

int GiD_WriteNurbsSurfaceVector( int id, int n, int num_comp, double* v )
{
#ifdef ENABLE_HDF5
  if( G_ResultFile->m_post_mode ==GiD_PostHDF5){
    return GiD_WriteNurbsSurfaceVector_HDF5( G_ResultFile->m_hdf5_file, id,n,num_comp,v);
  }
#endif
  
  return _GiDfiles_WriteNurbsSurfaceVector( G_ResultFile, id, n, num_comp, v );
}

int GiD_fWriteNurbsSurfaceVector( GiD_FILE fd, int id, int n, int num_comp, double* v )
{
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}
  
#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteNurbsSurfaceVector_HDF5( my_hdf5_file, id, n, num_comp, v );
  }
#endif

  return _GiDfiles_WriteNurbsSurfaceVector( File, id, n, num_comp, v );
}

/* User defined properties/atributes */
/* User defined properties defined inside Mesh or Result blocks
   hdf5: stored as properties/attributes (Name, value) of the current Mesh/N or Result/N folder
   ASCII / raw binary: stored as comments
     # Name: value
   Define the macro COMPASSIS_USER_ATTRIBUTES_FORMAT
   to have it like compassis wants:
     # ResultUserDefined \"%s\" \"%s\"      or
     # ResultUserDefined \"%s\" %s
*/

int GiD_WriteMeshUserAttribute( GP_CONST char *Name, GP_CONST char *Value ) {
  CPostFile * meshFile = _GiDfiles_GetMeshFile();

#ifdef ENABLE_HDF5
  if ( meshFile->m_post_mode == GiD_PostHDF5 ) {
    return GiD_WriteMeshUserAttribute_HDF5( meshFile->m_hdf5_file,  Name, Value );
  }
#endif

  return _GiDFiles_WriteMeshUserAttribute( meshFile, Name, Value );
}

int GiD_fWriteMeshUserAttribute( GiD_FILE fd, GP_CONST char *Name, GP_CONST char *Value ) {
  CPostFile *meshFile = NULL;
  if ( (  meshFile  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( meshFile->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = meshFile->m_hdf5_file;
    return GiD_WriteMeshUserAttribute_HDF5( my_hdf5_file, Name, Value );
  }
#endif

  return _GiDFiles_WriteMeshUserAttribute( meshFile, Name, Value );
  // return 1;
}

int GiD_WriteResultUserAttribute( GP_CONST char *Name, GP_CONST char *Value ) {
#ifdef ENABLE_HDF5
  if ( G_ResultFile->m_post_mode == GiD_PostHDF5 ) {
    return GiD_WriteResultUserAttribute_HDF5( G_ResultFile->m_hdf5_file,  Name, Value );
  }
#endif
  return _GiDFiles_WriteResultUserAttribute( G_ResultFile, Name, Value );
}

int GiD_fWriteResultUserAttribute( GiD_FILE fd, GP_CONST char *Name, GP_CONST char *Value ) {
  CPostFile *File = NULL;
  if ( (  File  = FD2FILE(  fd)) == NULL) { return -8;}

#ifdef ENABLE_HDF5
  if ( File->m_post_mode == GiD_PostHDF5 ) {
    CurrentHdf5WriteData *my_hdf5_file = File->m_hdf5_file;
    return GiD_WriteResultUserAttribute_HDF5( my_hdf5_file, Name, Value );
  }
#endif

  return _GiDFiles_WriteResultUserAttribute( File, Name, Value );
}

#ifndef WIN32
#pragma GCC diagnostic pop
#endif // WIN32
