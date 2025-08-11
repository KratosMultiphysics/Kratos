/* -*- mode: c++ -*-
 *
 *  gidpostInt.h --
 *
 *    This file declare some internal clases used in the
 *    implementation of the interface gidpost.h
 *
 */

#ifndef __GIDPOSTINT__
#define __GIDPOSTINT__

#include <stdarg.h>
#include <stdio.h>
#include "zlib.h"

#include "gidpost.h"

#ifdef _MSC_VER
#if _MSC_VER >= 1400
/* this is Visual C++ 2005 */
#define strdup _strdup
#elif _MSC_VER >= 1310
/* this is Visual c++ .NET 2003 */
#define strdup _strdup
#elif _MSC_VER > 1300
/* this is Visual C++ .NET 2002 */
#define strdup _strdup
#endif
#endif

#define LINE_SIZE 8192

GP_CONST char * GetResultTypeName(GiD_ResultType type, int dim);
void GetResultTypeMinMaxValues(GiD_ResultType type, int *min, int *max);

typedef struct _CPostFile CPostFile;

typedef struct _CBufferValues CBufferValues;

typedef struct _CurrentHdf5WriteData CurrentHdf5WriteData;

typedef enum {
  POST_UNDEFINED,
  POST_S0,               /* TOP level */
  POST_MESH_S0,          /* MESH header */
  POST_MESHGROUP_S0,     /* MESHGROUP TOP level */
  POST_MESH_COORD_INSIDE, /* inside a Coordinate block */
  POST_MESH_COORD_DONE,   /* after a Coordinate block but inside a MESH */
  POST_MESH_ELEM,        /* inside an Element block */
  POST_GAUSS_S0,         /* GAUSS point block: implicit */
  POST_GAUSS_GIVEN,      /* GAUSS point block: explicit */
  POST_RANGE_S0,         /* RANGE table block */
  POST_RESULT_ONGROUP,   /* OnGroup result */
  POST_RESULT_DEPRECATED,   /* Deprecated result block */
  POST_RESULT_SINGLE,    /* Result block */
  POST_RESULT_GROUP,     /* Result group block */
  POST_RESULT_DESC,      /* Result description block */
  POST_RESULT_VALUES     /* writing values */
} post_state;

#define STACK_STATE_SIZE 10
#define LOCAL_AXES_MAX_LEN 500

struct _CPostFile
{
  int m_LastID;
  int m_connectivity;
  void *m_FILE;
  int m_fail;
  CBufferValues *buffer_values;
  int GP_number_check;
  int gauss_written;
  int flag_isgroup;
  int flag_begin_values;
  int has_mesh;
  int has_meshgroup;
  char local_axes_format[ LOCAL_AXES_MAX_LEN ];

  post_state level_mesh;
  post_state level_res;

  post_state stack_state[STACK_STATE_SIZE];
  int stack_pos;
  
  int (*ptr_Open)            (CPostFile* post_file, GP_CONST char *name);
  int (*ptr_Close)           (CPostFile* post_file);
  int (*ptr_Flush)           (CPostFile* post_file);
  int (*ptr_IsBinary)        (CPostFile* post_file);
  int (*ptr_WriteString)     (CPostFile* post_file, GP_CONST char * str);
  int (*ptr_BeginCoordinates)(CPostFile* post_file);
  int (*ptr_BeginElements)   (CPostFile* post_file);
  int (*ptr_BeginValues)     (CPostFile* post_file);
  int (*ptr_EndValues)       (CPostFile* post_file);
  int (*ptr_WriteInteger)    (CPostFile* post_file, int i,    int op);
  int (*ptr_WriteDouble)     (CPostFile* post_file, double x, int op);
  int (*ptr_WriteValuesVA)   (CPostFile* post_file, int id,   int n, va_list ap);
  int ( *ptr_WriteValues )( CPostFile* post_file, int id, int n, GP_CONST double * );
  int ( *ptr_WriteValuesNS )( CPostFile* post_file, int id, int n, GP_CONST double * );
  int (*ptr_WriteValuesNSV)  (CPostFile* post_file, int id,   int n, int num_comp, GP_CONST double* );
  int (*ptr_Write2D)         (CPostFile* post_file, double x, double y);
  int (*ptr_Write3D)         (CPostFile* post_file, double x, double y, double z);
  int ( *ptr_WriteElement )( CPostFile* post_file, int id, int n, GP_CONST int nid[] );
  int (*ptr_WritePostHeader) (CPostFile* post_file);
  int (*ptr_WritePostHeaderIGA) (CPostFile* post_file);

  // if m_post_mode == GiD_PostHDF5 then CurrentHdf5WriteData is valid and above not
  /* TODO add here G_PostMode and CurrentHdf5WriteData *my_hdf5_file */
  GiD_PostMode m_post_mode;
  GiD_PostEncoding m_encoding_filename; //to open the file converting filename from utf-8 to local or not
  CurrentHdf5WriteData *m_hdf5_file;
};

int CPostFile_Release         (CPostFile* post_file);
int CPostFile_Open            (CPostFile* post_file, GP_CONST char * name);
int CPostFile_Close           (CPostFile* post_file);
post_state CPostFile_TopState ( CPostFile* post_file );
int CPostFile_PushState       ( CPostFile* post_file, post_state s );
post_state CPostFile_PopState (CPostFile* post_file);
int CPostFile_Flush           (CPostFile* post_file);
int CPostFile_IsBinary        (CPostFile* post_file);
int CPostFile_WriteString     (CPostFile* post_file, GP_CONST char * str);
int CPostFile_BeginCoordinates(CPostFile* post_file);
int CPostFile_BeginElements   (CPostFile* post_file);
int CPostFile_BeginValues     (CPostFile* post_file);
int CPostFile_EndValues       (CPostFile* post_file);
int CPostFile_WriteInteger    (CPostFile* post_file, int i,    int op);
int CPostFile_WriteDouble     (CPostFile* post_file, double x, int op);
int CPostFile_WriteValuesVA   (CPostFile* post_file, int id,   int n, ...);
int CPostFile_WriteValues( CPostFile* post_file, int id, int n, GP_CONST double * );
int CPostFile_WriteValuesNS( CPostFile* post_file, int id, int n, GP_CONST double * );
int CPostFile_WriteValuesNSV( CPostFile* post_file, int id, int n, int num_comp, GP_CONST double * );
int CPostFile_Write2D         (CPostFile* post_file, double x, double y);
int CPostFile_Write3D         (CPostFile* post_file, double x, double y, double z);
int CPostFile_WriteElement( CPostFile* post_file, int id, int n, GP_CONST int nid[] );
int CPostFile_WritePostHeader (CPostFile* post_file);

void CPostFile_ResetLastID     (CPostFile* post_file);
void CPostFile_SetConnectivity (CPostFile* post_file, int nnode);
int CPostFile_GetConnectivity  (CPostFile* post_file);
int CPostFile_MatchConnectivity(CPostFile* post_file, int written);

void CPostFile_ResultGroupOnBegin       (CPostFile* post_file);
void CPostFile_ResultGroupOnNewType     (CPostFile* post_file, GiD_ResultType t);
int  CPostFile_ResultGroupIsEmpty       (CPostFile* post_file);
void CPostFile_ResultGroupOnBeginValues (CPostFile* post_file);
int  CPostFile_ResultGroupWriteValues   (CPostFile* post_file, GiD_ResultType t, int id, int n, ...);

CPostFile* CPostAscii_Create( void );
CPostFile* CPostAsciiZ_Create( void );
CPostFile* CPostBinary_Create( void );
CPostFile* CPostHdf5_Create( void );

#endif
