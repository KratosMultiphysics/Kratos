/* gidpost 2.1 */
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

GP_CONST char * GetResultTypeName(GiD_ResultType type, size_t s);
void GetResultTypeMinMaxValues(GiD_ResultType type, size_t *min, size_t *max);

typedef struct _CPostFile CPostFile;

typedef struct _CBufferValues CBufferValues;

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

  post_state level_mesh;
  post_state level_res;
  
  int (*ptr_Open)            (CPostFile* this, GP_CONST char *name);
  int (*ptr_Close)           (CPostFile* this);
  int (*ptr_Flush)           (CPostFile* this);
  int (*ptr_IsBinary)        (CPostFile* this);
  int (*ptr_WriteString)     (CPostFile* this, GP_CONST char * str);
  int (*ptr_BeginCoordinates)(CPostFile* this);
  int (*ptr_BeginElements)   (CPostFile* this);
  int (*ptr_BeginValues)     (CPostFile* this);
  int (*ptr_EndValues)       (CPostFile* this);
  int (*ptr_WriteInteger)    (CPostFile* this, int i,    int op);
  int (*ptr_WriteDouble)     (CPostFile* this, double x, int op);
  int (*ptr_WriteValuesVA)   (CPostFile* this, int id,   int n, va_list ap);
  int (*ptr_WriteValues)     (CPostFile* this, int id,   int n, double* );
  int (*ptr_Write2D)         (CPostFile* this, double x, double y);
  int (*ptr_Write3D)         (CPostFile* this, double x, double y, double z);
  int (*ptr_WriteElement)    (CPostFile* this, int id,   int n, int nid[]);
  int (*ptr_WritePostHeader) (CPostFile* this);
};

int CPostFile_Release         (CPostFile*);
int CPostFile_Open            (CPostFile* this, GP_CONST char * name);
int CPostFile_Close           (CPostFile* this);
int CPostFile_Flush           (CPostFile* this);
int CPostFile_IsBinary        (CPostFile* this);
int CPostFile_WriteString     (CPostFile* this, GP_CONST char * str);
int CPostFile_BeginCoordinates(CPostFile* this);
int CPostFile_BeginElements   (CPostFile* this);
int CPostFile_BeginValues     (CPostFile* this);
int CPostFile_EndValues       (CPostFile* this);
int CPostFile_WriteInteger    (CPostFile* this, int i,    int op);
int CPostFile_WriteDouble     (CPostFile* this, double x, int op);
int CPostFile_WriteValuesVA   (CPostFile* this, int id,   int n, ...);
int CPostFile_WriteValues     (CPostFile* this, int id,   int n, double*);
int CPostFile_Write2D         (CPostFile* this, double x, double y);
int CPostFile_Write3D         (CPostFile* this, double x, double y, double z);
int CPostFile_WriteElement    (CPostFile* this, int id,   int n, int nid[]);
int CPostFile_WritePostHeader (CPostFile* this);

void CPostFile_ResetLastID     (CPostFile* this);
void CPostFile_SetConnectivity (CPostFile* this, int nnode);
int CPostFile_GetConnectivity  (CPostFile* this);
int CPostFile_MatchConnectivity(CPostFile* this, int written);

void CPostFile_ResultGroupOnBegin       (CPostFile* this);
void CPostFile_ResultGroupOnNewType     (CPostFile* this, GiD_ResultType t);
int  CPostFile_ResultGroupIsEmpty       (CPostFile* this);
void CPostFile_ResultGroupOnBeginValues (CPostFile* this);
int  CPostFile_ResultGroupWriteValues   (CPostFile* this, GiD_ResultType t, int id, int n, ...);

CPostFile *CPostAscii_Create();
CPostFile *CPostAsciiZ_Create();
CPostFile *CPostBinary_Create();

#endif
