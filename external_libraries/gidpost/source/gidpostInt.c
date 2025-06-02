/* gidpost */
/*
 *  gidpostInt.cc --
 *
 *    This file implemente the internal class declared in gidpostInt.h
 *
 */

#include <assert.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "gidpostInt.h"
#ifdef _WIN32
#include <windows.h>
#endif
#include "gidpostHDF5.h"  // for GiD_FlushPostFile_HDF5

static int ByteOrderCheck = 0x91d;  /* Magic number */

/* format to print reals to string, customizable only once at the begin*/
static char format_real[100]={"%.9g"};

GIDPOST_API int GiD_PostSetFormatReal(GP_CONST char * f){
  strcpy(format_real,f);
  return 0;
}

GIDPOST_API GP_CONST char* GiD_PostGetFormatReal() {
  return format_real;
}

/* special format to not truncate time steps converted to string, not customizable but centralize its use*/
GIDPOST_API GP_CONST char* GiD_PostGetFormatStep(){
  static const char format_step[]={"%.16g"};
  return format_step;
}

struct _CBufferValues
{
  double * buffer_values;
  int last_value;
  int values_size_min;
  int values_size_max;
  int buffer_values_size;
  GiD_ResultType *buffer_types;
  int buffer_types_size;
  int size_types;
  int next_type;
};

CBufferValues* CBufferValues_Create( void );
void CBufferValues_Destroy(CBufferValues* this);

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
  {"LocalAxes",               {3,  0, 0, 0}},
  {"ComplexScalar",           {2,  0, 0, 0}},
  {"ComplexVector",           {4,  6, 0, 0}},
  {"ComplexMatrix",           {6, 12, 0, 0}}
  /* Vector2 */
  /* Vector3 */
};

GP_CONST char * GetResultTypeName(GiD_ResultType type, size_t s)
{
  static char buffer[255];
  char * ptr;
  int i;
  
  strcpy(buffer, _ResultTypeInfo[(int)type].str);
  ptr = &(buffer[0]) + strlen(buffer);
  if (s) {
    for (i = 0; _ResultTypeInfo[(int)type].dims[i]; i++) {
      if (s == _ResultTypeInfo[(int)type].dims[i])
        break;
    }
    if (!(_ResultTypeInfo[(int)type].dims[i]))
      printf("Invalid dimension %u for type %s\n", ( unsigned int)s, buffer);
    else
      sprintf(ptr, ":%u", ( unsigned int)s);
  }
  return buffer;
}

void GetResultTypeMinMaxValues(GiD_ResultType type, size_t *min, size_t *max)
{
  int i = 0;
  
  *min = _ResultTypeInfo[(int)type].dims[0];
  for (i = 1; _ResultTypeInfo[(int)type].dims[i]; i++)
    ;
  *max = _ResultTypeInfo[(int)type].dims[i-1];
}

/* ---------------------------------------------------------------------------
 *
 *  Post Files classes implementation : gidpostInt.h
 *
 * ---------------------------------------------------------------------------
 */

static
CPostFile *CPostFile_Create( void )
{
  CPostFile *this = (CPostFile*)malloc(sizeof(CPostFile));
  assert(this);
  this->m_LastID = -1;
  this->m_FILE = NULL;
  this->m_fail = 0;
  this->buffer_values = CBufferValues_Create();
  this->GP_number_check = 0;
  this->gauss_written = 0;
  this->flag_isgroup = 0;
  this->flag_begin_values = 0;
  this->has_mesh = 0;
  this->has_meshgroup = 0;
  this->level_mesh = POST_UNDEFINED;
  this->level_res  = POST_UNDEFINED;
  this->stack_pos = -1;
  this->local_axes_format[ 0 ] = '\0';

  this->ptr_Open             = NULL;
  this->ptr_Close            = NULL;
  this->ptr_Flush            = NULL;
  this->ptr_IsBinary         = NULL;
  this->ptr_WriteString      = NULL;
  this->ptr_BeginCoordinates = NULL;
  this->ptr_BeginElements    = NULL;
  this->ptr_BeginValues      = NULL;
  this->ptr_EndValues        = NULL;
  this->ptr_WriteInteger     = NULL;
  this->ptr_WriteDouble      = NULL;
  this->ptr_WriteValuesVA    = NULL;
  this->ptr_WriteValues      = NULL;  
  this->ptr_WriteValuesNS    = NULL;    
  this->ptr_WriteValuesNSV   = NULL;    
  this->ptr_Write2D          = NULL;
  this->ptr_Write3D          = NULL;
  this->ptr_WriteElement     = NULL;
  this->ptr_WritePostHeader  = NULL;
  this->ptr_WritePostHeaderIGA=NULL;

  this->m_post_mode = GiD_PostUndefined;
  this->m_hdf5_file = NULL;

  return this;
}

int CPostFile_Release(CPostFile* this)
{
  int ret;
  
  assert(this);
  
  ret = CPostFile_Close(this);
  if (this->buffer_values)  {
    CBufferValues_Destroy(this->buffer_values);
    this->buffer_values = NULL;
  }
  free(this);
  return ret;
}

int CPostFile_Open(CPostFile* this, GP_CONST char * str)
{
  assert(this);
  assert(this->ptr_Open);
  return (*this->ptr_Open)(this, str);
}

int CPostFile_Close(CPostFile* this)
{
  assert(this);
  if ( this->m_post_mode == GiD_PostHDF5 ) {
    // already closed before with GiD_ClosePostResultFile_HDF5()
    return 0;
  }
  assert(this->ptr_Close);
  return (*this->ptr_Close)(this);
}

post_state CPostFile_TopState( CPostFile* this )
{
  return ( this->stack_pos < 0 ) ? POST_UNDEFINED : this->stack_state[this->stack_pos];
}

int CPostFile_PushState( CPostFile* this, post_state s )
{
  assert( this->stack_pos < STACK_STATE_SIZE - 1 );
  if ( this->stack_pos < STACK_STATE_SIZE - 1 )
    {
      this->stack_state[++this->stack_pos] = s;
      return 0;
    }
  else
    {
      return -1;
    }
}

post_state CPostFile_PopState(CPostFile* this)
{
  post_state top;
  if ( this->stack_pos < 0)
    {
      return POST_UNDEFINED;
    }
  top = this->stack_state[this->stack_pos--];
  return top;
}

int CPostFile_Flush(CPostFile* this)
{
  int res=0;
  assert(this);  
#ifdef ENABLE_HDF5
  if(this->m_post_mode==GiD_PostHDF5) {
    res=GiD_FlushPostFile_HDF5(this->m_hdf5_file);
  } else {
    assert(this->ptr_Flush);
    res=(*this->ptr_Flush)(this);
  }
#else //ENABLE_HDF5
  assert(this->ptr_Flush);
  res=(*this->ptr_Flush)(this);
#endif //ENABLE_HDF5
  return res;
}

int CPostFile_IsBinary(CPostFile* this)
{
  assert(this);
  assert(this->ptr_IsBinary);
  return (*this->ptr_IsBinary)(this);
}

int CPostFile_WriteString(CPostFile* this, GP_CONST char * str)
{
  assert(this);
  assert(this->ptr_WriteString);
  return (*this->ptr_WriteString)(this, str);
}

int CPostFile_BeginCoordinates(CPostFile* this)
{
  assert(this);
  assert(this->ptr_BeginCoordinates);
  
  return (*this->ptr_BeginCoordinates)(this);
}

int CPostFile_BeginElements(CPostFile* this)
{
  assert(this);
  assert(this->ptr_BeginElements);
  
  return (*this->ptr_BeginElements)(this);
}

int CPostFile_BeginValues(CPostFile* this)
{
  assert(this);
  assert(this->ptr_BeginValues);
  
  return (*this->ptr_BeginValues)(this);
}

int CPostFile_EndValues(CPostFile* this)
{
  assert(this);
  assert(this->ptr_EndValues);
  
  return (*this->ptr_EndValues)(this);
}

int CPostFile_WriteInteger(CPostFile* this, int i, int op)
{
  assert(this);
  assert(this->ptr_WriteInteger);
  
  return (*this->ptr_WriteInteger)(this, i, op);
}

int CPostFile_WriteDouble(CPostFile* this, double x, int op)
{
  assert(this);
  assert(this->ptr_WriteDouble);
  
  return (*this->ptr_WriteDouble)(this, x, op);
}

int CPostFile_WriteValuesVA(CPostFile* this, int id, int num_comp, ...)
{
  va_list ap;
  int ret;
  
  assert(this);
  assert(this->ptr_WriteValuesVA);

  va_start(ap, num_comp);
  ret = (*this->ptr_WriteValuesVA)(this, id, num_comp, ap);
  va_end(ap);
  
  return ret;
}

int CPostFile_WriteValues( CPostFile *this, int id, int n, GP_CONST double *values ) {
  assert(this);
  assert(this->ptr_WriteValues);
  
  return (*this->ptr_WriteValues)(this, id, n, values);
}

int CPostFile_WriteValuesNS( CPostFile *this, int id, int n, GP_CONST double *values ) {
  assert(this);
  assert(this->ptr_WriteValuesNS);
  
  return (*this->ptr_WriteValuesNS)(this, id, n, values);
}

int CPostFile_WriteValuesNSV( CPostFile *this, int id, int n, int num_comp, GP_CONST double *values ) {
  assert(this);
  assert(this->ptr_WriteValuesNSV);
  
  return (*this->ptr_WriteValuesNSV)(this, id, n, num_comp, values);
}

int CPostFile_Write2D(CPostFile *this, double x, double y )
{
  assert(this);
  assert(this->ptr_Write2D);
  
  return (*this->ptr_Write2D)(this, x, y);
}

int CPostFile_Write3D(CPostFile *this, double x, double y, double z )
{
  assert(this);
  assert(this->ptr_Write3D);
  
  return (*this->ptr_Write3D)(this, x, y, z);
}

int CPostFile_WriteElement( CPostFile *this, int id, int n, GP_CONST int nid[] ) {
  assert(this);
  assert(this->ptr_WriteElement);
  
  return (*this->ptr_WriteElement)(this, id, n, nid);
}

int CPostFile_WritePostHeader(CPostFile *this)
{
  assert(this);
  assert(this->ptr_WritePostHeader);
  
  return (*this->ptr_WritePostHeader)(this);
}

void CPostFile_ResetLastID(CPostFile* this)
{
  assert(this);
  
  this->m_LastID = -1;
}

void CPostFile_SetConnectivity(CPostFile* this, int nnode)
{
  assert(this);
  
  this->m_connectivity = nnode;
}

int CPostFile_GetConnectivity(CPostFile* this)
{
  assert(this);
  
  return this->m_connectivity;
}

int CPostFile_MatchConnectivity(CPostFile* this, int written)
{
  int conn = CPostFile_GetConnectivity(this);
  return (written == conn) ? 0 : (((written-1) == conn) ? 1 : 2);
}

/*
 *  class CPostAscii
 */

static
void CPostAscii_Init(CPostFile*);

CPostFile* CPostAscii_Create()
{
  CPostFile *this;

  this = CPostFile_Create();
  CPostAscii_Init(this);
  return this;
}

static int CPostAscii_Close(CPostFile *this);

static
int CPostAscii_Open(CPostFile *this, GP_CONST char *name )
{
#ifdef _WIN32
  wchar_t wname[MAX_PATH];
#endif
  assert(this);
  assert(this->ptr_Close);
  
  CPostAscii_Close(this);
#ifdef _WIN32  
  /* convert from utf-8 */
  MultiByteToWideChar(CP_UTF8,0,name,-1,wname,MAX_PATH);
  this->m_FILE=_wfopen(wname,L"w");  
#else
  this->m_FILE = fopen(name, "w");
#endif
  return this->m_FILE == NULL;
}

static
int CPostAscii_Close(CPostFile *this)
{
  assert(this);

  if ( this->m_FILE ) {
    this->m_fail = fclose(this->m_FILE);
    this->m_FILE = NULL;
  } else
    this->m_fail = 1;
  return this->m_fail;
}

static
int CPostAscii_Flush(CPostFile* this)
{
  return this->m_FILE ? fflush(this->m_FILE) : 1;
}

static
int CPostAscii_IsBinary(CPostFile* this)
{
  return 0;
}

static
int CPostAscii_WriteString(CPostFile* this, GP_CONST char * str)
{
  assert(this);
  
  fprintf(this->m_FILE, "%s\n", str );
  return 0;
}

static
int CPostAscii_BeginCoordinates(CPostFile* this)
{
  return CPostAscii_WriteString(this, "Coordinates");
}

static
int CPostAscii_BeginElements(CPostFile* this)
{
  return CPostAscii_WriteString(this, "Elements");
}

static
int CPostAscii_BeginValues(CPostFile* this)
{
  return CPostAscii_WriteString(this, "Values");
}

static
int CPostAscii_EndValues(CPostFile* this)
{
  return CPostAscii_WriteString(this, "End Values");
}

static
int CPostAscii_WriteInteger(CPostFile* this, int i, int op)
{
  assert(this);

  if ( op == 0 ) {
    fprintf( this->m_FILE, "%d", i );
  } else {
    fprintf( this->m_FILE, " %d", i );
    if ( op == 2 ) {
      fprintf( this->m_FILE, "\n" );
    }
  }
  return 0;
}

static
int CPostAscii_WriteDouble(CPostFile* this, double x, int op)
{
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  assert(op!=0);
  fprintf(this->m_FILE,local_format, x);
  if (op==2) {
    fprintf(this->m_FILE, "\n");
  }
  return 0;
}

static
int CPostAscii_WriteValuesVA(CPostFile* this, int id, int num_comp, va_list ap)
{
  int i;
  double value;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  assert(this->m_FILE);
  
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(this->m_FILE, "%d", id);
  
  for (i = 0; i < num_comp; i++) {
    value = va_arg(ap, double);
    fprintf(this->m_FILE, local_format, value);
  }
  fprintf(this->m_FILE, "\n");
  this->m_LastID = id;
  
  return 0;
}

static int CPostAscii_WriteValues( CPostFile *this, int id, int n, GP_CONST double *buffer ) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(this->m_FILE, "%d", id);
  for ( i = 0; i < n; i++ )
    fprintf(this->m_FILE, local_format, buffer[i]);
  fprintf(this->m_FILE, "\n");
  this->m_LastID = id;
  return 0;
}

static int CPostAscii_WriteValuesNS( CPostFile *this, int id, int n, GP_CONST double *buffer ) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s\n",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  
  //write surface Id
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(this->m_FILE, "%d\n", id);

  //write values
  for ( i = 0; i < n; i++ ) {
    if ( buffer[i] == GP_UNKNOWN) fprintf(this->m_FILE, "NR\n");
    else fprintf(this->m_FILE, local_format, buffer[i]);
  }
  this->m_LastID = id;
  return 0;
}

static int CPostAscii_WriteValuesNSV( CPostFile *this, int id, int n, int num_comp, GP_CONST double *buffer ) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s ",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  
  //write surface Id
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(this->m_FILE, "%d\n", id);

  //write values
  for ( i = 0; i < n*num_comp; i++ ) {
    if ( buffer[i] == GP_UNKNOWN) fprintf(this->m_FILE, "NR ");
    else fprintf(this->m_FILE,local_format,buffer[i]);
    if ( (i+1) % num_comp == 0) fprintf(this->m_FILE, "\n");
  }
  this->m_LastID = id;
  return 0;
}

static
int CPostAscii_Write2D(CPostFile *this, double x, double y )
{
  static char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format, x, y);
  return CPostAscii_WriteString(this, line);
}

static
int CPostAscii_Write3D(CPostFile *this, double x, double y, double z )
{
  static char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format,x,y,z);
  return CPostAscii_WriteString(this, line);
}

static int CPostAscii_WriteElement( CPostFile *this, int id, int n, GP_CONST int nid[] ) {
  int i;

  assert(this);

  fprintf(this->m_FILE, "%d", id);
  for (i = 0; i < n; i++) {
    fprintf(this->m_FILE, " %d", nid[i]);
  }
  fprintf(this->m_FILE, "\n");
  
  return 0;
}

static
int CPostAscii_WritePostHeader(CPostFile *this)
{
  return CPostAscii_WriteString(this, "GiD Post Results File 1.2");
}

static
void CPostAscii_Init(CPostFile* this)
{
  assert(this);
  
  this->ptr_Open             = &CPostAscii_Open;
  this->ptr_Close            = &CPostAscii_Close;
  this->ptr_Flush            = &CPostAscii_Flush;
  this->ptr_IsBinary         = &CPostAscii_IsBinary; 
  this->ptr_WriteString      = &CPostAscii_WriteString;
  this->ptr_BeginCoordinates = &CPostAscii_BeginCoordinates;
  this->ptr_BeginElements    = &CPostAscii_BeginElements;
  this->ptr_BeginValues      = &CPostAscii_BeginValues;
  this->ptr_EndValues        = &CPostAscii_EndValues;
  this->ptr_WriteInteger     = &CPostAscii_WriteInteger;
  this->ptr_WriteDouble      = &CPostAscii_WriteDouble;
  this->ptr_WriteValuesVA    = &CPostAscii_WriteValuesVA;
  this->ptr_WriteValues      = &CPostAscii_WriteValues;
  this->ptr_WriteValuesNS    = &CPostAscii_WriteValuesNS;
  this->ptr_WriteValuesNSV   = &CPostAscii_WriteValuesNSV;
  this->ptr_Write2D          = &CPostAscii_Write2D;
  this->ptr_Write3D          = &CPostAscii_Write3D;
  this->ptr_WriteElement     = &CPostAscii_WriteElement;
  this->ptr_WritePostHeader  = &CPostAscii_WritePostHeader;
}

/*
 *  class CPostAsciiZ --
 */

static
void CPostAsciiZ_Init(CPostFile*);

CPostFile* CPostAsciiZ_Create()
{
  CPostFile *this;

  this = CPostFile_Create();
  CPostAsciiZ_Init(this);
  return this;
}

static int CPostAsciiZ_Close(CPostFile *this);

static
int CPostAsciiZ_Open(CPostFile *this, GP_CONST char * name)
{

#ifdef _WIN32
  wchar_t wname[MAX_PATH];
#endif
  assert(this);
  assert(this->ptr_Close);
  
  CPostAsciiZ_Close(this);
#ifdef _WIN32  
  /* convert from utf-8 */
  MultiByteToWideChar(CP_UTF8,0,name,-1,wname,MAX_PATH);  
  this->m_FILE = gzopen_w(wname, "w1");  
#else  
  this->m_FILE = gzopen(name, "w1");  
#endif
  return this->m_FILE == NULL;
}

static
int CPostAsciiZ_Close(CPostFile *this)
{
  assert(this);
  
  if ( this->m_FILE ) {
    this->m_fail = gzclose(this->m_FILE);
    this->m_FILE = NULL;
  } else
    this->m_fail = 1;
  return this->m_fail;
}

static
int CPostAsciiZ_Flush(CPostFile *this)
{
  assert(this);
  
  return this->m_FILE ? gzflush(this->m_FILE, Z_FINISH) : 1;
}

static
int CPostAsciiZ_IsBinary(CPostFile* this)
{
  return 0;
}

static
int CPostAsciiZ_WriteString(CPostFile* this, GP_CONST char * str )
{
  assert(this);
  
  gzprintf(this->m_FILE, "%s\n", str );
  return 0;
}

static
int CPostAsciiZ_BeginCoordinates(CPostFile* this)
{
  return CPostAsciiZ_WriteString(this, "Coordinates");
}

static
int CPostAsciiZ_BeginElements(CPostFile* this)
{
  return CPostAsciiZ_WriteString(this, "Elements");
}

static
int CPostAsciiZ_BeginValues(CPostFile* this)
{
  return CPostAsciiZ_WriteString(this, "Values");
}

static
int CPostAsciiZ_EndValues(CPostFile* this)
{
  return CPostAsciiZ_WriteString(this, "End Values");
}

static
int CPostAsciiZ_WriteInteger(CPostFile* this, int i, int op)
{
  assert(this);
    
  if (op==1) {
    gzprintf(this->m_FILE, " ");
    gzprintf(this->m_FILE, "%d", i);
  } else {
    gzprintf(this->m_FILE, "%d", i);
    if (op==2) {
      gzprintf(this->m_FILE, "\n");
    }
  }
  
  return 0;
}

static
int CPostAsciiZ_WriteDouble(CPostFile* this, double x, int op)
{
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  assert(op!=0);

  gzprintf(this->m_FILE,local_format,x);
  if (op==2) {
    gzprintf(this->m_FILE, "\n");
  }

  return 0;
}

static
int CPostAsciiZ_WriteValuesVA(CPostFile* this, int id, int num_comp, va_list ap)
{
  int i;
  double value;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzprintf(this->m_FILE, "%d", id);
  for (i = 0; i < num_comp; i++) {
    value = va_arg(ap, double);
    gzprintf(this->m_FILE,local_format,value);
  }
  gzprintf(this->m_FILE, "\n");
  this->m_LastID = id;

  return 0;
}

static int CPostAsciiZ_WriteValues( CPostFile *this, int id, int n, GP_CONST double *buffer ) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzprintf(this->m_FILE, "%d", id);
  for ( i = 0; i < n; i++ )
    gzprintf(this->m_FILE,local_format,buffer[i]);
  gzprintf(this->m_FILE, "\n");
  this->m_LastID = id;

  return 0;
}

static int CPostAsciiZ_WriteValuesNS( CPostFile *this, int id, int n, GP_CONST double *buffer ) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s\n",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  
  //write surface Id
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(this->m_FILE, "%d\n", id);

  //write values
  for ( i = 0; i < n; i++ ) {
    if ( buffer[i] == GP_UNKNOWN) fprintf(this->m_FILE, "NR\n");
    else fprintf(this->m_FILE, local_format,buffer[i]);
  }
  this->m_LastID = id;
  return 0;
}

static int CPostAsciiZ_WriteValuesNSV( CPostFile *this, int id, int n, int num_comp, GP_CONST double *buffer ) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s ",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(this);
  
  //write surface Id
  if (this->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(this->m_FILE, "%d\n", id);

  //write values
  for ( i = 0; i < n*num_comp; i++ ) {
    if ( buffer[i] == GP_UNKNOWN) fprintf(this->m_FILE, "NR ");
    else fprintf(this->m_FILE,local_format,buffer[i]);
    if ( (i+1) % num_comp == 0) fprintf(this->m_FILE, "\n");
  }
  this->m_LastID = id;
  return 0;
}

static
int CPostAsciiZ_Write2D(CPostFile *this, double x, double y )
{
  static char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format,x,y);
  return CPostAsciiZ_WriteString(this, line);
}

static
int CPostAsciiZ_Write3D(CPostFile *this, double x, double y, double z )
{
  char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format,x,y,z);
  return CPostAsciiZ_WriteString(this, line);
}

static int CPostAsciiZ_WriteElement( CPostFile *this, int id, int n, GP_CONST int nid[] ) {
  int i;

  assert(this);

  gzprintf(this->m_FILE, "%d", id);
  for (i = 0; i < n; i++) {
    gzprintf(this->m_FILE, " %d", nid[i]);
  }
  gzprintf(this->m_FILE, "\n");

  return 0;
}
 
static
int CPostAsciiZ_WritePostHeader(CPostFile *this)
{
  return CPostAsciiZ_WriteString(this, "GiD Post Results File 1.2");
}


static
void CPostAsciiZ_Init(CPostFile* this)
{
  assert(this);
  
  this->ptr_Open             = &CPostAsciiZ_Open;
  this->ptr_Close            = &CPostAsciiZ_Close;
  this->ptr_Flush            = &CPostAsciiZ_Flush;
  this->ptr_IsBinary         = &CPostAsciiZ_IsBinary; 
  this->ptr_WriteString      = &CPostAsciiZ_WriteString;
  this->ptr_BeginCoordinates = &CPostAsciiZ_BeginCoordinates;
  this->ptr_BeginElements    = &CPostAsciiZ_BeginElements;
  this->ptr_BeginValues      = &CPostAsciiZ_BeginValues;
  this->ptr_EndValues        = &CPostAsciiZ_EndValues;
  this->ptr_Write2D          = &CPostAsciiZ_Write2D;
  this->ptr_Write3D          = &CPostAsciiZ_Write3D;
  this->ptr_WriteInteger     = &CPostAsciiZ_WriteInteger;
  this->ptr_WriteDouble      = &CPostAsciiZ_WriteDouble;
  this->ptr_WriteValuesVA    = &CPostAsciiZ_WriteValuesVA;
  this->ptr_WriteValues      = &CPostAsciiZ_WriteValues;
  this->ptr_WriteValuesNS    = &CPostAsciiZ_WriteValuesNS;
  this->ptr_WriteValuesNSV   = &CPostAsciiZ_WriteValuesNSV;
  this->ptr_WriteElement     = &CPostAsciiZ_WriteElement;
  this->ptr_WritePostHeader  = &CPostAsciiZ_WritePostHeader;
}

/*
 *  class CPostBinary --
 */

static
void CPostBinary_Init(CPostFile*);

CPostFile* CPostBinary_Create()
{
  CPostFile *this;

  this = CPostFile_Create();
  CPostBinary_Init(this);
  return this;
}

CPostFile *CPostHdf5_Create() {
  return CPostFile_Create();
}

static int CPostBinary_Close(CPostFile *this);

static
int CPostBinary_Open(CPostFile *this, GP_CONST char * name)
{

#ifdef _WIN32
  wchar_t wname[MAX_PATH];
#endif  
  assert(this);
  assert(this->ptr_Close);
  
  CPostBinary_Close(this);
#ifdef _WIN32  
  /* convert from utf-8 */
  MultiByteToWideChar(CP_UTF8,0,name,-1,wname,MAX_PATH);  
  this->m_FILE = gzopen_w(wname, "wb1");  
#else  
  this->m_FILE = gzopen(name, "wb1");  
#endif
  
  /* escribir el numero magico */
  if (this->m_FILE) {
    gzwrite(this->m_FILE, &ByteOrderCheck, sizeof(ByteOrderCheck));
  }

  return this->m_FILE == NULL;
}

static
int CPostBinary_Close(CPostFile *this)
{
  assert(this);
  
  if ( this->m_FILE ) {
    this->m_fail = gzclose(this->m_FILE);
    this->m_FILE = NULL;
  } else
    this->m_fail = 1;
  
  return this->m_fail;
}

static
int CPostBinary_Flush(CPostFile *this)
{
  assert(this);
  
  return this->m_FILE ? gzflush(this->m_FILE, Z_FULL_FLUSH) : 1;
}

static
int CPostBinary_IsBinary(CPostFile* this)
{
  return 1;
}

static
int CPostBinary_WriteString(CPostFile *this, GP_CONST char * str)
{
  int size, written = 0, tam;

  assert(this);

  this->m_fail = 1;
  if (this->m_FILE) {
    GP_CONST char *buf = "\0";
    
    if (str) buf = str;
    tam = (int)strlen(buf) + 1; /* incluido el \0 */
    written = gzwrite(this->m_FILE, &tam, sizeof(int));
    size = ( int)sizeof(char) * tam;
    written += gzwrite(this->m_FILE, (void*)buf, ( unsigned int)size);
    if ( written == size+(int)(sizeof(int)) )
      this->m_fail = 0;
  }
  return this->m_fail;
}

static
int CPostBinary_BeginCoordinates(CPostFile* this)
{
  return CPostBinary_WriteString(this, "Coordinates -1 Indexed"); 
}

static
int CPostBinary_BeginElements(CPostFile* this)
{
  return CPostBinary_WriteString(this, "Elements -1 Indexed"); 
}

static
int CPostBinary_BeginValues(CPostFile* this)
{
  return CPostBinary_WriteString(this, "Values -1 Indexed");
}

static
int CPostBinary_EndValues(CPostFile* this)
{
  int idxend = -1;

  assert(this);
  assert(this->m_FILE);
  
  if ( gzwrite(this->m_FILE, &idxend, sizeof(int)) != sizeof(int) )
    return 1;
  return CPostBinary_WriteString(this, "End Values");
}

static
int CPostBinary_WriteInteger(CPostFile* this, int i, int op)
{
  assert(this);
  assert(this->m_FILE);
  
  gzwrite(this->m_FILE, &i, sizeof(i));
  return 0;
}

static
int CPostBinary_WriteDouble(CPostFile* this, double x, int op)
{
  float v = (float)(x);
  
  assert(this);
  assert(this->m_FILE);
  
  gzwrite(this->m_FILE, &v, sizeof(v)); 
  return 0;
}

static
int CPostBinary_WriteValuesVA(CPostFile* this, int id, int num_comp, va_list ap)
{
  int i;
  float value;

  assert(this);
  assert(this->m_FILE);

  if ( this->m_LastID != id ) { /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzwrite(this->m_FILE, &id, sizeof(id));
    this->m_LastID = id;
  }
  for ( i = 0; i < num_comp; i++ ) {
    value = (float)va_arg(ap, double);
    gzwrite(this->m_FILE, &value, sizeof(value));
  }
  return 0;  
} 

static int CPostBinary_WriteValues( CPostFile *this, int id, int n, GP_CONST double *buffer ) {
  int i;
  float value;

  assert(this);
  assert(this->m_FILE);

  if ( this->m_LastID != id ) { /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzwrite(this->m_FILE, &id, sizeof(id));
    this->m_LastID = id;
  }
  for ( i = 0; i < n; i++ ) {
    value = (float) buffer[i];
    gzwrite(this->m_FILE, &value, sizeof(value));
  }
  return 0;  
} 

static int CPostBinary_WriteValuesNS( CPostFile *this, int id, int n, GP_CONST double *buffer ) {
  int i;
  float value;
  assert(this);
  //write surface Id
  gzwrite(this->m_FILE, &id, sizeof(id));
  //write amount of values (else without the geometry is not possible to read it)
  gzwrite(this->m_FILE, &n, sizeof(n));
  //write values
  for ( i = 0; i < n; i++ ) {
    value = (float) buffer[i];
    gzwrite(this->m_FILE, &value, sizeof(value));
  }
  return 0;
}

static int CPostBinary_WriteValuesNSV( CPostFile *this, int id, int n, int num_comp, GP_CONST double *buffer ) {
  int i;
  float value;
  int num_values=n*num_comp;
  assert(this);
  
  //write surface Id
  gzwrite(this->m_FILE, &id, sizeof(id));
  //write amount of values (else without the geometry is not possible to read it)
  gzwrite(this->m_FILE, &num_values, sizeof(num_values));
  //write values
  for ( i = 0; i < num_values; i++ ) {
    value = (float) buffer[i];
    gzwrite(this->m_FILE, &value, sizeof(value));
  }
  this->m_LastID = id;
  return 0;
} 

static
int CPostBinary_Write2D(CPostFile* this, double x, double y )
{
  float values[2];

  assert(this);
  assert(this->m_FILE);
  
  values[0] = (float)(x);
  values[1] = (float)(y);
  gzwrite(this->m_FILE, values, sizeof(float)*2);
  
  return 0;
}

static
int CPostBinary_Write3D(CPostFile* this, double x, double y, double z )
{
  float values[3];
  
  assert(this);
  assert(this->m_FILE);
  
  values[0] = (float)(x);
  values[1] = (float)(y);
  values[2] = (float)(z);
  gzwrite(this->m_FILE, values, sizeof(float)*3);
  
  return 0;
}

static int CPostBinary_WriteElement( CPostFile *this, int id, int n, GP_CONST int nid[] ) {
  int i;
  
  assert(this);
  assert(this->m_FILE);
  
  gzwrite(this->m_FILE, &id, sizeof(id));
  for ( i = 0; i < n; i++ ) {
    gzwrite(this->m_FILE, nid+i, sizeof(int));
  }
  /* verify connectivity */
  switch ( CPostFile_MatchConnectivity(this, n) ) {
    case 0:
      /* match exactly */
      i = 0;
      /* so write material 0 == no material */
      gzwrite(this->m_FILE, &i, sizeof(int));
    case 1:
      /* match connectivity and material
       * the material was already written */
      break;
    case 2:
      /* connectivity MISSMATCH */
      return 1;
  }
  return 0;
} 

static
int CPostBinary_WritePostHeader(CPostFile* this)
{
  return CPostBinary_WriteString(this, "GiDPostEx1.2");
}


static
void CPostBinary_Init(CPostFile* this)
{
  assert(this);
  
  this->ptr_Open             = &CPostBinary_Open;
  this->ptr_Close            = &CPostBinary_Close;
  this->ptr_Flush            = &CPostBinary_Flush;
  this->ptr_IsBinary         = &CPostBinary_IsBinary; 
  this->ptr_WriteString      = &CPostBinary_WriteString;
  this->ptr_BeginCoordinates = &CPostBinary_BeginCoordinates;
  this->ptr_BeginElements    = &CPostBinary_BeginElements;
  this->ptr_BeginValues      = &CPostBinary_BeginValues;
  this->ptr_EndValues        = &CPostBinary_EndValues;
  this->ptr_Write2D          = &CPostBinary_Write2D;
  this->ptr_Write3D          = &CPostBinary_Write3D;
  this->ptr_WriteInteger     = &CPostBinary_WriteInteger;
  this->ptr_WriteDouble      = &CPostBinary_WriteDouble;
  this->ptr_WriteValuesVA    = &CPostBinary_WriteValuesVA;
  this->ptr_WriteValues      = &CPostBinary_WriteValues;
  this->ptr_WriteValuesNS    = &CPostBinary_WriteValuesNS;
  this->ptr_WriteValuesNSV   = &CPostBinary_WriteValuesNSV;
  this->ptr_WriteElement     = &CPostBinary_WriteElement;
  this->ptr_WritePostHeader  = &CPostBinary_WritePostHeader;
}

/* CBufferValues implementation */

static 
void CBufferValues_Init(CBufferValues* this)
{
  assert(this);
  
  this->buffer_values = NULL;
  this->last_value = -1;
  this->values_size_min = 0;
  this->values_size_max = 0;
  this->buffer_values_size = 0;
  this->buffer_types = (GiD_ResultType*)malloc(sizeof(GiD_ResultType)*10);
  this->buffer_types_size = 10;
  this->size_types = 0;
  this->next_type = 0;
}

CBufferValues* CBufferValues_Create( void )
{
  CBufferValues* this = (CBufferValues*)malloc(sizeof(CBufferValues));
  CBufferValues_Init(this);

  return this;
}

void CBufferValues_Destroy(CBufferValues* this)
{
  assert(this);
  
  if (this->buffer_values) {
    free(this->buffer_values);
    this->buffer_values = NULL;
  }
  if (this->buffer_types) {
    free(this->buffer_types);
    this->buffer_types = NULL;
  }
  free(this);
}

void CPostFile_ResultGroupOnBegin(CPostFile* this) 
{
  CBufferValues *buffer;
  
  assert(this);
  assert(this->buffer_values);

  buffer = this->buffer_values;
  buffer->size_types = 0;
  buffer->next_type = 0;
  buffer->values_size_min = 0;
  buffer->values_size_max = 0;
}

void CPostFile_ResultGroupOnNewType(CPostFile* this, GiD_ResultType res_type) 
{
  CBufferValues *buffer;
  size_t min, max;
  
  assert(this);
  assert(this->buffer_values);

  buffer = this->buffer_values;
  
  if (buffer->size_types==buffer->buffer_types_size) {
    buffer->buffer_types_size += 10;
    buffer->buffer_types = (GiD_ResultType*)realloc(buffer->buffer_types,
                                                    sizeof(GiD_ResultType)*( size_t)buffer->buffer_types_size);
  }
  buffer->buffer_types[buffer->size_types++] = res_type;
  GetResultTypeMinMaxValues( res_type, &min, &max);
  buffer->values_size_min += (int)min;
  buffer->values_size_max += (int)max;
}

void CPostFile_ResultGroupOnBeginValues(CPostFile* this)
{
  CBufferValues *buffer;
  
  assert(this);
  assert(this->buffer_values);

  buffer = this->buffer_values;

  if (!buffer->buffer_values) {
    buffer->buffer_values_size=buffer->values_size_max;
    buffer->buffer_values = (double*)malloc(( size_t)buffer->buffer_values_size*sizeof(double));
  } else if (buffer->values_size_max > buffer->buffer_values_size) {
    buffer->buffer_values_size=buffer->values_size_max;
    buffer->buffer_values = (double*)realloc(buffer->buffer_values, ( size_t)buffer->buffer_values_size*sizeof(double));
  }
  buffer->last_value = -1;
}

int CPostFile_ResultGroupFlushValues(CPostFile* this, int id) 
{
  CBufferValues *buffer;
  
  assert(this);
  assert(this->buffer_values);

  buffer = this->buffer_values;

  assert(buffer->last_value>=buffer->values_size_min-1 ||
         buffer->last_value<=buffer->values_size_max-1);
  if (CPostFile_WriteValues(this, id, buffer->last_value+1, buffer->buffer_values)) {
    /* could not write buffer values */
    buffer->last_value = -1;
    return 1;
  }
  /* prepare for next group of values */
  buffer->last_value = -1;
  return 0;
}

int CBufferValues_OnWriteType(CBufferValues* buffer, GiD_ResultType res_type)
{
  assert(buffer);
  
  if (res_type == buffer->buffer_types[buffer->next_type++]) {
    if (buffer->next_type == buffer->size_types) {
      buffer->next_type = 0;
      return 1;
    }
    return 0;
  } else {
    printf("error expected '%s' instead of '%s'\n",
           GetResultTypeName(buffer->buffer_types[buffer->next_type-1], 0),
           GetResultTypeName( res_type, 0));
    return -1;
  }
}

int CPostFile_ResultGroupWriteValues(CPostFile* this, GiD_ResultType res_type,
                                     int id, int n_comp, ...)
{
  CBufferValues *buffer;
  va_list ap;
  int i;
  int flush;
  
  assert(this);
  assert(this->buffer_values);

  buffer = this->buffer_values;

  assert(buffer->last_value + n_comp < buffer->values_size_max);
  
  flush = CBufferValues_OnWriteType(buffer, res_type);
  if (flush==-1)
    return -1;
  va_start(ap, n_comp);
  for (i = 0; i < n_comp; i++)
    buffer->buffer_values[++buffer->last_value] = va_arg(ap, double);
  va_end(ap);
  
  if (flush)
    return CPostFile_ResultGroupFlushValues(this, id);
  return 0;  
}

int CPostFile_ResultGroupIsEmpty(CPostFile* this)
{
  CBufferValues *buffer;
  
  assert(this);
  assert(this->buffer_values);

  buffer = this->buffer_values;

  
  return !buffer->buffer_values || buffer->last_value==-1;
}
