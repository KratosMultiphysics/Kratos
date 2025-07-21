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

static int ByteOrderCheck=0x91d;  /* Magic number */

/* format to print reals to string, customizable only once at the begin*/
static char G_format_real[100]={"%.9g"};

GIDPOST_API int GiD_PostSetFormatReal(GP_CONST char* f){
  strcpy(G_format_real,f);
  return 0;
}

GIDPOST_API GP_CONST char *GiD_PostGetFormatReal() {
  return G_format_real;
}

/* special format to not truncate time steps converted to string, not customizable but centralize its use*/
GIDPOST_API GP_CONST char *GiD_PostGetFormatStep(){
  static const char format_step[] ={ "%.16g"};
  return format_step;
}

struct _CBufferValues
{
  double* buffer_values;
  int last_value;
  int values_size_min;
  int values_size_max;
  int buffer_values_size;
  GiD_ResultType* buffer_types;
  int buffer_types_size;
  int size_types;
  int next_type;
};

CBufferValues* CBufferValues_Create(void);
void CBufferValues_Destroy(CBufferValues* post_file);

#define NMAX_DIMS 4
typedef struct {
  GP_CONST char* str;
  int dims[NMAX_DIMS];
} SResultTypeInfo;

static SResultTypeInfo _ResultTypeInfo[]={
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

GP_CONST char* GetResultTypeName(GiD_ResultType type, int dim)
{
  static char buffer[255];
  char* ptr;
  int i;

  strcpy(buffer,_ResultTypeInfo[(int)type].str);
  ptr=&(buffer[0]) + strlen(buffer);
  if(dim) {
    for(i=0; _ResultTypeInfo[(int)type].dims[i]; i++) {
      if ( dim == _ResultTypeInfo[(int)type].dims[i])
        break;
    }
    if(!(_ResultTypeInfo[(int)type].dims[i]))
      printf("Invalid dimension %d for type %s\n", dim,buffer);
    else
      sprintf(ptr,":%d", dim);
  }
  return buffer;
}

void GetResultTypeMinMaxValues(GiD_ResultType type, int *min, int *max)
{
  int i=0;

  *min=_ResultTypeInfo[(int)type].dims[0];
  for(i=1; _ResultTypeInfo[(int)type].dims[i]; i++)
    ;
  *max=_ResultTypeInfo[(int)type].dims[i-1];
}

/* ---------------------------------------------------------------------------
 *
 *  Post Files classes implementation : gidpostInt.h
 *
 * ---------------------------------------------------------------------------
 */

static
CPostFile* CPostFile_Create(void)
{
  CPostFile* post_file=(CPostFile*)malloc(sizeof(CPostFile));
  assert(post_file);
  post_file->m_LastID=-1;
  post_file->m_FILE=NULL;
  post_file->m_fail=0;
  post_file->buffer_values=CBufferValues_Create();
  post_file->GP_number_check=0;
  post_file->gauss_written=0;
  post_file->flag_isgroup=0;
  post_file->flag_begin_values=0;
  post_file->has_mesh=0;
  post_file->has_meshgroup=0;
  post_file->level_mesh=POST_UNDEFINED;
  post_file->level_res=POST_UNDEFINED;
  post_file->stack_pos=-1;
  post_file->local_axes_format[0]='\0';

  post_file->ptr_Open=NULL;
  post_file->ptr_Close=NULL;
  post_file->ptr_Flush=NULL;
  post_file->ptr_IsBinary=NULL;
  post_file->ptr_WriteString=NULL;
  post_file->ptr_BeginCoordinates=NULL;
  post_file->ptr_BeginElements=NULL;
  post_file->ptr_BeginValues=NULL;
  post_file->ptr_EndValues=NULL;
  post_file->ptr_WriteInteger=NULL;
  post_file->ptr_WriteDouble=NULL;
  post_file->ptr_WriteValuesVA=NULL;
  post_file->ptr_WriteValues=NULL;
  post_file->ptr_WriteValuesNS=NULL;
  post_file->ptr_WriteValuesNSV=NULL;
  post_file->ptr_Write2D=NULL;
  post_file->ptr_Write3D=NULL;
  post_file->ptr_WriteElement=NULL;
  post_file->ptr_WritePostHeader=NULL;
  post_file->ptr_WritePostHeaderIGA=NULL;

  post_file->m_post_mode=GiD_PostUndefined;
  post_file->m_encoding_filename=GiD_PostEncodeLocale;
  post_file->m_hdf5_file=NULL;

  return post_file;
}

int CPostFile_Release(CPostFile* post_file)
{
  int ret;

  assert(post_file);

  ret=CPostFile_Close(post_file);
  if(post_file->buffer_values)  {
    CBufferValues_Destroy(post_file->buffer_values);
    post_file->buffer_values=NULL;
  }
  free(post_file);
  return ret;
}

int CPostFile_Open(CPostFile* post_file,GP_CONST char* str)
{
  assert(post_file);
  assert(post_file->ptr_Open);
  return (*post_file->ptr_Open)(post_file,str);
}

int CPostFile_Close(CPostFile* post_file)
{
  assert(post_file);
  if(post_file->m_post_mode == GiD_PostHDF5) {
    // already closed before with GiD_ClosePostResultFile_HDF5()
    return 0;
  }
  assert(post_file->ptr_Close);
  return (*post_file->ptr_Close)(post_file);
}

post_state CPostFile_TopState(CPostFile* post_file)
{
  return (post_file->stack_pos < 0)?POST_UNDEFINED:post_file->stack_state[post_file->stack_pos];
}

int CPostFile_PushState(CPostFile* post_file,post_state s)
{
  assert(post_file->stack_pos < STACK_STATE_SIZE - 1);
  if(post_file->stack_pos < STACK_STATE_SIZE - 1)
  {
    post_file->stack_state[++post_file->stack_pos]=s;
    return 0;
  } else
  {
    return -1;
  }
}

post_state CPostFile_PopState(CPostFile* post_file)
{
  post_state top;
  if(post_file->stack_pos < 0)
  {
    return POST_UNDEFINED;
  }
  top=post_file->stack_state[post_file->stack_pos--];
  return top;
}

int CPostFile_Flush(CPostFile* post_file)
{
  int res=0;
  assert(post_file);
#ifdef ENABLE_HDF5
  if(post_file->m_post_mode==GiD_PostHDF5) {
    res=GiD_FlushPostFile_HDF5(post_file->m_hdf5_file);
  } else {
    assert(post_file->ptr_Flush);
    res=(*post_file->ptr_Flush)(post_file);
  }
#else //ENABLE_HDF5
  assert(post_file->ptr_Flush);
  res=(*post_file->ptr_Flush)(post_file);
#endif //ENABLE_HDF5
  return res;
}

int CPostFile_IsBinary(CPostFile* post_file)
{
  assert(post_file);
  assert(post_file->ptr_IsBinary);
  return (*post_file->ptr_IsBinary)(post_file);
}

int CPostFile_WriteString(CPostFile* post_file,GP_CONST char* str)
{
  assert(post_file);
  assert(post_file->ptr_WriteString);
  return (*post_file->ptr_WriteString)(post_file,str);
}

int CPostFile_BeginCoordinates(CPostFile* post_file)
{
  assert(post_file);
  assert(post_file->ptr_BeginCoordinates);

  return (*post_file->ptr_BeginCoordinates)(post_file);
}

int CPostFile_BeginElements(CPostFile* post_file)
{
  assert(post_file);
  assert(post_file->ptr_BeginElements);

  return (*post_file->ptr_BeginElements)(post_file);
}

int CPostFile_BeginValues(CPostFile* post_file)
{
  assert(post_file);
  assert(post_file->ptr_BeginValues);

  return (*post_file->ptr_BeginValues)(post_file);
}

int CPostFile_EndValues(CPostFile* post_file)
{
  assert(post_file);
  assert(post_file->ptr_EndValues);

  return (*post_file->ptr_EndValues)(post_file);
}

int CPostFile_WriteInteger(CPostFile* post_file,int i,int op)
{
  assert(post_file);
  assert(post_file->ptr_WriteInteger);

  return (*post_file->ptr_WriteInteger)(post_file,i,op);
}

int CPostFile_WriteDouble(CPostFile* post_file,double x,int op)
{
  assert(post_file);
  assert(post_file->ptr_WriteDouble);

  return (*post_file->ptr_WriteDouble)(post_file,x,op);
}

int CPostFile_WriteValuesVA(CPostFile* post_file,int id,int num_comp,...)
{
  va_list ap;
  int ret;

  assert(post_file);
  assert(post_file->ptr_WriteValuesVA);

  va_start(ap,num_comp);
  ret=(*post_file->ptr_WriteValuesVA)(post_file,id,num_comp,ap);
  va_end(ap);

  return ret;
}

int CPostFile_WriteValues(CPostFile* post_file,int id,int n,GP_CONST double* values) {
  assert(post_file);
  assert(post_file->ptr_WriteValues);

  return (*post_file->ptr_WriteValues)(post_file,id,n,values);
}

int CPostFile_WriteValuesNS(CPostFile* post_file,int id,int n,GP_CONST double* values) {
  assert(post_file);
  assert(post_file->ptr_WriteValuesNS);

  return (*post_file->ptr_WriteValuesNS)(post_file,id,n,values);
}

int CPostFile_WriteValuesNSV(CPostFile* post_file,int id,int n,int num_comp,GP_CONST double* values) {
  assert(post_file);
  assert(post_file->ptr_WriteValuesNSV);

  return (*post_file->ptr_WriteValuesNSV)(post_file,id,n,num_comp,values);
}

int CPostFile_Write2D(CPostFile* post_file,double x,double y)
{
  assert(post_file);
  assert(post_file->ptr_Write2D);

  return (*post_file->ptr_Write2D)(post_file,x,y);
}

int CPostFile_Write3D(CPostFile* post_file,double x,double y,double z)
{
  assert(post_file);
  assert(post_file->ptr_Write3D);

  return (*post_file->ptr_Write3D)(post_file,x,y,z);
}

int CPostFile_WriteElement(CPostFile* post_file,int id,int n,GP_CONST int nid[]) {
  assert(post_file);
  assert(post_file->ptr_WriteElement);

  return (*post_file->ptr_WriteElement)(post_file,id,n,nid);
}

int CPostFile_WritePostHeader(CPostFile* post_file)
{
  assert(post_file);
  assert(post_file->ptr_WritePostHeader);

  return (*post_file->ptr_WritePostHeader)(post_file);
}

void CPostFile_ResetLastID(CPostFile* post_file)
{
  assert(post_file);

  post_file->m_LastID=-1;
}

void CPostFile_SetConnectivity(CPostFile* post_file,int nnode)
{
  assert(post_file);

  post_file->m_connectivity=nnode;
}

int CPostFile_GetConnectivity(CPostFile* post_file)
{
  assert(post_file);

  return post_file->m_connectivity;
}

int CPostFile_MatchConnectivity(CPostFile* post_file,int written)
{
  int conn=CPostFile_GetConnectivity(post_file);
  return (written == conn)?0:(((written-1) == conn)?1:2);
}

/*
 *  class CPostAscii
 */

static
void CPostAscii_Init(CPostFile*);

CPostFile* CPostAscii_Create()
{
  CPostFile* post_file;

  post_file=CPostFile_Create();
  CPostAscii_Init(post_file);
  return post_file;
}

static int CPostAscii_Close(CPostFile* post_file);

static
int CPostAscii_Open(CPostFile* post_file,GP_CONST char* name) {
  assert(post_file);
  assert(post_file->ptr_Close);

  CPostAscii_Close(post_file);
  if(post_file->m_encoding_filename == GiD_PostEncodeUTF8) {
    /* convert from utf-8 */
#ifdef _WIN32
    wchar_t wname[MAX_PATH];
    MultiByteToWideChar(CP_UTF8,0,name,-1,wname,MAX_PATH);
    post_file->m_FILE=_wfopen(wname,L"w");
#else
    post_file->m_FILE=fopen(name,"w");
#endif
  } else {
    assert(post_file->m_encoding_filename == GiD_PostEncodeLocale);
    post_file->m_FILE=fopen(name,"w");
  }
  return post_file->m_FILE == NULL;
}

static
int CPostAscii_Close(CPostFile* post_file)
{
  assert(post_file);

  if(post_file->m_FILE) {
    post_file->m_fail=fclose(post_file->m_FILE);
    post_file->m_FILE=NULL;
  } else
    post_file->m_fail=1;
  return post_file->m_fail;
}

static
int CPostAscii_Flush(CPostFile* post_file)
{
  return post_file->m_FILE?fflush(post_file->m_FILE):1;
}

static
int CPostAscii_IsBinary(CPostFile* post_file)
{
  return 0;
}

static
int CPostAscii_WriteString(CPostFile* post_file,GP_CONST char* str)
{
  assert(post_file);

  fprintf(post_file->m_FILE,"%s\n",str);
  return 0;
}

static
int CPostAscii_BeginCoordinates(CPostFile* post_file)
{
  return CPostAscii_WriteString(post_file,"Coordinates");
}

static
int CPostAscii_BeginElements(CPostFile* post_file)
{
  return CPostAscii_WriteString(post_file,"Elements");
}

static
int CPostAscii_BeginValues(CPostFile* post_file)
{
  return CPostAscii_WriteString(post_file,"Values");
}

static
int CPostAscii_EndValues(CPostFile* post_file)
{
  return CPostAscii_WriteString(post_file,"End Values");
}

static
int CPostAscii_WriteInteger(CPostFile* post_file,int i,int op)
{
  assert(post_file);

  if(op == 0) {
    fprintf(post_file->m_FILE,"%d",i);
  } else {
    fprintf(post_file->m_FILE," %d",i);
    if(op == 2) {
      fprintf(post_file->m_FILE,"\n");
    }
  }
  return 0;
}

static
int CPostAscii_WriteDouble(CPostFile* post_file,double x,int op)
{
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);
  assert(op!=0);
  fprintf(post_file->m_FILE,local_format,x);
  if(op==2) {
    fprintf(post_file->m_FILE,"\n");
  }
  return 0;
}

static
int CPostAscii_WriteValuesVA(CPostFile* post_file,int id,int num_comp,va_list ap)
{
  int i;
  double value;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);
  assert(post_file->m_FILE);

  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(post_file->m_FILE,"%d",id);

  for(i=0; i < num_comp; i++) {
    value=va_arg(ap,double);
    fprintf(post_file->m_FILE,local_format,value);
  }
  fprintf(post_file->m_FILE,"\n");
  post_file->m_LastID=id;

  return 0;
}

static int CPostAscii_WriteValues(CPostFile* post_file,int id,int n,GP_CONST double* buffer) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);

  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(post_file->m_FILE,"%d",id);
  for(i=0; i < n; i++)
    fprintf(post_file->m_FILE,local_format,buffer[i]);
  fprintf(post_file->m_FILE,"\n");
  post_file->m_LastID=id;
  return 0;
}

static int CPostAscii_WriteValuesNS(CPostFile* post_file,int id,int n,GP_CONST double* buffer) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s\n",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);

  //write surface Id
  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(post_file->m_FILE,"%d\n",id);

  //write values
  for(i=0; i < n; i++) {
    if(buffer[i] == GP_UNKNOWN) fprintf(post_file->m_FILE,"NR\n");
    else fprintf(post_file->m_FILE,local_format,buffer[i]);
  }
  post_file->m_LastID=id;
  return 0;
}

static int CPostAscii_WriteValuesNSV(CPostFile* post_file,int id,int n,int num_comp,GP_CONST double* buffer) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s ",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);

  //write surface Id
  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(post_file->m_FILE,"%d\n",id);

  //write values
  for(i=0; i < n*num_comp; i++) {
    if(buffer[i] == GP_UNKNOWN) fprintf(post_file->m_FILE,"NR ");
    else fprintf(post_file->m_FILE,local_format,buffer[i]);
    if((i+1) % num_comp == 0) fprintf(post_file->m_FILE,"\n");
  }
  post_file->m_LastID=id;
  return 0;
}

static
int CPostAscii_Write2D(CPostFile* post_file,double x,double y)
{
  static char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format,x,y);
  return CPostAscii_WriteString(post_file,line);
}

static
int CPostAscii_Write3D(CPostFile* post_file,double x,double y,double z)
{
  static char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format,x,y,z);
  return CPostAscii_WriteString(post_file,line);
}

static int CPostAscii_WriteElement(CPostFile* post_file,int id,int n,GP_CONST int nid[]) {
  int i;

  assert(post_file);

  fprintf(post_file->m_FILE,"%d",id);
  for(i=0; i < n; i++) {
    fprintf(post_file->m_FILE," %d",nid[i]);
  }
  fprintf(post_file->m_FILE,"\n");

  return 0;
}

static
int CPostAscii_WritePostHeader(CPostFile* post_file)
{
  return CPostAscii_WriteString(post_file,"GiD Post Results File 1.2");
}

static
void CPostAscii_Init(CPostFile* post_file)
{
  assert(post_file);

  post_file->ptr_Open=&CPostAscii_Open;
  post_file->ptr_Close=&CPostAscii_Close;
  post_file->ptr_Flush=&CPostAscii_Flush;
  post_file->ptr_IsBinary=&CPostAscii_IsBinary;
  post_file->ptr_WriteString=&CPostAscii_WriteString;
  post_file->ptr_BeginCoordinates=&CPostAscii_BeginCoordinates;
  post_file->ptr_BeginElements=&CPostAscii_BeginElements;
  post_file->ptr_BeginValues=&CPostAscii_BeginValues;
  post_file->ptr_EndValues=&CPostAscii_EndValues;
  post_file->ptr_WriteInteger=&CPostAscii_WriteInteger;
  post_file->ptr_WriteDouble=&CPostAscii_WriteDouble;
  post_file->ptr_WriteValuesVA=&CPostAscii_WriteValuesVA;
  post_file->ptr_WriteValues=&CPostAscii_WriteValues;
  post_file->ptr_WriteValuesNS=&CPostAscii_WriteValuesNS;
  post_file->ptr_WriteValuesNSV=&CPostAscii_WriteValuesNSV;
  post_file->ptr_Write2D=&CPostAscii_Write2D;
  post_file->ptr_Write3D=&CPostAscii_Write3D;
  post_file->ptr_WriteElement=&CPostAscii_WriteElement;
  post_file->ptr_WritePostHeader=&CPostAscii_WritePostHeader;
}

/*
 *  class CPostAsciiZ --
 */

static
void CPostAsciiZ_Init(CPostFile*);

CPostFile* CPostAsciiZ_Create()
{
  CPostFile* post_file;

  post_file=CPostFile_Create();
  CPostAsciiZ_Init(post_file);
  return post_file;
}

static int CPostAsciiZ_Close(CPostFile* post_file);

static
int CPostAsciiZ_Open(CPostFile* post_file,GP_CONST char* name) {
  assert(post_file);
  assert(post_file->ptr_Close);

  CPostAsciiZ_Close(post_file);
  if(post_file->m_encoding_filename == GiD_PostEncodeUTF8) {
    /* convert from utf-8 */
#ifdef _WIN32
    wchar_t wname[MAX_PATH];
    MultiByteToWideChar(CP_UTF8,0,name,-1,wname,MAX_PATH);
    post_file->m_FILE=gzopen_w(wname,"w1");
#else
    post_file->m_FILE=gzopen(name,"w1");
#endif
  } else {
    assert(post_file->m_encoding_filename == GiD_PostEncodeLocale);
    post_file->m_FILE=gzopen(name,"w1");
  }
  return post_file->m_FILE == NULL;
}

static
int CPostAsciiZ_Close(CPostFile* post_file)
{
  assert(post_file);

  if(post_file->m_FILE) {
    post_file->m_fail=gzclose(post_file->m_FILE);
    post_file->m_FILE=NULL;
  } else
    post_file->m_fail=1;
  return post_file->m_fail;
}

static
int CPostAsciiZ_Flush(CPostFile* post_file)
{
  assert(post_file);

  return post_file->m_FILE?gzflush(post_file->m_FILE,Z_FINISH):1;
}

static
int CPostAsciiZ_IsBinary(CPostFile* post_file)
{
  return 0;
}

static
int CPostAsciiZ_WriteString(CPostFile* post_file,GP_CONST char* str)
{
  assert(post_file);

  gzprintf(post_file->m_FILE,"%s\n",str);
  return 0;
}

static
int CPostAsciiZ_BeginCoordinates(CPostFile* post_file)
{
  return CPostAsciiZ_WriteString(post_file,"Coordinates");
}

static
int CPostAsciiZ_BeginElements(CPostFile* post_file)
{
  return CPostAsciiZ_WriteString(post_file,"Elements");
}

static
int CPostAsciiZ_BeginValues(CPostFile* post_file)
{
  return CPostAsciiZ_WriteString(post_file,"Values");
}

static
int CPostAsciiZ_EndValues(CPostFile* post_file)
{
  return CPostAsciiZ_WriteString(post_file,"End Values");
}

static
int CPostAsciiZ_WriteInteger(CPostFile* post_file,int i,int op)
{
  assert(post_file);

  if(op==1) {
    gzprintf(post_file->m_FILE," ");
    gzprintf(post_file->m_FILE,"%d",i);
  } else {
    gzprintf(post_file->m_FILE,"%d",i);
    if(op==2) {
      gzprintf(post_file->m_FILE,"\n");
    }
  }

  return 0;
}

static
int CPostAsciiZ_WriteDouble(CPostFile* post_file,double x,int op)
{
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);
  assert(op!=0);

  gzprintf(post_file->m_FILE,local_format,x);
  if(op==2) {
    gzprintf(post_file->m_FILE,"\n");
  }

  return 0;
}

static
int CPostAsciiZ_WriteValuesVA(CPostFile* post_file,int id,int num_comp,va_list ap)
{
  int i;
  double value;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);

  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzprintf(post_file->m_FILE,"%d",id);
  for(i=0; i < num_comp; i++) {
    value=va_arg(ap,double);
    gzprintf(post_file->m_FILE,local_format,value);
  }
  gzprintf(post_file->m_FILE,"\n");
  post_file->m_LastID=id;

  return 0;
}

static int CPostAsciiZ_WriteValues(CPostFile* post_file,int id,int n,GP_CONST double* buffer) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);

  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzprintf(post_file->m_FILE,"%d",id);
  for(i=0; i < n; i++)
    gzprintf(post_file->m_FILE,local_format,buffer[i]);
  gzprintf(post_file->m_FILE,"\n");
  post_file->m_LastID=id;

  return 0;
}

static int CPostAsciiZ_WriteValuesNS(CPostFile* post_file,int id,int n,GP_CONST double* buffer) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s\n",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);

  //write surface Id
  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(post_file->m_FILE,"%d\n",id);

  //write values
  for(i=0; i < n; i++) {
    if(buffer[i] == GP_UNKNOWN) fprintf(post_file->m_FILE,"NR\n");
    else fprintf(post_file->m_FILE,local_format,buffer[i]);
  }
  post_file->m_LastID=id;
  return 0;
}

static int CPostAsciiZ_WriteValuesNSV(CPostFile* post_file,int id,int n,int num_comp,GP_CONST double* buffer) {
  int i;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s ",GiD_PostGetFormatReal());
    create_format=0;
  }
  assert(post_file);

  //write surface Id
  if(post_file->m_LastID != id) /* Hey this is only useful when writing values for gauss points in elements !!!! */
    fprintf(post_file->m_FILE,"%d\n",id);

  //write values
  for(i=0; i < n*num_comp; i++) {
    if(buffer[i] == GP_UNKNOWN) fprintf(post_file->m_FILE,"NR ");
    else fprintf(post_file->m_FILE,local_format,buffer[i]);
    if((i+1) % num_comp == 0) fprintf(post_file->m_FILE,"\n");
  }
  post_file->m_LastID=id;
  return 0;
}

static
int CPostAsciiZ_Write2D(CPostFile* post_file,double x,double y)
{
  static char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format,x,y);
  return CPostAsciiZ_WriteString(post_file,line);
}

static
int CPostAsciiZ_Write3D(CPostFile* post_file,double x,double y,double z)
{
  char line[256];
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s %s %s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
  sprintf(line,local_format,x,y,z);
  return CPostAsciiZ_WriteString(post_file,line);
}

static int CPostAsciiZ_WriteElement(CPostFile* post_file,int id,int n,GP_CONST int nid[]) {
  int i;

  assert(post_file);

  gzprintf(post_file->m_FILE,"%d",id);
  for(i=0; i < n; i++) {
    gzprintf(post_file->m_FILE," %d",nid[i]);
  }
  gzprintf(post_file->m_FILE,"\n");

  return 0;
}

static
int CPostAsciiZ_WritePostHeader(CPostFile* post_file)
{
  return CPostAsciiZ_WriteString(post_file,"GiD Post Results File 1.2");
}


static
void CPostAsciiZ_Init(CPostFile* post_file)
{
  assert(post_file);

  post_file->ptr_Open=&CPostAsciiZ_Open;
  post_file->ptr_Close=&CPostAsciiZ_Close;
  post_file->ptr_Flush=&CPostAsciiZ_Flush;
  post_file->ptr_IsBinary=&CPostAsciiZ_IsBinary;
  post_file->ptr_WriteString=&CPostAsciiZ_WriteString;
  post_file->ptr_BeginCoordinates=&CPostAsciiZ_BeginCoordinates;
  post_file->ptr_BeginElements=&CPostAsciiZ_BeginElements;
  post_file->ptr_BeginValues=&CPostAsciiZ_BeginValues;
  post_file->ptr_EndValues=&CPostAsciiZ_EndValues;
  post_file->ptr_Write2D=&CPostAsciiZ_Write2D;
  post_file->ptr_Write3D=&CPostAsciiZ_Write3D;
  post_file->ptr_WriteInteger=&CPostAsciiZ_WriteInteger;
  post_file->ptr_WriteDouble=&CPostAsciiZ_WriteDouble;
  post_file->ptr_WriteValuesVA=&CPostAsciiZ_WriteValuesVA;
  post_file->ptr_WriteValues=&CPostAsciiZ_WriteValues;
  post_file->ptr_WriteValuesNS=&CPostAsciiZ_WriteValuesNS;
  post_file->ptr_WriteValuesNSV=&CPostAsciiZ_WriteValuesNSV;
  post_file->ptr_WriteElement=&CPostAsciiZ_WriteElement;
  post_file->ptr_WritePostHeader=&CPostAsciiZ_WritePostHeader;
}

/*
 *  class CPostBinary --
 */

static
void CPostBinary_Init(CPostFile*);

CPostFile* CPostBinary_Create()
{
  CPostFile* post_file;

  post_file=CPostFile_Create();
  CPostBinary_Init(post_file);
  return post_file;
}

CPostFile* CPostHdf5_Create() {
  return CPostFile_Create();
}

static int CPostBinary_Close(CPostFile* post_file);

static
int CPostBinary_Open(CPostFile* post_file,GP_CONST char* name) {
  assert(post_file);
  assert(post_file->ptr_Close);

  CPostBinary_Close(post_file);
  if(post_file->m_encoding_filename == GiD_PostEncodeUTF8) {
    /* convert from utf-8 */
#ifdef _WIN32
    wchar_t wname[MAX_PATH];
    MultiByteToWideChar(CP_UTF8,0,name,-1,wname,MAX_PATH);
    post_file->m_FILE=gzopen_w(wname,"wb1");
#else
    post_file->m_FILE=gzopen(name,"wb1");
#endif
  } else {
    assert(post_file->m_encoding_filename == GiD_PostEncodeLocale);
    post_file->m_FILE=gzopen(name,"wb1");
  }

  /* escribir el numero magico */
  if(post_file->m_FILE) {
    gzwrite(post_file->m_FILE,&ByteOrderCheck,sizeof(ByteOrderCheck));
  }

  return post_file->m_FILE == NULL;
}

static
int CPostBinary_Close(CPostFile* post_file)
{
  assert(post_file);

  if(post_file->m_FILE) {
    post_file->m_fail=gzclose(post_file->m_FILE);
    post_file->m_FILE=NULL;
  } else
    post_file->m_fail=1;

  return post_file->m_fail;
}

static
int CPostBinary_Flush(CPostFile* post_file)
{
  assert(post_file);

  return post_file->m_FILE?gzflush(post_file->m_FILE,Z_FULL_FLUSH):1;
}

static
int CPostBinary_IsBinary(CPostFile* post_file)
{
  return 1;
}

static
int CPostBinary_WriteString(CPostFile* post_file,GP_CONST char* str)
{
  int size,written=0,tam;

  assert(post_file);

  post_file->m_fail=1;
  if(post_file->m_FILE) {
    GP_CONST char* buf="\0";

    if(str) buf=str;
    tam=(int)strlen(buf) + 1; /* incluido el \0 */
    written=gzwrite(post_file->m_FILE,&tam,sizeof(int));
    size=(int)sizeof(char) * tam;
    written+=gzwrite(post_file->m_FILE,(void*)buf,(unsigned int)size);
    if(written == size+(int)(sizeof(int)))
      post_file->m_fail=0;
  }
  return post_file->m_fail;
}

static
int CPostBinary_BeginCoordinates(CPostFile* post_file)
{
  return CPostBinary_WriteString(post_file,"Coordinates -1 Indexed");
}

static
int CPostBinary_BeginElements(CPostFile* post_file)
{
  return CPostBinary_WriteString(post_file,"Elements -1 Indexed");
}

static
int CPostBinary_BeginValues(CPostFile* post_file)
{
  return CPostBinary_WriteString(post_file,"Values -1 Indexed");
}

static
int CPostBinary_EndValues(CPostFile* post_file)
{
  int idxend=-1;

  assert(post_file);
  assert(post_file->m_FILE);

  if(gzwrite(post_file->m_FILE,&idxend,sizeof(int)) != sizeof(int))
    return 1;
  return CPostBinary_WriteString(post_file,"End Values");
}

static
int CPostBinary_WriteInteger(CPostFile* post_file,int i,int op)
{
  assert(post_file);
  assert(post_file->m_FILE);

  gzwrite(post_file->m_FILE,&i,sizeof(i));
  return 0;
}

static
int CPostBinary_WriteDouble(CPostFile* post_file,double x,int op)
{
  float v=(float)(x);

  assert(post_file);
  assert(post_file->m_FILE);

  gzwrite(post_file->m_FILE,&v,sizeof(v));
  return 0;
}

static
int CPostBinary_WriteValuesVA(CPostFile* post_file,int id,int num_comp,va_list ap)
{
  int i;
  float value;

  assert(post_file);
  assert(post_file->m_FILE);

  if(post_file->m_LastID != id) { /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzwrite(post_file->m_FILE,&id,sizeof(id));
    post_file->m_LastID=id;
  }
  for(i=0; i < num_comp; i++) {
    value=(float)va_arg(ap,double);
    gzwrite(post_file->m_FILE,&value,sizeof(value));
  }
  return 0;
}

static int CPostBinary_WriteValues(CPostFile* post_file,int id,int n,GP_CONST double* buffer) {
  int i;
  float value;

  assert(post_file);
  assert(post_file->m_FILE);

  if(post_file->m_LastID != id) { /* Hey this is only useful when writing values for gauss points in elements !!!! */
    gzwrite(post_file->m_FILE,&id,sizeof(id));
    post_file->m_LastID=id;
  }
  for(i=0; i < n; i++) {
    value=(float)buffer[i];
    gzwrite(post_file->m_FILE,&value,sizeof(value));
  }
  return 0;
}

static int CPostBinary_WriteValuesNS(CPostFile* post_file,int id,int n,GP_CONST double* buffer) {
  int i;
  float value;
  assert(post_file);
  //write surface Id
  gzwrite(post_file->m_FILE,&id,sizeof(id));
  //write amount of values (else without the geometry is not possible to read it)
  gzwrite(post_file->m_FILE,&n,sizeof(n));
  //write values
  for(i=0; i < n; i++) {
    value=(float)buffer[i];
    gzwrite(post_file->m_FILE,&value,sizeof(value));
  }
  return 0;
}

static int CPostBinary_WriteValuesNSV(CPostFile* post_file,int id,int n,int num_comp,GP_CONST double* buffer) {
  int i;
  float value;
  int num_values=n*num_comp;
  assert(post_file);

  //write surface Id
  gzwrite(post_file->m_FILE,&id,sizeof(id));
  //write amount of values (else without the geometry is not possible to read it)
  gzwrite(post_file->m_FILE,&num_values,sizeof(num_values));
  //write values
  for(i=0; i < num_values; i++) {
    value=(float)buffer[i];
    gzwrite(post_file->m_FILE,&value,sizeof(value));
  }
  post_file->m_LastID=id;
  return 0;
}

static
int CPostBinary_Write2D(CPostFile* post_file,double x,double y)
{
  float values[2];

  assert(post_file);
  assert(post_file->m_FILE);

  values[0]=(float)(x);
  values[1]=(float)(y);
  gzwrite(post_file->m_FILE,values,sizeof(float)*2);

  return 0;
}

static
int CPostBinary_Write3D(CPostFile* post_file,double x,double y,double z)
{
  float values[3];

  assert(post_file);
  assert(post_file->m_FILE);

  values[0]=(float)(x);
  values[1]=(float)(y);
  values[2]=(float)(z);
  gzwrite(post_file->m_FILE,values,sizeof(float)*3);

  return 0;
}

static int CPostBinary_WriteElement(CPostFile* post_file,int id,int n,GP_CONST int nid[]) {
  int i;

  assert(post_file);
  assert(post_file->m_FILE);

  gzwrite(post_file->m_FILE,&id,sizeof(id));
  for(i=0; i < n; i++) {
    gzwrite(post_file->m_FILE,nid+i,sizeof(int));
  }
  /* verify connectivity */
  switch(CPostFile_MatchConnectivity(post_file,n)) {
  case 0:
    /* match exactly */
    i=0;
    /* so write material 0 == no material */
    gzwrite(post_file->m_FILE,&i,sizeof(int));
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
int CPostBinary_WritePostHeader(CPostFile* post_file)
{
  return CPostBinary_WriteString(post_file,"GiDPostEx1.2");
}


static
void CPostBinary_Init(CPostFile* post_file)
{
  assert(post_file);

  post_file->ptr_Open=&CPostBinary_Open;
  post_file->ptr_Close=&CPostBinary_Close;
  post_file->ptr_Flush=&CPostBinary_Flush;
  post_file->ptr_IsBinary=&CPostBinary_IsBinary;
  post_file->ptr_WriteString=&CPostBinary_WriteString;
  post_file->ptr_BeginCoordinates=&CPostBinary_BeginCoordinates;
  post_file->ptr_BeginElements=&CPostBinary_BeginElements;
  post_file->ptr_BeginValues=&CPostBinary_BeginValues;
  post_file->ptr_EndValues=&CPostBinary_EndValues;
  post_file->ptr_Write2D=&CPostBinary_Write2D;
  post_file->ptr_Write3D=&CPostBinary_Write3D;
  post_file->ptr_WriteInteger=&CPostBinary_WriteInteger;
  post_file->ptr_WriteDouble=&CPostBinary_WriteDouble;
  post_file->ptr_WriteValuesVA=&CPostBinary_WriteValuesVA;
  post_file->ptr_WriteValues=&CPostBinary_WriteValues;
  post_file->ptr_WriteValuesNS=&CPostBinary_WriteValuesNS;
  post_file->ptr_WriteValuesNSV=&CPostBinary_WriteValuesNSV;
  post_file->ptr_WriteElement=&CPostBinary_WriteElement;
  post_file->ptr_WritePostHeader=&CPostBinary_WritePostHeader;
}

/* CBufferValues implementation */

static
void CBufferValues_Init(CBufferValues* post_file)
{
  assert(post_file);

  post_file->buffer_values=NULL;
  post_file->last_value=-1;
  post_file->values_size_min=0;
  post_file->values_size_max=0;
  post_file->buffer_values_size=0;
  post_file->buffer_types=(GiD_ResultType*)malloc(sizeof(GiD_ResultType)*10);
  post_file->buffer_types_size=10;
  post_file->size_types=0;
  post_file->next_type=0;
}

CBufferValues* CBufferValues_Create(void)
{
  CBufferValues* post_file=(CBufferValues*)malloc(sizeof(CBufferValues));
  CBufferValues_Init(post_file);

  return post_file;
}

void CBufferValues_Destroy(CBufferValues* post_file)
{
  assert(post_file);

  if(post_file->buffer_values) {
    free(post_file->buffer_values);
    post_file->buffer_values=NULL;
  }
  if(post_file->buffer_types) {
    free(post_file->buffer_types);
    post_file->buffer_types=NULL;
  }
  free(post_file);
}

void CPostFile_ResultGroupOnBegin(CPostFile* post_file)
{
  CBufferValues* buffer;

  assert(post_file);
  assert(post_file->buffer_values);

  buffer=post_file->buffer_values;
  buffer->size_types=0;
  buffer->next_type=0;
  buffer->values_size_min=0;
  buffer->values_size_max=0;
}

void CPostFile_ResultGroupOnNewType(CPostFile* post_file,GiD_ResultType res_type)
{
  CBufferValues *buffer = NULL;
  int min = 0, max = 0;

  assert(post_file);
  assert(post_file->buffer_values);

  buffer=post_file->buffer_values;

  if(buffer->size_types==buffer->buffer_types_size) {
    buffer->buffer_types_size+=10;
    buffer->buffer_types=(GiD_ResultType*)realloc(buffer->buffer_types,
      sizeof(GiD_ResultType)*(size_t)buffer->buffer_types_size);
  }
  buffer->buffer_types[buffer->size_types++]=res_type;
  GetResultTypeMinMaxValues(res_type,&min,&max);
  buffer->values_size_min += min;
  buffer->values_size_max += max;
}

void CPostFile_ResultGroupOnBeginValues(CPostFile* post_file)
{
  CBufferValues* buffer;

  assert(post_file);
  assert(post_file->buffer_values);

  buffer=post_file->buffer_values;

  if(!buffer->buffer_values) {
    buffer->buffer_values_size=buffer->values_size_max;
    buffer->buffer_values=(double*)malloc((size_t)buffer->buffer_values_size*sizeof(double));
  } else if(buffer->values_size_max > buffer->buffer_values_size) {
    buffer->buffer_values_size=buffer->values_size_max;
    buffer->buffer_values=(double*)realloc(buffer->buffer_values,(size_t)buffer->buffer_values_size*sizeof(double));
  }
  buffer->last_value=-1;
}

int CPostFile_ResultGroupFlushValues(CPostFile* post_file,int id)
{
  CBufferValues* buffer;

  assert(post_file);
  assert(post_file->buffer_values);

  buffer=post_file->buffer_values;

  assert(buffer->last_value>=buffer->values_size_min-1 ||
    buffer->last_value<=buffer->values_size_max-1);
  if(CPostFile_WriteValues(post_file,id,buffer->last_value+1,buffer->buffer_values)) {
    /* could not write buffer values */
    buffer->last_value=-1;
    return 1;
  }
  /* prepare for next group of values */
  buffer->last_value=-1;
  return 0;
}

int CBufferValues_OnWriteType(CBufferValues* buffer,GiD_ResultType res_type)
{
  assert(buffer);

  if(res_type == buffer->buffer_types[buffer->next_type++]) {
    if(buffer->next_type == buffer->size_types) {
      buffer->next_type=0;
      return 1;
    }
    return 0;
  } else {
    printf("error expected '%s' instead of '%s'\n",
      GetResultTypeName(buffer->buffer_types[buffer->next_type-1],0),
      GetResultTypeName(res_type,0));
    return -1;
  }
}

int CPostFile_ResultGroupWriteValues(CPostFile* post_file,GiD_ResultType res_type,
  int id,int n_comp,...)
{
  CBufferValues* buffer;
  va_list ap;
  int i;
  int flush;

  assert(post_file);
  assert(post_file->buffer_values);

  buffer=post_file->buffer_values;

  assert(buffer->last_value + n_comp < buffer->values_size_max);

  flush=CBufferValues_OnWriteType(buffer,res_type);
  if(flush==-1)
    return -1;
  va_start(ap,n_comp);
  for(i=0; i < n_comp; i++)
    buffer->buffer_values[++buffer->last_value]=va_arg(ap,double);
  va_end(ap);

  if(flush)
    return CPostFile_ResultGroupFlushValues(post_file,id);
  return 0;
}

int CPostFile_ResultGroupIsEmpty(CPostFile* post_file)
{
  CBufferValues* buffer;

  assert(post_file);
  assert(post_file->buffer_values);

  buffer=post_file->buffer_values;


  return !buffer->buffer_values || buffer->last_value==-1;
}
