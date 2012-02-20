/* gidpost 1.7 */
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

static int ByteOrderCheck = 0x91d;

// int GetResultTypeMinValues(GiD_ResultType type);
// int GetResultTypeMaxValues(GiD_ResultType type);

/* ---------------------------------------------------------------------------
 *
 *  Post Files classes implementation : gidpostInt.h
 *
 * ---------------------------------------------------------------------------
 */

int CPostFile::fail = 0;

CPostFile::CPostFile()
{
  _LastID() = -1;
}

CPostFile::~CPostFile()
{
}

int CPostFile::BeginCoordinates()
{
  return WriteString("Coordinates");
}

int CPostFile::BeginElements()
{
  return WriteString("Elements");
}

int CPostFile::BeginValues()
{
  return WriteString("Values");
}

int CPostFile::EndValues()
{
  return WriteString("End Values");
}

int CPostFile::Write2D( double x, double y )
{
  char line[256];

  sprintf( line, "%g %g", x, y );
  return WriteString(line);
}

int CPostFile::Write3D( double x, double y, double z )
{
  char line[256];

  sprintf( line, "%g %g %g", x, y, z );
  return WriteString(line);
}

int CPostFile::WritePostHeader()
{
  return WriteString("GiD Post Results File 1.0");
}

int CPostFile::MatchConnectivity( int written )
{
  return (written == _Connectivity()) ? 0 : (((written-1) == _Connectivity()) ? 1 : 2); 
}

/*
 *  class CPostAscii
 */

CPostAscii::CPostAscii()
{
  m_file = NULL;
}

CPostAscii::~CPostAscii()
{
  Close();
}

int CPostAscii::Open( GP_CONST char * name )
{
  Close();
  m_file = fopen(name, "w");
  return m_file == NULL;
}

int CPostAscii::Close()
{
  if ( m_file ) {
    CPostFile::fail = fclose( m_file );
    m_file = NULL;
  } else
    CPostFile::fail = 1;
  return CPostFile::fail;
}

int CPostAscii::Flush()
{
  return m_file ? fflush(m_file) : 1;
}

int CPostAscii::IsBinary()
{
  return 0;
}

int CPostAscii::WriteString( GP_CONST char * str )
{
  fprintf(m_file, "%s\n", str );
  return 0;
}

int CPostAscii::WriteInteger(int i, int op)
{
  if (op==1) {
    fprintf(m_file, " ");
    fprintf(m_file, "%d", i);
  } else {
    fprintf(m_file, "%d", i);
    if (op==2) {
      fprintf(m_file, "\n");
    }
  }
  return 0;
}

int CPostAscii::WriteDouble(double x, int op)
{
  assert(op!=0);
  fprintf(m_file, " %g", x);
  if (op==2) {
    fprintf(m_file, "\n");
  }
  return 0;
}

int CPostAscii::WriteValues(int id, int n, ... )
{
  va_list ap;
  int i;
  double value;

  if ( _LastID() != id )
    fprintf(m_file, "%d", id);
  va_start(ap, n);
  for ( i = 0; i < n; i++ ) {
    value = va_arg(ap, double);
    fprintf(m_file, " %g", value);
  }
  fprintf(m_file, "\n");
  va_end(ap);
  _LastID() = id;
  return 0;
}

int CPostAscii::WriteValues( int id, int n, double * buffer)
{
  int i;

  if ( _LastID() != id )
    fprintf(m_file, "%d", id);
  for ( i = 0; i < n; i++ )
    fprintf(m_file, " %g", buffer[i]);
  fprintf(m_file, "\n");
  _LastID() = id;
  return 0;
}

int CPostAscii::WriteElement( int id, int n, int nid[] )
{
  int i;

  fprintf(m_file, "%d", id);
  for ( i = 0; i < n; i++ ) {
    fprintf(m_file, " %d", nid[i]);
  }
  fprintf(m_file, "\n");
  return 0;
}

/*
 *  class CPostAsciiZ --
 */

CPostAsciiZ::CPostAsciiZ()
{
  m_file = NULL;
}

CPostAsciiZ::~CPostAsciiZ()
{
  Close();
}

int CPostAsciiZ::Open( GP_CONST char * name )
{
  Close();
  m_file = gzopen(name, "w1");  
  return m_file == NULL;
}

int CPostAsciiZ::Close()
{
  return m_file ? gzclose(m_file) : 1;
}

int CPostAsciiZ::Flush()
{
  return m_file ? gzflush(m_file,Z_FINISH) : 1;
}

int CPostAsciiZ::IsBinary()
{
  return 0;
}

int CPostAsciiZ::WriteString( GP_CONST char * str )
{
  gzprintf( m_file, "%s\n", str );
  return 0;
}

int CPostAsciiZ::WriteInteger(int i, int op)
{
  if (op==1) {
    gzprintf(m_file, " ");
    gzprintf(m_file, "%d", i);
  } else {
    gzprintf(m_file, "%d", i);
    if (op==2) {
      gzprintf(m_file, "\n");
    }
  }
  return 0;
}

int CPostAsciiZ::WriteDouble(double x, int op)
{
  assert(op!=0);
  gzprintf(m_file, " %g", x);
  if (op==2) {
    gzprintf(m_file, "\n");
  }
  return 0;
}

int CPostAsciiZ::WriteValues( int id, int n, ... )
{
  va_list ap;
  int i;
  double value;

  if ( _LastID() != id )
    gzprintf(m_file, "%d", id);
  va_start(ap, n);
  for ( i = 0; i < n; i++ ) {
    value = va_arg(ap, double);
    gzprintf(m_file, " %g", value);
  }
  gzprintf(m_file, "\n");
  va_end(ap);
  _LastID() = id;
  return 0;
}
 
int CPostAsciiZ::WriteValues(int id, int n, double * buffer)
{
  int i;

  if ( _LastID() != id )
    gzprintf(m_file, "%d", id);
  for (i = 0; i < n; i++)
    gzprintf(m_file, " %g", buffer[i]);
  gzprintf(m_file, "\n");
  _LastID() = id;
  return 0;
}
 
int CPostAsciiZ::WriteElement( int id, int n, int nid[] )
{
  int i;

  gzprintf(m_file, "%d", id);
  for ( i = 0; i < n; i++ ) {
    gzprintf(m_file, " %d", nid[i]);
  }
  gzprintf(m_file, "\n");
  return 0;
}
 
/*
 *  class CPostBinary --
 */

CPostBinary::CPostBinary()
{
  m_file = NULL;
}

CPostBinary::~CPostBinary()
{
  Close();
}

int CPostBinary::Open( GP_CONST char * name )
{
  Close();
  m_file = gzopen(name, "wb1");
  /* escribir el numero magico */
  gzwrite(m_file, &ByteOrderCheck, sizeof(ByteOrderCheck));
  return m_file == NULL;
}

int CPostBinary::Close()
{
  if (m_file) {
    CPostFile::fail = gzclose(m_file);
    m_file = NULL;
  } else
    CPostFile::fail = 1;
  return CPostFile::fail;
}

int CPostBinary::Flush()
{
  return m_file ? gzflush(m_file,Z_FULL_FLUSH) : 1;
}

int CPostBinary::IsBinary()
{
  return 1;
}

int CPostBinary::WriteString( GP_CONST char * str )
{
  int _fail = 1;
  int size, written;

  if ( m_file ) {
    GP_CONST char *buf = "\0";
    if (str) buf = str;
    int tam = strlen(buf) + 1; // incluido el \0
    written = gzwrite(m_file, &tam, sizeof(int));
    size = sizeof(char) * tam;
    written += gzwrite(m_file, (void*)buf, size);
    if ( written == size+int(sizeof(int)) )
      _fail = 0;
  }
  return _fail;
}

int CPostBinary::BeginCoordinates()
{
  return WriteString( "Coordinates -1 Indexed"); 
}

int CPostBinary::BeginElements()
{
  return WriteString( "Elements -1 Indexed"); 
}

int CPostBinary::BeginValues()
{
  return WriteString("Values -1 Indexed");
}

int CPostBinary::EndValues()
{
  int idxend = -1;

  if ( gzwrite(m_file, &idxend, sizeof(int)) != sizeof(int) )
    return 1;
  return WriteString("End Values");
}

int CPostBinary::WriteInteger(int i, int)
{
  gzwrite(m_file, &i, sizeof(i));
  return 0;
}

int CPostBinary::WriteDouble(double x, int)
{
  float v = float(x);
  
  gzwrite(m_file, &v, sizeof(v)); 
  return 0;
}

int CPostBinary::WriteValues( int id, int n, ... )
{
  va_list ap;
  int i;
  float value;

  if ( _LastID() != id ) {
    gzwrite(m_file, &id, sizeof(id));
    _LastID() = id;
  }
  va_start(ap, n);
  for ( i = 0; i < n; i++ ) {
    value = (float)va_arg(ap, double);
    gzwrite(m_file, &value, sizeof(value));
  }
  va_end(ap);
  return 0;  
} 

int CPostBinary::WriteValues(int id, int n, double * buffer)
{
  int i;
  float value;

  if ( _LastID() != id ) {
    gzwrite(m_file, &id, sizeof(id));
    _LastID() = id;
  }
  for ( i = 0; i < n; i++ ) {
    value = (float) buffer[i];
    gzwrite(m_file, &value, sizeof(value));
  }
  return 0;  
} 

int CPostBinary::Write2D( double x, double y )
{
  float values[2];
  
  values[0] = float(x);
  values[1] = float(y);
  gzwrite(m_file, values, sizeof(float)*2);
  return 0;
}

int CPostBinary::Write3D( double x, double y, double z )
{
  float values[3];
  
  values[0] = float(x);
  values[1] = float(y);
  values[2] = float(z);
  gzwrite(m_file, values, sizeof(float)*3);
  return 0;
}

int CPostBinary::WriteElement( int id, int n, int nid[] )
{
  int i;
  
  gzwrite(m_file, &id, sizeof(id));
  for ( i = 0; i < n; i++ ) {
    gzwrite(m_file, nid+i, sizeof(int));
  }
  /* verify connectivity */
  switch ( MatchConnectivity(n) ) {
    case 0:
      /* match exactly */
      i = 1;
      /* so write material 1 */
      gzwrite(m_file, &i, sizeof(int));
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

int CPostBinary::WritePostHeader()
{
  return WriteString( "GiDPostEx1.1");
}

CBufferValues::CBufferValues()
{
  file = NULL;
  buffer_values = NULL;
  last_value = -1;
  values_size_min = values_size_max = buffer_values_size = 0;
  buffer_types = (GiD_ResultType*)malloc(sizeof(GiD_ResultType)*10);
  buffer_types_size = 10;
  size_types = 0;
  next_type = 0;
}

CBufferValues::~CBufferValues()
{
  if (buffer_values) {
    free(buffer_values);
  }
  if (buffer_types) {
    free(buffer_types);
  }
}

void CBufferValues::OnNewType(GiD_ResultType t) 
{
  size_t min, max;
  
  if (size_types==buffer_types_size) {
    buffer_types_size += 10;
    buffer_types = (GiD_ResultType*)realloc(buffer_types,
                                            sizeof(GiD_ResultType)*buffer_types_size);
  }
  buffer_types[size_types++] = t;
  GetResultTypeMinMaxValues(t, min, max);
  values_size_min += min;
  values_size_max += max;
}

void CBufferValues::OnBeginValues()
{
  if (!buffer_values)
    buffer_values = (double*)malloc((buffer_values_size=values_size_max)*sizeof(double));
  else if (values_size_max > buffer_values_size)
    buffer_values = (double*)realloc(buffer_values, (buffer_values_size=values_size_max)*sizeof(double));
  last_value = -1;
}

int CBufferValues::FlushValues(int id) 
{
  assert(file);
  assert(last_value>=values_size_min-1 || last_value<=values_size_max-1);
  if (file->WriteValues(id, last_value+1, buffer_values)) {
    /* could not write buffer values */
    last_value = -1;
    return 1;
  }
  /* prepare for next group of values */
  last_value = -1;
  return 0;
}

int CBufferValues::WriteValues(GiD_ResultType t, int id, int n, ...)
{
  va_list ap;
  int i;
  int flush;
  
  assert(last_value+n<values_size_max);

  flush = OnWriteType(t);
  if (flush==-1)
    return -1;
  va_start(ap, n);
  for (i = 0; i < n; i++)
    buffer_values[++last_value] = va_arg(ap, double);
  va_end(ap);
  
  if (flush)
    return FlushValues(id);
  return 0;  
}

int CBufferValues::OnWriteType(GiD_ResultType t) 
{
  if (t==buffer_types[next_type++]) {
    if (next_type == size_types) {
      next_type = 0;
      return 1;
    }
    return 0;
  } else {
    printf("error expected '%s' instead of '%s'\n",
           GetResultTypeName(buffer_types[next_type-1]),
           GetResultTypeName(t));
    return -1;
  }
}
