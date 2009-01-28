/* gidpost 1.7 */
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

#include <stdio.h>
#include "zlib.h"

#include "gidpost.h"

#define LINE_SIZE 8192

/*GP_CONST*/ char * GetResultTypeName(GiD_ResultType type, size_t s = 0);
void GetResultTypeMinMaxValues(GiD_ResultType type, size_t &min, size_t &max);

class CPostFile
{
public:
  static int fail;
  CPostFile();
  virtual ~CPostFile();
  virtual int Open( GP_CONST char * name ) = 0;
  virtual int Close() = 0;
  virtual int Flush() = 0;
  virtual int IsBinary() = 0;
  virtual int WriteString( GP_CONST char * str ) = 0;
  virtual int BeginCoordinates();
  virtual int BeginElements();
  virtual int BeginValues();
  virtual int EndValues();
  virtual int WriteValues( int id, int n, ... ) = 0;
  virtual int WriteValues(int id, int n, double *) = 0;
  virtual int Write2D( double x, double y );
  virtual int Write3D( double x, double y, double z );
  virtual int WriteElement( int id, int n, int nid[] ) = 0;
  virtual int WritePostHeader();
  void ResetLastID() {
    _LastID() = -1;
  }
  void SetConnectivity( int nnode ) {
    _Connectivity() = nnode;
  }
  int GetConnectivity() {
    return _Connectivity();
  }
  int MatchConnectivity( int written );
protected:
  int & _LastID() {
    return m_LastID;
  }
  int & _Connectivity() {
    return m_connectivity;
  }
private:
  int m_LastID;
  int m_connectivity;
};

class CPostAscii : public CPostFile
{
public:
  CPostAscii();
  virtual ~CPostAscii();
  virtual int Open( GP_CONST char * name );
  virtual int Close();
  virtual int Flush();
  virtual int IsBinary();
  virtual int WriteString( GP_CONST char * str );
  virtual int WriteValues( int id, int n, ... );  
  virtual int WriteValues(int id, int n, double *);
  virtual int WriteElement( int id, int n, int nid[] );
protected:
private:
  FILE * m_file;
};

class CPostAsciiZ :  public CPostFile
{
public:
  CPostAsciiZ();
  virtual ~CPostAsciiZ();
  virtual int Open( GP_CONST char * name );
  virtual int Close();
  virtual int Flush();
  virtual int IsBinary();
  virtual int WriteString( GP_CONST char * str );
  virtual int WriteValues( int id, int n, ... );  
  virtual int WriteValues(int id, int n, double *);
  virtual int WriteElement( int id, int n, int nid[] );
protected:
private:
  gzFile m_file;
};

class CPostBinary : public CPostFile
{
public:
  CPostBinary();
  virtual ~CPostBinary();
  virtual int Open( GP_CONST char * name );
  virtual int Close();
  virtual int Flush();
  virtual int IsBinary();
  virtual int WriteString( GP_CONST char * str );
  virtual int BeginCoordinates();
  virtual int BeginElements();
  virtual int BeginValues();
  virtual int EndValues();
  virtual int WriteValues( int id, int n, ... );  
  virtual int WriteValues(int id, int n, double *);
  virtual int Write2D( double x, double y );
  virtual int Write3D( double x, double y, double z );
  virtual int WriteElement( int id, int n, int nid[] );
  virtual int WritePostHeader();
protected:
private:
  gzFile m_file;
};

class CBufferValues {
public:
  CBufferValues();
  ~CBufferValues();
  void OnBeginResultGroup() 
    {
      size_types = 0;
      next_type = 0;
      values_size_min = 0;
      values_size_max = 0;
    }
  int NumberOfResults() 
    {
      return size_types;
    }
  void OnNewType(GiD_ResultType t);
  void OnBeginValues();
  double * BufferValues() {
    return buffer_values;
  }
  int IsEmpty() {
    return !buffer_values || last_value==-1;
  }
  int FlushValues(int id);
  int WriteValues(GiD_ResultType t, int id, int n, ... );
  void SetFile(CPostFile * f) {
    file = f;
  }
protected:
  int OnWriteType(GiD_ResultType t);
private:
  double * buffer_values;
  int last_value;
  int values_size_min;
  int values_size_max;
  int buffer_values_size;
  CPostFile * file;
  GiD_ResultType* buffer_types;
  int buffer_types_size;
  int size_types;
  int next_type;
};

#endif
