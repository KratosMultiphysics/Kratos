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
#include <stdarg.h>

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

#include "gidpost.h"

int SetupComponents( int CompNumber, GP_CONST char * CompBuffer[], ... )
{
  va_list ap;
  char * comp;
  int i;

  va_start(ap, CompBuffer);
  for (i = 0; i < CompNumber; i++) {
    comp = va_arg(ap,char*);
    if ( !comp || !*comp )
      break;
    CompBuffer[i] = comp;
  }
  va_end(ap);
  return i;
}

int GiD_BeginScalarResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where, GP_CONST char * GaussPointsName,
		          GP_CONST char * RangeTable, GP_CONST char * Comp1)
{
  int ncomp = 0;
  
  if ( Comp1 && *Comp1 )
    ncomp = 1;
  
  return GiD_BeginResult(Result, Analysis, step, GiD_Scalar, Where, GaussPointsName, RangeTable, ncomp, &Comp1);
}

int GiD_fBeginScalarResult(GiD_FILE fd,
		           GP_CONST char * Result, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where, GP_CONST char * GaussPointsName,
		           GP_CONST char * RangeTable, GP_CONST char * Comp1)
{
  int ncomp = 0;
  
  if ( Comp1 && *Comp1 )
    ncomp = 1;
  return GiD_fBeginResult(fd, Result, Analysis, step, GiD_Scalar, Where, GaussPointsName, RangeTable, ncomp, &Comp1);
}

int GiD_BeginVectorResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		          GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);

  return GiD_BeginResult(Result, Analysis, step, GiD_Vector, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBeginVectorResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		           GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		           GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];
  
  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);
  return GiD_fBeginResult(fd, Result, Analysis, step, GiD_Vector, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_Begin2DMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		         GiD_ResultLocation Where,
		         GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		         GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);

  return GiD_BeginResult(Result, Analysis, step, GiD_Matrix, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBegin2DMatResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);
  return GiD_fBeginResult(fd,Result, Analysis, step, GiD_Matrix, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_Begin3DMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		         GiD_ResultLocation Where,
		         GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		         GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		         GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[6];

  ncomp = SetupComponents(6, CompBuffer, Comp1, Comp2, Comp3, Comp4, Comp5, Comp6);

  return GiD_BeginResult(Result, Analysis, step, GiD_Matrix, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBegin3DMatResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		          GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[6];

  ncomp = SetupComponents(6, CompBuffer, Comp1, Comp2, Comp3, Comp4, Comp5, Comp6);
  return GiD_fBeginResult(fd,Result, Analysis, step, GiD_Matrix, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_BeginPDMMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		          GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);

  return GiD_BeginResult(Result, Analysis, step, GiD_PlainDeformationMatrix,
		         Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBeginPDMMatResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		           GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		           GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);
  return GiD_fBeginResult(fd, Result, Analysis, step, GiD_PlainDeformationMatrix,
		          Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_BeginMainMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		           GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
		           GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
		           GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
		           GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[12];

  ncomp = SetupComponents(12, CompBuffer,
		          Comp1, Comp2, Comp3, Comp4,
		          Comp5, Comp6, Comp7, Comp8,
		          Comp9, Comp10, Comp11, Comp12);

  return GiD_BeginResult(Result, Analysis, step, GiD_MainMatrix,
		         Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBeginMainMatResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		            GiD_ResultLocation Where,
		            GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		            GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
		            GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
		            GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
		            GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[12];

  ncomp = SetupComponents(12, CompBuffer,
		          Comp1, Comp2, Comp3, Comp4,
		          Comp5, Comp6, Comp7, Comp8,
		          Comp9, Comp10, Comp11, Comp12);

  return GiD_fBeginResult(fd, Result, Analysis, step, GiD_MainMatrix,
		          Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_BeginLAResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		      GiD_ResultLocation Where,
		      GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		      GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);

  return GiD_BeginResult(Result, Analysis, step, GiD_LocalAxes,
		         Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBeginLAResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		       GiD_ResultLocation Where,
		       GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		       GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);

  return GiD_fBeginResult(fd, Result, Analysis, step, GiD_LocalAxes,
		         Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_BeginComplexScalarResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Re, GP_CONST char * Im)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[2];

  ncomp = SetupComponents(2, CompBuffer, Re, Im);

  return GiD_BeginResult(Result, Analysis, step, GiD_ComplexScalar, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBeginComplexScalarResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		           GP_CONST char * Re, GP_CONST char * Im)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[2];
  
  ncomp = SetupComponents(2, CompBuffer, Re, Im);
  return GiD_fBeginResult(fd, Result, Analysis, step, GiD_ComplexScalar, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_BeginComplexVectorResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Rex, GP_CONST char * Imx, GP_CONST char * Rey,
		          GP_CONST char * Imy, GP_CONST char * Rez, GP_CONST char * Imz)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[6];

  ncomp = SetupComponents(6, CompBuffer, Rex, Imx, Rey, Imy, Rez, Imz);

  return GiD_BeginResult(Result, Analysis, step, GiD_ComplexVector, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}

int GiD_fBeginComplexVectorResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		           GP_CONST char * Rex, GP_CONST char * Imx, GP_CONST char * Rey,
		           GP_CONST char * Imy, GP_CONST char * Rez, GP_CONST char * Imz)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[6];
  
  ncomp = SetupComponents(6, CompBuffer, Rex, Imx, Rey, Imy, Rez, Imz);
  return GiD_fBeginResult(fd, Result, Analysis, step, GiD_ComplexVector, Where, GaussPointsName, RangeTable, ncomp, CompBuffer);
}


int GiD_ScalarComp(GP_CONST char * Comp1)
{
  assert(Comp1 && *Comp1);

  if (Comp1 && *Comp1) {
    return GiD_ResultComponents(1, &Comp1);
  }
  return 0;
}

int GiD_fScalarComp(GiD_FILE fd, GP_CONST char *Comp1)
{
  assert(Comp1 && *Comp1);

  if (Comp1 && *Comp1) {
    return GiD_fResultComponents(fd, 1, &Comp1);
  }
  return 0;
}

int GiD_VectorComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3, GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_fVectorComp(GiD_FILE fd, GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3, GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);
  return GiD_fResultComponents(fd, ncomp, CompBuffer);  
}

int GiD_2DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_f2DMatrixComp(GiD_FILE fd, GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);
  return GiD_fResultComponents(fd, ncomp, CompBuffer);  
}

int GiD_3DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		     GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[6];

  ncomp = SetupComponents(6, CompBuffer, Comp1, Comp2, Comp3, Comp4, Comp5, Comp6);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_f3DMatrixComp(GiD_FILE fd, GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		      GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[6];

  ncomp = SetupComponents(6, CompBuffer, Comp1, Comp2, Comp3, Comp4, Comp5, Comp6);
  return GiD_fResultComponents(fd, ncomp, CompBuffer);  
}

int GiD_PDMComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3, GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_fPDMComp(GiD_FILE fd,
		 GP_CONST char *Comp1,
		 GP_CONST char *Comp2,
		 GP_CONST char *Comp3,
		 GP_CONST char *Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);
  return GiD_fResultComponents(fd, ncomp, CompBuffer);  
}

int GiD_MainMatrixComp(GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
		       GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
		       GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
		       GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[12];

  ncomp = SetupComponents(12, CompBuffer,
		          Comp1, Comp2, Comp3, Comp4, Comp5, Comp6,
		          Comp7, Comp8, Comp9, Comp10, Comp11, Comp12);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_fMainMatrixComp(GiD_FILE fd,
		        GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
		        GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
		        GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
		        GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[12];

  ncomp = SetupComponents(12, CompBuffer,
		          Comp1, Comp2, Comp3, Comp4, Comp5, Comp6,
		          Comp7, Comp8, Comp9, Comp10, Comp11, Comp12);
  return GiD_fResultComponents(fd, ncomp, CompBuffer);  
}

int GiD_LAComponents(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_fLAComponents(GiD_FILE fd,
		      GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);
  return GiD_fResultComponents(fd, ncomp, CompBuffer);  
}

int GiD_ComplexScalarComp(GP_CONST char * Re,GP_CONST char * Im)
{
  assert(Re && *Re);
  assert(Im && *Im);
  if (Re && *Re && Im && *Im) {   
    GP_CONST char* compv[2];
    compv[0]=Re;
    compv[1]=Im;
    return GiD_ResultComponents(2, compv);
  }
  return 0;
}

int GiD_fComplexScalarComp(GiD_FILE fd, GP_CONST char *Re, GP_CONST char *Im)
{
  assert(Re && *Re);
  assert(Im && *Im);
  if (Re && *Re && Im && *Im) {
    GP_CONST char* compv[2];
    compv[0]=Re;
    compv[1]=Im;
    return GiD_fResultComponents(fd, 2,compv);
  }
  return 0;
}

int GiD_ComplexVectorComp(GP_CONST char * Rex,GP_CONST char * Imx,
  GP_CONST char * Rey,GP_CONST char * Imy,GP_CONST char * Rez,GP_CONST char * Imz)
{
  assert(Rex && *Rex);
  assert(Imx && *Imx);
  if (Rex && *Rex && Imx && *Imx && Rey && *Rey && Imy && *Imy && Rez && *Rez && Imz && *Imz) {   
    GP_CONST char* compv[6];
    compv[0]=Rex;
    compv[1]=Imx;
    compv[2]=Rey;
    compv[3]=Imy;
    compv[4]=Rez;
    compv[5]=Imz;
    return GiD_ResultComponents(6, compv);
  }
  return 0;
}

int GiD_fComplexVectorComp(GiD_FILE fd, GP_CONST char * Rex,GP_CONST char * Imx,
  GP_CONST char * Rey,GP_CONST char * Imy,GP_CONST char * Rez,GP_CONST char * Imz)
{
  assert(Rex && *Rex);
  assert(Imx && *Imx);
  if (Rex && *Rex && Imx && *Imx && Rey && *Rey && Imy && *Imy && Rez && *Rez && Imz && *Imz) {   
    GP_CONST char* compv[6];
    compv[0]=Rex;
    compv[1]=Imx;
    compv[2]=Rey;
    compv[3]=Imy;
    compv[4]=Rez;
    compv[5]=Imz;
    return GiD_fResultComponents(fd, 6, compv);
  }
  return 0;
}

#ifndef WIN32
#pragma GCC diagnostic pop
#endif // WIN32
