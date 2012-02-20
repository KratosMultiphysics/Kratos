/* gidpost 1.7 */

#include <assert.h>
#include <stdarg.h>
#include "cfortran/cfortran.h"
#include "gidpost.h"

#include "gidpostfor.h"

#ifdef fcallsc
#undef fcallsc
#endif
#define fcallsc(UN,LN) append_fcallsc(_,_,UN,LN)

#include "gidpostfor.h"

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

int GiD_ScalarComp(GP_CONST char * Comp1)
{
  assert(Comp1 && *Comp1);

  if (Comp1 && *Comp1) {
    return GiD_ResultComponents(1, &Comp1);
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

int GiD_2DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_3DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		     GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[6];

  ncomp = SetupComponents(6, CompBuffer, Comp1, Comp2, Comp3, Comp4, Comp5, Comp6);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}

int GiD_PDMComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3, GP_CONST char * Comp4)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[4];

  ncomp = SetupComponents(4, CompBuffer, Comp1, Comp2, Comp3, Comp4);
  return GiD_ResultComponents(ncomp, CompBuffer);  
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

int GiD_LAComponents(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3)
{
  int ncomp = 0;
  GP_CONST char * CompBuffer[3];

  ncomp = SetupComponents(3, CompBuffer, Comp1, Comp2, Comp3);
  return GiD_ResultComponents(ncomp, CompBuffer);  
}
