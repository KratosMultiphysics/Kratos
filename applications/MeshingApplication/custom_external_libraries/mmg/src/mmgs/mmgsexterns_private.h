#ifndef MMGSEXTERNS_H
#define MMGSEXTERNS_H

#include "libmmgtypes.h"
#include "mmgcommon_private.h"

#ifndef MMG_EXTERN
#define MMG_EXTERN extern
#define MMG_ASSIGN_NULL
#endif

FUNCTION_POINTER ( double (*MMG5_calelt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) );
FUNCTION_POINTER ( int    (*MMGS_defsiz)(MMG5_pMesh mesh,MMG5_pSol met) );
FUNCTION_POINTER ( int    (*MMGS_gradsiz)(MMG5_pMesh mesh,MMG5_pSol met) );
FUNCTION_POINTER ( int    (*MMGS_gradsizreq)(MMG5_pMesh mesh,MMG5_pSol met) );
FUNCTION_POINTER ( int    (*intmet)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,double s) );
FUNCTION_POINTER ( int    (*movintpt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *list,int ilist) );
FUNCTION_POINTER ( int    (*movridpt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *list,int ilist) );

#undef MMG_EXTERN
#undef MMG_ASSIGN_NULL

#endif
