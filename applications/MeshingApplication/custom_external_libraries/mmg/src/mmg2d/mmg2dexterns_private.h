#ifndef MMG2DEXTERNS_H
#define MMG2DEXTERNS_H

#include "libmmgtypes.h"
#include "mmgcommon_private.h"

#ifndef MMG_EXTERN
#define MMG_EXTERN extern
#define MMG_ASSIGN_NULL
#endif

FUNCTION_POINTER ( int   (*MMG2D_defsiz)(MMG5_pMesh ,MMG5_pSol ) );
FUNCTION_POINTER ( int   (*MMG2D_intmet)(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int8_t ,MMG5_int ,double ) );
FUNCTION_POINTER ( double(*MMG2D_lencurv)(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int ) );
FUNCTION_POINTER ( int   (*MMG2D_gradsizreq)(MMG5_pMesh ,MMG5_pSol ) );
FUNCTION_POINTER ( double(*MMG2D_caltri)(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria ) );
FUNCTION_POINTER ( int   (*MMG2D_gradsiz)(MMG5_pMesh ,MMG5_pSol ) );

#undef MMG_EXTERN
#undef MMG_ASSIGN_NULL

#endif
