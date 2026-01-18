#ifndef MMG3DEXTERNS_H
#define MMG3DEXTERNS_H

#include "libmmgtypes.h"
#include "PRoctree_3d_private.h"
#include "mmgcommon_private.h"

#ifndef MMG_EXTERN
#define MMG_EXTERN extern
#define MMG_ASSIGN_NULL
#endif

FUNCTION_POINTER ( double (*MMG5_lenedg)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra ) );
FUNCTION_POINTER ( double (*MMG5_lenedgspl)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra ) );
FUNCTION_POINTER ( double (*MMG5_caltet)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt) );
FUNCTION_POINTER ( double (*MMG5_caltri)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) );
FUNCTION_POINTER ( int    (*MMG3D_defsiz)(MMG5_pMesh ,MMG5_pSol ) );
FUNCTION_POINTER ( int    (*MMG3D_gradsiz)(MMG5_pMesh ,MMG5_pSol ) );
FUNCTION_POINTER ( int    (*MMG3D_gradsizreq)(MMG5_pMesh ,MMG5_pSol ) );
FUNCTION_POINTER ( int    (*MMG5_intmet)(MMG5_pMesh,MMG5_pSol,MMG5_int,int8_t,MMG5_int, double) );
FUNCTION_POINTER ( int    (*MMG5_interp4bar)(MMG5_pMesh,MMG5_pSol,MMG5_int,MMG5_int,double *) );
FUNCTION_POINTER ( int    (*MMG5_movintpt)(MMG5_pMesh ,MMG5_pSol, MMG3D_pPROctree ,int64_t *, int , int ) );
FUNCTION_POINTER ( int    (*MMG5_movbdyregpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int64_t*, int, MMG5_int*, int, int ,int) );
FUNCTION_POINTER ( int    (*MMG5_movbdyrefpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int64_t*, int, MMG5_int*, int ,int) );
FUNCTION_POINTER ( int    (*MMG5_movbdynompt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int64_t*, int, MMG5_int*, int ,int) );
FUNCTION_POINTER ( int    (*MMG5_movbdyridpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int64_t*, int, MMG5_int*, int ,int) );
FUNCTION_POINTER ( int    (*MMG5_cavity)(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int ,int64_t *,int ,double) );
FUNCTION_POINTER ( int    (*MMG3D_PROctreein)(MMG5_pMesh ,MMG5_pSol ,MMG3D_pPROctree ,MMG5_int,double ) );
FUNCTION_POINTER ( int    (*MMG3D_resetRef)(MMG5_pMesh) );
FUNCTION_POINTER ( int    (*MMG3D_setref)(MMG5_pMesh,MMG5_pSol) );
FUNCTION_POINTER ( int    (*MMG3D_snpval)(MMG5_pMesh,MMG5_pSol) );
FUNCTION_POINTER ( int    (*MMG3D_cuttet)(MMG5_pMesh,MMG5_pSol,MMG5_pSol) );

#undef MMG_EXTERN
#undef MMG_ASSIGN_NULL

#endif
