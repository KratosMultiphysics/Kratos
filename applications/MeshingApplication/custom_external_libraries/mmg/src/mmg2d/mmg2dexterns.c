#include "libmmgtypes.h"
#include "mmg2d_export.h"

#define MMG_EXTERN
#define MMG_ASSIGN_NULL =NULL

#include  "mmg2dexterns_private.h"

LIBMMG2D_EXPORT int    (*MMG2D_doSol)(MMG5_pMesh ,MMG5_pSol )=NULL;
