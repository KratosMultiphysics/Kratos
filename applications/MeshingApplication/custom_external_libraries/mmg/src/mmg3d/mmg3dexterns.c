#include "libmmg3d_private.h"
#include "mmg3d_export.h"


#define MMG_EXTERN
#define MMG_ASSIGN_NULL =NULL

#include  "mmg3dexterns_private.h"

LIBMMG3D_EXPORT double (*MMG3D_lenedgCoor)(double *ca,double *cb,double *sa,double *sb)=NULL;
LIBMMG3D_EXPORT int    (*MMG3D_doSol)(MMG5_pMesh mesh,MMG5_pSol met)=NULL;
