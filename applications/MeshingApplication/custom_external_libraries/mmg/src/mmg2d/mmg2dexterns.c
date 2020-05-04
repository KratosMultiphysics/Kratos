#include "mmgexterns.c"
#include "mmg2d.h"

int    (*MMG2D_defsiz)(MMG5_pMesh ,MMG5_pSol );
int    (*MMG2D_intmet)(MMG5_pMesh ,MMG5_pSol ,int ,char ,int ,double );
double (*MMG2D_lencurv)(MMG5_pMesh ,MMG5_pSol ,int ,int );
int    (*MMG2D_gradsizreq)(MMG5_pMesh ,MMG5_pSol );
double (*MMG2D_caltri)(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria );
int    (*MMG2D_gradsiz)(MMG5_pMesh ,MMG5_pSol );

