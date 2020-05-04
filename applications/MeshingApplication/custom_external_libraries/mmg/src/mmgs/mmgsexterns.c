#include "mmgexterns.c"
#include "mmgs.h"

int    (*movintpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    (*MMGS_defsiz)(MMG5_pMesh mesh,MMG5_pSol met);
double (*MMG5_calelt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
int    (*MMGS_gradsiz)(MMG5_pMesh mesh,MMG5_pSol met);
int    (*MMGS_gradsizreq)(MMG5_pMesh mesh,MMG5_pSol met);
int    (*intmet)(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    (*movridpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
