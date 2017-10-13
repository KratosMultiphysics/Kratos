#include "u_Types.h"
#include "u_qualityMetrics.h"

#define  CRITERIA_OPT_MIN 0
#define	 CRITERIA_OPT_AVG 1
#define CRITERIA_OPT_AVG_AND_MIN 2

void iterativeMeshSmooth(TMesh* aMesh, int mxIter, int subI, double minExpectedQ = 10000) ;
void laplacianMeshSmooth(TMesh* aMesh);
void wightedlaplacianMeshSmooth(TMesh* aMesh);
void localMeshSmooth(int i, int thId  ,TObject* destObject);
void setSmoothParams(double gbm , int mxIter , int subI);