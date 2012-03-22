#include <math.h>
#include "u_Types.h"
#include "Math3D.h"

double calidadxArea(float4 v1,float4 v2,float4 v3);
double calidadxArea(TVertex* v1,TVertex* v2,TVertex* v3);
double getmaxEdgeLength(TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3);
double vrelaxQuality(TVertex* v0, TVertex* v1,TVertex* v2,TVertex* v3);
double vrelaxQuality(TVertex* vertexes[]);
float relaxQuality(TObject* o) ;
void CalcAng(float4 coor0,float4 coor1,float4 coor2,float4 coor3 , double* d, double *c, double *vol);
double diedralAngle(float4 v0, float4 v1, float4 v2, float4 v3);
