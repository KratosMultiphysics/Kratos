#include "u_Types.h"
#include "Math3D.h"
#include <cmath>
#include <stdio.h>

class    TetQuality 
{
public :
	double DieAveMax,DieAveMin, FaceAveMax,FaceAveMin ,
			fVolMax,fVolMin ,
			fDieMax,fDieMin, fFaceMax,fFaceMin, vol ,fMinQuality , fMaxQuality;
	int fNumChanges,fiter, nonPositive,fNumElements ;
	double fOldVMin , fOldVMax ;
	TMesh* aMesh;
	TetQuality(TMesh* volMesh);

	void getQuality();

	void refresh() ;

	void print();

};

void exportMetrics(std::string outfilename , TMesh* m,double time);
