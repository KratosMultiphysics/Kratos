
#include "Math3D.h"
#include <math.h>
#include <stdio.h>
#include "u_ShowMetrics.h"
#include "u_TetraFunctions.h"
#include "u_qualityMetrics.h"
#include <iostream>
#include <fstream>

using namespace std;

void exportMetrics(std::string outfilename , TMesh* m,double time)
{

	ofstream myfile;
	myfile.open (outfilename.c_str());
	myfile<<"Req time "<<time<<"\n";
	for (int i=0; i<m->elements->Count();i++)
	{
		//-------------------------------
		TTetra *t = (TTetra*)(m->elements->elementAt(i));
		myfile<<i<<"\t"<<vrelaxQuality(t->vertexes)<<"\t"<<diedralAngle(t->vertexes)<<"\n";
	}
	myfile.close();
}

TetQuality::TetQuality(TMesh* volMesh)
{ 
	this->aMesh = volMesh;
}

void TetQuality::getQuality()
{
	nonPositive = 0;
	fFaceMax = fVolMax = fDieMax = -100000000;
	fFaceMin = fVolMin = fDieMin = 100000000;
	FaceAveMin = 0;
	DieAveMin = 0;
	DieAveMax= 0;
	fMinQuality =  100000000;
	fMaxQuality = -100000000;

	for (int i=0 ; i< aMesh->elements->Count();i++)
	{
		TTetra* t = (TTetra*)(aMesh->elements->elementAt(i));
		if (t == NULL) 
			continue;
		t->update();
		this->fDieMax = Max(fDieMax,t->fMaxDiedralAngle);
		this->fDieMin = Min(fDieMin,t->fMinDiedralAngle);
		this->fFaceMax = Max(fFaceMax,t->fFaceAngle);
		this->fFaceMin = Min(fFaceMin,t->fFaceAngle);
		this->fVolMax = Max(fVolMax,t->fVolume);
		this->fVolMin = Min(fVolMin,t->fVolume);										
		double q = vrelaxQuality(t->vertexes);
		this->fMinQuality = Min(fMinQuality,q);	
		this->fMaxQuality = Max(fMaxQuality,q);	
		FaceAveMin = FaceAveMin + t->getminEdgeLength() + t->getmaxEdgeLength();
		DieAveMin = DieAveMin + t->fMinDiedralAngle;
		if (t->fVolume<-0) nonPositive++;
	}
	FaceAveMin =  FaceAveMin /(2*aMesh->elements->Count());

	DieAveMin = DieAveMin / aMesh->elements->Count();


}

void TetQuality::refresh() 
{ 
	this->getQuality(); 
}

void TetQuality::print()
{
	std :: cout <<"----------------------------------------"<<"\n";
	std :: cout << "Mesh quality: " <<"\n";
	std :: cout << "Diedral Angle. Min :" <<this->fDieMin<< " Max :" <<this->fDieMax <<"\n";	
	std :: cout << "Face Angle. Min :" <<this->fFaceMin << " Max :" <<this->fFaceMax <<"\n";
	std :: cout << "Volume. Min :" <<this->fVolMin<< " Max :" <<this->fVolMax <<"\n";	
	std :: cout << "Quality. Min :" <<this->fMinQuality<< " Max :" <<this->fMaxQuality <<"\n";	
	std :: cout << "Average diedral Angle :" <<this->DieAveMin <<"\n";
	std :: cout << "Num negative elements :" <<this->nonPositive <<"\n";
	int ne = this->aMesh->elements->Count();
	int nv = this->aMesh->vertexes->Count();
	int nt = 0;
	if (aMesh->fFaces )
		nt = aMesh->fFaces->Count();
	std :: cout << "#Elements :" << ne << " #Vertexes: "<< nv << " #Faces: " << nt <<"\n";
	
	

}
