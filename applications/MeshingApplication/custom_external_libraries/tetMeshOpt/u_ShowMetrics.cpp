
#include "Math3D.h"
#include <math.h>
#include <stdio.h>
#include "u_ShowMetrics.h"
#include "u_TetraFunctions.h"



TetQuality::TetQuality(TMesh* volMesh)
			{ 
				 this->aMesh = volMesh;
		    }
        void TetQuality::getQuality()
		{
			nonPositive = 0;
			fFaceMax = fVolMax = fDieMax = -100000000;
			fFaceMin = fVolMin = fDieMin = 100000000;
			  for (int i=0 ; i< aMesh->elements->Count();i++)
			  {
				    TTetra* t = (TTetra*)(aMesh->elements->elementAt(i));
					t->update();
					this->fDieMax = Max(fDieMax,t->fDiedralAngle);
					this->fDieMin = Min(fDieMin,t->fDiedralAngle);
					this->fFaceMax = Max(fFaceMax,t->fFaceAngle);
					this->fFaceMin = Min(fFaceMin,t->fFaceAngle);
					this->fVolMax = Max(fVolMax,t->fVolume);
					this->fVolMin = Min(fVolMin,t->fVolume);					
					this->FaceAveMin = Min(FaceAveMin,t->getminEdgeLength());
					FaceAveMin = FaceAveMin + t->getminEdgeLength() + t->getmaxEdgeLength();
					if (t->fVolume<-0) nonPositive++;
			  }
			  FaceAveMin =  FaceAveMin /(2*aMesh->elements->Count());


		}
		void TetQuality::refresh() { this->getQuality(); }

		void TetQuality::print()
		{
		    std :: 	cout <<"----------------------------------------"<<"\n";
		    std :: 	cout << "Mesh quality: " <<"\n";
			std :: cout << "Min diedral Angle :" <<this->fDieMin <<"\n";
			std :: cout << "Max diedral Angle :" <<this->fDieMax <<"\n";
			std :: cout << "Min Face Angle :" <<this->fFaceMin <<"\n";
			std :: cout << "Max Face Angle :" <<this->fFaceMax <<"\n";
			std :: cout << "Min Volume :" <<this->fVolMin <<"\n";
			std :: cout << "Max Volume :" <<this->fVolMax <<"\n";
			std :: cout << "Num negative elements :" <<this->nonPositive <<"\n";
			std :: cout << "Num elements :" <<this->aMesh->elements->Count() <<"\n";
		}