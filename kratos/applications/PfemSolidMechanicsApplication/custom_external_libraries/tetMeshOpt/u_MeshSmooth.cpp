#include "u_Types.h"
#include "u_TetraFunctions.h"
#include "u_qualityMetrics.h"
#include "u_MeshSmooth.h"
//-------------------------------
//----- wightedlaplacianMeshSmooth  :INCOMPLETE
void wightedlaplacianMeshSmooth(TMesh* aMesh)
{
	  int i,j;
  TList<TVertex*> *vList;
  TList<TObject*> *eList;
  TTetra *t;
  float4 centerPos ;
  bool hasNegative;
  double minq, q;
  int nchanges = 0;

   for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 1;
	   tr->vertexes[1]->fixed = 1;
	   tr->vertexes[2]->fixed = 1;
   }
   vList = new TList<TVertex*>();   
   
   for (i=0; i<aMesh->vertexes->Count();i++)
   {
	   TVertex *v = aMesh->vertexes->elementAt(i);
	   if (v->fixed == 1) continue;
	   eList = v->elementsList;
	   vList->Clear();
	   v->getVertexNeighboursByElem(vList);	   
	   if (vList->Count() == 0) continue;
	   //calculo la calidad minima del set
		minq = 1000000;
		double totQ = 0;
		hasNegative = false;
		centerPos = Float4(0.0);
		// Calculo la calidad inicial del cluster
		for (j=0 ; j<eList->Count() ; j++)
		{
			t = (TTetra*)(eList->elementAt(j)) ;
			if (t->getVolume() < 0) hasNegative = true;
			q = diedralAngle(t->vertexes);
			totQ +=q;			
		}
		
   }
   
   delete vList;
   if (nchanges == 0) return;
   for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 0;
	   tr->vertexes[1]->fixed = 0;
	   tr->vertexes[2]->fixed = 0;
   }
   aMesh->updateIndexes(0);
}

//-------------------------------
//----- laplacianMeshSmooth
void laplacianMeshSmooth(TMesh* aMesh) 
{
  int i,j;
  TList<TVertex*> *vList;
  TList<TObject*> *eList;
  TTetra *t;
  float4 centerPos ,origPos;
  bool hasNegative,hasNegative2;
  double minq,minq2, q;
  int nchanges = 0;

   for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 1;
	   tr->vertexes[1]->fixed = 1;
	   tr->vertexes[2]->fixed = 1;
   }
   vList = new TList<TVertex*>();   
   
   for (i=0; i<aMesh->vertexes->Count();i++)
   {
	   TVertex *v = aMesh->vertexes->elementAt(i);
	   if (v->fixed == 1) continue;
	   eList = v->elementsList;
	   vList->Clear();
	   v->getVertexNeighboursByElem(vList);	   
	   if (vList->Count() == 0) continue;
	   //calculo la calidad minima del set
		minq = 1000000;
		
		hasNegative = false;
		centerPos = Float4(0.0);
		for (j=0 ; j<vList->Count() ; j++)
		{
			TVertex *v2 = (TVertex*)(vList->elementAt(j)) ;
			centerPos =centerPos + v2->fPos;
		}
		// Calculo la calidad inicial del cluster
		for (j=0 ; j<eList->Count() ; j++)
		{
			t = (TTetra*)(eList->elementAt(j)) ;
			if (t->getVolume() < 0) hasNegative = true;

			q = diedralAngle(t->vertexes);
			minq =Min(minq,q );			
		}
		//Promedio la posicion
		centerPos = centerPos * (1.0f / vList->Count());
		origPos = v->fPos;

		v->fPos = centerPos;
		minq2 = 1000000;
		hasNegative2 = false;
		// Calculo nuevamente la calidad del cluster
		for (j=0 ; j<eList->Count() ; j++)
		{
			t = (TTetra*)(eList->elementAt(j)) ;
			if (t->getVolume() < 0) hasNegative2 = true;

			q = diedralAngle(t->vertexes);
			minq2 =Min(minq2,q );
					
		}
		 // Si mejoro!
		// Acepto 3 casos : NegNeg - PosPos - NegPos 
	     if ( (minq2>minq) && (( hasNegative && hasNegative2 ) || (!hasNegative && !hasNegative2)  || (hasNegative && !hasNegative2)) )
							{
								v->fPos = v->fPos;
								nchanges++;
								   
							   }
		 else
		 {
			 v->fPos = origPos;
		 }
		
   }
   
   delete vList;
   if (nchanges == 0) return;
   for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 0;
	   tr->vertexes[1]->fixed = 0;
	   tr->vertexes[2]->fixed = 0;
   }
   aMesh->updateIndexes(0);

}

double gbminExpectedQ ;
int gbmxIter;
int gbsubI;

void setSmoothParams(double gbm , int mxIter , int subI)
{
	gbminExpectedQ =gbm;
	gbmxIter = mxIter;
	gbsubI = subI;

}


void localMeshSmooth(int i, int thId  ,TObject* destObject)
{
  int j, nproposals,iter,criteria;
  TList<TObject*> *vList;
  TTetra *t;
  float4 initialPos ,proposed;
  bool hasNegative;
  double minq,avgq ,minq2, avgq2 ,q,minVol,radius;

	   TVertex *v = (TVertex*)destObject;
	   if (v->fixed == 1) return;
	   vList = v->elementsList;
	   if (vList == NULL) return ;
	   if (vList->Count() == 0) return;
	   //calculo la calidad minima del set
		minq = 1000000;
		avgq = 0;
		radius = 0;
		hasNegative = false;
		for (j=0 ; j<vList->Count() ; j++)
		{
			t = (TTetra*)(vList->elementAt(j)) ;
			if (t->getVolume() < 0) hasNegative = true;

			q = diedralAngle(t->vertexes);
			minq =Min(minq,q );
			avgq += q;
			radius += (t->getmaxEdgeLength() + t->getminEdgeLength())*0.5;
		}

	     avgq = avgq /vList->Count();
		 radius /=2*vList->Count();
		 if (minq > gbminExpectedQ) return;

	   criteria =CRITERIA_OPT_MIN;
		//pruebo variando en un entorno
		initialPos =v->fPos;
		double randX, randY, randZ;
		for (iter = 0 ; iter<gbmxIter ; iter++)
		{
			for (nproposals = 0 ;nproposals<gbsubI ;nproposals++)
			{
				randX =  (double)(rand()*1.0/RAND_MAX);
				randY =  (double)(rand()*1.0/RAND_MAX);
				randZ =  (double)(rand()*1.0/RAND_MAX);
				proposed =initialPos + Float4( (randX-0.5)*radius,(randY-0.5)*radius,(randZ-0.5)*radius);
				v->fPos =  proposed;
				minq2 = 1000000;
				avgq2 = 0;
				minVol = 100000000;
				for (j=0 ; j<vList->Count() ; j++)
				{	
					t = (TTetra*)(vList->elementAt(j)) ;
					q = diedralAngle(t->vertexes);
					minq2 =Min(minq2,q );
					avgq2 += q;
					minVol = Min(minVol, tetraVolume(t->vertexes[0]->fPos,t->vertexes[1]->fPos,t->vertexes[2]->fPos,t->vertexes[3]->fPos) );
				}

				if (minVol<0) continue;
				
				avgq2 = avgq2 /vList->Count();


       //si mejora lo acepto
			if (criteria == CRITERIA_OPT_MIN)
			{
                               if (minq2>minq)
							   {
								   initialPos = v->fPos;
									minq = minq2;
							   }
			}
			else if (criteria == CRITERIA_OPT_AVG)
			{
                               if (avgq2>avgq) 
							   {
                                 initialPos = v->fPos;
                                 avgq = avgq2;
							   }
			}
			else if (criteria == CRITERIA_OPT_AVG_AND_MIN)
			{
                               if ((minq2>minq) && (avgq2>avgq )   )
							   {
                                  initialPos = v->fPos;
									avgq = avgq2;
									minq = minq2;
							   }
			}

       

			}
		}
		v->fPos = initialPos;
}



void iterativeMeshSmooth(TMesh* aMesh, int mxIter, int subI, double minExpectedQ) 
{
  int i;  

   for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 1;
	   tr->vertexes[1]->fixed = 1;
	   tr->vertexes[2]->fixed = 1;
   }
   gbmxIter = mxIter;
   gbsubI = subI;
   gbminExpectedQ = minExpectedQ;
   for (i=0; i<aMesh->vertexes->Count();i++)
   {
	   localMeshSmooth(i,0,aMesh->vertexes->elementAt(i));
   }

   for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 0;
	   tr->vertexes[1]->fixed = 0;
	   tr->vertexes[2]->fixed = 0;
   }
   aMesh->updateIndexes(0);

}