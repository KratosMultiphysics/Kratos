#include "stdafx.h"
#include "u_Types.h"
#include "u_qualityMetrics.h"

bool sortByID(TObject* i, TObject* j)
{
	return ((TValuedObject*)(i))->id < ((TValuedObject*)(j))->id ;
}

TList<TObject*>* TVertex::getvNeigh(int depth , int mode, TList<TObject*>* toL)
{
   TList<TObject*>* storeList ;
   
   if (toL != NULL)
        storeList = toL;   
   else
   {
     if (neighV == NULL) 
       neighV = new TList<TObject*>();
	 neighV->Clear();
     storeList = neighV;
   }
   if (mode = 0)
   {
	   for (int i=0;i< neighTr->Count() ; i++)
	   {
		   TTriangle *t = (TTriangle*)( neighTr->elementAt(i));
		   for (int j= 0; j<3 ;j++)
		   {
			   if (t->vertexes[j] == this) continue;
			   if (storeList->indexOf( t->vertexes[j])<0)
				   storeList->Add(t->vertexes[j]);
		   }
	   }
   }
   else
   {
	   for (int i=0;i< elementsList->Count() ; i++)
	   {
		   TTetra *t = (TTetra*)( elementsList->elementAt(i));
		   if (t == NULL) continue;
		   if (t->isdestroyed) continue;
		   for (int j= 0; j<4 ;j++)
		   {
			   if (t->vertexes[j] == this) continue;
			   if (storeList->indexOf( t->vertexes[j])<0)
				   storeList->Add(t->vertexes[j]);
		   }
	   }
   }
 return storeList;  
}

TList<TObject*> *TVertex::getElemNeighbours(TList<TObject*> *toL )
			 {
               int j,k;
               TVertex  *v2;
			   TList<TObject*> *storeList ;


               if (toL == NULL)
			   {
				   if (neighV = NULL)  neighV = new TList<TObject*>();
				   storeList = neighV;
			   }
			   else
				   storeList= toL;

               storeList->Clear();			   
			   
			   // Compute neighbours
			   for (j = 0 ; j<elementsList->Count() ; j++)
			   {
                   TTetra *t  = (TTetra*)(elementsList->elementAt(j));
				   if (t == NULL ) continue;
				   if (t->isdestroyed) continue;
		           for (k = 0 ; k<4 ; k++)
				   {
					   v2 =  t->vertexes[k];
					   if (v2 == this) continue;                     
					   if (storeList->indexOf(v2)<0)
						   storeList->Add(v2);
				   }
			   }
               return storeList;
}

void TTetra::update()
			{
				double f, c, v;
				CalcAng(vertexes[0]->fPos, vertexes[1]->fPos, vertexes[2]->fPos , vertexes[3]->fPos,&f,&c,&v);
				fDiedralAngle = f;
				fFaceAngle = c;
				fVolume = v;
				//fVolume =tetraVolume( ;			
			}



