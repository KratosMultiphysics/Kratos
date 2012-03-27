#include "u_ParallelFunctions.h"
#include "u_elementCluster.h"
#include "u_ProcessTime.h"
#include "u_qualityMetrics.h"
#include <cstdio>

bool* _faces;


void localProcessI(int i, int thId  ,TObject* destObject)
{
  int j ,k;
  TVertex* v0, *v1, *v2;
  TTetra *t,*t2;
  bool res;

    t = (TTetra*)( destObject);
	if (t ==  NULL ) return ;

    for (j = 0 ; j<4; j++)
	{
      _faces[i*4+j] = false;
      v0 = t->vertexes[ TTetraFaces[j*3] ];
      v1 = t->vertexes[ TTetraFaces[j*3+1] ];
      v2 = t->vertexes[ TTetraFaces[j*3+2] ];
      res = true;
	  for (k = 0 ; k<v0->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v0->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
      if (!res) continue;
      for (k = 0 ; k<v1->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v1->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
      if (!res) continue;
      for (k = 0 ; k<v2->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v2->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
       if (res)       
        _faces[i*4+j] = true;
	}
}

void fastGetSurfaceTriangles(TMesh* aModel)
{
	
  int i,j , numT;
  TVertex*v0,*v1,*v2;
  TTetra *t;
  TTriangle *tr;

  numT = aModel->elements->Count();
  _faces = new bool[numT * 4];
  for (i = 0 ; i<numT*4 ; i++)
	   _faces[i] = false;

  for (int i=0;i<aModel->fFaces->Count();i++)
  {
	  TTriangle *tr = (TTriangle*)(aModel->fFaces->structure[i]);
	  delete tr;
  }
  
  aModel->fFaces->Clear();
  
  startProcess((char*)("parallelPart"));
    //parallel_for(blocked_range<size_t>(0,numT), TParallelIterator(localProcessI,aModel->elements) );    
     	  #pragma omp parallel for
              for (i = 0; i<numT ; i++)
					localProcessI(i,0,aModel->elements->elementAt(i));
		  
	 
  endProcess((char*)("parallelPart"));

  startProcess((char*)("createPart"));
  for (i = 0 ; i<numT ; i++)
  {
	t = (TTetra*)( aModel->elements->elementAt(i));
    for (j = 0 ; j<4 ;j++)
	{
      if (!_faces[i*4+j]) continue;
      v0 = t->vertexes[ TTetraFaces[j*3] ];
      v1 = t->vertexes[ TTetraFaces[j*3+1] ];
      v2 = t->vertexes[ TTetraFaces[j*3+2] ];
      tr = new TTriangle( v0,v1,v2);
	  tr->calcNormal();
	  aModel->addTriangle(tr);
	}
  }
  endProcess((char*)("createPart"));
  delete _faces;
  
}


void lpEvaluateCluster(int i, int thId  ,TObject* destObject)
{
	TElementsCluster* resC =  (TElementsCluster*)destObject;
	 resC->generateSubMesh();
     resC->evaluateSet();
     resC->updateMesh(true);
}


void assignVertexesAvoidingVisited(TList<TVertex*> *vs, TList<TObject*> *vRes ,int iter, int maxAssignment)
{
  int i,j, nAssigned;
  TList<TObject*>  *lneigh;
  bool wasVisited;

  lneigh = new TList<TObject*>();
  /*
  for (i = 0 ; i<vs->Count() ; i++)
  {
	  ((TVertex*)(vs->elementAt(i)))->marked = false;
  }
  */
  nAssigned = 0;
  for (i = 0 ; i<vs->Count() ; i++)
  {
     TVertex *v = vs->elementAt(i);
	 if (v == NULL ) continue;
	 // Ya fue visitado en una iteracion previa
	 if (v->visited>0) continue;
	 // Ya lo marco otro thread!!
	 if (v->flag == iter) continue;	    
	 if (v->isdestroyed) continue;

	 lneigh->Clear();
	 v->getElemNeighbours(lneigh);    
	 if (lneigh->Count() == 0 ) continue;
     
     wasVisited = false;
     // Chequeo de no ver los vecinos
	 for (j = 0 ; j<lneigh->Count() ; j++)
	 {
		 TVertex* v2 = (TVertex*)(lneigh->elementAt(j));
		 if (v2 == NULL ) continue;
		 if (v2->flag == iter)
		 {
           wasVisited = true;
           break ;
		 }
	 }
     //
     if (wasVisited)  continue;

	 v->visited = iter;
	 v->flag = iter;
	 vRes->Add(v);
     
     for (j = 0 ; j<lneigh->Count() ; j++)
	 {
		 TVertex* v2 = (TVertex*)(lneigh->elementAt(j));
		 if (v2 == NULL ) continue;
		 v2->flag = iter;		 
	 }
	 
     nAssigned ++;
     if (vRes->Count()  >= maxAssignment) break;         
  }
}


void ParallelEvaluateClusterByNode(TMesh *aMesh , TVertexesEvaluator fc)
{  
  int	iv ,i ,nsimCh;
  TList<TObject*> *inspectedElements, *vRes;
  TList<TVertex*> *vertexesCopy;
  TList<TObject*> * resultedClusters; 
  TVertex *inspVertex ; 
  
  nsimCh = 2048;

  resultedClusters = new TList<TObject*>();
  //setLength( resultedClusters ,nsimCh);

 for (i = 0 ; i<nsimCh ; i++)
 {
    TElementsCluster* e = new TElementsCluster(aMesh,vrelaxQuality) ;
	 resultedClusters->Add( (TObject*)(e));
 }

 vertexesCopy = new TList<TVertex*>();
 vRes = new TList<TObject*>();
 
 vertexesCopy->Assign(aMesh->vertexes);
 for (i = 0 ;i<vertexesCopy->Count() ; i++)
 {
	 vertexesCopy->elementAt(i)->visited = 0;
	 vertexesCopy->elementAt(i)->flag = 0;
	 vertexesCopy->elementAt(i)->isdestroyed = false;
 }
 TStringList* st = new TStringList();
 startProcess((char*)("evaluateClustersInParallel"));
 char* iter = new char[0];
 TParallelIterator* pi = new TParallelIterator();
 
  
 ////-- Facil de Paralelizar
 for (iv = 0  ; iv<=Max(20.0f,(float)(vertexesCopy->Count() / nsimCh)-1 ) ; iv++) 
 {
     st->Add(iter,intToStr(iv));
	 vRes->Clear();
     //Distribuyo la carga
     //assignVertexes(vertexesCopy,vRes,iv+1,nsimCh);
	  startProcess((char*)("assignVertexesAvoidingVisited"));
        assignVertexesAvoidingVisited(vertexesCopy,vRes,iv+1,nsimCh-1);
	  endProcess((char*)("assignVertexesAvoidingVisited"));
      if (vRes->Count() == 0 ) break;
      
	  startProcess((char*)("clearVars"));
    //--Limpio las variables
	  // por cada vertice, tengo un cluster
	  for (i = 0 ; i<vRes->Count() ; i++)
	  {
		 inspVertex =(TVertex*)(vRes->elementAt(i));
         inspectedElements = inspVertex->elementsList;
		 inspVertex->elementsList->Pack();
		 TElementsCluster* resC = (TElementsCluster*)(resultedClusters->elementAt(i));
		 resC->inspectedElements->Assign( inspectedElements) ;         
         resC->doRemoveElements = false;
         resC->testCenter = false;
         resC->perturbCenter = 0;
         resC->checkMaxLength = true;        
	  }
	  endProcess((char*)("clearVars"));

	  startProcess((char*)("parallelPart"));
	  //for (int k = 0; k<vRes->Count() ; k++)
	 //	  lpEvaluateCluster(k,0,resultedClusters->elementAt(k));
  	     pi->forloop(0,vRes->Count()-1,resultedClusters, lpEvaluateCluster);      
      endProcess((char*)("parallelPart"));
      /* end of parallel section */
      for (i = 0 ; i<vRes->Count() ; i++)
	  {
		  TElementsCluster* ec = (TElementsCluster*)( resultedClusters->elementAt(i));
		  TVertex* _v = (TVertex*)(vRes->elementAt(i));
		  st->Add(intToStr(_v->id), intToStr(ec->newElements->Count()), intToStr(ec->inspectedElements->Count()));	  
		  ec->genElements();		  
	  }
     //-- Muchos elementos para remover. Limpio las estructuras
    // if (aMesh->elementsToRemove->Count()>100000 )
	// aMesh->updateRefs();
	
 }
  endProcess((char*)("evaluateClustersInParallel"));
  
  aMesh->updateRefs();
  
}
