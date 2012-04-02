#include "u_ParallelFunctions.h"
#include "u_elementCluster.h"
#include "u_ProcessTime.h"
#include "u_qualityMetrics.h"
#include <cstdio>

bool* _faces;

//#define KRATOS
void parallelFor(int from, int to,  TList<TObject*>* elements,TStForLoopElement call)
{
	#if !defined(KRATOS)
	parallel_for(blocked_range<size_t>(from, to), TParallelIterator(call,elements ) );
    #else if
		#pragma omp parallel for
		for (int i=from ; i<=to ; i++)
			call(i,0,elements->elementAt(i));
		   
		#pragma omp barrier
    #endif
}
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

void fastGetSurfaceTriangles(TMesh* aMesh)
{
	
  int i,j , numT;
	TVertex*v0,*v1,*v2;
	TTetra *t;
	TTriangle *tr;

	numT = aMesh->elements->Count();
	_faces = new bool[numT * 4];


#pragma omp parallel for
	for (i = 0 ; i<numT*4 ; i++)


		_faces[i] = false;
	startProcess((char*)("fastGetSurfaceTriangles : delete old triangles"));
	for (int i=0;i<aMesh->fFaces->Count();i++)
	{







		TTriangle *tr = (TTriangle*)(aMesh->fFaces->structure[i]);
		// Si tengo el pool activado, lo guardo ahi 
		if (aMesh->memPool == NULL ) 
			delete tr;
		else
			aMesh->memPool->releaseInstance(tr);
	}

	endProcess((char*)("fastGetSurfaceTriangles : delete old triangles"));




	aMesh->fFaces->Clear();

	startProcess((char*)("fastGetSurfaceTriangles : parallelPart"));
	//parallel_for(blocked_range<size_t>(0,numT), TParallelIterator(localProcessI,aMesh->elements) );    
#pragma omp parallel for
	for (i = 0; i<numT ; i++)
		localProcessI(i,0,aMesh->elements->elementAt(i));


	endProcess((char*)("fastGetSurfaceTriangles : parallelPart"));

	startProcess((char*)("fastGetSurfaceTriangles : createPart"));
	for (i = 0 ; i<numT ; i++)
	{
		t = (TTetra*)( aMesh->elements->elementAt(i));
		for (j = 0 ; j<4 ;j++)
		{
			if (!_faces[i*4+j]) continue;
			v0 = t->vertexes[ TTetraFaces[j*3] ];
			v1 = t->vertexes[ TTetraFaces[j*3+1] ];
			v2 = t->vertexes[ TTetraFaces[j*3+2] ];
			if (aMesh->memPool == NULL ) 
				tr = new TTriangle( v0,v1,v2);
			else
				tr = aMesh->memPool->getTriangleInstance(v0,v1,v2);
			tr->calcNormal();
			aMesh->addTriangle(tr);
		}
	}
	endProcess((char*)("fastGetSurfaceTriangles : createPart"));
	delete _faces;
  
}
// Simple.. 
void lpEvaluateClusterByFace(int i, int thId  ,TObject* destObject)
{
	TElementsCluster* aCluster =  (TElementsCluster*)destObject;	
	TList<TObject*>* vl ;	
	TVertex *v0, *v1, *v2, *vi;
	TTetra *t0 , *t1;
	// obtengo nuevamente los vecinos del nodo
	vi = aCluster->inspVertex;
	vl = vi->elementsList;
	// Para todos las caras de los elementos q contienen este nodo
	for (int ie = 0; ie< vl->Count() ; ie ++)
	{
		 t0 =(TTetra*)( vl->elementAt(ie));
		 // Si es el nodo con el ID mas chico del elemento, sigo
		 if ( (vi->id <= t0->vertexes[0]->id) &&
			  (vi->id <= t0->vertexes[1]->id) &&
			  (vi->id <= t0->vertexes[2]->id) &&
			  (vi->id <= t0->vertexes[3]->id) )
		 {			 
			 // Para todas las caras de los elementos
			 for (int iface = 0 ; iface<4 ; iface++)
			 {
				v0 = t0->vertexes[TTetraFaces[iface*3]];
				v1 = t0->vertexes[TTetraFaces[iface*3+1]];
				v2 = t0->vertexes[TTetraFaces[iface*3+2]];
				//Si esta cara contiene al nodo, sigo
				if (!((v0 == vi) || (v1 == vi) || (v2 == vi))) continue;
				//Obtengo el Tetra vecino
				t1 = t0->getTetraNeighbour(0,v0,v1,v2,NULL);
				if (t1 == NULL ) continue;
				if (t0->id>t1->id) continue;
				// Evaluo
				aCluster->inspectedElements->Clear();
				aCluster->inspectedElements->Add(t0);
				aCluster->inspectedElements->Add(t1);
				aCluster->generateSubMesh();				
				if (aCluster->evaluateSet()>0) 
				{
					aCluster->updateMesh(true);
				}
			 }
		 }
	}
		 	
}
// Para ese vertice evalua con todas las aristas vecinas
void lpEvaluateClusterByEdge(int i, int thId  ,TObject* destObject)
{
	TElementsCluster* aCluster =  (TElementsCluster*)destObject;
	TList<TObject*>* vl = new TList<TObject*>();
	
	TTetra *t;
	TVertex *v0,*v1;
	int j,k;

	v0 = aCluster->inspVertex;
	// obtengo nuevamente los vecinos de orden 1
	vl->Clear();
	v0->getVertexNeighboursByElem(vl,1);
	
	//  Recorro los vertices vecinos
	for ( j = 0 ; j<vl->Count() ; j++)
	{
			v1 = (TVertex*)(vl->elementAt(j));
			if (v0->id>=v1->id) continue;
			if (v1->elementsList == NULL) continue;
			
			aCluster->inspectedElements->Clear();
			//--- Veo los vecinos de arista de ambos vertices
			//Para el vertice 0
			for (k = 0 ; k<v0->elementsList->Count() ; k++)
			{
				t = (TTetra*)(v0->elementsList->elementAt(k));
				if (t == NULL) continue;
				if (t->isdestroyed) continue;
				if  (!t->hasEdge(v0,v1)) continue;
				if (aCluster->inspectedElements->indexOf(t)>=0) continue;
				// InnerFlag
				//if (t->flag == innerFlag) continue;
				//t->flag = innerFlag;
				aCluster->inspectedElements->Add(t);
			}
			//Para el vertice 1
			for (k = 0 ; k<v1->elementsList->Count() ; k++)
			{
				t = (TTetra*)(v1->elementsList->elementAt(k));
				if (t == NULL) continue;
				if (t->isdestroyed) continue;
				if  (!t->hasEdge(v0,v1)) continue;
				if (aCluster->inspectedElements->indexOf(t)>=0) continue;
				//if (t->flag == innerFlag) continue;
				//t->flag = innerFlag;
				aCluster->inspectedElements->Add(t);
			}
						
			aCluster->generateSubMesh( );
			if (aCluster->evaluateSet()>0) 
			{
				aCluster->updateMesh(true);	
				break;
			}
	}
}


void lpEvaluateClusterByNode(int i, int thId  ,TObject* destObject)
{
	TElementsCluster* aCluster =  (TElementsCluster*)destObject;
	TVertex *v0;
	v0 = aCluster->inspVertex;
	v0->elementsList->Pack();
	aCluster->inspectedElements->Clear();
	aCluster->inspectedElements->Assign(v0->elementsList);
	// obtengo nuevamente los vecinos de orden 1	
	aCluster->generateSubMesh();
    aCluster->evaluateSet();
    if (aCluster->evaluateSet()>0) 
	{
		aCluster->updateMesh(true) ;			
	}
}


void assignVertexesAvoidingVisited(TList<TVertex*> *vs, TList<TObject*> *vRes ,int iter, int maxAssignment)
{
  int i,j, nAssigned;
  TList<TObject*>  *lneigh;
  bool wasVisited;

  lneigh = new TList<TObject*>();
    
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
	 v->getVertexNeighboursByElem(lneigh);    	 
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
  delete lneigh;
}


void ParallelEvaluateClusterByNode(TMesh *aMesh , TVertexesEvaluator fc)
{  
  int	iv ,i ,nsimCh;
	TList<TObject*> *vRes;
	TList<TVertex*> *vertexesCopy;
	TList<TObject*> * resultedClusters; 
	// Tama�o maximo de procesos simultaneos
	nsimCh = 2048;
	//----------------------------------------
	// Initialization part!
	startProcess((char*)("Initialization"));
	resultedClusters = new TList<TObject*>();	
	for (i = 0 ; i<nsimCh ; i++)
	{
		TElementsCluster* e = new TElementsCluster(aMesh,vrelaxQuality) ;
		resultedClusters->Add( (TObject*)(e));
	}

	vertexesCopy = new TList<TVertex*>();
	vRes = new TList<TObject*>();

	vertexesCopy->Assign(aMesh->vertexes);
	double minQ;
	for (i = 0 ;i<vertexesCopy->Count() ; i++)
	{
		TVertex *v = vertexesCopy->elementAt(i);
		v->visited = 0;
		v->flag = 0;
		v->isdestroyed = false;
		minQ =  50000;
		v->elementsList->Pack();
		for (int j = 0; j<v->elementsList->Count();j++)
		{
			TTetra *t = (TTetra*)(v->elementsList->elementAt(j));
			t->calidad = diedralAngle(t->vertexes);
			minQ = Min( minQ , t->calidad);
		}
		v->calidad = minQ ; 
	}
	
	vertexesCopy->Sort(sortByQuality);
	endProcess((char*)("Initialization"));
	//----------------------------------------	
	//----------------------------------------
	////-- Facil de Paralelizar
	startProcess((char*)("evaluateClustersInParallel"));
	for (iv = 0  ; iv<=Max(20.0f,(float)(vertexesCopy->Count() / nsimCh)-1 ) ; iv++) 
	{
		vRes->Clear();
		//Distribuyo la carga por Arista		
		startProcess((char*)("assignVertexesAvoidingVisited"));
				assignVertexesAvoidingVisited(vertexesCopy,vRes,iv+1,nsimCh-1);  		
		endProcess((char*)("assignVertexesAvoidingVisited"));
		if (vRes->Count() == 0 ) break;

		startProcess((char*)("clearVars"));
		//--Limpio las variables
		// por cada vertice, tengo un cluster
		for (i = 0 ; i<vRes->Count() ; i++)
		{
			TElementsCluster* resC = (TElementsCluster*)(resultedClusters->elementAt(i));			
			resC->inspVertex = (TVertex*)(vRes->elementAt(i));
			resC->doRemoveElements = false;
			resC->testCenter = false;
			resC->perturbCenter = 0;
			resC->checkMaxLength = true;        
		}
		endProcess((char*)("clearVars"));

		startProcess((char*)("parallelPart"));
			//pi->forloop(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByEdge);      
		    parallelFor(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByNode);    
		    //for (int vi = 0 ; vi<vRes->Count() ; vi++)
			 //  lpEvaluateClusterByEdge(vi,0,resultedClusters->elementAt(vi));
		endProcess((char*)("parallelPart"));
		/* end of parallel section */

		startProcess((char*)("Generating new elements"));
		for (i = 0 ; i<vRes->Count() ; i++)
		{
			TElementsCluster* ec = (TElementsCluster*)( resultedClusters->elementAt(i));			
			ec->genElements();		  
		}

		if (aMesh->elementsToRemove->Count()>ELEMENTS_TO_FORCE_UPDATE )
				aMesh->updateRefs();
		endProcess((char*)("Generating new elements"));
	}
	endProcess((char*)("evaluateClustersInParallel"));	

	startProcess((char*)("updateRefs"));
	aMesh->updateRefs();
	aMesh->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
	endProcess((char*)("updateRefs"));
}

void ParallelEvaluateClusterByEdge(TMesh *aMesh , TVertexesEvaluator fc)
{  
	int	iv ,i ,nsimCh;
	TList<TObject*> *vRes;
	TList<TVertex*> *vertexesCopy;
	TList<TObject*> * resultedClusters; 
	// Tama�o maximo de procesos simultaneos
	nsimCh = 2048;
	//----------------------------------------
	// Initialization part!
	startProcess((char*)("Initialization"));
	resultedClusters = new TList<TObject*>();	
	for (i = 0 ; i<nsimCh ; i++)
	{
		TElementsCluster* e = new TElementsCluster(aMesh,vrelaxQuality) ;
		resultedClusters->Add( (TObject*)(e));
	}

	vertexesCopy = new TList<TVertex*>();
	vRes = new TList<TObject*>();

	vertexesCopy->Assign(aMesh->vertexes);
	double minQ;
	for (i = 0 ;i<vertexesCopy->Count() ; i++)
	{
		TVertex *v = vertexesCopy->elementAt(i);
		v->visited = 0;
		v->flag = 0;
		v->isdestroyed = false;
		minQ =  50000;
		v->elementsList->Pack();
		for (int j = 0; j<v->elementsList->Count();j++)
		{
			TTetra *t = (TTetra*)(v->elementsList->elementAt(j));
			t->calidad = diedralAngle(t->vertexes);
			minQ = Min( minQ , t->calidad);
		}
		v->calidad = minQ ; 
	}
	
	vertexesCopy->Sort(sortByQuality);
	endProcess((char*)("Initialization"));
	//----------------------------------------	
	//----------------------------------------
	////-- Facil de Paralelizar
	startProcess((char*)("evaluateClustersInParallel"));
	for (iv = 0  ; iv<=Max(20.0f,(float)(vertexesCopy->Count() / nsimCh)-1 ) ; iv++) 
	{
		vRes->Clear();
		//Distribuyo la carga por Arista		
		startProcess((char*)("assignVertexesAvoidingVisited"));
				assignVertexesAvoidingVisited(vertexesCopy,vRes,iv+1,nsimCh-1);  		
		endProcess((char*)("assignVertexesAvoidingVisited"));
		if (vRes->Count() == 0 ) break;

		startProcess((char*)("clearVars"));
		//--Limpio las variables
		// por cada vertice, tengo un cluster
		for (i = 0 ; i<vRes->Count() ; i++)
		{
			TElementsCluster* resC = (TElementsCluster*)(resultedClusters->elementAt(i));			
			resC->inspVertex = (TVertex*)(vRes->elementAt(i));
			resC->doRemoveElements = false;
			resC->testCenter = false;
			resC->perturbCenter = 0;
			resC->checkMaxLength = true;        
		}
		endProcess((char*)("clearVars"));

		startProcess((char*)("parallelPart"));
			//pi->forloop(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByEdge);      
		    parallelFor(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByEdge);    
		    //for (int vi = 0 ; vi<vRes->Count() ; vi++)
			 //  lpEvaluateClusterByEdge(vi,0,resultedClusters->elementAt(vi));
		endProcess((char*)("parallelPart"));
		/* end of parallel section */

		startProcess((char*)("Generating new elements"));
		for (i = 0 ; i<vRes->Count() ; i++)
		{
			TElementsCluster* ec = (TElementsCluster*)( resultedClusters->elementAt(i));			
			ec->genElements();		  
		}

		if (aMesh->elementsToRemove->Count()>ELEMENTS_TO_FORCE_UPDATE )
				aMesh->updateRefs();
		endProcess((char*)("Generating new elements"));
	}
	endProcess((char*)("evaluateClustersInParallel"));	

	startProcess((char*)("updateRefs"));
	aMesh->updateRefs();
	aMesh->updateIndexes( KEEP_ORIG_IDS);
	endProcess((char*)("updateRefs"));

}

void ParallelEvaluateClusterByFace(TMesh *aMesh , TVertexesEvaluator fc)
{  
	int	iv ,i ,nsimCh;
	TList<TObject*> *vRes;
	TList<TVertex*> *vertexesCopy;
	TList<TObject*> * resultedClusters; 
	// Tama�o maximo de procesos simultaneos
	nsimCh = 2048;
	//----------------------------------------
	// Initialization part!
	startProcess((char*)("Initialization"));
	resultedClusters = new TList<TObject*>();	
	for (i = 0 ; i<nsimCh ; i++)
	{
		TElementsCluster* e = new TElementsCluster(aMesh,vrelaxQuality) ;
		resultedClusters->Add( (TObject*)(e));
	}

	vertexesCopy = new TList<TVertex*>();
	vRes = new TList<TObject*>();

	vertexesCopy->Assign(aMesh->vertexes);
	double minQ;
	for (i = 0 ;i<vertexesCopy->Count() ; i++)
	{
		TVertex *v = vertexesCopy->elementAt(i);
		v->visited = 0;
		v->flag = 0;
		v->isdestroyed = false;
		minQ =  50000;
		v->elementsList->Pack();
		for (int j = 0; j<v->elementsList->Count();j++)
		{
			TTetra *t = (TTetra*)(v->elementsList->elementAt(j));
			t->calidad = diedralAngle(t->vertexes);
			minQ = Min( minQ , t->calidad);
		}
		v->calidad = minQ ; 
	}
	
	vertexesCopy->Sort(sortByQuality);
	endProcess((char*)("Initialization"));
	//----------------------------------------	
	//----------------------------------------
	////-- Facil de Paralelizar
	startProcess((char*)("evaluateClustersInParallel"));
	for (iv = 0  ; iv<=Max(20.0f,(float)(vertexesCopy->Count() / nsimCh)-1 ) ; iv++) 
	{
		vRes->Clear();
		//Distribuyo la carga por Arista		
		startProcess((char*)("assignVertexesAvoidingVisited"));
			assignVertexesAvoidingVisited(vertexesCopy,vRes,iv+1,nsimCh-1);  		
		endProcess((char*)("assignVertexesAvoidingVisited"));
		if (vRes->Count() == 0 ) break;

		startProcess((char*)("clearVars"));
		//--Limpio las variables
		// por cada vertice, tengo un cluster
		for (i = 0 ; i<vRes->Count() ; i++)
		{
			TElementsCluster* resC = (TElementsCluster*)(resultedClusters->elementAt(i));			
			resC->inspVertex = (TVertex*)(vRes->elementAt(i));
			resC->doRemoveElements = false;
			resC->testCenter = false;
			resC->perturbCenter = 0;
			resC->checkMaxLength = true;        
		}
		endProcess((char*)("clearVars"));

		startProcess((char*)("parallelPart"));
			//pi->forloop(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByEdge);      
		    parallelFor(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByFace);    
		    //for (int vi = 0 ; vi<vRes->Count() ; vi++)
			 //  lpEvaluateClusterByFace(vi,0,resultedClusters->elementAt(vi));
		endProcess((char*)("parallelPart"));
		/* end of parallel section */

		startProcess((char*)("Generating new elements"));
		for (i = 0 ; i<vRes->Count() ; i++)
		{
			TElementsCluster* ec = (TElementsCluster*)( resultedClusters->elementAt(i));			
			ec->genElements();		  
		}

		if (aMesh->elementsToRemove->Count()>ELEMENTS_TO_FORCE_UPDATE )
				aMesh->updateRefs();
		endProcess((char*)("Generating new elements"));
	}
	endProcess((char*)("evaluateClustersInParallel"));	

	startProcess((char*)("updateRefs"));
	aMesh->updateRefs();
		aMesh->updateIndexes( KEEP_ORIG_IDS);
	endProcess((char*)("updateRefs"));

}