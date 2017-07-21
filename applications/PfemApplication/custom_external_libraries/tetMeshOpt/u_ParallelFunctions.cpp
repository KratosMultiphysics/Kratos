#include "u_ParallelFunctions.h"
#include "u_elementCluster.h"
#include "u_ProcessTime.h"
#include "u_qualityMetrics.h"
//#include "cl_VolumeMeshSmooth.h"

#include <cstdio>
#if !defined(KRATOS)	
#include <tbb\mutex.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// #endif

#if !defined(KRATOS)	
#include "cl_utils.h"
#endif

bool* _faces;
int innerID = 0;
int innerNumThreads =12;
TList<TObject*> * resultedClusters = NULL;  	

void setNumThreads(int nt)
{
	innerNumThreads = nt;
#ifdef _OPENMP
  omp_set_num_threads(nt);
#endif 

}

double minExpectedAngle;
#if !defined(KRATOS)	
tbb::mutex innerMutex;
#endif

void enterCritical()
{
#if !defined(KRATOS)	
	innerMutex.lock();
#else   
#pragma omp barrier
#endif
}

void leaveCritical()
{
#if !defined(KRATOS)	
	innerMutex.unlock();
#else   
#pragma omp barrier
#endif
}

int incrementalID()
{
#if !defined(KRATOS)	
	innerMutex.lock();
	innerID++;
	innerMutex.unlock();
#else   
#pragma omp barrier
	innerID++;
#endif

	return innerID;
}

void parallelFor(int from, int to,  TList<TObject*>* elements,TStForLoopElement call)
{
#if !defined(KRATOS)
	if (innerNumThreads == 1)
	{
		for (int i=from ;i<=to ; i++)
			call(i,0,elements->structure[i]);
	}
	else
		parallel_for(blocked_range<size_t>(from, to), TParallelIterator(call,elements ) );
#else
#pragma omp parallel for
	for (int i=from ; i<=to ; i++)
		call(i,0,elements->elementAt(i));

#pragma omp barrier
#endif
}

void parallelFor(int from, int to,  TList<TVertex*>* elements,TStForLoopElement call)
{
	TList<TObject*>* lo = new TList<TObject*>();
	for (int i=0; i<elements->Count();i++)
		lo->Add(elements->structure[i]);

#if !defined(KRATOS)
	if (innerNumThreads == 1)
		for (int i=from ;i<=to ; i++)
			call(i,0,lo->structure[i]);
	else
		parallel_for(blocked_range<size_t>(from, to), TParallelIterator(call,lo ) );
#else
#pragma omp parallel for
	for (int i=from ; i<=to ; i++)
		call(i,0,elements->elementAt(i));

#pragma omp barrier
#endif

	delete lo;
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
	#if !defined(KRATOS)	
	startProcess((char*)("fastGetSurfaceTriangles : delete old triangles"));
	#endif
	for (int i=0;i<aMesh->fFaces->Count();i++)
	{
		TTriangle *tr = (TTriangle*)(aMesh->fFaces->structure[i]);
		// Si tengo el pool activado, lo guardo ahi 
		if (aMesh->memPool == NULL ) 
			delete tr;
		else
			aMesh->memPool->releaseInstance(tr);
	}
	#if !defined(KRATOS)	
	endProcess((char*)("fastGetSurfaceTriangles : delete old triangles"));
	#endif

	aMesh->fFaces->Clear();
	#if !defined(KRATOS)	
	startProcess((char*)("fastGetSurfaceTriangles : parallelPart"));
	#endif
	//parallel_for(blocked_range<size_t>(0,numT), TParallelIterator(localProcessI,aMesh->elements) );    
#pragma omp parallel for
	for (i = 0; i<numT ; i++)
		localProcessI(i,0,aMesh->elements->elementAt(i));
	#if !defined(KRATOS)	
	endProcess((char*)("fastGetSurfaceTriangles : parallelPart"));

	startProcess((char*)("fastGetSurfaceTriangles : createPart"));
	#endif
	
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
	#if !defined(KRATOS)	
	endProcess((char*)("fastGetSurfaceTriangles : createPart"));
	#endif
	
	delete _faces;

}
// Simple.. 
void lpEvaluateClusterByFace(int i, int thId  ,TObject* destObject)
{
	TElementsCluster* aCluster =  (TElementsCluster*)destObject;	
	TList<TObject*>* vl ;	
	TVertex *v0, *v1, *v2, *vi;
	TTetra *t0 , *t1;
	bool wasCreated = false;
	// obtengo nuevamente los vecinos del nodo
	vi = aCluster->inspVertex;
	vl = vi->elementsList;
	
	// Para todos las caras de los elementos q contienen este nodo
	for (int ie = 0; (ie< vl->Count()) && (!wasCreated) ; ie ++)
	{
		t0 =(TTetra*)( vl->elementAt(ie));
		if (t0->fMinDiedralAngle>minExpectedAngle) continue;
		for (int iface = 0 ; (iface<4) && (!wasCreated); iface++)
		{
			// Para todas las caras de los elementos
			v0 = t0->vertexes[TTetraFaces[iface*3]];
			v1 = t0->vertexes[TTetraFaces[iface*3+1]];
			v2 = t0->vertexes[TTetraFaces[iface*3+2]];
			//Si esta cara contiene al nodo, sigo
			if (!((v0 == vi) || (v1 == vi) || (v2 == vi))) continue;
			// Si es el nodo con el ID mas chico del elemento, sigo
			if ( (vi->getID()< v0->getID()) &&
				 (vi->getID() < v1->getID()) &&
				 (vi->getID() < v2->getID())  )
			{	
				//Obtengo el Tetra vecino
				t1 = t0->getTetraNeighbour(0,v0,v1,v2,NULL);
				if (t1 == NULL ) continue;
				if (t1->fMinDiedralAngle>minExpectedAngle) continue;
				if (t0->getID()>t1->getID()) continue;
				// Evaluo
				aCluster->inspectedElements->Clear();
				aCluster->inspectedElements->Add(t0);
				aCluster->inspectedElements->Add(t1);
				//if (aCluster->getMinAngle() > minExpectedAngle) continue;
				aCluster->generateSubMesh();				

				if (aCluster->evaluateSet()>0) 
				{
					aCluster->updateMesh(true);
					wasCreated = true;
					break;
				}
			}
		}
	}

}

void getElementsByEdgeFaster(TVertex *v0, TVertex*v1 , TList<TObject*>* inspectedElements)
{
	TTetra *t, *tn, *ti;
	int k;
	TVertex *nv0,*nv1,*nv2,*fv0,*fv1,*fv2;
	//--- Obtengo el primer elemento que tenga esta arista
	ti = NULL;
	t = NULL;

	//Para el vertice 0
	for (k = 0 ; k<v0->elementsList->Count() ; k++)
	{
		t = (TTetra*)(v0->elementsList->elementAt(k));
		if (t == NULL) continue;
		if (t->isdestroyed) continue;
		if  (!t->hasEdge(v0,v1)) continue;	
		ti = t;
		break;				
	}

	inspectedElements->Add(ti);
	// Itero en las 2 direcciones
	for (int iface = 0; iface<3 ; iface++)
	{
		fv0 = ti->vertexes[TTetraFaces[iface*3]];
		fv1 = ti->vertexes[TTetraFaces[iface*3+1]];
		fv2 = ti->vertexes[TTetraFaces[iface*3+2]];
		// Contiene la arista
		if ( ((fv0 == v0 )&& (fv1 ==v1)) || ((fv1 == v0 )&& (fv2 ==v1)) || ((fv2 == v0 )&& (fv0 ==v1)) )
		{
			// Guardo la cara de inicio
			nv0 = fv0 ; nv1 = fv1; nv2 = fv2;			
		}
		else
			continue;

		// Recorro los vecinos por cara
		while (true)
		{				
			tn = t->getNeighByFace(nv0,nv1,nv2);
			//Di la vuelta
			if (tn == ti) {iface = iface+5; break;}
			//Es una cara de superficie
			if (tn == NULL) break;
			// Sigo iterando
			t = tn;
			inspectedElements->Add(t);
			// Busco el vecino
			for (int k=0; k<4;k++)
			{
				fv0 = t->vertexes[TTetraFaces[k*3]];
				fv1 = t->vertexes[TTetraFaces[k*3+1]];
				fv2 = t->vertexes[TTetraFaces[k*3+2]];
				// Contiene la arista
				if ( ((fv0 == v0 )&& (fv1 ==v1)) || ((fv1 == v0 )&& (fv2 ==v1)) || ((fv2 == v0 )&& (fv0 ==v1)) )
				{
					//Por esta cara ya vine
					if( ((fv0 == nv0) || (fv0 == nv1) || (fv0 == nv2) ) &&
						((fv1 == nv0) || (fv1 == nv1) || (fv1 == nv2) ) &&
						((fv2 == nv0) || (fv2 == nv1) || (fv2 == nv2) ) )
						continue;
					//es la cara opuesta
					else
					{
						nv0 = fv0 ;
						nv1 = fv1 ; 
						nv2 = fv2;
						break;
					}
				}
			}
		}
	}
}

double getElementsByEdge(TVertex *v0, TVertex*v1 , TList<TObject*>* inspectedElements)
{
	TTetra *t;
	int k;
	double minQ = 500000;
	//--- Veo los vecinos de arista de ambos vertices
	//Para el vertice 0
	for (k = 0 ; k<v0->elementsList->Count() ; k++)
	{
		t = (TTetra*)(v0->elementsList->elementAt(k));
		if (t == NULL) continue;
		if (t->isdestroyed) continue;
		if  (!t->hasEdge(v0,v1)) continue;				
		minQ = Min( t->fMinDiedralAngle , minQ);
		if (inspectedElements->contains(t)) continue;
		inspectedElements->Add(t);
	}	
	return minQ;
}
// Para ese vertice evalua con todas las aristas vecinas
void lpEvaluateClusterByEdge(int i, int thId  ,TObject* destObject)
{
	TElementsCluster* aCluster =  (TElementsCluster*)destObject;
	TList<TVertex*>* vl = new TList<TVertex*>();
	TVertex *v0,*v1;
	int j;
	
	v0 = aCluster->inspVertex;
	// obtengo nuevamente los vecinos de orden 1	
	v0->getVertexNeighboursByElem(vl,1,true);
	//vl->Sort(sortByID);
	//  Recorro los vertices vecinos
	for ( j = 0 ; j<vl->Count() ; j++)
	{
		v1 = (TVertex*)(vl->elementAt(j));
		//if (v0->id>=v1->id) continue;
		if (v1->elementsList == NULL) continue;
		//int innerFlag = rand();

		aCluster->inspectedElements->Clear();			
		double mq = getElementsByEdge(v0,v1 , aCluster->inspectedElements);

		if (mq > minExpectedAngle) continue;
		aCluster->generateSubMesh( );			

		if (aCluster->evaluateSet()>0) 
		{
			aCluster->updateMesh(true);	
			break;
		}
	}

	delete vl;
	
}

double _gbminExpectedQ ;
int _gbmxIter;
int _gbsubI;


void lpMeshSmooth(int i, int thId  ,TObject* destObject)
{
  int j, nproposals,iter;
  TList<TObject*> *vList;
  TTetra *t;
  float4 initialPos ,proposed;
  bool hasNegative;
  double minq,avgq ,minq2, avgq2 ,q,minVol,radius;
       TElementsCluster* aCluster =  (TElementsCluster*)destObject;
	   TVertex *v = aCluster->inspVertex;
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
		 if (minq > _gbminExpectedQ) return;

	   	//pruebo variando en un entorno
		initialPos =v->fPos;
		double randX, randY, randZ;
		for (iter = 0 ; iter<_gbmxIter ; iter++)
		{
			for (nproposals = 0 ;nproposals<_gbsubI ;nproposals++)
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
                               if (minq2>minq)
							   {
								   initialPos = v->fPos;
									minq = minq2;
							   }

			}
		}
		v->fPos = initialPos;
}


void lpEvaluateClusterByNode(int i, int thId  ,TObject* destObject)
{
	TElementsCluster* aCluster =  (TElementsCluster*)destObject;
	TVertex *v0;
	v0 = aCluster->inspVertex;
	v0->elementsList->Pack();
	
	if (v0->elementsList->Count() == 0) return;
	aCluster->inspectedElements->Clear();
	aCluster->inspectedElements->Assign(v0->elementsList);
	// obtengo nuevamente los vecinos de orden 1	
	
	if (aCluster->getMinAngle() > minExpectedAngle) return;
	aCluster->generateSubMesh();

	aCluster->evaluateSet();
	if (aCluster->evaluateSet()>0) 
	{
		aCluster->updateMesh(true) ;			
	}
	
}

bool validateAssignment(TList<TVertex*> *vs,TList<TObject*> *vRes)
{
	TList<TVertex*>  *lneigh;
	lneigh = new TList<TVertex*>();
	for (int i = 0 ; i<vs->Count() ; i++)
	{
		TVertex *v = (TVertex*)(vs->elementAt(i));
		v->innerFlag = 0;
	}


	for (int i=0; i<vRes->Count(); i++)
	{
		TVertex *v = (TVertex*)(vRes->elementAt(i));

		lneigh->Clear();
		v->getVertexNeighboursByElem(lneigh);    	 
		for (int j=0; j<lneigh->Count() ;j++)
		{
			TVertex *v2 = (TVertex*)(lneigh->elementAt(j));
			if (v2->innerFlag == 1) return false;
			v2->innerFlag = 1;
		}
	}
	return true;
}

void assignVertexesAvoidingVisited(TList<TVertex*> *vs, TList<TObject*> *vRes ,int iter, int maxAssignment, double minExpQ)
{
	int i,j, nAssigned;
	TList<TVertex*>  *lneigh;
	bool wasVisited;

	lneigh = new TList<TVertex*>();

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
		if (v->calidad>minExpQ) break;

		lneigh->Clear();
		v->getVertexNeighboursByElem(lneigh);    	 
		//lneigh = v->neighV;
		if (lneigh == NULL) continue;
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

void lpComputeQuality(int i, int thId  ,TObject* destObject)
{
	TTetra *t = (TTetra*)(destObject);
	if (t == NULL) return;
	t->calidad = vrelaxQuality(t->vertexes);
	t->fMinDiedralAngle = diedralAngle(t->vertexes);
}

void lpCalculateVertexVertexNeigh(int i, int thId  ,TObject* destObject)
{
	TVertex *v = (TVertex*)(destObject);
	if (v == NULL) return ;
	if ( v->neighV == NULL) 
		v->neighV = new TList<TVertex*>();
	v->neighV->Clear();
	v->getVertexNeighboursByElem(NULL,1);
}

double ParallelEvaluateCluster(TMesh *aMesh , TVertexesEvaluator fc, int mode, bool sort)
{
	int	iv ,i ,nsimCh;
	TList<TObject*> *vRes;
	TList<TVertex*> *vertexesCopy;
    
	// Tama�o maximo de procesos simultaneos
	nsimCh = ASSIGNMENT_SIZE;
	//----------------------------------------
	// Initialization part!
	startProcess((char*)("Initialization"));

	// Initialize variables
	vRes = new TList<TObject*>();
	if (resultedClusters == NULL)
	{
		std :: cout <<"*********** Creating clusters**************** \n";
		resultedClusters = new TList<TObject*>();	
		for (i = 0 ; i<ASSIGNMENT_SIZE ; i++)
		{
			TElementsCluster* e = new TElementsCluster(aMesh,vrelaxQuality) ;
			resultedClusters->Add( (TObject*)(e));
		}
	}
	else
	{
		for (i = 0 ; i<ASSIGNMENT_SIZE ; i++)
		{
			TElementsCluster* e = (TElementsCluster*)(  resultedClusters->elementAt(i)); 
			e->originalMesh = aMesh;
		}
	}
	

	vertexesCopy = aMesh->vertexes;

	int numClusters, numEvaluatedClusters;
	numClusters = numEvaluatedClusters = 0;

	if (sort)
		{
			parallelFor(0,aMesh->elements->Count()-1,aMesh->elements,lpComputeQuality);
			double minQ;

			for (i = 0 ;i<(int)(vertexesCopy->Count()) ; i++)
			{
				TVertex *v = vertexesCopy->elementAt(i);
				v->visited = 0;
				v->flag = 0;
				v->isdestroyed = false;
				minQ =  50000;
				v->elementsList->Pack();			
				for (int j = 0; j<(int)(v->elementsList->Count());j++)
				{
					TTetra *t = (TTetra*)(v->elementsList->elementAt(j));				
					minQ = Min( minQ , t->fMinDiedralAngle);
				}
				v->calidad = minQ ; 
			}

			vertexesCopy->Sort(sortByQuality);
		}
		else
		{
			for (i = 0 ;i<(int)(vertexesCopy->Count()) ; i++)
			{
				TVertex *v = vertexesCopy->structure[i];
				v->visited = 0;
				v->flag = 0;
				v->isdestroyed = false;			
				v->elementsList->Pack();	
			}
		}


	endProcess((char*)("Initialization"));
	//----------------------------------------	
	//----------------------------------------
	////-- Facil de Paralelizar
	startProcess((char*)("evaluateClustersInParallel"));
	for (iv = 0  ; iv<=Max(20.0f,(float)(vertexesCopy->Count() / nsimCh)-1 ) ; iv++) 
	{
		vRes->Clear();		
		//parallelFor(0,aMesh->vertexes->Count()-1,aMesh->vertexes,lpCalculateVertexVertexNeigh);

		//Distribuyo la carga por Arista		
		startProcess((char*)("assignVertexesAvoidingVisited"));
		assignVertexesAvoidingVisited(vertexesCopy,vRes,iv+1,nsimCh-1 ,minExpectedAngle );  		
		endProcess((char*)("assignVertexesAvoidingVisited"));
		if (vRes->Count() == 0 ) break;

		/*if (!validateAssignment(vertexesCopy, vRes))			
		{
		std::cout << "Invalid Assignment"<< "\n";
		break;
		}
		*/
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
			numClusters ++ ;
		}
		endProcess((char*)("clearVars"));

		startProcess((char*)("parallelPart"));
		if (mode ==0) 
			parallelFor(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByNode);    
		else
			if (mode ==1) 
				parallelFor(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByEdge);    
			else
				if (mode ==2) 
				parallelFor(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByFace);    
		        if (mode ==3)
					parallelFor(0,vRes->Count()-1,resultedClusters, lpMeshSmooth);    
                #if !defined(KRATOS)	
				else
				if (mode ==4)					
					clBestVertexPos(vRes,_gbminExpectedQ,_gbmxIter,_gbsubI);
				#endif
		//pi->forloop(0,vRes->Count()-1,resultedClusters, lpEvaluateClusterByEdge);      
		endProcess((char*)("parallelPart"));
		/* end of parallel section */

		startProcess((char*)("Generating new elements"));
		for (i = 0 ; i<vRes->Count() ; i++)
		{
			TElementsCluster* ec = (TElementsCluster*)( resultedClusters->elementAt(i));			
			if (ec->newElements->Count()>0)
				 numEvaluatedClusters++; 
			ec->genElements();	
			
		}

		if (aMesh->elementsToRemove->Count()>ELEMENTS_TO_FORCE_UPDATE )
			aMesh->updateRefs();
		endProcess((char*)("Generating new elements"));
	}
	endProcess((char*)("evaluateClustersInParallel"));	

	startProcess((char*)("updateRefs"));
	aMesh->updateRefs();	
	endProcess((char*)("updateRefs"));

	// Clear variables;
	delete vRes;
	/*
	for (int i = 0; i<resultedClusters->Count(); i++)
	{
		TElementsCluster* resC = (TElementsCluster*)(resultedClusters->elementAt(i));					
	}

	for (int i = 0; i<resultedClusters->Count(); i++)
	{
		TElementsCluster* resC = (TElementsCluster*)(resultedClusters->elementAt(i));			
		resultedClusters->setElementAt(i, NULL);
		delete resC ;
	}
	
	delete resultedClusters ;	
	*/
	return numEvaluatedClusters;
}

void preparePool()
{
		resultedClusters = new TList<TObject*>();	
		for (int i = 0 ; i<ASSIGNMENT_SIZE ; i++)
		{
			TElementsCluster* e = new TElementsCluster(NULL,vrelaxQuality) ;
			resultedClusters->Add( (TObject*)(e));
		}

}
void clearPool()
{
	for (int i = 0; i<ASSIGNMENT_SIZE; i++)
	{
		TElementsCluster* resC = (TElementsCluster*)(resultedClusters->elementAt(i));			
		resultedClusters->setElementAt(i, NULL);
		delete resC ;
	}
	
	delete resultedClusters ;	

	resultedClusters = NULL;
	
}

double ParallelEvaluateClusterByNode(TMesh *aMesh , TVertexesEvaluator fc, double minExpAngle )
{  
	minExpectedAngle = minExpAngle;
	return ParallelEvaluateCluster(aMesh,fc,0,true);
	//aMesh->updateIndexes( GENERATE_SURFACE | KEEP_ORIG_IDS);
}

double ParallelEvaluateClusterByEdge(TMesh *aMesh , TVertexesEvaluator fc, double minExpAngle )
{  
	minExpectedAngle = minExpAngle;
	return ParallelEvaluateCluster(aMesh,fc,1,true);
	//aMesh->updateIndexes( GENERATE_SURFACE | KEEP_ORIG_IDS);
}

double ParallelEvaluateClusterByFace(TMesh *aMesh , TVertexesEvaluator fc, double minExpAngle )
{  
	minExpectedAngle = minExpAngle;
	return ParallelEvaluateCluster(aMesh,fc,2,true);
	//aMesh->updateIndexes( GENERATE_SURFACE | KEEP_ORIG_IDS);	
}

double ParallelSmoothMesh(TMesh *aMesh , TVertexesEvaluator fc, double minExpAngle )
{  
	minExpectedAngle = minExpAngle;
	_gbminExpectedQ = minExpAngle;
	_gbmxIter = 10;
	_gbsubI = 10;	
	int i;

	for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 1;
	   tr->vertexes[1]->fixed = 1;
	   tr->vertexes[2]->fixed = 1;
   }

	double d = ParallelEvaluateCluster(aMesh,fc,3,true);

	for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 0;
	   tr->vertexes[1]->fixed = 0;
	   tr->vertexes[2]->fixed = 0;
   }
    return d;
	//aMesh->updateIndexes( GENERATE_SURFACE | KEEP_ORIG_IDS);	
}

#if !defined(KRATOS)	
double ParallelSmoothMeshInGPU(TMesh *aMesh , TVertexesEvaluator fc, double minExpAngle )
{  
	minExpectedAngle = minExpAngle;
	_gbminExpectedQ = minExpAngle;
	_gbmxIter = 10;
	_gbsubI = 10;	
	int i;

	for (i = 0 ; i<aMesh->vertexes->Count() ; i++)
   {
	   TVertex *v = (TVertex*)( aMesh->vertexes->elementAt(i));
	   v->fixed = 0;
   }

	for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 1;
	   tr->vertexes[1]->fixed = 1;
	   tr->vertexes[2]->fixed = 1;
   }
	clPrepareMesh(aMesh);

	double d = ParallelEvaluateCluster(aMesh,fc,4,true);

	clReleaseMesh(aMesh);

	for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
   {
	   TTriangle *tr = (TTriangle*)( aMesh->fFaces->elementAt(i));
	   tr->vertexes[0]->fixed = 0;
	   tr->vertexes[1]->fixed = 0;
	   tr->vertexes[2]->fixed = 0;
   }    

	return d;
	//aMesh->updateIndexes( GENERATE_SURFACE | KEEP_ORIG_IDS);	
}
#endif
