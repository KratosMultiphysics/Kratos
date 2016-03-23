#include "u_Types.h"
#include "u_elementCluster.h"
#include "u_qualityMetrics.h"
#include "u_tools.h"
#include "u_ProcessTime.h"
#include "u_TetraFunctions.h"


//-----------------------------------------------------------------------------
///  calcula los triangulos que forman la superficie
/// con la regla que si una cara pertenece a un unico elemento es un triangulo superficial
void dgetSurfaceTriangles(TList<TObject*>* orig , TList<TObject*>* res , TList<TObject*>* surfVertexes)
{
	int i,j  ;
	TTetra *t;
	TList<TObject*>* tempL , *tr;

	tempL = new TList<TObject*>();

	for (i = 0 ; i<orig->Count() ; i++)
	{
		t = (TTetra*)(orig->elementAt(i));
		if (t == NULL) continue;

		tr = t->getSurfaceTriangle(false, tempL,NULL);
		if (tr!=NULL)
		{
			for (j = 0 ;j<tr->Count() ; j++)
			{
				res->Add(tr->elementAt(j));
				if (surfVertexes->indexOf(tr->elementAt(j))<0)
					surfVertexes->Add(tr->elementAt(j));
			}
			tempL->Clear();
		}
	}
	delete tempL;

}


double evaluteVertexChange(TVertex * vOrig , TVertex* vNew ,
	TList<TObject*> *velList, double origVolume , double newVolume)
{
	int i,cnt;
	double result;

	// Reemplazo el viejo vertice
	origVolume = 0;
	newVolume = 0;

	for (i = 0 ; i<velList->Count()/ 4 ; i++)
	{
		TVertex* _v0 = (TVertex*)(velList->elementAt(4*i));
		TVertex* _v1 = (TVertex*)(velList->elementAt(4*i+1));
		TVertex* _v2 = (TVertex*)(velList->elementAt(4*i+2));
		TVertex* _v3 = (TVertex*)(velList->elementAt(4*i+3));     
		origVolume =origVolume+ tetraVolume( _v0->fPos,_v1->fPos,_v2->fPos,_v3->fPos);
	}

	for (i = 0 ; i<velList->Count() ; i++)
	{
		if (velList->elementAt(i) == vOrig) 
			velList->setElementAt(i, vNew);
	}
	//---------------------------------
	for (i = 0 ; i<velList->Count()/ 4 ; i++)
	{
		cnt = 0;
		TVertex* _v0 = (TVertex*)(velList->elementAt(4*i));
		TVertex* _v1 = (TVertex*)(velList->elementAt(4*i+1));
		TVertex* _v2 = (TVertex*)(velList->elementAt(4*i+2));
		TVertex* _v3 = (TVertex*)(velList->elementAt(4*i+3));     
		newVolume =newVolume+ tetraVolume( _v0->fPos,_v1->fPos,_v2->fPos,_v3->fPos);
		// Valido el elemento
		if (_v0 == vNew) cnt++;
		if (_v1 == vNew) cnt++;
		if (_v2 == vNew) cnt++;
		if (_v3 == vNew) cnt++;
		// Si es invalido, lo borro
		if (cnt>1)
		{
			velList->setElementAt(4*i, NULL);
			velList->setElementAt(4*i+1, NULL);
			velList->setElementAt(4*i+2, NULL);
			velList->setElementAt(4*i+3, NULL);
		}
	}
	// elimino invalidos
	velList->Pack();
	// Calculo la minima calidad
	if (velList->Count() == 0 )
		result = -1000000;
	else
	{
		result = 1000;
		for (i = 0 ; i<velList->Count()/ 4 ; i++)
		{
			cnt = 0;
			TVertex* _v0 = (TVertex*)(velList->elementAt(4*i));
			TVertex* _v1 = (TVertex*)(velList->elementAt(4*i+1));
			TVertex* _v2 = (TVertex*)(velList->elementAt(4*i+2));
			TVertex* _v3 = (TVertex*)(velList->elementAt(4*i+3));
			result = Min ( result, vrelaxQuality(_v0,_v1,_v2,_v3) );
		}
	}
	return result;
}


void innercollapseEdges(TVolumeMesh*  aMesh , int  maxSteps,
	TList<TObject*>*  vertexList , double minQuality , int removeSurface ) 
{
	TList<TObject*> *elementsList,   *elemVertexes;
	TList<TVertex*> *vNeigh;
	int i,j,k,nneg;
	TVertex *vOrig, *vNew;
	TTetra *t ;
	double oV = 0 ,nV = 0,nQ = 0,oQ = 0;
	TElementsCluster *aCluster  ;

	//Elijo los malos elementos
	aMesh->selectedElements->Clear();
	//--------------------------------------------------------
	elementsList = new TList<TObject*>();

	for (i=0; i<vertexList->Count() ; i++)
	{
		vOrig = (TVertex*)(vertexList->elementAt(i));
		vOrig->flag = 0;
		vOrig->fmetrica = 500000;
		nneg = 1;
		for (j = 0 ; j<vOrig->elementsList->Count() ; j++)
		{
			t =  (TTetra*)( vOrig->elementsList->elementAt(j));
			vOrig->fmetrica = Min(vOrig->fmetrica, vrelaxQuality(t->vertexes));
			if (t->fVolume<0) nneg++;
		}
		vOrig->fmetrica  = vOrig->fmetrica * nneg;
		vOrig->isSurface = false;
	}

	for (i = 0 ; i<aMesh->fFaces->Count() ; i++)
	{
		TTriangle *tr = (TTriangle*)(aMesh->fFaces->elementAt(i));
		tr->vertexes[0]->isSurface = true;
		tr->vertexes[1]->isSurface = true;
		tr->vertexes[2]->isSurface = true;
	}

	//  vertexList.sort(compareMetrics);
	aCluster = new TElementsCluster(aMesh,vrelaxQuality);
	vNeigh = new TList<TVertex*>();
	elemVertexes = new TList<TObject*>();
	for (i=0 ; i<vertexList->Count() ; i++)
	{
		vOrig = (TVertex*)( vertexList->elementAt(i));
		if (vOrig->flag >0) continue;
		// 1 = Remuevo los vertices de superficie
		if ((removeSurface==1) && (! vOrig->isSurface)) continue;
		// 2 = Remuevo los vertices internos
		if ((removeSurface==2) && (vOrig->isSurface)) continue;
		//Verifico que no lo haya visitado
		vNeigh->Clear();
		vOrig->getVertexNeighboursByElem(vNeigh, 1);
		for (j = 0 ; j<vNeigh->Count() ; j++)
		{
			vNew = (TVertex*)(vNeigh->elementAt(j));
			// 1 = Remuevo los vertices de superficie
			if ((removeSurface==1) && (! vNew->isSurface) ) continue;
			// 2 = Remuevo los vertices internos
			if ((removeSurface==2) && (vNew->isSurface)) continue;

			if (vNew == vOrig) continue;
			if (vNew->flag >0) continue;
			if (vNew->isdestroyed) continue;
			if (vNew->getID()<vOrig->getID())continue;
			if (vNew->elementsList == NULL)  continue;
			// Reemplazo un vertice por otro que tambien cumpla las mismas propiedades
			if (vNew->isSurface && !vOrig->isSurface) continue;
			if (!vNew->isSurface && vOrig->isSurface) continue;

			elementsList->Clear();
			elementsList->Assign(vOrig->elementsList);
			elementsList->Assign(vNew->elementsList,laOr);
			elemVertexes->Clear();
			for (k = 0 ; k<elementsList->Count() ; k++)
			{
				t = (TTetra*)( elementsList->elementAt(k));
				elemVertexes->Add(t->vertexes[0]);
				elemVertexes->Add(t->vertexes[1]);
				elemVertexes->Add(t->vertexes[2]);
				elemVertexes->Add(t->vertexes[3]);
			}
			aCluster->inspectedElements->Assign(elementsList);
			aCluster->copyEL->Assign(elementsList);
			oQ = aCluster->getMinQuality();
			if (oQ>minQuality) continue;

			//2.Evaluo si el resultado es mejor
			nQ = evaluteVertexChange(vOrig,vNew,elemVertexes,oV,nV);
			if (oQ>=nQ)  continue;
			if ((oV>0) && (nV<0) ) continue;
			if (Min (fabs(double(nV)/oV) ,fabs(double(oV)/nV))<0.95) continue;

			aCluster->goodTetraList->Assign(elemVertexes);
			//3.Actualizo la malla
			if (aCluster->updateMesh(true))
			{
				aCluster->genElements(false);
				aMesh->vertexes->Extract(vOrig);
				vOrig->isdestroyed = true;
				vOrig->flag = 1;                                        
				vNew->flag = 1;
				break;
			}
		}
	}
	aMesh->updateRefs();

	aMesh->updateIndexes(0);
}


TElementsCluster::TElementsCluster(TMesh* aMesh, TVertexesEvaluator functor  )
{
	originalMesh = aMesh;
	elements = new TList<TObject*>();
	surfaceT = new TList<TVertex*>();
	copyEL = new TList<TObject*>();
	vertexes = new TList<TVertex*>();
	newElements=new TList<TObject*>();   
	tempTetraList = new TList<TObject*>();
	goodTetraList = new TList<TObject*>();   
	elements2=new TList<TObject*>();
	inspectedElements = new TList<TObject*>();
	tempL = new TList<TObject*>();

	fc = functor;
	vC =NULL;

	//avgEdgeLength = aMesh.AvgEdgeLength;
};

double TElementsCluster::getMinQuality()
{
	int i;
	TTetra* t;
	minQuality = 5000000;

	for (i = 0 ;i< inspectedElements->Count() ; i++)
	{
		t = (TTetra*)( inspectedElements->elementAt(i));		
		t->calidad = fc(t->vertexes);
		minQuality = Min(minQuality, t->calidad );
	}
	return minQuality;
}

double TElementsCluster::getMinAngle()
{
	int i;
	TTetra* t;
	minDAngle = 5000000;

	for (i = 0 ;i< inspectedElements->Count() ; i++)
	{
		t = (TTetra*)( inspectedElements->elementAt(i));		
		t->fMinDiedralAngle = diedralAngle(t->vertexes);
		minDAngle = Min(minDAngle,t->fMinDiedralAngle);		
	}
	return minDAngle;
}

TElementsCluster::~TElementsCluster()
{
	fc =NULL;
	vC =NULL;
	delete(  goodTetraList);
	delete(  surfaceT);
	delete(  vertexes);
	delete(  copyEL);
	delete(  elements2);
	delete(  newElements);
	delete(  tempL );
    delete (elements );
	delete (tempTetraList);
	delete (inspectedElements);
}



int TElementsCluster::evaluateSet() 
{  
	#if !defined(KRATOS)	
		startProcess((char*)("evaluateSet"));
	#endif
	double	tempQuality,actualQuality ;
	int i,j;
	TVertex *v0,*v1,*v2,*v3;
	float4 cn ;

	vC = NULL;

	if (testCenter) 
	{
		cn = 0;
		for (i = 0 ; i<vertexes->Count() ; i++)
		{
			TVertex *_v = (TVertex*)(vertexes->elementAt(i));		   
			cn = cn+  _v->fPos;
		}
		cn = cn * (1.0f/vertexes->Count());
		cn = cn ;
		for (i = 0 ; i<=0; i++)
		{
			vC = new TVertex(cn);
			vC->elementsList = new TList<TObject*>();
			vertexes->Add(vC);
		};
	};
	//--- Evaluate Set
	goodMinQuality =minQuality;
	vertexes->Sort(sortByID);
	
	//---------------------------------------------------
	// Para todos los vertices
	
	bool isBetter = false;
	for (i = 0 ; i<(int)(vertexes->structure.size()) ; i++)
	{
		actualQuality = 50000;
		TVertex *_v = vertexes->structure[i];
		bool clusterWasImproved = true;

		// Para todas las caras
		for (j = 0 ; j<(int)(surfaceT->structure.size() /3) ; j++)
		{
			v0 = surfaceT->structure[3*j];
			v1 = surfaceT->structure[3*j+1];
			v2 = surfaceT->structure[3*j+2];
			v3 =  _v;
			if ((v3 == v1) || (v3 == v2) || (v3 == v0)) continue;
			tempQuality =vrelaxQuality(v0,v1,v2,v3);
			if (tempQuality<goodMinQuality)
			//if (!improvedTetra(v0,v1,v2,v3,this,false) )
			{
				clusterWasImproved = false;
				break;
			}
			
			// Guardo los vertices
			tempTetraList->Add(v0);
			tempTetraList->Add(v1);
			tempTetraList->Add(v2);
			tempTetraList->Add(v3);
			
			actualQuality = Min(actualQuality,tempQuality);
		}
		// if not improvedClusterByBands(tempTetraList,copyEl,self)		
		if (clusterWasImproved)
		{				
			if (!improvedCluster(tempTetraList,this, false))
			{
					tempTetraList->Clear();
					continue;
			}
			goodMinQuality = actualQuality;
			goodTetraList->Assign(tempTetraList);
			// First found
			isBetter = true;
		}
		tempTetraList->Clear();
		if (isBetter) break;
	}
	elements->Clear();
	int result = goodTetraList->Count();
	numChanges = result;
    #if !defined(KRATOS)	
	endProcess((char*)("evaluateSet"));
	#endif
	
	return result;
}

void  TElementsCluster::genElements(bool vertexIsNew)
{
	int i ;  


	if ((!testCenter) && (vC != NULL)) 
		delete vC;

	if (newElements->Count()>0)
	{
		for (i = 0 ;i<newElements->Count() ; i++)
			originalMesh->elementsToAdd->Add(newElements->elementAt(i));

		for (i = 0 ;i<copyEL->Count() ; i++)
			originalMesh->elementsToRemove->Add(copyEL->elementAt(i));

		if (testCenter && vertexIsNew)
			originalMesh->vertexes->Add(vC);

		newElements->Clear();
		copyEL->Clear();
	}

}


bool TElementsCluster::updateMesh(bool checkIfInvalid )
{
	#if !defined(KRATOS)	
	startProcess("updateMesh");
	#endif
	int i ;
	TVertex *v0,*v1,*v2,*v3;
	TTetra *t ;
	bool /*centerUsed,*/  invalidCOnfig , result;

//	centerUsed = false;
	result = false;
	newElements->Clear();
	if (goodTetraList->Count()>0)
	{
		copyEL->Assign(inspectedElements);
	

		// Extraigo los elementos anteriores
		for (i = 0 ; i<copyEL->Count() ; i++)
		{
			TTetra *_t = (TTetra*)(copyEL->elementAt(i));      			
			_t->isdestroyed = true;
			_t->removeVertexRef();
		}
		//Añado los vertices centrales
//		if ((vC!= NULL) && (goodTetraList->indexOf(vC)>0) )
//			centerUsed = true;


		// Creo los nuevos elementos
		for (i = 0 ; i<goodTetraList->Count()  ; i = i+4)
		{

			v0 = (TVertex*)(goodTetraList->elementAt(i));
			v1 = (TVertex*)(goodTetraList->elementAt(i+1));
			v2 = (TVertex*)(goodTetraList->elementAt(i+2));
			v3 = (TVertex*)(goodTetraList->elementAt(i+3));

			t = new TTetra(NULL,v0,v1,v2,v3,true);          
			newElements->Add(t);
		}

		if (checkIfInvalid)
		{
			invalidCOnfig = false;
			
			tempL->Clear();
			for (i = 0 ; i<newElements->Count() ; i++)
			{
				TTetra* _t = (TTetra*)( newElements->elementAt(i));			
				if (_t->getNeighboursByFace(1 , tempL)->Count()>4)
				{
					//---  Invalid!!
					invalidCOnfig = true;
					break;
				}
			}

			if (invalidCOnfig )
			{
				//remuevo los elementos nuevos
				for (i = 0 ; i<newElements->Count();i++)
				{
					TTetra *_t = (TTetra *)(newElements->elementAt(i));					
					if (_t != NULL)
					{
					  _t->removeVertexRef();
					  newElements->setElementAt(i,NULL);
					  delete _t;
					}
				}
				newElements->Clear();

				// recupero los viejos
				for (i = 0 ; i<copyEL->Count();i++)
				{
					TTetra* _t = (TTetra*)( copyEL->elementAt(i));					
					_t->isdestroyed = false;
					_t->updateVertexRef();
				}
				if (vC!= NULL)
				{ 
					delete vC;
					vC = NULL;
				}
				newElements->Clear();
				copyEL->Clear();
			}  
			else
				result = true;
		}

		goodTetraList->Clear();
	}
	else
		result = true;
	//genElements();
	#if !defined(KRATOS)	
	endProcess("updateMesh");
	#endif
	return result;

} ;

bool TElementsCluster::generateSubMesh(double minExpectedQuality , double minExpectedAngle)
{
	int i;
	TTetra* t;

	elements->Clear();
	surfaceT->Clear();
	copyEL->Clear();
	vertexes->Clear();
	tempTetraList->Clear();
	goodTetraList->Clear();
	tempL->Clear();
	vC = NULL;

	for (i = 0 ; i< inspectedElements->Count() ; i++)
	{
		TTetra* _t = (TTetra*)( inspectedElements->elementAt(i));
		if (_t == NULL) continue;
		if (! _t->isdestroyed )
			elements->Add(_t);
	}

	//-- Asocio a mano los vecinos
	minQuality = 5000000;
	minDAngle = 50000000;
	for (i = 0 ; i< elements->Count() ;i++)
	{
		TTetra* _t =t = (TTetra*)( elements->elementAt(i));
		if (_t == NULL) continue;		
		_t->calidad = fc(_t->vertexes);       
		_t->fMinDiedralAngle = diedralAngle(_t->vertexes);
		minQuality = Min(minQuality,_t->calidad );
		minDAngle = Min(minDAngle,_t->fMinDiedralAngle);
	}
	//---------------------------------------
	innerGetElementsSurface(elements,surfaceT, vertexes);
	// dgetSurfaceTriangles(elements, surfaceT, vertexes);
	//--------------------------------------
	return true;
}


//************************************************
// Optimizacion por CARA
//************************************************
void evaluateClusterByFace(TMesh *aMesh , double minExpectedQuality,TVertexesEvaluator fc)
{
	TList<TObject*> *inspectedElements , *elements, *nFaces;
	int i,iv,j;
	double minQuality ,   meshMinQ;
	TElementsCluster *aCluster;

	meshMinQ = 50000;
	for (i = 0 ; i<aMesh->elements->Count() ; i++)
	{
		TTetra* _t = (TTetra*)(aMesh->elements->elementAt(i));
		_t->isdestroyed = false;
		_t->calidad = fc(_t->vertexes);
		meshMinQ = Min(meshMinQ , _t->calidad);
	} //Creacion de variables
	aMesh->elementsToRemove->Clear();
	aMesh->selectedElements->Clear();
	aCluster = new TElementsCluster(aMesh,fc);
	elements = aMesh->elements;

	nFaces = new TList<TObject*>();
	#if !defined(KRATOS)	
	startProcess((char*)("evaluateClusterByFace"));
	#endif
	inspectedElements = new TList<TObject*>();
	for (iv = 0 ; iv<elements->Count()  ; iv++)
	{
		TTetra *_t = (TTetra*)(elements->elementAt(iv)); 
		if (_t == NULL ) continue;
		if (_t->isdestroyed) continue;
		nFaces->Clear();
		_t->getNeighboursByFace(1,nFaces);

		if (nFaces->Count() == 0) continue;

		for (j = 0  ; j<nFaces->Count() ; j++) 
		{
			TTetra *_t2 = (TTetra*)(nFaces->elementAt(j));
			if (_t2 == NULL) continue;
			if (_t2->isdestroyed) continue;

			inspectedElements->Clear();
			inspectedElements->Add(_t);
			inspectedElements->Add(_t2);

			aCluster->inspectedElements->Assign(inspectedElements);

			minQuality =aCluster->getMinQuality();
			if (minQuality>minExpectedQuality) continue;
			//---------------------
			//--Limpio las variables
			//Cluster 1

			aCluster->doRemoveElements = false;

			aCluster->testCenter = false;
			aCluster->perturbCenter = 0.0;
			aCluster->checkMaxLength =  false;
			aCluster->generateSubMesh( );

			if (aCluster->evaluateSet()>0)
			{
				if (aCluster->updateMesh(true))
				{
				  aCluster->genElements();
				  _t->isdestroyed = true;
				  _t2->isdestroyed = true;
				}
				break;
			}
		}
	}

	aMesh->updateRefs();
	#if !defined(KRATOS)	
	endProcess((char*)("evaluateClusterByFace"));
	#endif


	delete aCluster;
	//aMesh.updateIndexes(0);
}

//************************************************
// Optimizacion por NODO
//************************************************
void evaluateClusterByNode(TMesh *aMesh , double minExpectedQuality,TVertexesEvaluator fc)
{
	//declaracion de variables
	
	TList<TVertex*> *vertexesCopy;
	int i,iv;
	TVertex* inspVertex;
	double minQuality ,   meshMinQ;
	TElementsCluster *aCluster  ;

	//cuerpo
	meshMinQ = 50000;
	for (i = 0 ;i< aMesh->elements->Count(); i++)
	{
		TTetra *t = (TTetra*)(aMesh->elements->elementAt(i));
		if (t==NULL) continue;
		t->calidad =vrelaxQuality( t->vertexes);
		meshMinQ = Min(meshMinQ,  t->calidad);
	}
	//Creacion de variables

	aCluster = new TElementsCluster(aMesh,fc);
	aCluster->doRemoveElements = false;
	aCluster->testCenter = false;
	aCluster->perturbCenter = 0.0f;
	aCluster->checkMaxLength = false;

	vertexesCopy = new TList<TVertex*>();
	vertexesCopy->Assign(aMesh->vertexes);
	#if !defined(KRATOS)	
	startProcess((char*)("evaluateClusterByNode"));
	#endif
	
	////-- Facil de Paralelizar
	for (iv = 0 ; iv<vertexesCopy->Count() ; iv++)
	{
		inspVertex =(TVertex*)(vertexesCopy->elementAt(iv));
		if (inspVertex->elementsList== NULL)  continue;
		inspVertex->elementsList->Pack();
		if (inspVertex->elementsList->Count() == 0) continue;
		//--Limpio las variables
		aCluster->inspectedElements->Assign( inspVertex->elementsList) ;		

		minQuality =aCluster->getMinQuality();
		if (minQuality>minExpectedQuality) continue;
		//---------------------
		//--Limpio las variables          
		//Cluster 1
		//1.Genero la malla de superficie
		aCluster->generateSubMesh( );

		aCluster->testCenter = false;

		//2.Evaluo si el resultado es mejor
		if (aCluster->evaluateSet()>0) 
		{
			//3.Actualizo la malla
			if (aCluster->updateMesh(true))
				aCluster->genElements();
		}

		if ( aMesh->elementsToRemove->Count() > ELEMENTS_TO_FORCE_UPDATE)
			((TVolumeMesh*)(aMesh))->updateRefs();

	} 
	((TVolumeMesh*)(aMesh))->updateRefs();
	#if !defined(KRATOS)	
	endProcess((char*)("evaluateClusterByNode"));
	#endif
	
	delete (vertexesCopy);
	delete (aCluster);

}

//************************************************
// Optimizacion por ARISTA
//************************************************
void evaluateClusterByEdge(TMesh *aMesh , double minExpectedQuality,TVertexesEvaluator fc)
{

	int i,j,k ;
	TList<TObject*> *inspElements ; //, elements,vertexes,surfaceT,copyEL,vl: TList;
	TList<TVertex*> *vl;
	TVertex *v0, *v1;
//	double meshMinQ;
	TTetra *t;
	TElementsCluster *aCluster  ;

//	meshMinQ = 50000;
	aMesh->elementsToRemove->Clear();
	aMesh->selectedElements->Clear();

	aCluster = new TElementsCluster(aMesh,fc);
	aCluster->doRemoveElements = false;
	aCluster->testCenter = false;
	aCluster->perturbCenter = 0.0;
	aCluster->checkMaxLength =  false;
			
	inspElements = new TList<TObject*>();
	vl = new TList<TVertex*>();
	#if !defined(KRATOS)	
	startProcess((char*)("evaluateClusterByEdge"));
	#endif
	//--------------------------------------------
	for (i = 0 ; i<aMesh->vertexes->Count() ; i++)
	{
		v0 =  aMesh->vertexes->elementAt(i);
		vl->Clear();
		v0->getVertexNeighboursByElem(vl,1,true);

		if (vl == NULL) continue;
		if (vl->Count() == 0 ) continue;

		//  Recorro los vertices vecinos
		for ( j = 0 ; j<vl->Count() ; j++)
		{
			v1 = (TVertex*)(vl->elementAt(j));
			if (v0->getID()>=v1->getID()) continue;
			if (v1->elementsList == NULL) continue;
			inspElements->Clear();
			//--- Veo los vecinos de arista de ambos vertices
			//Para el vertice 0
			for (k = 0 ; k<v0->elementsList->Count() ; k++)
			{
				t = (TTetra*)(v0->elementsList->elementAt(k));
				if (t == NULL) continue;
				if (t->isdestroyed) continue;
				if (!t->hasEdge(v0,v1)) continue;
				if (t->innerFlag == i*1000+j) continue;
				t->innerFlag = i*1000+j;
				inspElements->Add(t);
			}
			//Para el vertice 1
			for (k = 0 ; k<v1->elementsList->Count() ; k++)
			{
				t = (TTetra*)(v1->elementsList->elementAt(k));
				if (t == NULL) continue;
				if (t->isdestroyed) continue;
				if (!t->hasEdge(v0,v1)) continue;
				if (t->innerFlag == i*1000+j) continue;
				t->innerFlag = i*1000+j;
				inspElements->Add(t);
			}

			//-------------------
			//-- Asigno la superficie
			#if !defined(KRATOS)	
			startProcess("generateSubMesh");
			#endif
			aCluster->inspectedElements->Assign( inspElements);
			aCluster->generateSubMesh( );
			#if !defined(KRATOS)	
			endProcess("generateSubMesh");
			#endif
			if (aCluster->evaluateSet()>0) 
			{
				if (aCluster->updateMesh(true))
					aCluster->genElements();
				break;
			}

		}

	}
	#if !defined(KRATOS)	
	endProcess((char*)("evaluateClusterByEdge"));
	startProcess((char*)("updateRefs"));
	#endif
	aMesh->updateRefs();
	#if !defined(KRATOS)	
	endProcess((char*)("updateRefs"));
	#endif
	delete aCluster;

}


//
// Metodo de control para verificar si mejora un cluster. Teni}o en cuenta la minima CALIDAD
//

bool improvedTetra(TVertex *v0,TVertex *v1,TVertex *v2,TVertex *v3,	              
	               TElementsCluster* cl,
	                bool maxEdgeLengthConstrain )
{
	double goodMinQuality,  tempQuality ;

	//Datos ya calculados	
	goodMinQuality = cl->goodMinQuality;	
	
	tempQuality = vrelaxQuality(v0,v1,v2,v3);
	//Si la calidad es peor, termino
	//if (tempQuality > goodMinQuality)
	if (tempQuality<goodMinQuality)
				return false;
	return true;
}
//
// Metodo de control para verificar si mejora un cluster. Teni}o en cuenta la minima CALIDAD
//
bool improvedCluster(TList<TObject*>* c1,	
	TElementsCluster* cl,
	bool maxEdgeLengthConstrain )
{

	double avgEdgeLength , goodMinQuality, origDAngle ;
	TVertex* v0, *v1 , *v2, *v3  ;
	int j;
	double tempQuality , maxEdgeLength ,mindAngle ;

	//Datos ya calculados
	avgEdgeLength = cl->avgEdgeLength ;
	goodMinQuality = cl->goodMinQuality;
	origDAngle = cl->minDAngle;
	

	for (j = 0 ; j<(c1->Count() / 4) ;j++)
	{
		v0 =  (TVertex*)(c1->elementAt(4*j));
		v1 =  (TVertex*)(c1->elementAt(4*j+1));
		v2 =  (TVertex*)(c1->elementAt(4*j+2));
		v3 =  (TVertex*)(c1->elementAt(4*j+3));

		tempQuality = vrelaxQuality(v0,v1,v2,v3);
		mindAngle = minDiedralAngle(v0->fPos,v1->fPos,v2->fPos,v3->fPos);
		//Si la calidad es peor, termino
		if ( ((tempQuality>0) && ( mindAngle>origDAngle) ) || 
			 ( (tempQuality<0) &&(tempQuality>goodMinQuality*ACCEPTANCE_TOLERANCE) ) )
		{
			maxEdgeLength = getmaxEdgeLength(v0,v1,v2,v3);
			if ((maxEdgeLengthConstrain) && (maxEdgeLength/avgEdgeLength>4))
			{
				return false;
			}
		}
		else
			return false;
	}

	return true;
}

double testTetraSplit4(TTetra *t, TVertex* v, TList<TObject*>* lRes , TVertexesEvaluator qualityFunction)
{
	TVertex *vs[4];
	int i;
	double q;
	q = 5000000;
	for (i = 0 ; i<4 ; i++)
	{
		vs[0] = t->vertexes[ TTetraFaces[i*3] ];
		vs[1] = t->vertexes[ TTetraFaces[i*3+1] ];
		vs[2] = t->vertexes[ TTetraFaces[i*3+2] ];
		vs[3] = v;
		q = Min(q , qualityFunction(vs) );
		if (lRes != NULL )
		{
			lRes->Add(vs[0]);
			lRes->Add(vs[1]);
			lRes->Add(vs[2]);
			lRes->Add(vs[3]);
		}
	}
	return q;
}



int vertexTetraReInsertion(TMesh *am ,TList<TVertex*>* vertexesList )
{
	int i,j,res;
	TVertex *v;
	TTetra *t, *tNew ;
	double minMeshQ,oldQ , newQ;
	TList<TObject*>* lRes;
	TVolumeMesh *aMesh = (TVolumeMesh*)(am);

	lRes = new TList<TObject*>();
	minMeshQ = 50000000;
	res = 0;
	for (i = 0 ; i<aMesh->elements->Count() ; i++)
	{
		t = (TTetra*)(aMesh->elements->structure[i]);
		minMeshQ = Min(minMeshQ ,diedralAngle(t->vertexes) );
	}


	for ( i = 0 ; i<vertexesList->Count() ; i++)
	{
		v = (TVertex*)(vertexesList->elementAt(i));
		t = aMesh->isPointInside(v->fPos);
		// Ningun tetra lo contiene
		if (t == NULL)  continue;      
		oldQ = minMeshQ;

		lRes->Clear();
		//     newQ := testTetraSplit6(aMesh,t,v,lRes, qualityFunction,t2);

		newQ = testTetraSplit4(t,v,lRes,diedralAngle);
		// No empeora la calidad actual
		if (newQ>=oldQ)  
		{
			//Pude volver a insertar el nodo          
			vertexesList->setElementAt(i, NULL); 
			res ++ ;
			
			// Create new Tetra!!
			for (j = 0 ; j<lRes->Count();j=j+4)
			{
				tNew = new TTetra(NULL, (TVertex*)(lRes->structure[j]),
					(TVertex*)(lRes->structure[j+1]),
					(TVertex*)(lRes->structure[j+2]),
					(TVertex*)(lRes->structure[j+3]));
				aMesh->elements->Add(tNew);
			}

			//Remove old TTetra
			aMesh->elements->Extract(t);
			t->removeVertexRef();
			delete t;
		}
	}

	vertexesList->Pack();
	return res;
}


