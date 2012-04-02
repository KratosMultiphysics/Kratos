#include "u_Types.h"
#include "u_tools.h"
#include "u_TetraFunctions.h"
#include "u_qualityMetrics.h"
#include "u_ParallelFunctions.h"

bool swapTetra(TVertex* v0,TVertex*  v1,TVertex*  v2,TVertex*  v3,TVertex* v4)
{

	double q1,q2 , nq1,nq2 , vl1,vl2,vl3,vl4 ;

	vl1 =tetraVolume(v0->fPos,v1->fPos,v2->fPos,v3->fPos);
	vl2 =tetraVolume(v0->fPos,v4->fPos,v1->fPos,v3->fPos);
	vl3 =tetraVolume(v1->fPos,v2->fPos,v4->fPos,v3->fPos);
	vl4 =tetraVolume(v0->fPos,v4->fPos,v2->fPos,v3->fPos);
	if ((vl1+vl2) != (vl3+vl4))  return false;

	q1 = vrelaxQuality(v0,v1,v2,v3);
	q2 = vrelaxQuality(v0,v4,v1,v3);

	nq1 = vrelaxQuality(v1,v2,v4,v3);
	nq2 = vrelaxQuality(v0,v4,v2,v3);

	return  Min(q1,q2)<Min(nq1,nq2);
}


void igetSurfaceTriangles(TList<TObject*> *elements , TList<TObject*> *res, TList<TObject*> *surfVertexes  )
{
	int i,j;
	TList<TObject*>* tr;
	TList<TObject*>* _tempL;    

	_tempL =new TList<TObject*>();

	for (i = 0 ;i< elements->Count(); i++)
	{
		TTetra* t = (TTetra*)(elements->elementAt(i));
		if (t==NULL) continue;

		tr = t->getSurfaceTriangle(false, _tempL,NULL);

		if (tr!=NULL)
		{
			for ( j = 0 ; j<tr->Count();j++)
			{
				TVertex* v =(TVertex*)(tr->elementAt(j)); 
				res->Add(v);

				if (surfVertexes!=NULL)
				{
					if (surfVertexes->indexOf(v)<0)
						surfVertexes->Add(v);
				}
			}
		}
		_tempL->Clear();
	}

}

bool innerGetElementsSurface(TList<TObject*>* elements,TList<TObject*>* surfaceT,TList<TObject*>* vertexes)
{
	int i,j ,k;
	bool res;
	TVertex *v0, *v1, *v2;
	TTetra *t,*t2;  
	bool faceVisited;


	for (i = 0; i<elements->Count() ; i++)
	{
		t =(TTetra*)( elements->elementAt(i));
		if ( t == NULL) continue;
		for (j = 0 ;j<4 ;j++)
		{
			v0 = t->vertexes[ TTetraFaces[j*3] ];
			v1 = t->vertexes[ TTetraFaces[j*3+1]];
			v2 = t->vertexes[ TTetraFaces[j*3+2]];
			faceVisited = false;
			// Esta cara ya fue visitada
			for (k = 0 ; k<surfaceT->Count()/3 ; k = k+1)
			{
				if( ((surfaceT->elementAt(3*k) == v0)   ||  (surfaceT->elementAt(3*k) == v1)   || (surfaceT->elementAt(3*k) == v2)) &&
					((surfaceT->elementAt(3*k+1) == v0) ||  (surfaceT->elementAt(3*k+1) == v1) || (surfaceT->elementAt(3*k+1) == v2)) &&
					((surfaceT->elementAt(3*k+2) == v0) ||  (surfaceT->elementAt(3*k+2) == v1) || (surfaceT->elementAt(3*k+2) == v2)) )
				{
					faceVisited = true;
					break;
				}
			}
			if (faceVisited) continue;


			res = true;
			for (k = 0 ; k<v0->elementsList->Count() ; k++)
			{
				t2 = (TTetra*)(v0->elementsList->elementAt(k));
				if (t2 == t) continue;
				if (elements->indexOf(t2)<0) continue;
				if (t2->hasFace(v0,v1,v2)) { res = false; break; };
			}

			if (!res) continue;
			for (k = 0 ; k<v1->elementsList->Count(); k++)
			{
				t2 = (TTetra*)(v1->elementsList->elementAt(k));
				if (t2 == t) continue;
				if (elements->indexOf(t2)<0) continue;
				if (t2->hasFace(v0,v1,v2)) { res = false; break; };
			}

			if (!res) continue;
			for (k = 0; k<v2->elementsList->Count() ; k++)
			{
				t2 = (TTetra*)(v2->elementsList->elementAt(k));
				if (t2 == t) continue;
				if (elements->indexOf(t2)<0) continue;
				if (t2->hasFace(v0,v1,v2)) { res = false; break; };
			}

			if (res)
			{
				surfaceT->Add(v0);
				surfaceT->Add(v1);
				surfaceT->Add(v2);
				if (vertexes->indexOf(v0)<0)     vertexes->Add(v0);
				if (vertexes->indexOf(v1)<0)     vertexes->Add(v1);
				if (vertexes->indexOf(v2)<0)     vertexes->Add(v2);
			}
		}
	}
	return NULL;
}
// -----------------------------
// Remove vertexes that do not have any related element 
void TVolumeMesh::removeFreeVertexes()
{
	for (int i=0; i<vertexes->Count() ; i++)
	{ 
		if (vertexes->structure[i]->elementsList->Count() == 0)
			vertexes->setElementAt(i,NULL);
	}

	vertexes->Pack();
}

TVolumeMesh::TVolumeMesh():
TMesh()
{
	elements = new TList<TObject*>();
	elementsToAdd = new TList<TObject*>();
	elementsToRemove= new TList<TObject*>();
}

TVolumeMesh::~TVolumeMesh()
{
	delete vertexesToRemove;

	for (int i=0 ; i<elements->Count() ; i++)
	{
		TTetra *t = (TTetra *)(elements->structure[i]);
		delete t;
	}
	delete elements;	
	delete elementsToAdd;
	delete elementsToRemove;
	for (int i=0 ; i<vertexes->Count() ; i++)
	{
	    TVertex *v =vertexes->structure[i];
		delete v;
	}
	delete vertexes;
	
	if (fFaces != NULL)
	{
		for (int i=0 ; i<fFaces->Count() ; i++)
		{
			TTriangle *tr =(TTriangle *)(fFaces->structure[i]);
			delete tr;
		}
		fFaces->Clear();
		delete fFaces;		
	}
	fFaces = NULL;
	vertexes = NULL;
	elements = NULL;
	vertexesToRemove = NULL;
	selectedElements = NULL;
}
void TVolumeMesh::updateIndexes(int flag )
{
	numVertices = vertexes->Count();
	numElements = elements->Count();

	fMassCenter = 0.0f;
	//- Limpio la estructura
	for (int i=0; i< numVertices; i++)
	{
		TVertex* v= vertexes->elementAt(i);
		if (v==NULL) continue;  
		if ( (flag & KEEP_ORIG_IDS) == 0)
			v->id = i;
		fMassCenter = fMassCenter +v->fPos;
		if (v->elementsList!=NULL)
			v->elementsList->Clear();
		else
			v->elementsList = new TList<TObject*>();
	}

	// Regeenero la lista de elementos vecinos
	for (int i=0 ; i<elements->Count(); i++)
	{
		TTetra* t =  (TTetra*)this->elements->elementAt(i);
		if (t == NULL ) continue;
		if ((flag & KEEP_ORIG_IDS) == 0) 
			t->id = i;
		for (int k=0;k<4;k++)
		{
			if ( t->vertexes[k]->elementsList !=NULL) 
			{
				t->vertexes[k]->elementsList->Add(t);
			}
		}
	}
	// Remuevo vertices que no pertenecen a ningun elemento
	vertexesToRemove->Clear();
	for (int i=0; i<vertexes->Count();i++)
	{
		TVertex *_v = vertexes->elementAt(i);
		if (_v->elementsList->Count()== 0 )
			vertexesToRemove->Add(_v);
	}
	/*
	for i:=0 to vertexes.count-1 do
	begin
	if (vertexes[i].elementsList=nil) or
	(vertexes[i].elementslist.count = 0) then
	vertexes[i] := nil;
	end;
	vertexes.Pack();
	*/
	numVertices = vertexes->Count();
	//--- obtengo los triangulos de la superficie
	if (flag & GENERATE_SURFACE)
		//getSurfaceTriangles();
		fastGetSurfaceTriangles(this);
}

void TVolumeMesh::updateRefs()
{
	int i;

	elements->Pack();
	elementsToAdd->Pack();
	elementsToRemove->Pack();
	//std :: cout << "elements to update " << elements->Count() << " " << elementsToAdd->Count() << " " << elementsToRemove->Count() << "\n";
	for (i = 0 ; i< elements->Count() ;i++)
	{
		TTetra* _t = (TTetra*)(elements->elementAt(i));		
		_t->flag = 1;
	}

	//-- remuevo duplicados
	for (i = 0 ; i< elementsToRemove->Count() ; i++)
	{
		TTetra* _nt = (TTetra*)(elementsToRemove->elementAt(i));				
		_nt->flag = 2;
	}
	//elementsToRemove->Pack();*/
	// Guardo cuantos elementos tengo
	int origNumElements = elements->Count();

	for (i = 0 ;i<elementsToAdd->Count() ; i++)
	{
		TTetra* _nt = (TTetra*)(elementsToAdd->elementAt(i));		
		this->elements->Add(_nt);
	}
	//std :: cout << "elements to update :Paso 1 " << elements->Count() << " " << elementsToAdd->Count() << " " << elementsToRemove->Count()<< "\n";
	for (i = 0 ; i< elements->Count() ; i++)
	{
		TTetra* _nt = (TTetra*)(elements->elementAt(i));		
		if ( _nt->flag==2)
		{
			elements->setElementAt(i,NULL);
			//delete _nt;
		}
	}


	int origNumElementsToRemove = elementsToRemove->Count();
	int origNumElementsToAdd = elementsToAdd->Count();
	//std :: cout << "elements to update :Paso 2 " << elements->Count() << " " << elementsToAdd->Count() << " " << elementsToRemove->Count()<< "\n";
	elementsToAdd->Clear();
	elementsToRemove->Clear();	
	elements->Pack(); 	
	//std :: cout << "elements to update :Paso 3 " << elements->Count() << " " << elementsToAdd->Count() << " " << elementsToRemove->Count()<< "\n";
	if (origNumElements +origNumElementsToAdd - origNumElementsToRemove != elements->Count()) 	
		std :: cout << " Invalid update Ref Configuration : "  << "\n";

	vertexes->Pack();
	fFaces->Pack();

}
/*
void TVolumeMesh::updateRefs()
{
int i,ind ;

//-----
for (i = 0 ; i< elements->Count() ; i++)
((TTetra*)(elements->elementAt(i)))->flag = 1;

for (i = 0 ; i< elementsToRemove->Count() ; i++)
((TTetra*)(elementsToRemove->elementAt(i)))->flag = 2;

for (i = 0 ; i< elementsToAdd->Count() ; i++)
{
TTetra* _t =(TTetra*)(elementsToAdd->elementAt(i));
if (_t->flag != 2)
elements->Add(_t);
}

for (i = 0 ; i< elements->Count() ; i++)
{
TTetra* _t =(TTetra*)(elements->elementAt(i));
if (_t->flag == 2)
elements->setElementAt(i,NULL);
}

elementsToRemove->Clear();
elements->Pack();
elementsToAdd->Clear();
}
*/
TTetra* TVolumeMesh::isPointInside( float4 pos )
{
	int i;

	for ( i=0 ; i<elements->Count() ; i++)
	{
		TTetra *t = (TTetra*)(elements->structure[i]);
		if (t->isInside(pos)) 
			return t;
	}
	return NULL;
}

void TVolumeMesh::validate(bool showMessages)
{

	int i,j,nv,ne,ntv, result;
	TTriangle *tr;
	TList<TObject*>* nfList;

	result = 0;
	nfList = new TList<TObject*>();
	//---------
	for (i=0 ; i<fFaces->Count() ; i++)
	{
		tr = (TTriangle*)(fFaces->elementAt(i));
		for (j=0 ; j<3 ; j++)
		{
			tr->vertexes[j]->flag = -1;
		}
	}
	//- Limpio la estructura
	ntv = 0;

	i = 0;
	while (i<fFaces->Count())
	{
		tr = (TTriangle*) (fFaces->elementAt(i));
		tr->calcEdges();
		i++;
	}
	nv = 0;
	for (i = 0 ; i<vertexes->Count() ; i++ )
	{
		TVertex* _v = vertexes->elementAt(i);
	}
	// Control de Elementos
	if (ntv > 0)
	{
		if (showMessages)
			printMessage( (char *)("Triangulos con vecinos por arista "));
	}

	if (nv > 0) 
	{
		if (showMessages)
			printMessage( (char *)("Vertices sueltos "));
	}
	ne = 0;
	for (i = 0 ; i< elements->Count() ; i++)
	{
		TTetra *_t = (TTetra*)(elements->elementAt(i));
		_t->id = i;	 
		nfList->Clear();
		_t->getNeighboursByFace(1,nfList);
		if (nfList->Count()>4)
		{
			ne++;
		}
	}
	result +=ne;
	if (ne > 0)
	{
		if (showMessages) 
			printMessage((char *)("Elementos con mas de 4 vecinos por cara "));
	}

	if (result == 0)
	{
		if (showMessages) 
			printMessage((char *)("No se detectaron errores"));
	}

}

void TVolumeMesh::getSurfaceTriangles()
{
	int j;
	TList<TObject*>* res;
	TTriangle* t;

	this->fFaces->Clear();
	res = new TList<TObject*>();

	igetSurfaceTriangles(this->elements ,res,NULL);
	for (j = 0 ; j<res->Count(); j=j+3)
	{
		t = new TTriangle( (TVertex*)res->elementAt(j),  (TVertex*)res->elementAt(j+1),  (TVertex*)res->elementAt(j+2));

		t->calcNormal();
		this->fFaces->Add(t);       
	}
}

bool swapVolumeMesh(TVolumeMesh* aMesh)
{
	int i,k;
	TList<TObject*> *l, *l2;
	TTetra *t0,*t1;
	TTriangle *tr1;
	TVertex *v0,* v1,*v2,*v3,*v4;
	// aMesh.updateIndexes(0);
	l   = new TList<TObject*>();
	l2   = new TList<TObject*>();
	// l2 := TList.create; 
	// Para todos los elementos de superficie
	for (i = 0 ; i<aMesh->fFaces->Count();i++)
	{
		tr1 =(TTriangle*)( aMesh->fFaces->elementAt(i));
		tr1->calcEdges();
		tr1->id = i;
	}

	for (i = 0 ; i< aMesh->fFaces->Count() ; i++)
	{
		tr1 =(TTriangle*)( aMesh->fFaces->elementAt(i));
		v0 = tr1->vertexes[0];
		v1 = tr1->vertexes[1];
		v2 = tr1->vertexes[2];
		l->Clear();
		l->Assign(v0->elementsList,laAppend);
		l->Assign(v1->elementsList,laOr);    
		l->Assign(v2->elementsList,laOr);
		t0 = nil;
		t1 = nil;

		for (k=0; k<l->Count();k++)
		{
			TTetra* _t = (TTetra*)(l->elementAt(k));
			if (_t->hasFace(v0,v1,v2) )
			{
				t0 = _t;
				v3 = t0->oppositeVertex(v0,v1,v2);
				break;           
			}
		}

		for (k = 0 ; k< l->Count() ; k++)
		{
			TTetra* _t = (TTetra*)(l->elementAt(k));
			if ((_t!=t0) &&  (_t->hasFace(v0,v1,v3)))
			{
				l2->Clear();
				_t->getSurfaceTriangle(false,l2,NULL);
				if (l2->Count() == 0) continue;
				t1 = _t;
				break;
			}
		}
		//No encontro ningun elemento vecino
		if (t0 == NULL)  continue;
		if (t1 == NULL) continue;
		if (v3 == NULL)  continue;
		v4 = t1->oppositeVertex(v0,v1,v3);
		if (v4 == NULL) continue;

		// verifico si mejora
		if (swapTetra(v0,v1,v2,v3,v4) )
		{
			t0->removeVertexRef();
			t1->removeVertexRef();

			t0->Create(NULL,v1,v2,v4,v3);
			t1->Create(NULL,v0,v4,v2,v3);
			aMesh->selectedElements->Add(t0);
			aMesh->selectedElements->Add(t1);
			//break;
		}

	}
	aMesh->updateIndexes();
	return true;
}
