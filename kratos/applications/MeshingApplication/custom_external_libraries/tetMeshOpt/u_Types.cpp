#include "u_Types.h"
#include "u_qualityMetrics.h"

bool sortByID(TObject* i, TObject* j)
{
	return ((TValuedObject*)(i))->id < ((TValuedObject*)(j))->id ;
}
/// Class TValuedObject
TValuedObject::TValuedObject(void) 
{};
TValuedObject::~TValuedObject(void)
{};

float4 TValuedObject::getCenter() 
{ 	
	return Float4(0.0f,0.0f,0.0f); 
}

/// Class TVertex
TVertex::TVertex(float x,float y, float z)
{ 
	fPos.x= x; fPos.y = y; fPos.z = z; 
	neighTr = new  TList<TObject*>();
	elementsList = new TList<TObject*>();			 
}

TVertex::TVertex(float4 v)
{ 
	fPos.x= v.x; fPos.y = v.y; fPos.z = v.z; 
	neighTr = new  TList<TObject*>();
	elementsList = new TList<TObject*>();			 
}

TVertex::~TVertex()
{
	neighTr->Clear();
	elementsList->Clear();
	delete neighTr;
	delete elementsList;

}

void TVertex::calcNormal()
{
};

float4 TVertex::pos()
{ 
	return fPos; 
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
	if (mode == 0)
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
		if (neighV == NULL)  neighV = new TList<TObject*>();
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

// Class TTriangle

TTriangle::TTriangle(TVertex* v0,TVertex* v1,TVertex* v2)
{
	vertexes[0] = v0;
	vertexes[1] = v1;
	vertexes[2] = v2;
}
void TTriangle::calcNormal()
{
	this->normal =Normal(vertexes[0]->pos(),vertexes[1]->pos(),vertexes[2]->pos());
}

void TTriangle::calcEdges()
{
}

bool TTriangle::verticesEncomun(TTriangle* otherT,TVertex* &v1,TVertex* &v2)
{
	int i,j;
	v1=NULL;
	v2=NULL;

	for (i=0;i<2;i++) for (j=0;j<2;j++) 
	{

		if ( vertexes[i] == otherT->vertexes[i])
		{
			if (v1==NULL) 
				v1 = vertexes[i] ;
			else
				if (v1!=vertexes[i] ) 
					v2 = vertexes[i] ;
		}
	}

	return ((v1!=NULL) && (v2!=NULL));
}

//Class TElement
TElement::TElement() {}
TElement::TElement(TVertex v0,TVertex v1,TVertex v2,TVertex v3){};
BoundBox TElement::CalcBound() 
{ 
	BoundBox b = calcBound(vertexes[0]->fPos, vertexes[1]->fPos); 
	return  b;
};

// Class TTetra
void TTetra::update()
{
	double f, c, v;
	CalcAng(vertexes[0]->fPos, vertexes[1]->fPos, vertexes[2]->fPos , vertexes[3]->fPos,&f,&c,&v);
	fDiedralAngle = f;
	fFaceAngle = c;
	fVolume = v;
	//fVolume =tetraVolume( ;			
}


double TTetra::getVolume() 
{ 
	return tetraVolume(vertexes[0]->fPos, vertexes[1]->fPos, vertexes[2]->fPos , vertexes[3]->fPos);
}

double TTetra::getPerimeter() 
{ 
	return 0;
}

double TTetra::getQuadPerimeter() 
{
	return 0;
}

int TTetra::isInside(float4 ps )
{	
	int i ;
	TVertex *v0, *v1, *v2;
	
	//The iteration method
	for (i = 0 ; i<4 ; i++)
	{
		v0 = vertexes[TTetraFaces[3*i]];
		v1 = vertexes[TTetraFaces[3*i+1]];
		v2 = vertexes[TTetraFaces[3*i+2]];
		if (tetraVolume(v0->fPos,v1->fPos,v2->fPos,ps)<0)
			return 0;
	}
	return 1;
}

double TTetra::getminEdgeLength() 
{
	double distances[6] = {distanceD(vertexes[0]->fPos , vertexes[1]->fPos) , 
		distanceD(vertexes[0]->fPos , vertexes[2]->fPos) ,
		distanceD(vertexes[0]->fPos , vertexes[3]->fPos) ,
		distanceD(vertexes[1]->fPos , vertexes[2]->fPos) ,
		distanceD(vertexes[1]->fPos , vertexes[3]->fPos) ,
		distanceD(vertexes[2]->fPos , vertexes[3]->fPos) };
	return MinValue(distances, 6);
}

double TTetra::getmaxEdgeLength() 
{
	double distances[6] = {distanceD(vertexes[0]->fPos , vertexes[1]->fPos) , 
		distanceD(vertexes[0]->fPos , vertexes[2]->fPos) ,
		distanceD(vertexes[0]->fPos , vertexes[3]->fPos) ,
		distanceD(vertexes[1]->fPos , vertexes[2]->fPos) ,
		distanceD(vertexes[1]->fPos , vertexes[3]->fPos) ,
		distanceD(vertexes[2]->fPos , vertexes[3]->fPos) };
	return MaxValue(distances, 6);
}

double TTetra::getSurface(){return 0;}
void TTetra::copyFrom(TValuedObject e ){}
TTetra::TTetra(TObject* owner) { this->isdestroyed = false;  return ; }

void TTetra::Create(TObject* owner, TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3, bool publish )
{
	this->rColor = Float4(1.0f);
	this->isdestroyed = false; 
	this->initialized = false;
	this->NeighBourByFace = NULL;
	this->Neighbours = NULL;
	vertexes[0] = v0;
	vertexes[1] = v1;
	vertexes[2] = v2;
	vertexes[3] = v3;
}

TTetra::TTetra(TObject* owner, TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3, bool publish ) 
{
	this->Create(owner,v0,v1,v2,v3,publish);
	if (publish)
		updateVertexRef();
} 

TTetra::~TTetra() 
{
	removeVertexRef();
	//delete normals;
	//delete vertexes;
	if (Neighbours)
		delete Neighbours;
	if (NeighBourByFace)
		delete NeighBourByFace;
}
// Obtener el vecino que comparte estos 3 vertices  
TTetra* TTetra::getTetraNeighbour(int faceI, TVertex* v0,TVertex* v1,TVertex* v2, TList<TObject*>*  tl)
{
	int i;
	TTetra* t2;
	//----------------
	for (i=0;i <v0->elementsList->Count();i++)
	{
		t2 = (TTetra*)(v0->elementsList->elementAt(i));
		if (t2==this) continue;
		if (t2->isdestroyed) continue;

		if (t2->hasFace(v0,v1,v2))  { return t2; }

	}

	for (i=0;i <v1->elementsList->Count();i++)
	{
		t2 = (TTetra*)(v1->elementsList->elementAt(i));
		if (t2==this) continue;
		if (t2->isdestroyed) continue;

		if (t2->hasFace(v0,v1,v2))  { return t2; }

	}

	for (i=0;i <v2->elementsList->Count();i++)
	{
		t2 = (TTetra*)(v2->elementsList->elementAt(i));
		if (t2==this) continue;
		if (t2->isdestroyed) continue;

		if (t2->hasFace(v0,v1,v2))  { return t2; }				   
	}

	return nil;
}



// Obtener triangulos que sean de la superficie
TList<TObject*>* TTetra::getSurfaceTriangle(bool bmpMode , TList<TObject*>* res, TList<TObject*>* ts)
{ 
	int i;

	TVertex* v0;
	TVertex* v1;
	TVertex* v2;  

	TList<TObject*>*  result;				

	if (res !=NULL)
		result =res;
	else
		result = NULL ;

	for ( i=0;i<4;i++)
	{
		if (i==0) 
		{
			v0 = vertexes[0];   v1 = vertexes[1];     v2 = vertexes[2];
		} 
		else
			if (i==1)
			{
				v0 = vertexes[0];   v1 = vertexes[2];     v2 = vertexes[3];
			}
			else
				if (i==2)
				{
					v0 = vertexes[0];  v1 = vertexes[3];     v2 = vertexes[1];
				}
				else
				{
					v0 = vertexes[1];   v1 = vertexes[3];     v2 = vertexes[2];
				}

				if (getTetraNeighbour(i,v0,v1,v2,result)== NULL)
				{
					if (result==NULL)
						result = new TList<TObject*>();
					result->Add( v0);
					result->Add( v1);
					result->Add( v2);
				}
	}
	return result;
} 

// Verificar si comparte esta cara
bool TTetra::hasFace(TVertex* v0,TVertex*  v1,TVertex* v2) 			
{
	return ((v0 == vertexes[0]) || (v0 == vertexes[1]) || (v0 == vertexes[2]) || (v0 == vertexes[3])) &&
		((v1 == vertexes[0]) || (v1 == vertexes[1]) || (v1 == vertexes[2]) || (v1 == vertexes[3]))&&
		((v2 == vertexes[0]) || (v2 == vertexes[1]) || (v2 == vertexes[2]) || (v2 == vertexes[3]));
}

TVertex* TTetra::oppositeVertex(TVertex *v0,TVertex *v1,TVertex *v2)
{  
	if ((vertexes[0]!=v0) && (vertexes[0]!=v1) && (vertexes[0]!=v2))  return  vertexes[0];
	else if ((vertexes[1]!=v0) && (vertexes[1]!=v1) && (vertexes[1]!=v2))  return vertexes[1];
	else if ((vertexes[2]!=v0) && (vertexes[2]!=v1) && (vertexes[2]!=v2))  return  vertexes[2];
	else if ((vertexes[3]!=v0) && (vertexes[3]!=v1) && (vertexes[3]!=v2))  return  vertexes[3];
	else return NULL;
}
bool TTetra::hasEdge(TVertex *v0,TVertex *v1)
{   
	return ((v0 == vertexes[0]) || (v0 == vertexes[1]) || (v0 == vertexes[2]) || (v0 == vertexes[3])) &&
		((v1 == vertexes[0]) || (v1 == vertexes[1]) || (v1 == vertexes[2]) || (v1 == vertexes[3])) ;				
}
bool TTetra::hasFace(TTetra* t1 )
{
	int i ; 
	TVertex *v0, *v1,*v2;
	bool result;
	for (i=0 ; i<4 ; i++)
	{
		v0 = t1->vertexes[TTetraFaces[3*i]];
		v1 = t1->vertexes[TTetraFaces[3*i+1]];
		v2 = t1->vertexes[TTetraFaces[3*i+2]];
		result = hasFace(v0,v1,v2);
		if (result) return true;
	}
	return false;
}

int TTetra::hasVertex(TTetra t1 ){ return 0;} 
bool TTetra::hasVertex(TVertex v0 ) { return 0;}
TTetra* TTetra::getTetraNeighbour(int faceI ,TVertex v0,TVertex v1,TVertex v2, TList<TObject> *tl) { return NULL;}


void TTetra::replaceTriangle(TTriangle oldTr,TTriangle newTr ) {}
void TTetra::replaceVertex(TVertex oldV, TVertex newV){}


bool TTetra::isInvalid() {return 0;}

TList<TObject*>* TTetra::getNeighboursByFace(int depth ,TList<TObject*>* nFL ) 
{
	int i,j;
	TTetra* t2;
	TList<TObject*>* result ;

	if (nFL )
	{
		result = nFL;
	}
	else
	{
		if (this->NeighBourByFace == NULL) NeighBourByFace = new TList<TObject*>();
		result = NeighBourByFace;
	}
	result->Clear();
	if (this->isdestroyed )  return result;

	//----------------
	for (i = 0 ; i<4 ; i++)
	{
		for (j = 0 ; j< vertexes[i]->elementsList->Count() ; j++)
		{
			t2 =  (TTetra*)(vertexes[i]->elementsList->elementAt(j));
			if (t2 == NULL) continue;
			if  (t2== this) continue;
			if  (t2->isdestroyed) continue;

			if (!hasFace(t2)) continue;
			if ( result->indexOf(t2)<0 ) 
				result->Add(t2);

		}
	}
	return result;
}

TList<TObject>* TTetra::getNeighbours(int depth) {return NULL;}

void TTetra::removeVertexRef()
{
	clearVertexRef();
	for (int i=0 ; i<4 ; i++)
		if (vertexes[i]->elementsList != NULL )
			vertexes[i]->elementsList->Pack();
}
void TTetra::updateVertexRef()
{
	for (int i = 0 ; i<4 ; i++)
	{				   
		if (vertexes[i]->elementsList == NULL)     
			vertexes[i]->elementsList = new TList<TObject*>();
		if ( vertexes[i]->elementsList->indexOf(this)<0 )
			vertexes[i]->elementsList->Add(this);

	}
}


void TTetra::clearVertexRef()
{
	for (int i=0; i<4 ; i++) 
	{
		if (vertexes[i]->elementsList != NULL)
		{					 
			TList<TObject*>* _lv = vertexes[i]->elementsList;
			for (int j = 0 ; j<_lv->Count()  ; j++)
			{
				if (_lv->elementAt(j) == this)
					_lv->setElementAt(j,NULL);
			}
		}
	}

}

/// CLASS TElementsPool
TElementsPool::TElementsPool()
{
	availableElements = new TList<TTetra*>();
}

TTetra* TElementsPool::getInstance()
{
	if (availableElements->Count() == 0)
		return new TTetra(NULL);
	else
	{
		TTetra *t = availableElements->structure.back();
		availableElements->structure.pop_back();
		return t;
	}
}

TTetra* TElementsPool::getInstance(TVertex *v0,TVertex *v1,TVertex *v2,TVertex *v3)
{
	if (availableElements->Count() == 0)
		return new TTetra(NULL,v0,v1,v2,v3);
	else
	{
		TTetra *t = availableElements->structure.back();
		availableElements->structure.pop_back();
		t->vertexes[0] = v0;
		t->vertexes[1] = v1;
		t->vertexes[2] = v2;
		t->vertexes[3] = v3;
		t->updateVertexRef();
		return t;
	}
}


void TElementsPool::releaseInstance(TTetra *t)
{
	availableElements->structure.push_back(t);
}

void TElementsPool::releaseMemmory()
{
	for (int i= 0 ; i<availableElements->Count() ; i++)
	{
		TTetra *t = availableElements->structure[i];
		delete t;
	}
	availableElements->Clear();
}



