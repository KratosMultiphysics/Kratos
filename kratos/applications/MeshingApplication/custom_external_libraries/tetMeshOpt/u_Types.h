#pragma once

#include "u_delphiClasses.h"
#include "Math3D.h"

//------ MENUS
#define SURFACE 1
#define VOLUME 2
#define COLORIZE 3
#define UPDATE_IDS  0
#define GENERATE_SURFACE  1
#define KEEP_ORIG_IDS  2

typedef char* String;
// Solution using a typedef: Define a pointer to a function which is taking
// two floats and returns a float
typedef float(*fncEvaluation)(TObject*);

const  int TTetraFaces[] = {0,1,2, // Face 0
	0,2,3,// Face 1
	0,3,1,// Face 2
	1,3,2 }; // Face 3]

bool sortByID(TObject* i, TObject* j);

struct  object
{
	virtual ~object() = 0;
};

class TValuedObject : public TObject
{
public:
	double	fmetrica , sortvalue ,calidad;
	bool marked, changed ,locked,isdestroyed;
	int id , groupID,visited,flag , innerFlag;
	void* linkSet ;
	int _color ;
	float4 rColor;
	object* userData;

	TValuedObject(void) ;
	~TValuedObject(void);
	float4 getCenter() ;
};

class TNeighboured : public TValuedObject{};

class TVertex: public TNeighboured
{
public:   
	float4 fPos ;
	int fixed ;
	TList<TObject*>* neighTr;
	TList<TObject*>* elementsList;
	TList<TObject*>* neighV;
	bool isSurface;
	float4 normal ;

	TVertex(float x,float y, float z);
	TVertex(float4 v);
	~TVertex();
	void calcNormal();
	float4 pos();
	TList<TObject*> *getVertexNeighboursByElem(TList<TObject*> *toL = NULL , int depth = 1);			 
	TList<TObject*>* getVertexNeighboursByTriangle(TList<TObject*>* toL = NULL , int depth = 1);
};

class TVertex3D : public TVertex
{

};

typedef double(*TVertexesEvaluator)(TVertex*[]);

class TPolygon: public TNeighboured
{
private : 

public :
	TVertex* vertexes[3];
	float4 normal;
};


class TTriangle: public TPolygon
{
public :    
	TTriangle(TVertex* v0,TVertex* v1,TVertex* v2);

	void calcNormal();

	void calcEdges();

	bool verticesEncomun(TTriangle* otherT,TVertex* &v1,TVertex* &v2);
};



class TElement : public TNeighboured
{
public:      
	TVertex* vertexes[4];
	float4 normals[4];
	bool initialized ;
	TElement();
	TElement(TVertex v0,TVertex v1,TVertex v2,TVertex v3);
	virtual BoundBox CalcBound() ;	
};


class TTetra : public TElement
{
public : 
	double fVolume,fDiedralAngle, fFaceAngle,fSurface ;
	int ts[6];

	TTriangle* triangles[4];
	TList<TObject*>* NeighBourByFace;
	TList<TValuedObject>* Neighbours;

	int materialID;
	float4 velocity;

	double getVolume() ;

	double getPerimeter() ;

	double getQuadPerimeter() ;

	double getminEdgeLength() ;

	double getmaxEdgeLength() ;

	double getSurface();
	void copyFrom(TValuedObject e );
	TTetra(TObject* owner);

	void Create(TObject* owner, TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3, bool publish = true);

	TTetra(TObject* owner, TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3, bool publish = true) ;

	~TTetra() ;
	// Obtener el vecino que comparte estos 3 vertices  
	TTetra* getTetraNeighbour(int faceI, TVertex* v0,TVertex* v1,TVertex* v2, TList<TObject*>*  tl);

	// Obtener triangulos que sean de la superficie
	TList<TObject*>* getSurfaceTriangle(bool bmpMode , TList<TObject*>* res, TList<TObject*>* ts);

	// Verificar si comparte esta cara
	bool hasFace(TVertex* v0,TVertex*  v1,TVertex* v2) 			;

	TVertex* oppositeVertex(TVertex *v0,TVertex *v1,TVertex *v2);

	bool hasEdge(TVertex *v0,TVertex *v1);

	bool hasFace(TTetra* t1 );

	int hasVertex(TTetra t1 );
	bool hasVertex(TVertex v0 );
	TTetra* getTetraNeighbour(int faceI ,TVertex v0,TVertex v1,TVertex v2, TList<TObject> *tl) ;
	void coincidenteFace(TTetra t1,  TList<TObject> *tl);

	void replaceTriangle(TTriangle oldTr,TTriangle newTr ) ;
	void replaceVertex(TVertex oldV, TVertex newV);
	int isInside(float4 ps );

	bool isInvalid() ;

	TList<TObject*>* getNeighboursByFace(int depth ,TList<TObject*>* nFL= NULL , bool ommitLowerIds = false) ;

	TList<TObject>* getNeighbours(int depth) ;

	void removeVertexRef();

	void updateVertexRef();

	void clearVertexRef();

	void update();

	/*
	property Surface get getSurface set fSurface ;
	property Volume : double read getVolume write fVolume ;
	property Perimeter : double read getPerimeter ;
	property QuadPerimeter : double read getQuadPerimeter ;
	property diedralAngle : double read fDiedralAngle;
	property FaceAngle : double read fFaceAngle;
	property MinEdgeLength : double read getminEdgeLength;
	property MaxEdgeLength : double read getmaxEdgeLength;

	end;
	*/
};

class TElementsPool
{
public : 
	TList<TTetra*>* availableElements;
	TList<TTriangle*>* availableTriangles;
	TList<TVertex*>* availableVertexes;
	
	TTetra* getTetraInstance();
	TTetra* getTetraInstance(TVertex *v0,TVertex *v1,TVertex *v2,TVertex *v3);

	TTriangle* getTriangleInstance(TVertex *v0,TVertex *v1,TVertex *v2);
	TVertex* getVertexInstance(float4 fpos);

	TElementsPool();
	~TElementsPool();
	void releaseInstance(TTetra *t);
	void releaseInstance(TTriangle *tr);
	void releaseInstance(TVertex *v);
	void releaseMemmory();
};

class TMesh
{
public:
	TList<TVertex*>* vertexes ;
	TList<TPolygon*>* fFaces;
	TList<TObject*>* patchList;
	TList<TObject*>* elements;
	float4 fMassCenter;
	float4 scale ;

	TList<TObject*>* selectedElements;
	TList<TObject*>* elementsToAdd;
	TList<TObject*>* elementsToRemove;
	TList<TVertex*>* vertexesToRemove;
	TElementsPool* memPool;

	TMesh(void) 
	{
		patchList = NULL;
		selectedElements = new TList<TObject*>();
		vertexes = new TList<TVertex*>();
		fFaces = new TList<TPolygon*>();
		vertexesToRemove = new TList<TVertex*>();

		scale = Float4(1.0,1.0,1.0);
		memPool = new TElementsPool();

	};
	~TMesh(void)
	{
		delete memPool;
		delete selectedElements;
		delete vertexes;
		delete fFaces;
		delete selectedElements;
		delete vertexesToRemove;

	};

	void addTriangle(TPolygon* tr)
	{
		fFaces->Add(tr);
		tr->id = fFaces->Count();
	};

	void addVertex(TVertex* v)
	{
		vertexes->Add(v);

	};

	float4 size()
	{
		float4 minS, maxS ;
		minS = Float4(10000.0,10000.0,10000.0);
		maxS = Float4(-10000.0,-10000.0,-10000.0);
		for (int j=0 ; j<this->vertexes->Count();j++)
		{
			TVertex* v = this->vertexes->elementAt(j);
			minS.x = (float) Min(v->fPos.x, minS.x);
			minS.y = (float) Min(v->fPos.y, minS.y);
			minS.z = (float) Min(v->fPos.z, minS.z);

			maxS.x = (float) Max(v->fPos.x, maxS.x);
			maxS.y = (float) Max(v->fPos.y, maxS.y);
			maxS.z = (float) Max(v->fPos.z, maxS.z);
		}
		return maxS - minS ;
	}

	virtual void updateRefs() { };

	TVertex* innerBinSearch( int key, int imin, int imax)
	{

		// test if array is empty
		if (imax < imin)
			// set is empty, so return value showing not found
			return NULL;
		else
		{
			// calculate midpoint to cut set in half

			int imid = (imin + imax) / 2;
			if (imid>= vertexes->Count() )
				return NULL;
			TVertex * v = vertexes->elementAt(imid);

			if (v == NULL )
				return NULL;		 

			// three-way comparison
			if (v->id > key)
				// key is in lower subset
				return innerBinSearch( key, imin, imid-1);
			else if (v->id  < key)
				// key is in upper subset
				return innerBinSearch( key, imid+1, imax);
			else
				// key has been found
				return v;
		}
	}

	TVertex* findVertexById(int id)
	{ 				  
		return innerBinSearch(id, 0, vertexes->Count());
	}

	void normalizeMesh(float normSize )
	{
		float4 sz;
		float minSz;
		sz = this->size();
		//sz = getMeshSize(aMesh);
		minSz =(float)Min( (float)Min((float)normSize/sz.x,(float)normSize/sz.y),(float)normSize/sz.z);
		this->setScale(minSz,minSz,minSz);
		// mn := getMeshCenter(aMesh);
		// getMeshMinMAx(aMesh,mn,mx);
		// aMesh.position :=-mn;

	}

	void setScale(float nx,float ny, float nz)
	{
		float4 nw = Float4(nx,ny,nz);
		for (int i=0; i< this->vertexes->Count() ; i++)
		{
			TVertex* v = vertexes->elementAt(i);
			v->fPos.x *= nw.x/scale.x;
			v->fPos.y *= nw.y/scale.y;
			v->fPos.z *= nw.z/scale.z;

		}
		scale = nw; 

	}

	virtual void updateIndexes(int flag){};

	void removeTriangle(TPolygon* tr)
	{
		//       lremovedElements.add(tr);
		fFaces->Extract(tr);
		tr->isdestroyed = true;
	};

	void removeVertex(TVertex* v)
	{
		vertexes->Extract(v);
		v->isdestroyed = true;
	};
};






