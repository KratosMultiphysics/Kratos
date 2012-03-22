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

struct  /* __declspec(novtable) */ object
{
   virtual ~object() = 0 { /* empty */ }
};

class TValuedObject : public TObject
{
public:
  double	fmetrica , sortvalue ,calidad;
  bool marked, changed ,locked,isdestroyed;
  int id , groupID,visited,flag;
  void* linkSet ;
  int _color ;
  float4 rColor;
  object* userData;

  TValuedObject(void) {};
  ~TValuedObject(void){};

  float4 getCenter() { float4 f;
	                   return f; }
};

class TNeighboured : public TValuedObject{};
class TVertex;
class TTriangle;

class TPolygon: public TNeighboured
{
private : 
	

public :
	TVertex* vertexes[3];
	float4 normal;
};


class TVertex: public TNeighboured
{
public:   
	         float4 fPos ;
             int fixed ;
			 TList<TObject*>* neighTr;
			 TList<TObject*>* elementsList;
			 TList<TObject*>* neighV;
             bool isSurface;
			 

			 TVertex(float x,float y, float z)
			 { 
				  fPos.x= x; fPos.y = y; fPos.z = z; 
				  neighTr = new  TList<TObject*>();
			      elementsList = new TList<TObject*>();			 
			 }

			 TVertex(float4 v)
			 { 
				  fPos.x= v.x; fPos.y = v.y; fPos.z = v.z; 
				  neighTr = new  TList<TObject*>();
			      elementsList = new TList<TObject*>();			 
			 }

			 

             
             float4 normal ;
			 void calcNormal(){};
			 float4 pos(){ return fPos; }
			 TList<TObject*> *getElemNeighbours(TList<TObject*> *toL = NULL);			 
			 TList<TObject*>* getvNeigh(int depth =1 , int mode = 0, TList<TObject*>* toL = NULL);
};

class TVertex3D : public TVertex
{

};

typedef double(*TVertexesEvaluator)(TVertex*[]);

class TTriangle: public TPolygon
{
public :    
	TTriangle(TVertex* v0,TVertex* v1,TVertex* v2)
	{
		 vertexes[0] = v0;
		 vertexes[1] = v1;
		 vertexes[2] = v2;
	}
	void calcNormal()
	{
       this->normal =Normal(vertexes[0]->pos(),vertexes[1]->pos(),vertexes[2]->pos());
	}

	void calcEdges()
	{
	}

	bool verticesEncomun(TTriangle* otherT,TVertex* &v1,TVertex* &v2)
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

	TMesh(void) 
	{
	  patchList = NULL;
	  selectedElements = new TList<TObject*>();
	  vertexes = new TList<TVertex*>();
	  fFaces = new TList<TPolygon*>();
	  scale = Float4(1.0,1.0,1.0);
	
	};
   ~TMesh(void){};

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

class TElement : public TNeighboured
{
public:      
	TVertex* vertexes[4];
	float4 normals[4];
	bool initialized ;
	TElement() {}
	TElement(TVertex v0,TVertex v1,TVertex v2,TVertex v3){};
	virtual BoundBox CalcBound() 
	{ 
		BoundBox b = calcBound(vertexes[0]->fPos, vertexes[1]->fPos); 
		return  b;
	};
};


class TTetra : public TElement
{
public : 
            double fVolume,fDiedralAngle, fFaceAngle,fSurface ;
			double getVolume() 
			{ 
				 return tetraVolume(vertexes[0]->fPos, vertexes[1]->fPos, vertexes[2]->fPos , vertexes[3]->fPos);
			}

            double getPerimeter() { return 0;}
            double getQuadPerimeter() {return 0;}
            double getminEdgeLength() {return 0;}
            double getmaxEdgeLength() {return 0;}

            int ts[6];

            TTriangle* triangles[4];
			TList<TObject*>* NeighBourByFace;
			TList<TValuedObject>* Neighbours;
			
            int materialID;
            float4 velocity;

          
			double getSurface(){return 0;}
			void copyFrom(TValuedObject e ){}
			TTetra(TObject* owner) { this->isdestroyed = false;  return ; }

			void Create(TObject* owner, TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3, bool publish = true)
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

			TTetra(TObject* owner, TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3, bool publish = true) 
			{
				this->Create(owner,v0,v1,v2,v3,publish);
				if (publish)
					updateVertexRef();
			} 
			
			~TTetra() 
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
			TTetra* getTetraNeighbour(int faceI, TVertex* v0,TVertex* v1,TVertex* v2, TList<TObject*>*  tl)
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
			TList<TObject*>* getSurfaceTriangle(bool bmpMode , TList<TObject*>* res, TList<TObject*>* ts)
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
			bool hasFace(TVertex* v0,TVertex*  v1,TVertex* v2) 			
			{
            return ((v0 == vertexes[0]) || (v0 == vertexes[1]) || (v0 == vertexes[2]) || (v0 == vertexes[3])) &&
                   ((v1 == vertexes[0]) || (v1 == vertexes[1]) || (v1 == vertexes[2]) || (v1 == vertexes[3]))&&
                   ((v2 == vertexes[0]) || (v2 == vertexes[1]) || (v2 == vertexes[2]) || (v2 == vertexes[3]));
			}

			TVertex* oppositeVertex(TVertex *v0,TVertex *v1,TVertex *v2)
			{  
				  if ((vertexes[0]!=v0) && (vertexes[0]!=v1) && (vertexes[0]!=v2))  return  vertexes[0];
					else if ((vertexes[1]!=v0) && (vertexes[1]!=v1) && (vertexes[1]!=v2))  return vertexes[1];
					else if ((vertexes[2]!=v0) && (vertexes[2]!=v1) && (vertexes[2]!=v2))  return  vertexes[2];
					else if ((vertexes[3]!=v0) && (vertexes[3]!=v1) && (vertexes[3]!=v2))  return  vertexes[3];
					else return NULL;
			}
			bool hasEdge(TVertex *v0,TVertex *v1)
			{   
				return ((v0 == vertexes[0]) || (v0 == vertexes[1]) || (v0 == vertexes[2]) || (v0 == vertexes[3])) &&
                          ((v1 == vertexes[0]) || (v1 == vertexes[1]) || (v1 == vertexes[2]) || (v1 == vertexes[3])) ;				
			}
			bool hasFace(TTetra* t1 )
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
	
			int hasVertex(TTetra t1 ){ return 0;} 
			bool hasVertex(TVertex v0 ) { return 0;}
			TTetra* getTetraNeighbour(int faceI ,TVertex v0,TVertex v1,TVertex v2, TList<TObject> *tl) { return NULL;}
            void coincidenteFace(TTetra t1,  TList<TObject> *tl);
           
			void replaceTriangle(TTriangle oldTr,TTriangle newTr ) {}
			void replaceVertex(TVertex oldV, TVertex newV){}
			int isInside(float4 ps ){ return 0;}
            
			bool isInvalid() {return 0;}

			TList<TObject*>* getNeighboursByFace(int depth ,TList<TObject*>* nFL= NULL ) 
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

			TList<TObject>* getNeighbours(int depth) {return NULL;}
            
			void removeVertexRef()
			{
				clearVertexRef();
				for (int i=0 ; i<4 ; i++)
				 if (vertexes[i]->elementsList != NULL )
				    vertexes[i]->elementsList->Pack();
			}
            void updateVertexRef()
			{
			   for (int i = 0 ; i<4 ; i++)
               {				   
				   if (vertexes[i]->elementsList == NULL)     
                      vertexes[i]->elementsList = new TList<TObject*>();
                   if ( vertexes[i]->elementsList->indexOf(this)<0 )
                       vertexes[i]->elementsList->Add(this);
				   
			   }
			}
			
			
			void clearVertexRef()
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





