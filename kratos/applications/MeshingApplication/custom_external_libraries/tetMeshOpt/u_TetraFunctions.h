#include "u_Types.h"
#include "u_delphiClasses.h"

//---------------------------------------------------------------------------------------------------------
class TVolumeMesh : public TMesh
{
public :    	
	int numVertices ; 
	int numElements ; 
    TVolumeMesh();	
	void getSurfaceTriangles();
	void updateRefs();
	void validate(bool showMessages);	
	void updateIndexes(int flag = 1);	
	void removeFreeVertexes();
};

//---------------------------------------------------------------------------------------------------------
void igetSurfaceTriangles(TList<TObject*> *elements , TList<TObject*> *res, TList<TObject*> *surfVertexes  );
bool innerGetElementsSurface(TList<TObject*>* elements,TList<TObject*>* surfaceT,TList<TObject*>* vertexes);
bool swapTetra(TVertex* v0,TVertex*  v1,TVertex*  v2,TVertex*  v3,TVertex* v4);
bool swapVolumeMesh(TVolumeMesh* aMesh);

