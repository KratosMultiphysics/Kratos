#include <math.h>
#include "Math3D.h"
#include "u_Types.h"

#define ACCEPTANCE_TOLERANCE 1.0
#define ELEMENTS_TO_FORCE_UPDATE 100000

const float cotas[] ={0.0f,0.001f, 0.05f, 1.0f, 5000000000.0f};

class TElementsCluster
{
public :
	TMesh* originalMesh;
	double goodMinQuality ,   avgEdgeLength ,    minDAngle ;
	TList<TObject*>* vertexes ;
	TList<TObject*>* surfaceT;
	TList<TObject*>* elements ;
	TList<TObject*>* copyEL ;
	//TList<TObject*>* origElements ;
	TVertex* vC , *inspVertex;

	TList<TObject*>* tempTetraList;
	TList<TObject*>* goodTetraList;
	TList<TObject*>*  elements2;
	TList<TObject*>* tempL;
	TList<TObject*>* newElements ;
	TVertexesEvaluator fc;
	int numChanges ;
	double minQuality ;
	///----- Variables to store parameters
	TList<TObject*>* inspectedElements;
	bool doRemoveElements ,    doTestSwap ,   testCenter;
	float perturbCenter ;
	bool checkMaxLength ;
	float minExpectedquality ;

	// Constructor	   
	void computeMetrics();

	int evaluateSet();
	void  genElements(bool vertexIsNew = true);
	bool updateMesh(bool checkIfInvalid );
	bool generateSubMesh(double minExpectedQuality =500000, double minExpectedAngle = 50000);
	TElementsCluster(TMesh* aMesh, TVertexesEvaluator functor  );
	double getMinQuality();
	~TElementsCluster();

} ;



void evaluateClusterByNode(TMesh *aMesh , double minExpectedQuality,TVertexesEvaluator fc) ;
void evaluateClusterByFace(TMesh *aMesh , double minExpectedQuality);
void evaluateClusterByFace(TMesh *aMesh , double minExpectedQuality,TVertexesEvaluator fc);

void evaluateClusterByEdge(TMesh *aMesh , double minExpectedQuality,TVertexesEvaluator fc);

bool improvedCluster(TList<TObject*>* c1,
	TList<TObject*>* cRef ,
	TElementsCluster* cl,
	bool maxEdgeLengthConstrain );

bool improvedTetra(TVertex *v0,TVertex *v1,TVertex *v2,TVertex *v3,
	                TList<TObject*>* cRef ,
	               TElementsCluster* cl,
	                bool maxEdgeLengthConstrain );
double testTetraSplit4(TTetra *t, TVertex* v, TList<TObject*>* lRes , TVertexesEvaluator qualityFunction) ;
int vertexTetraReInsertion(TMesh *am ,TList<TVertex*>* vertexesList ) ;

/*


void innerSwapSurface(TList<TObject*>* surfaceT , TList<TObject*>* elementsList)
{

TVertex* v0, *v1, *v2, *v3,*v4,*v5,*sV0,*sV1,*sV2,*sV3;
int i,j,k,l,cnt  ;
TTetra*t0,*t1;
bool isSurface ;
double vl1,vl2;
float4 vc;

//Visito todos los triangulos de la superficie
for (i = 0 ; i<(surfaceT->Count() / 3) ; i++)
for (j = i+1 ; j<(surfaceT->Count() /3) ; j++)
{
//Triangle 0
v0 =(TVertex*)(surfaceT->elementAt(3*i));
v1 =(TVertex*)(surfaceT->elementAt(3*i+1));
v2 =(TVertex*)(surfaceT->elementAt(3*i+2));
// verifico si es un triangulo de Superficie
isSurface = false;
for (l = 0 ; l< elementsList->Count() ;l++)
{
t0 = (TTetra*)(elementsList->elementAt(l));
for (k = 0 ; k<=3; k++)
{
if (t0->getTetraNeighbour(k,v0,v1,v2,NULL)==NULL)
isSurface = true;
}
}
if (!isSurface) continue ;
//Triangle 1
v3 =(TVertex*)(surfaceT->elementAt(3*j+0));
v4 =(TVertex*)(surfaceT->elementAt(3*j+1));
v5 =(TVertex*)(surfaceT->elementAt(3*j+2));
// Triangulo replicado
if (((v0= v3) || (v0 = v4) || (v0 = v5)) &&
((v1= v3) || (v1 = v4) || (v1 = v5)) &&
((v2= v3) || (v2 = v4) || (v2 = v5))) continue;

// verifico si es un triangulo de Superficie
isSurface = false;
for (l = 0 ; l<elementsList->Count();l++)
{
t1 = (TTetra*)(elementsList->elementAt(l));
for (k = 0 ; k<4; k++)
{
if (t1->getTetraNeighbour(k,v3,v4,v5,NULL)==NULL)
isSurface = true;
}
}
if (!isSurface) continue ;

sV0 = NULL ; sV1 = NULL;
// 1st Rot : v3 & V4
if ((v0==v4) && (v1 ==v3)) { sV0 = v2 ; sV3 = v5; sV1 = v0; sV2 = v1; }
if ((v1==v4) && (v2 ==v3)) { sV0 = v0 ; sV3 = v5; sV1 = v1; sV2 = v2; }
if ((v2==v4) && (v0 ==v3)) { sV0 = v1 ; sV3 = v5; sV1 = v2; sV2 = v0; }
// 1st Rot : v4 & V5
if ((v0==v5) && (v1 ==v4)) { sV0 = v2 ; sV3 = v3; sV1 = v0; sV2 = v1; }
if ((v1==v5) && (v2 ==v4)) { sV0 = v0 ; sV3 = v3; sV1 = v1; sV2 = v2; }
if ((v2==v5) && (v0 ==v4)) { sV0 = v1 ; sV3 = v3; sV1 = v2; sV2 = v0; }
// 1st Rot : v5 & V3
if ((v0==v3) && (v1 ==v5)) { sV0 = v2 ; sV3 = v4; sV1 = v0; sV2 = v1; }
if ((v1==v3) && (v2 ==v5)) { sV0 = v0 ; sV3 = v4; sV1 = v1; sV2 = v2; }
if ((v2==v3) && (v0 ==v5)) { sV0 = v1 ; sV3 = v4; sV1 = v2; sV2 = v0; }
// Sino comparte arista
if ((sV0 == NULL) && (sV1== NULL )) continue;
//Sino mejora la calidad
if (Min(calidadxArea(v0,v1,v2),calidadxArea(v3,v4,v5)) >
Min(calidadxArea(sV0,sV1,sV3),calidadxArea(sV0,sV3,sV2))) continue;
// Medir la perdida de volumen

vc = t0->oppositeVertex(v0,v1,v2)->pos();
vl1 =( tetraVolume(v0->pos(),v1->pos(),v2->pos(),vc)) +( tetraVolume(v3->pos(),v4->pos(),v5->pos(),vc));
vl2 =( tetraVolume(sV0->pos(),sV1->pos(),sV3->pos(),vc)) +( tetraVolume(sV0->pos(),sV3->pos(),sV2->pos(),vc));

if ((vl2!=0) && (vl1/vl2<0.95))  continue;
if ((vl1!=0) && (vl2/vl1<0.95))  continue;

// reemplazo los elementos
surfaceT->setElementAt(3*i, sV0);
surfaceT->setElementAt(3*i+1, sV1);
surfaceT->setElementAt(3*i+2, sV3);

surfaceT->setElementAt(3*j, sV0);
surfaceT->setElementAt(3*j+1, sV3);
surfaceT->setElementAt(3*j+2, sV2);
return ;
}

}
*/

