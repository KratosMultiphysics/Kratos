#include "u_delphiClasses.h"
#include "u_Types.h"

typedef void(*TStForLoopElement)(int , int , TObject*);
//-----------------------------------------
// This class encapsulates how the parallel for is invoked : TBB , OpenMP
class TParallelIterator
{
	TStForLoopElement innercall;
	TList<TObject*>* elements;
public :

	TParallelIterator()
	{
	}

	TParallelIterator(TStForLoopElement ic,  TList<TObject*>* el) : innercall(ic), elements(el)
	{  }


	void forloop(int from, int to, TList<TObject*>* elements,TStForLoopElement call)
	{
        #pragma omp parallel for
		for (int i=from ; i<=to ; i++)
			call(i,0,elements->elementAt(i));
		//parallel_for(blocked_range<size_t>(from,to), this );
		#pragma omp barrier
	}
};



void assignVertexesAvoidingVisited(TList<TObject*> *vs, TList<TObject*> *vRes ,int iter, int maxAssignment);
void fastGetSurfaceTriangles(TMesh* aModel);
void ParallelEvaluateClusterByNode(TMesh *aMesh , TVertexesEvaluator fc);
void ParallelEvaluateClusterByEdge(TMesh *aMesh , TVertexesEvaluator fc);