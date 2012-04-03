// PrSpain.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <omp.h>
#include "u_Types.h"
#include "u_MeshLoaders.h"
#include "u_TetraFunctions.h"
#include "u_elementCluster.h"
#include "u_ParallelFunctions.h"
#include "u_ProcessTime.h"
#include "u_ShowMetrics.h"
#include "u_qualityMetrics.h"

using namespace std;

float fu(float x)
{
    return (x*x) * ( (1.0-x)*(1.0-x) );
}

float df(float x)
{
    return 2.0*x*( (1.0-x)*(1.0-x) ) - 2.0*(x*x)*(1.0-x);
}

float d2f(float x)
{
    return 2.0*( (1.0-x)*(1.0-x) ) - 8.0*x*(1.0-x) + 2.0*(x*x);
}

float d3f(float x)
{
    return -12.0*(1.0-x) + 12.0*x;
}   


void MoveNodes(TMesh* amesh, float Dt,int substeps, float h)
{
  float dt_small = Dt / float(substeps);
    
  for (int i=0; i<amesh->vertexes->Count(); i++)
  {
    TVertex *v = amesh->vertexes->elementAt(i);
	if (v->fixed) continue;
	 
      for (int step = 0; step<substeps; step++)
	  {
        float x = v->fPos.x;
        float y = v->fPos.y;
        float ht = 1.0;
        float vx = 100.0 * ht * fu(x) * df(y);
        float vy = -100.0 * ht * fu(y) * df(x);
        v->fPos.x = x + vx*dt_small;
        v->fPos.y = y + vy*dt_small;
	  }
  }
}

int main (int argc, char* argv[])
{

	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/cavityMesh_0_00_600k.vwm"));   
	//TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/f1_tetgen_mesh.vwm"));   
	//TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/kratosProof.vwm"));   
	//TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/cavityMesh_0_00.vwm"));   
	omp_set_num_threads(1);
	//m->normalizeMesh(1000);
	m->updateIndexes(GENERATE_SURFACE);
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	

	for (int i=0 ; i<m->vertexes->Count();i++)
	{
		m->vertexes->elementAt(i)->fixed = false;
	}

	for (int i=0 ; i<m->fFaces->Count() ; i++)
	{
		TTriangle *tr = (TTriangle*)( m->fFaces->elementAt(i));
		tr->vertexes[0]->fixed = 1;
		tr->vertexes[1]->fixed = 1;
		tr->vertexes[2]->fixed = 1;		
	}

	
	//swapVolumeMesh((TVolumeMesh*)(m));
	//startProcess("optimize by node");
	cout <<".............................................." <<"\n"; 
	cout <<"...Preparing optimization" <<"\n"; 
	TetQuality *qt = new TetQuality(m);
	/*qt->refresh();   qt->print();
	evaluateClusterByFace(m,0.5,vrelaxQuality);
	qt->refresh();   qt->print();

	m->validate(true);
	*/
	//cout <<"...Optimizing by Node" <<"\n"; 
	//evaluateClusterByNode( (TVolumeMesh*)(m),0.5,vrelaxQuality);
	cout <<".............................................." <<"\n"; 
	cout <<"...Initial Mesh quality" <<"\n"; 
	qt->refresh();   qt->print();
	cout <<".............................................." <<"\n"; 
	
	startProcess("Mesh evaluation");
	stopTimers();
	//ParallelEvaluateClusterByNode(m,vrelaxQuality);   
	for (int i=0 ;i<40 ; i++)
	{
		int substeps = 10;
		float dt = 0.01;
		MoveNodes(m,dt,substeps,1);
		qt->refresh();   qt->print();

		for (int j=0; j<3;j++)
		{
			cout <<"--------- ITERATION "<< j <<"--------------------"<<"\n"; 
			
			cout <<"...Optimizing by Node" <<"\n"; 
			//evaluateClusterByNode( (TVolumeMesh*)(m),5000000,vrelaxQuality);
			cout <<"...Optimizing by Edge" <<"\n"; 
			//evaluateClusterByEdge(m,50000000000, vrelaxQuality);
			cout <<"...Optimizing by Face" <<"\n"; 
			//evaluateClusterByFace(m,50000000, vrelaxQuality);
			cout <<"...Parallel Optimizing by Node" <<"\n"; 
			ParallelEvaluateClusterByNode(m,vrelaxQuality);
			cout <<"...Parallel Optimizing by Face" <<"\n"; 
			//ParallelEvaluateClusterByFace(m,vrelaxQuality);
			cout <<"...Parallel Optimizing by Edge" <<"\n"; 
			ParallelEvaluateClusterByEdge(m,vrelaxQuality);
						
		}
		
		qt->refresh();   qt->print();
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
		m->validate(true);
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

		if (qt->nonPositive>0) 
		{
			std :: cout<< " Reached iteration " << i <<"\n";
			break;
		}
		/*
		
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
		qt->refresh();   qt->print();
		m->validate(true);
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
		int vToR =m->vertexesToRemove->Count();
		int ri = vertexTetraReInsertion(m , m->vertexesToRemove);
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
		std :: cout<< " Reinsert vertexes " << ri << " of " <<  vToR <<"\n";
		qt->refresh();   qt->print();
		m->validate(true);
		*/
	}
	startTimers();
	endProcess("Mesh evaluation");
	    showProcessTime();
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
		qt->refresh();   qt->print();
		m->validate(true);
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

	//m->removeFreeVertexes();
	delete m;

	char cw;
	std::cout <<"...waiting to finish" <<"\n";         
	cin >> cw;
	
	return 1;
	/*

	cout <<"...Optimizing by Edge" <<"\n"; 
	evaluateClusterByEdge( (TVolumeMesh*)(m),50000,vrelaxQuality);
	m->updateIndexes(0);
	qt->refresh();   qt->print();
	m->validate(true);

	cout <<"...Optimizing by Face" <<"\n"; 
	evaluateClusterByFace( (TVolumeMesh*)(m),50000,vrelaxQuality);
	m->updateIndexes(0);
	qt->refresh();   qt->print();
	m->validate(true);

	cout <<"...Optimizing by Edge" <<"\n"; 
	evaluateClusterByEdge( (TVolumeMesh*)(m),50000,vrelaxQuality);
	m->updateIndexes(0);
	qt->refresh();   qt->print();


	m->validate(true);
	*/
	//endProcess("optimize by node");
	showProcessTime();
	char c;
	std::cout <<"...Input 'S' to save" <<"\n";         
	cin >> c;
	if ((c == 's' ) || (c == 'S'))
	{
		TMeshLoader* ml2 = new TVMWLoader();
		std::string s("");
	    s = "D:/out_MeshFromKratos" + intToString(0)+".vwm";
		ml2->save(s , m);
		delete ml2;
	}


	////showProcessTime();

	return 0;
} 
