// PrSpain.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "u_Types.h"
#include "u_MeshLoaders.h"
#include "u_TetraFunctions.h"
#include "u_elementCluster.h"
#include "u_ParallelFunctions.h"
#include "u_ProcessTime.h"
#include "u_ShowMetrics.h"
#include "u_qualityMetrics.h"

using namespace std;

int main (int argc, char* argv[])
{

	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/kratosProof.vwm"));   
	//TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/out_MeshFromKratos0.vwm"));   
	
	//m->normalizeMesh(1000);
	m->updateIndexes(GENERATE_SURFACE);
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

	

	for (int i=0 ; i<m->vertexes->Count() ; i++)
	{
		int index = i+1;
		if (m->findVertexById(index) == NULL) 
			cout << "Error ID "<< index <<"\n";
	}

	
	//swapVolumeMesh((TVolumeMesh*)(m));
	//startProcess("optimize by node");
	cout <<"...Optimizing by Face" <<"\n"; 
	TetQuality *qt = new TetQuality(m);
	/*qt->refresh();   qt->print();
	evaluateClusterByFace(m,0.5,vrelaxQuality);
	qt->refresh();   qt->print();

	m->validate(true);
	*/
	//cout <<"...Optimizing by Node" <<"\n"; 
	//evaluateClusterByNode( (TVolumeMesh*)(m),0.5,vrelaxQuality);
	qt->refresh();   qt->print();
	cout <<"...Parallel Optimizing by Node" <<"\n"; 
	startProcess("Parallel evaluation");

	//ParallelEvaluateClusterByNode(m,vrelaxQuality);   
	for (int i=0 ;i<1 ; i++)
	{
		evaluateClusterByNode( (TVolumeMesh*)(m),5000000,vrelaxQuality);
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
		endProcess("Parallel evaluation");
		qt->refresh();   qt->print();
		m->validate(true);
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

		evaluateClusterByFace( (TVolumeMesh*)(m),5000000,vrelaxQuality);
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
		endProcess("Parallel evaluation");
		qt->refresh();   qt->print();
		m->validate(true);
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

    	evaluateClusterByEdge( (TVolumeMesh*)(m),5000000,vrelaxQuality);
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS) ;
		endProcess("Parallel evaluation");
		qt->refresh();   qt->print();
		m->validate(true);
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

	}

	m->removeFreeVertexes();
	showProcessTime();

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
