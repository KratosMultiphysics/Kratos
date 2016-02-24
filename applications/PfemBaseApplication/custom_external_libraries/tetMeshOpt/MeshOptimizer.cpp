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

#include "u_TetGenInterface.h"
#include "cl_VolumeMeshSmooth.h"
#include "u_MeshSmooth.h"

#include "cl_utils.h"

using namespace std;

#define APPLY_BY_NODE 2
#define APPLY_BY_EDGE 3
#define APPLY_BY_FACE 5
#define SMOOTH_MESH 7

ofstream myfile;

float fu(float x)
{
    return (x*x) * ( (1.0f-x)*(1.0f-x) );
}

float df(float x)
{
    return 2.0f*x*( (1.0f-x)*(1.0f-x) ) - 2.0f*(x*x)*(1.0f-x);
}

float d2f(float x)
{
    return 2.0f*( (1.0f-x)*(1.0f-x) ) - 8.0f*x*(1.0f-x) + 2.0f*(x*x);
}

float d3f(float x)
{
    return -12.0f*(1.0f-x) + 12.0f*x;
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

int  optimizeMesh(TVolumeMesh *m, int numIter, double minExpQ ,int flag)
{
	int evaluations = 0;
	for (int j=0; j<numIter;j++)
		{
			cout <<"--------- ITERATION "<< j <<"--------------------"<<"\n"; 			
			cout <<"...Parallel Optimizing by Node" <<"\n"; 
			if (flag % APPLY_BY_NODE == 0) 
			{
			 evaluations +=ParallelEvaluateClusterByNode(m,vrelaxQuality,minExpQ);
			 
			}
			
			cout <<"...Parallel Optimizing by Edge" <<"\n"; 
			if (flag % APPLY_BY_EDGE == 0) 
			{
			  evaluations += ParallelEvaluateClusterByEdge(m,vrelaxQuality,minExpQ);
			 
			}
			
			cout <<"...Parallel Optimizing by Face" <<"\n"; 
			if (flag % APPLY_BY_FACE == 0) 
			{
			  evaluations += ParallelEvaluateClusterByFace(m,vrelaxQuality,minExpQ);
			  
			}
			if (flag % SMOOTH_MESH == 0) 
			{
				cout <<"...Smoothing" <<"\n"; 
				//laplacianMeshSmooth(m);
				//iterativeMeshSmooth(m,25,25,Min(20.0,minExpQ));
				ParallelSmoothMesh(m,vrelaxQuality,minExpQ);
			}
		}
	return evaluations;
}

/// Optimize Mesh in Parallel
void benchmark0(std::string meshFileName,int numIter, double minExpQ ,int flag)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 0 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	TetQuality *qt = new TetQuality(m);
	
	qt->refresh();   qt->print();
	double origAve = qt->DieAveMin;
	double minD = qt->fDieMin;
	cout <<".............................................." <<"\n"; 	
	startProcess("Mesh evaluation in parallel :  benchmark0");
	stopTimers();
	double evaluations = 0.0;
	int numEvaluations = 0;
	numEvaluations = optimizeMesh(m,numIter,minExpQ ,flag);

	startTimers();			
	double t= endProcess("Mesh evaluation in parallel :  benchmark0");		
	qt->refresh();   qt->print();
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	m->validate(true);
	std :: cout<< " Required time " << t<<"\n";
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	/// File NumIter Smooth Time WorstQuality 
	myfile << meshFileName<<"\t" << evaluations<<" "<<numEvaluations<<" "<<m->elements->Count()<<"\t"<<minD <<"\t"<<origAve <<"\t"<<minExpQ <<"\t" <<numIter<<"\t" <<flag <<"\t" << t<< "\t"<< qt->fMinQuality <<"\t" <<qt->fDieMin <<"\t" <<qt->DieAveMin<<"\n";

	std::string s("");
	if (flag % SMOOTH_MESH == 0 )
	  s = "benchmark0_" + meshFileName+"_smoothed.vwm" ;
	else
	  s = "benchmark0_" + meshFileName+".vwm";
	ml->save(s , m);

	//exportMetrics(meshFileName+"_opt.cvs",m,t);
  
  delete qt;
  delete m;
  delete ml;

}


/// Optimize Mesh in Parallel
// Se reordenan las operaciones
void benchmark0Bis(std::string meshFileName,int numIter, double minExpQ ,int flag)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 0 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	
	
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	TetQuality *qt = new TetQuality(m);
	
	qt->refresh();   qt->print();
	double minD = qt->fDieMin;
	cout <<".............................................." <<"\n"; 	
	startProcess("Mesh evaluation in parallel :  benchmark0");
	stopTimers();

	for (int j=0; j<numIter;j++)
	{
			cout <<"--------- ITERATION "<< j <<"--------------------"<<"\n"; 			
			cout <<"...Parallel Optimizing by Node" <<"\n"; 
			if (flag % APPLY_BY_NODE == 0) 
			 ParallelEvaluateClusterByNode(m,vrelaxQuality,minExpQ);
	}	
	for (int j=0; j<numIter;j++)
	{
			cout <<"...Parallel Optimizing by Edge" <<"\n"; 
			if (flag % APPLY_BY_EDGE == 0) 
			  ParallelEvaluateClusterByEdge(m,vrelaxQuality,minExpQ);
	}		
	for (int j=0; j<numIter;j++)
	{
			cout <<"...Parallel Optimizing by Face" <<"\n"; 
			if (flag % APPLY_BY_FACE == 0) 
			  ParallelEvaluateClusterByFace(m,vrelaxQuality,minExpQ);
	}
	for (int j=0; j<numIter;j++)
	{
			if (flag % SMOOTH_MESH == 0) 
			{
				cout <<"...Smoothing" <<"\n"; 
				//laplacianMeshSmooth(m);
				iterativeMeshSmooth(m,25,25,Min(20.0,minD*3.0));
			}
		}
	startTimers();			
	double t= endProcess("Mesh evaluation in parallel :  benchmark0");		
	qt->refresh();   qt->print();
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	m->validate(true);
	std :: cout<< " Required time " << t<<"\n";
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	/// File NumIter Smooth Time WorstQuality 
	myfile << meshFileName<<"\t"<<minExpQ <<"\t" <<numIter<<"\t" <<flag <<"\t" << t<< "\t"<< qt->fMinQuality <<"\t" <<qt->fDieMin <<"\n";

	std::string s("");
	if (flag % SMOOTH_MESH == 0 )
	  s = "benchmark0_" + meshFileName+"_smoothed.vwm" ;
	else
	  s = "benchmark0_" + meshFileName+".vwm";
	//ml->save(s , m);

	//exportMetrics(meshFileName+"_opt.cvs",m,t);
  
  delete qt;
  delete m;
  delete ml;

}
/// Optimize Mesh - Serial
void benchmark1(std::string meshFileName , bool smooth)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 0 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	TetQuality *qt = new TetQuality(m);

	cout <<".............................................." <<"\n"; 	
	startProcess("Mesh evaluation serial :  benchmark1");
	//stopTimers();

	for (int j=0; j<3;j++)
		{
			cout <<"--------- ITERATION "<< j <<"--------------------"<<"\n"; 
			cout <<"...Optimizing by Node" <<"\n"; 
			//evaluateClusterByNode( (TVolumeMesh*)(m),5000000,vrelaxQuality);
			cout <<"...Optimizing by Edge" <<"\n"; 
			evaluateClusterByEdge(m,50000000000, vrelaxQuality);
			cout <<"...Optimizing by Face" <<"\n"; 
			evaluateClusterByFace(m,50000000, vrelaxQuality);	
			if (smooth)
			{
				cout <<"...Smoothing" <<"\n"; 
				laplacianMeshSmooth(m);
				//iterativeMeshSmooth(m,10,10);
			}
		}
	//startTimers();
	double t=endProcess("Mesh evaluation serial :  benchmark1");			
	qt->refresh();   qt->print();
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	m->validate(true);
	std :: cout<< " Required time " << t<<"\n";
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

  endProcess("Mesh evaluation F1 :  benchmark0");
  delete qt;
  delete m;
  delete ml;
}

/// Cavity movement
void benchmark4(std::string meshFileName, float dt, int maxIter)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 4 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	
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
	
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	TetQuality *qt = new TetQuality(m);

	cout <<".............................................." <<"\n"; 	
	startProcess("Mesh evaluation Cavity :  benchmark2");

	double elapsedTime = 0;
	int iterations = 0 ;

	ofstream myfile;
	myfile.open ("benchmark4.txt");
	
	for (int i=0 ;i<maxIter ; i++)
	{
		int substeps = 10;
		MoveNodes(m,dt,substeps,1);		
		
		qt->refresh(); 
		std :: cout<<i<< " Min diedral angle " << qt->fDieMin <<" Num negative elements "<<qt->nonPositive<<"\n";

		myfile<<i<<"\t"<<qt->fFaceMin<<"\t"<<qt->fDieMax<<"\t"<<qt->fDieMin<<"\t"<<qt->nonPositive<<"\n";
		if (i%10 == 0)
			ml->save(meshFileName+"it"+intToStr(i)+".vwm",m);
	}
	
	startTimers();			
	endProcess("Mesh evaluation Cavity :  benchmark2");	

	myfile.close();
  
  delete qt;
  delete m;
  delete ml;

}

/// Optimize Cavity
void benchmark2(std::string meshFileName, float dt, int maxIter, int numMejoras)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 2 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		

	TVolumeMesh* mOrig = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	mOrig->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		

	for (int i=0 ; i<m->vertexes->Count();i++)
	{
		m->vertexes->elementAt(i)->fixed = false;
		mOrig->vertexes->elementAt(i)->fixed = false;
	}

	for (int i=0 ; i<m->fFaces->Count() ; i++)
	{
		TTriangle *tr = (TTriangle*)( m->fFaces->elementAt(i));
		tr->vertexes[0]->fixed = 1;
		tr->vertexes[1]->fixed = 1;
		tr->vertexes[2]->fixed = 1;		
	}

	for (int i=0 ; i<mOrig->fFaces->Count() ; i++)
	{
		TTriangle *tr = (TTriangle*)( mOrig->fFaces->elementAt(i));
		tr->vertexes[0]->fixed = 1;
		tr->vertexes[1]->fixed = 1;
		tr->vertexes[2]->fixed = 1;		
	}
	
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	TetQuality *qt = new TetQuality(m);

	cout <<".............................................." <<"\n"; 	
	startProcess("Mesh evaluation Cavity :  benchmark2");

	double elapsedTime = 0;
	int iterations = 0 ;
	ofstream myfile;
	myfile.open ("benchmark2.txt");

	for (int i=0 ;i<maxIter ; i++)
	{
		int substeps = 10;
		MoveNodes(m,dt,substeps,1);	
		MoveNodes(mOrig,dt,substeps,1);	
		
		startProcess("Mesh Optimization");
		stopTimers();

		optimizeMesh(m,numMejoras,220, APPLY_BY_FACE * APPLY_BY_EDGE);

		startTimers();			
		double et = endProcess("Mesh Optimization");
		elapsedTime += et;	
		iterations ++;

		qt->refresh(); 
		std :: cout<<i<< " Min diedral angle " << qt->fDieMin <<" Num negative elements "<<qt->nonPositive<< " Time "<< et<<"\n";
		/*if ((numMejoras>0)&&(qt->nonPositive>0) )
		{
			std :: cout<< " Reached iteration " << i <<"\n";
			break;
		}
		*/
		myfile<<i<<"\t"<<et<<"\t"<<qt->fFaceMin<<"\t"<<qt->fDieMax<<"\t"<<qt->fDieMin<<"\t"<<qt->nonPositive<<"\n";		
		
		ml->save(meshFileName+"_opt_it"+intToStr(i)+".vwm",m);
		ml->save(meshFileName+"_orig_it"+intToStr(i)+".vwm",mOrig);

	}

	myfile.close();
	
	startTimers();			
	endProcess("Mesh evaluation Cavity :  benchmark2");	
	showProcessTime();
	std :: cout<< " Mean time required " << elapsedTime/iterations <<"\n";
	
	qt->refresh();   qt->print();
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	m->validate(true);
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

  
  delete qt;
  delete m;
  delete ml;

}

/// Optimize Cavity
void benchmark3(std::string meshFileName, float dt, int maxIter, int numMejoras)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 3 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		

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
	
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	TetQuality *qt = new TetQuality(m);

	cout <<".............................................." <<"\n"; 	
	startProcess("Mesh evaluation Cavity :  benchmark3");

	double elapsedTime = 0;
	int iterations = 0 ;
	
	for (int i=0 ;i<maxIter ; i++)
	{
		int substeps = 10;
		MoveNodes(m,dt,substeps,1);		
		
		startProcess("Mesh Optimization");
		stopTimers();
	
		for (int j=0; j<numMejoras;j++)
		{			
			evaluateClusterByEdge(m,5000000,vrelaxQuality);
			evaluateClusterByFace(m,5000000,vrelaxQuality);
		}

		startTimers();			
		double et = endProcess("Mesh Optimization");
		elapsedTime += et;	
		iterations ++;

		qt->refresh(); 
		std :: cout<<i<< " Min diedral angle " << qt->fDieMin <<" Num negative elements "<<qt->nonPositive<< " Time "<< et<<"\n";
		
		if ( (numMejoras>0)&&(qt->nonPositive>0) )
		{
			std :: cout<< " Reached iteration " << i <<"\n";
			break;
		}

	}
	
	startTimers();			
	endProcess("Mesh evaluation Cavity :  benchmark2");	

	std :: cout<< " Mean time required " << elapsedTime/iterations <<"\n";
	
	qt->refresh();   qt->print();
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	m->validate(true);
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

  
  delete qt;
  delete m;
  delete ml;

}


// Generate Delaunnay
void benchmark5(std::string meshFileName)
{

	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 5 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	TetQuality *qt = new TetQuality(m);
	qt->refresh(); qt->print();
	delete qt;
	
	TVolumeMesh* m2 = (TVolumeMesh*)(GenerateMesh(m));
	delete m;
	m = m2;

	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	qt = new TetQuality(m);
	qt->refresh(); qt->print();

	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
	std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
	delete m;
}

// Generate Delaunnay
void benchmark7(std::string meshFileName, int maxIter, float dt)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 7 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName + ".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	TetQuality *qt = new TetQuality(m);
	qt->refresh(); qt->print();
	delete qt;

    qt = new TetQuality(m);

	for (int i=0 ; i<m->vertexes->Count();i++)
		{
			m->vertexes->elementAt(i)->fixed = false;
		}

		for (int i=0 ; i<m->fFaces->Count() ; i++)
		{
			TTriangle *tr = (TTriangle*)( m->fFaces->elementAt(i));
			if (tr == NULL) continue;
			tr->vertexes[0]->fixed = true;
			tr->vertexes[1]->fixed = true;
			tr->vertexes[2]->fixed = true;		
		}
	for (int i=0 ;i<maxIter ; i++)
	{	
		int substeps = 10;		
		MoveNodes(m,dt,substeps,1);			
		qt->refresh(); qt->print();		
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";
		myfile<<i<<"\t"<<0<<"\t"<<qt->fFaceMin<<"\t"<<qt->fDieMax<<"\t"<<qt->fDieMin<<"\t"<<qt->nonPositive<<"\n";
	}
    /// First Part
	std :: cout<< " ........... FIRST PART ............"<<"\n";
	TVolumeMesh* meshes[50] ; 
	for (int i=0 ;i<maxIter ; i++)
	{	
		int substeps = 10;
		
		MoveNodes(m,dt,substeps,1);	
		startProcess("Delaunnay Optimization");		
		TVolumeMesh* m2 = (TVolumeMesh*)(GenerateMesh(m));
		double et = endProcess("Delaunnay Optimization");
		
		meshes[i] = m;

		m = m2;

		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
		
		for (int i=0 ; i<m->vertexes->Count();i++)
		{
			m->vertexes->elementAt(i)->fixed = false;
		}

		for (int i=0 ; i<m->fFaces->Count() ; i++)
		{
			TTriangle *tr = (TTriangle*)( m->fFaces->elementAt(i));
			if (tr == NULL) continue;
			tr->vertexes[0]->fixed = true;
			tr->vertexes[1]->fixed = true;
			tr->vertexes[2]->fixed = true;		
		}

		qt = new TetQuality(m);
		qt->refresh(); qt->print();		
		std :: cout<< " Number of faces" << m->fFaces->Count()<<"\n";

		myfile<<i<<"\t"<<et<<"\t"<<qt->fFaceMin<<"\t"<<qt->fDieMax<<"\t"<<qt->fDieMin<<"\t"<<qt->nonPositive<<"\n";
	}

	/// Second Part
	std :: cout<< " ........... SECOND PART ............"<<"\n";
	m = (TVolumeMesh*)(ml->load(meshFileName + ".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	qt = new TetQuality(m);
	for (int i=0 ;i<maxIter ; i++)
	{	
		int substeps = 10;
		
		MoveNodes(m,dt,substeps,1);	
		startProcess("Mesh Optimization");
		stopTimers();
        optimizeMesh(m,3,15 ,APPLY_BY_EDGE * APPLY_BY_FACE);
     	startTimers();			
	    double et= endProcess("Mesh Optimization");		
	
	    qt->refresh();   qt->print();
	    m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
		for (int i=0 ; i<m->vertexes->Count();i++)
		{
			m->vertexes->elementAt(i)->fixed = false;
		}

		for (int i=0 ; i<m->fFaces->Count() ; i++)
		{
			TTriangle *tr = (TTriangle*)( m->fFaces->elementAt(i));
			if (tr == NULL) continue;
			tr->vertexes[0]->fixed = true;
			tr->vertexes[1]->fixed = true;
			tr->vertexes[2]->fixed = true;		
		}
	    m->validate(true);
    	/// File NumIter Smooth Time WorstQuality 
	    myfile<<i<<"\t"<<et<<"\t"<<qt->fFaceMin<<"\t"<<qt->fDieMax<<"\t"<<qt->fDieMin<<"\t"<<qt->nonPositive<<"\n";
	}
    myfile.close();
	
	for (int i=0 ;i<maxIter ; i++)
	  delete meshes[i];
	
}


///-- Test predefined meshes
void benchmark8()
{
	std::string fns[6];
	fns[0] = "dragon.vwm";
	fns[1] = "staypuft.vwm";
	fns[2] = "sculpt10kv.vwm";
	fns[3] = "cavityMesh_0_00_600k.vwm";
	fns[4] = "fiat_cil.vwm";
	fns[5] = "f1_tetgen_mesh.vwm";

	TMeshLoader* ml = new TVMWLoader();
	
	for (int i = 0; i<1 ; i++)
	{
	std::cout<<"Mesh file: "<<fns[i];
	for (int j = 1; j<32 ; j= j*2)
	{
		setNumThreads(j);	
		TVolumeMesh* m = (TVolumeMesh*)(ml->load(fns[i]));   
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
		TetQuality *qt = new TetQuality(m);
		qt->refresh(); 
		
	
		startProcess("Mesh Optimization");
		stopTimers();
		optimizeMesh(m,3,220 ,APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE * SMOOTH_MESH);
     	startTimers();			
	    double et= endProcess("Mesh Optimization");		
		m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		

		TetQuality *qt2 = new TetQuality(m);
		qt2->refresh(); 
		/// File NumIter Smooth Time WorstQuality 
	    myfile<<j<<"\t"<<i<<"\t"<<et<<"\t"<<"\t"<<qt->fDieMin<<"\t"<<"\t"<<qt2->fDieMin<<"\t"<<qt->nonPositive<<"\n";

		delete m;
		delete qt ; 
		delete qt2;
	}
	}
	 myfile.close();
}

void benchmark9(std::string meshFileName)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 9 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName + ".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	TetQuality *qt = new TetQuality(m);
	qt->refresh(); qt->print();

	for (int i=0; i<m->vertexes->Count();i++)
		m->vertexes->elementAt(i)->id = i;
	startProcess("ParallelSmoothMesh"); 
	stopTimers();
	ParallelSmoothMesh(m,vrelaxQuality,200);
	startTimers();
	endProcess("ParallelSmoothMesh");
	cout <<".........optimized........................" <<"\n"; 
	qt->refresh(); qt->print();
	delete qt;
	delete m;

	TVolumeMesh* m2 = (TVolumeMesh*)(ml->load(meshFileName + ".vwm"));   
	m2->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	qt = new TetQuality(m2);
	qt->refresh(); qt->print();
	
	initGPU();
	initGPUProgram("smooth.cl","lpMeshSmoothKernel");
	startProcess("ParallelSmoothMeshInGPU");
	stopTimers();
	ParallelSmoothMeshInGPU(m2,vrelaxQuality,200);
	startTimers();
	endProcess("ParallelSmoothMeshInGPU");
	cout <<".........optimized by GPU........................" <<"\n"; 
	qt->refresh(); qt->print();

	showProcessTime();
	delete qt;
	delete m2;

}

// Smooth Mesh
void benchmark6(std::string meshFileName)
{
	cout <<".............................................." <<"\n"; 
	cout <<".........BENCHMARK 6 ........................." <<"\n"; 
	cout <<" Loading Mesh..." <<"\n"; 
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	TetQuality *qt0 = new TetQuality(m);	
	qt0->refresh();
	
	cout <<".............................................." <<"\n"; 
	cout <<"...Initial Mesh quality" <<"\n"; 
	cout <<"Initial Quality. Min diedral : "<< qt0->fDieMin<< " Num neg:"<<qt0->nonPositive <<"\n"; 
	
	// Metodo laplaciano
	TetQuality *qt = new TetQuality(m);
	startProcess("laplacianMeshSmooth");
	laplacianMeshSmooth(m);	
	double time0 = endProcess("laplacianMeshSmooth");
	qt->refresh();
	cout <<"Laplacian. Min diedral : "<< qt->fDieMin<< " Num neg:"<<qt->nonPositive <<" Time "<<time0<<"\n"; 

	//TVolumeMesh* m2 = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	//m2->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);	
	
	// Metodo iterativo
	/*startProcess("iterativeMeshSmooth0");
	iterativeMeshSmooth(m2,50,10);
	double time1 = endProcess("iterativeMeshSmooth0");
	TetQuality *qt2 = new TetQuality(m2);
	qt2->refresh();
	cout <<"Iterativ conf 0. Min diedral : "<< qt2->fDieMin<< " Num neg:"<<qt2->nonPositive<<" Time "<<time1 <<"\n"; 
	*/

	// Metodo iterativo
	TVolumeMesh* m3 = (TVolumeMesh*)(ml->load(meshFileName+".vwm"));   
	m3->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);		
	startProcess("iterativeMeshSmooth1");
	iterativeMeshSmooth(m3,50,50,2);
	double time2 = endProcess("iterativeMeshSmooth1");
	TetQuality *qt3 = new TetQuality(m3);
	qt3->refresh();
	cout <<"Iterativ conf 1. Min diedral : "<< qt3->fDieMin<< " Num neg:"<<qt3->nonPositive<<" Time "<<time2 <<"\n"; 


	delete m;
	delete qt;

	//delete m2;
	//delete qt2;

	delete m3;
	delete qt3;


}


int main (int argc, char* argv[])
{
	omp_set_num_threads(16);
	cout <<".............................................." <<"\n"; 
	cout <<"Num procs "<<omp_get_num_procs()<< " Num threads "<< omp_get_num_threads()<<"\n"; 
	cout <<".............................................." <<"\n"; 
	
	myfile.open ("benchmark_9.txt");
	myfile<<"benchmark9"<<"\n";		
	setNumThreads(1);
	benchmark9("dragon");
	myfile.close();
	return 0;
	///------------------------------	
	
	myfile.open ("benchmark_2.txt");
	myfile<<"benchmark2"<<"\n";		
	setNumThreads(32);
	benchmark2("cavityMesh_0_00_600k",0.001,100,3);
	myfile.close();
	return 0;
    
	
	myfile.open ("benchmark_8.txt");
	myfile<<"benchmark8"<<"\n";		
	setNumThreads(32);
	benchmark0("f1_tetgen_mesh",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE * SMOOTH_MESH );
	return 0;



	myfile.open ("benchmark_2.txt");
	myfile<<"benchmark7"<<"\n";
	
	setNumThreads(32);
	myfile<<" FILE cavityMesh_0_00_600k"<<"\n";
	benchmark7("cavityMesh_0_00_600k" ,9, 0.0025);
    

	return 0;
	
	setNumThreads(1);	
	myfile<<"benchmark0"<<"\n";

	setNumThreads(32);

	benchmark0("cavityMesh_0_00",3,220 , APPLY_BY_EDGE);	
	setNumThreads(32);
	benchmark0("cavityMesh_0_00",3,220 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE);
	benchmark0("cavityMesh_0_00",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE );
	
	setNumThreads(1);
	benchmark0("fiat_cil",3,220 , APPLY_BY_EDGE);
	setNumThreads(32);
	benchmark0("fiat_cil",3,220 , APPLY_BY_EDGE);
	benchmark0("fiat_cil",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE);
	benchmark0("fiat_cil",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE );

	setNumThreads(1);
	benchmark0("cavityMesh_0_00_600k",3,220 , APPLY_BY_EDGE);
	setNumThreads(32);
	benchmark0("cavityMesh_0_00_600k",3,220 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00_600k",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE);
	benchmark0("cavityMesh_0_00_600k",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE );
	

	setNumThreads(1);
	benchmark0("f1_tetgen_mesh",3,220 , APPLY_BY_EDGE);
	setNumThreads(32);
	benchmark0("f1_tetgen_mesh",3,220 , APPLY_BY_EDGE);
	benchmark0("f1_tetgen_mesh",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE);
	benchmark0("f1_tetgen_mesh",3,220 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE );
	
	myfile<<"Filtered remeshing"<<"\n" ;

	benchmark0("cavityMesh_0_00",3,5 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00",3,10 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00",3,20 , APPLY_BY_EDGE);	
	benchmark0("cavityMesh_0_00",3,50 , APPLY_BY_EDGE);


	benchmark0("fiat_cil",3,5 , APPLY_BY_EDGE);
	benchmark0("fiat_cil",3,10 , APPLY_BY_EDGE);
	benchmark0("fiat_cil",3,20 , APPLY_BY_EDGE);	
	benchmark0("fiat_cil",3,50 , APPLY_BY_EDGE);
	
	benchmark0("cavityMesh_0_00_600k",3,5 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00_600k",3,10 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00_600k",3,20 , APPLY_BY_EDGE);	
	benchmark0("cavityMesh_0_00_600k",3,50 , APPLY_BY_EDGE);

	benchmark0("f1_tetgen_mesh",3,1 , APPLY_BY_EDGE);
	benchmark0("f1_tetgen_mesh",3,3 , APPLY_BY_EDGE);
	benchmark0("f1_tetgen_mesh",3,5 , APPLY_BY_EDGE);	
	benchmark0("f1_tetgen_mesh",3,10 , APPLY_BY_EDGE);


	//setNumThreads(32);
	//benchmark0Bis("f1_tetgen_mesh",3,2 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE);
	//benchmark0("f1_tetgen_mesh",3,2 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE);
	
	/*benchmark0("cavityMesh_0_00_600k",3,5 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00_600k",3,10 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00_600k",3,20 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00_600k",3,30 , APPLY_BY_EDGE);
	benchmark0("cavityMesh_0_00_600k",3,50 , APPLY_BY_EDGE);

	benchmark0("cavityMesh_0_00_600k",3,5 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE * SMOOTH_MESH);
	benchmark0("cavityMesh_0_00_600k",3,10 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE * SMOOTH_MESH);

	benchmark0("f1_tetgen_mesh",3,2 , APPLY_BY_EDGE);
	benchmark0("f1_tetgen_mesh",3,5 , APPLY_BY_EDGE);
	benchmark0("f1_tetgen_mesh",3,10 , APPLY_BY_EDGE);
	benchmark0("f1_tetgen_mesh",3,20 , APPLY_BY_EDGE);
	benchmark0("f1_tetgen_mesh",3,50 , APPLY_BY_EDGE);

	benchmark0("f1_tetgen_mesh",3,2 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE * SMOOTH_MESH);
	benchmark0("f1_tetgen_mesh",3,5 , APPLY_BY_EDGE * APPLY_BY_FACE * APPLY_BY_NODE * SMOOTH_MESH);

	
	*/
	myfile.close();

//benchmark5("cavityMesh_0_00.vwm");
	//benchmark5("cavityMesh_0_00.vwm");
	//benchmark5("cavityMesh_0_00.vwm");
	//benchmark5("cavityMesh_0_00.vwm");

	//benchmark0("cavityMesh_0_00.vwm");
	//benchmark5("f1_tetgen_mesh.vwm");

	//benchmark5("cavityMesh_0_00.vwm");

	//benchmark2("cavityMesh_0_00.vwm",0.01,40,0);


	//benchmark3("cavityMesh_0_00.vwm",0.02,40);

	


	/*
	std :: cout<< " Loading mesh.............."<<"\n";
	TMeshLoader* ml = new TVMWLoader();
	TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/cavity3D.vwm"));   
	//TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/f1_tetgen_mesh.vwm"));   
	//TVolumeMesh* m = (TVolumeMesh*)(ml->load("d:/posDoc/kratosProof.vwm"));   
	TVolumeMesh* m = (TVolumeMesh*)(ml->load("f1_tetgen_mesh.vwm"));   
	
	
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
	omp_set_num_threads(4);
	cout <<".............................................." <<"\n"; 
	cout <<"Num procs "<<omp_get_num_procs()<< " Num threads "<< omp_get_num_threads()<<"\n"; 
	cout <<".............................................." <<"\n"; 
	cout <<"...Preparing optimization" <<"\n"; 
	TetQuality *qt = new TetQuality(m);
	
	//cout <<"...Optimizing by Node" <<"\n"; 
	//evaluateClusterByNode( (TVolumeMesh*)(m),0.5,vrelaxQuality);
	cout <<".............................................." <<"\n"; 
	cout <<"...Initial Mesh quality" <<"\n"; 
	qt->refresh();   qt->print();
	cout <<".............................................." <<"\n"; 
	
	startProcess("Mesh evaluation");
	stopTimers();

	int origNumFaces = m->fFaces->Count();
	//ParallelEvaluateClusterByNode(m,vrelaxQuality);   
	for (int i=0 ;i<1 ; i++)
	{
		int substeps = 10;
		float dt = 0.01;
		//MoveNodes(m,dt,substeps,1);
		

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
			//ParallelEvaluateClusterByNode(m,vrelaxQuality,220);
			
			cout <<"...Parallel Optimizing by Edge" <<"\n"; 
			ParallelEvaluateClusterByEdge(m,vrelaxQuality,220);
			
			cout <<"...Parallel Optimizing by Face" <<"\n"; 
			//ParallelEvaluateClusterByFace(m,vrelaxQuality,220);
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
		*/
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
	*/
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
	/*showProcessTime();
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
	*/
	return 0;
} 
