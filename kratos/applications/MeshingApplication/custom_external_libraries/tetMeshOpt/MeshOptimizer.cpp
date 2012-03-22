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
   //m->normalizeMesh(1000);
   m->updateIndexes(GENERATE_SURFACE);

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
     evaluateClusterByNode( (TVolumeMesh*)(m),5000000,vrelaxQuality);
     m->updateIndexes(GENERATE_SURFACE);
   endProcess("Parallel evaluation");
   qt->refresh();   qt->print();
   m->validate(true);

   cout <<"...Parallel Optimizing by Node" <<"\n"; 
   startProcess("Parallel evaluation");
     //ParallelEvaluateClusterByNode(m,vrelaxQuality);   
     evaluateClusterByNode( (TVolumeMesh*)(m),5000000,vrelaxQuality);
     m->updateIndexes(GENERATE_SURFACE);
   endProcess("Parallel evaluation");
   qt->refresh();   qt->print();
   m->validate(true);

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
	           ml2->save("D:/out_MeshC_Opt.vwm" , m);
			   delete ml2;
			}
   
   
   ////showProcessTime();
      
   return 0;
} 
