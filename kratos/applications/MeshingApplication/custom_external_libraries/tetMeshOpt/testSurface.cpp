/*
    Copyright 2005-2011 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#include <iostream>
#include <string>
#include <algorithm>
#include "u_ProcessTime.h"
#include "u_Types.h"
#include "u_TetraFunctions.h"
#include "u_MeshLoaders.h"


using namespace std;
#include <omp.h>

static const size_t N = 23;

	
typedef void(*TStForLoopElement)(int , int , TObject*);

bool* faces;


void localProcess(int i, int thId  ,TObject* destObject)
{
  int j ,k;
  TVertex* v0, *v1, *v2;
  TTetra *t,*t2;
  bool res;

    t = (TTetra*)( destObject);
	if (t ==  NULL ) return ;

    for (j = 0 ; j<4; j++)
	{
      faces[i*4+j] = false;
      v0 = t->vertexes[ TTetraFaces[j*3] ];
      v1 = t->vertexes[ TTetraFaces[j*3+1] ];
      v2 = t->vertexes[ TTetraFaces[j*3+2] ];
      res = true;
	  for (k = 0 ; k<v0->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v0->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
      if (!res) continue;
      for (k = 0 ; k<v1->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v1->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
      if (!res) continue;
      for (k = 0 ; k<v2->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v2->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
       if (res)       
        faces[i*4+j] = true;
	}
}

void serialGetSurfaceTriangles(TVolumeMesh* aModel)
{
  int i,j , numT;
  TVertex*v0,*v1,*v2;
  TTetra *t;
  TTriangle *tr;
  

  numT = aModel->numElements;
  faces = new bool[numT * 4];
  for (i = 0 ; i<numT*4 ; i++)
	   faces[i] = false;
  
  aModel->fFaces->Clear();
  
  startProcess("parallelPart");
    for (i = 0 ; i< numT ; i++)
	{
		TTetra *t = (TTetra*)(aModel->elements->elementAt(i));
		localProcess(i,0,t);
	}
   // parallel_for(blocked_range<size_t>(0,numT), TParallelIterator(localProcess,aModel->elements) );    
  endProcess("parallelPart");

  startProcess("createPart");
  for (i = 0 ; i<numT ; i++)
  {
	t = (TTetra*)( aModel->elements->elementAt(i));
    for (j = 0 ; j<4 ;j++)
	{
      if (!faces[i*4+j]) continue;
      v0 = t->vertexes[ TTetraFaces[j*3] ];
      v1 = t->vertexes[ TTetraFaces[j*3+1] ];
      v2 = t->vertexes[ TTetraFaces[j*3+2] ];
      tr = new TTriangle( v0,v1,v2);
	  tr->calcNormal();
	  aModel->addTriangle(tr);
	}
  }
  endProcess("createPart");
  delete faces;
}


void ilocalProcess(int i, int thId  ,TObject* destObject)
{
  int j ,k;
  TVertex* v0, *v1, *v2;
  TTetra *t,*t2;
  bool res;

    t = (TTetra*)( destObject);
	if (t ==  NULL ) return ;

    for (j = 0 ; j<4; j++)
	{
      faces[i*4+j] = false;
      v0 = t->vertexes[ TTetraFaces[j*3] ];
      v1 = t->vertexes[ TTetraFaces[j*3+1] ];
      v2 = t->vertexes[ TTetraFaces[j*3+2] ];
      res = true;
	  for (k = 0 ; k<v0->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v0->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
      if (!res) continue;
      for (k = 0 ; k<v1->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v1->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
      if (!res) continue;
      for (k = 0 ; k<v2->elementsList->Count(); k++)
	  {
		  t2 = (TTetra*)(v2->elementsList->elementAt(k));
          if (t2 == t) continue;           
		  if (t2->hasFace(v0,v1,v2) )
			  { 
				   res = false; 
				   break; 
		      }
	  }
       if (res)       
        faces[i*4+j] = true;
	}
}

void fastGetSurfaceTriangle(TMesh* aModel)
{
	
  int i,j , numT;
  TVertex*v0,*v1,*v2;
  TTetra *t;
  TTriangle *tr;

  numT = aModel->elements->Count();
  faces = new bool[numT * 4];
  for (i = 0 ; i<numT*4 ; i++)
	   faces[i] = false;
  
  aModel->fFaces->Clear();
  
  startProcess("parallelPart");
    //parallel_for(blocked_range<size_t>(0,numT), TParallelIterator(localProcessI,aModel->elements) );    
    #pragma omp parallel num_threads(4)  
   {  			
			    for (i = 0; i<numT ; i++)					
					ilocalProcess(i,0,aModel->elements->elementAt(i));
   }	  
  endProcess("parallelPart");

  startProcess("createPart");
  for (i = 0 ; i<numT ; i++)
  {
	t = (TTetra*)( aModel->elements->elementAt(i));
    for (j = 0 ; j<4 ;j++)
	{
      if (!faces[i*4+j]) continue;
      v0 = t->vertexes[ TTetraFaces[j*3] ];
      v1 = t->vertexes[ TTetraFaces[j*3+1] ];
      v2 = t->vertexes[ TTetraFaces[j*3+2] ];
      tr = new TTriangle( v0,v1,v2);
	  tr->calcNormal();
	  aModel->addTriangle(tr);
	}
  }
  endProcess("createPart");
  delete faces;
  
}
int main() 
{
   int nthreads = omp_get_num_threads();
      cout << "There are " << nthreads << " threads" << '\n';
   TMeshLoader* ml = new TVMWLoader();
   TVolumeMesh* m = (TVolumeMesh*)(ml->load("D:/posDOC/cube_sphere2.vwm"));
   //TVolumeMesh* m = (TVolumeMesh*)(ml->load("D:/posDOC/simple.vwm"));
 	m->updateIndexes(0);

   m->normalizeMesh(1000);
   m->fFaces->Clear();
   startProcess("Serial surface");
   serialGetSurfaceTriangles(m );
   endProcess("Serial surface");

   //TMeshLoader* mn = new TElementTetraLoader();
   //mn->save("d:/outC.NEIGH",m);
   //delete mn;

   m->fFaces->Clear();
   startProcess("Parallel surface");
   fastGetSurfaceTriangle(m );
   endProcess("Parallel surface");

   //mn = new TElementTetraLoader();
   //mn->save("d:/outCP.NEIGH",m);
   //delete mn;
   showProcessTime();
   cout <<"...Waiting for keyboard input to save" ;
   getchar();   
   
   TMeshLoader* ms = new TSurLoader();
     
   ms->save("D:/posDOC/test.sur" , m);
   delete ms;
   cout <<"...Save OK. Now exit" ;
   getchar();   
   return 0;
/*
  string str[N] = { string("a"), string("b") };
  for (size_t i = 2; i < N; ++i) str[i] = str[i-1]+str[i-2];
  string &to_scan = str[N-1]; 
  size_t num_elem = to_scan.size();

  size_t *max = new size_t[num_elem];
  size_t *pos = new size_t[num_elem];

  parallel_for(blocked_range<size_t>(0, num_elem ),
               SubStringFinder( to_scan, max, pos ) );

  for (size_t i = 0; i < num_elem; ++i)
    cout << " " << max[i] << "(" << pos[i] << ")" << endl;
  delete[] pos;
  delete[] max;
*/
  return 0;
}

