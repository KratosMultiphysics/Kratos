#include <iostream>
#include <string>
#include <cstdio>
#include <string>

#include "tetgen.h" // Defined tetgenio, tetrahedralize().
#include "u_TetGenInterface.h"
#include "u_TetraFunctions.h"

 TMesh* GenerateMesh(TMesh* inM )
        {
            tetgenio in, out;
			TVolumeMesh *aMesh = (TVolumeMesh*)(inM);
			
			// All indices start from 0.
            in.initialize();
            in.firstnumber = 0;
			
			in.numberofpoints = aMesh->vertexes->Count();
            in.pointlist = new REAL[in.numberofpoints * 3];
			
            //give the coordinates to the mesher            
            for (int i = 0; i < aMesh->vertexes->Count(); i++)
            {
                int base = (int)(i * 3);
				TVertex* v = (TVertex*)( aMesh->vertexes->elementAt(i));
				v->setID( i );
                in.pointlist[base] = v->fPos.x;
                in.pointlist[base + 1] = v->fPos.y;
                in.pointlist[base + 2] = v->fPos.z;
            }
			
            //prepare the list of faces
			in.numberoffacets = aMesh->fFaces->Count();
            in.facetlist = new tetgenio::facet[in.numberoffacets];
            in.facetmarkerlist = new int[in.numberoffacets];
            //give the surface connectivity to the mesher
			for (int i = 0; i < aMesh->fFaces->Count(); i++)
            {
				TTriangle* tr = (TTriangle*)(aMesh->fFaces->elementAt(i));
				
				tetgenio::facet* f = &in.facetlist[i];
                f->numberofpolygons = 1;
                f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
                f->numberofholes = 0;
                f->holelist = NULL;
                tetgenio::polygon* p = &f->polygonlist[0];
                p->numberofvertices = 3;
                p->vertexlist = new int[p->numberofvertices];
				p->vertexlist[0] = tr->vertexes[0]->getID();
                p->vertexlist[1] = tr->vertexes[1]->getID();
                p->vertexlist[2] = tr->vertexes[2]->getID();                
			}
            //give the hole list to the mesher
            in.numberofholes = 0;
            
            char* tetgen_options = (char*)("-pY");
            //perform meshing with tetgen
            tetrahedralize(tetgen_options, &in, &out);
			
			TVolumeMesh* outMesh = new TVolumeMesh();
			

			//generate  nodes            
            for(int i=0; i< out.numberofpoints; i++)
            {
                int base = i * 3;
                double x = out.pointlist[base];
                double y = out.pointlist[base+1];
                double z = out.pointlist[base+2];
				TVertex*v = new TVertex(x,y,z);
				v->setID( i );
				outMesh->vertexes->Add(v);
			}

			//generate new Elements    
            for(int i=0; i<out.numberoftetrahedra; i++)
            {
                int base=i*4;
				TTetra *t ;
				TVertex *v0 = outMesh->findVertexById( out.tetrahedronlist[base]); 
				TVertex *v1 = outMesh->findVertexById( out.tetrahedronlist[base+1]); 
				TVertex *v2 = outMesh->findVertexById( out.tetrahedronlist[base+2]); 
				TVertex *v3 = outMesh->findVertexById( out.tetrahedronlist[base+3]); 

				t = new TTetra(NULL, v0,v1,v2,v3);
				outMesh->elements->Add(t);                 
            }

			
            in.deinitialize();
            out.deinitialize();
            in.initialize(); //better deinitialize and initialize again
            out.initialize();

			return outMesh;


        }