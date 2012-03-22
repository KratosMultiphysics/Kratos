
#include <iostream>
#include <string>
#include "u_MeshLoaders.h"
#include "u_TetraFunctions.h"

	
	bool TSurLoader:: save(char* filename, TMesh* aMesh )
	{ 
		  int  i,nG;
          TTriangle *tr;

		TStringList* st = new TStringList();
    	 st->Add("*ELEMENT GROUPS");

		 char* s = new char[120];
		 
		if (aMesh->patchList== NULL)
			nG = 0 ;
		else
			nG = aMesh->patchList->Count();
		if (nG==0)
		{
			st->Add("1");
			sprintf(s, "%s%d%s","1 ", aMesh->fFaces->Count(), " Tri3");   
			st->Add(s);
		}
		else
		{
			sprintf(s, "%d", aMesh->patchList->Count());   
			st->Add(s);	
		}

		st->Add("*SHELL GROUPS");
		st->Add("1");
		st->Add("*INCIDENCE");
		for (i=0 ; i<aMesh->vertexes->Count() ; i++)
			aMesh->vertexes->elementAt(i)->id = i+1;
		for (i = 0 ; i< aMesh->fFaces->Count() ; i++)
		{
			tr = (TTriangle*)(  aMesh->fFaces->elementAt(i));
			sprintf(s, "%d%s%d%s%d",  tr->vertexes[0]->id," ",tr->vertexes[1]->id," ",
									tr->vertexes[2]->id," ");
			st->Add(s);
		}
        
		st->Add("*COORDINATES");	
		sprintf(s, "%d", aMesh->vertexes->Count());
		st->Add(s);
		for (i=0; i<aMesh->vertexes->Count() ;i++)
		{
			TVertex* v = aMesh->vertexes->elementAt(i);
			sprintf(s, "%d%s%f%s%f%s%f",v->id," ", v->fPos.x," ",v->fPos.y," ",v->fPos.z);
			st->Add(s);
		}
 
		st->saveToFile(filename);
		delete st;
 		return false;  
	}

	

	bool TVMWLoader::save(char* aMeshName, TMesh* aMesh )
	{ 
		//----------------------------
		TStringList* F;
        //S, S2 : String;
        int i;
		TVolumeMesh* m ;
        m = (TVolumeMesh*)(aMesh);
        F = new TStringList();
        F->Add("*COORDINATES");
		char* s = new char[120];
		sprintf(s, "%d", m->vertexes->Count());
		F->Add(s);
		for (i=0; i<m->vertexes->Count() ;i++)
		{
			TVertex* v = m->vertexes->elementAt(i);
			v->id = i+1;
			
			sprintf(s, "%d%s%f%s%f%s%f",v->id," ", v->fPos.x," ",v->fPos.y," ",v->fPos.z);
			F->Add(s);
		}

		F->Add("*ELEMENT GROUPS");
		F->Add( "1");		
		sprintf(s, "%s%d%s","1 ", m->elements->Count()," Tetra4");
        F->Add( s);
        F->Add( "*INCIDENCE");

		for (i=0; i<m->elements->Count();i++)
		{
			TTetra* t= (TTetra*)(m->elements->elementAt(i));
			sprintf(s, "%d%s%d%s%d%s%d",  t->vertexes[0]->id," ",t->vertexes[1]->id," ",
									t->vertexes[2]->id," ",t->vertexes[3]->id);
			F->Add(s);
		}
		F->Add( "*END");
		F->saveToFile(aMeshName);    

		delete F;

		return false;  
	}
	TMesh* TVMWLoader::load(char* aMeshName)
	{ 
	    FILE* fMesh;
		TVolumeMesh* m = new TVolumeMesh();
   
        fopen_s(&fMesh, aMeshName,"rb"); //xx = rb, wb, read and write binary, more 		
		
		char line[200];		
		int nCoords;		
		// Read first Line
        fread (line,14,1,fMesh); 
		// Read coords
		fscanf_s(fMesh, "%d", &nCoords);	
		for (int i=0; i<nCoords;i++)
		{
			int id;
			float x,y,z;
			fscanf_s(fMesh, "%d", &id);
			fscanf_s(fMesh, "%f", &x);
			fscanf_s(fMesh, "%f", &y);
			fscanf_s(fMesh, "%f", &z);			

			TVertex* v = new TVertex(x,y,z);
			v->id = id;
			m->addVertex(v);
		}

		// Read *ELEMENT GROUPS
		int elGroups, numElements;
		
		fread (line,17,1,fMesh); 
		fscanf_s(fMesh, "%d", &elGroups);	
		fscanf_s(fMesh, "%d%d", &elGroups,&numElements);	
		fread (line,9,1,fMesh); 
		fread (line,11,1,fMesh); 
		for (int i=0; i<numElements;i++)
		{
			int iv0,iv1,iv2,iv3;
			fscanf_s(fMesh, "%d", &iv0);
			fscanf_s(fMesh, "%d", &iv1);
			fscanf_s(fMesh, "%d", &iv2);			
			fscanf_s(fMesh, "%d", &iv3);	

			TVertex* v0 = m->findVertexById(iv0);
			TVertex* v1 = m->findVertexById(iv1);
			TVertex* v2 = m->findVertexById(iv2);
			TVertex* v3 = m->findVertexById(iv3);

			TTetra* t = new TTetra(NULL, v0,v1,v2,v3);
			m->elements->Add(t);
		}
		fclose(fMesh);

		return m; 
	
	}


	bool TElementTetraLoader::save(char* aMeshName, TMesh* aMesh )
	{
      TStringList* st;
      int i,j;
      TVertex *v;
      TTriangle* tr;
      TTetra* t;
      TVolumeMesh *m;
	  char* s;
      // s ,fdir, fn : string;

	// ....................
	st = new TStringList();
	m = (TVolumeMesh*)(aMesh);
	st->Add("NEIGH");
	s = intToStr(m->vertexes->Count());
	st->Add(s);

	for (i=0; i<m->elements->Count();i++)
	{
		((TTetra*)(m->elements->elementAt(i)))->id = i;
	}

	for (i=0; i<m->vertexes->Count();i++)
	{
		((TVertex*)(m->vertexes->elementAt(i)))->id = i;
	}

	for (i=0; i<m->fFaces->Count();i++)
	{
		((TTriangle*)(m->fFaces->elementAt(i)))->id = i;
	}

	
	// Tetra Vertex Neighbours
	s = new char[400];	
	char* s2 = new char[100];	
 
	for (i = 0 ; i<m->vertexes->Count() ; i++)
	{
		v = m->vertexes->elementAt(i);
		if (v->elementsList == NULL) continue;
		intToStr(v->id,s);
		strcat(s, " ");
				
		for (j = 0 ; j<v->elementsList->Count() ; j++)
		{
			t = (TTetra*)(v->elementsList->elementAt(j));
			if (!t) continue;
			intToStr(t->id,s2);
	        strcat(s, s2);
			strcat(s, " ");
		}
		st->Add(s);
	}
	
	//delete st;
	st->Add("FACES");
    intToStr(m->fFaces->Count(),s);
    st->Add(s);
    
	for (i=0 ; i<m->fFaces->Count();i++)
	{
		tr = (TTriangle*)(m->fFaces->elementAt(i));
        intToStr(tr->id,s);
		strcat(s, " ");
				
		for (j = 0 ; j<3 ; j++)
		{
			v = tr->vertexes[j];
			if (!v) continue;
			intToStr(v->id,s2);
	        strcat(s, s2);
			strcat(s, " ");
		}
		st->Add(s);
	}
	// VERTEXES
	st->Add("VERTEXES");
    intToStr(m->vertexes->Count(),s);
    st->Add(s);
    
	for (i=0 ; i<m->vertexes->Count();i++)
	{
		v = (TVertex*)(m->vertexes->elementAt(i));
        intToStr(v->id,s);
		strcat(s, " ");
				
		floatToStr(v->fPos.x,s2);	    strcat(s, s2); 		strcat(s, " ");
		floatToStr(v->fPos.y,s2);	    strcat(s, s2); 		strcat(s, " ");
		floatToStr(v->fPos.z,s2);	    strcat(s, s2); 
		
		st->Add(s);
	}
	// ELEMENTS
	st->Add("ELEMENTS");
	intToStr(m->elements->Count(),s);
    st->Add(s);
    
	for (i=0 ; i<m->elements->Count();i++)
	{
		t = (TTetra*)(m->elements->elementAt(i));
        intToStr(t->id,s);
		strcat(s, " ");
				
		for (j = 0 ; j<4 ; j++)
		{
			v = t->vertexes[j];
			if (!v) continue;
			intToStr(v->id,s2);
	        strcat(s, s2);
			strcat(s, " ");
		}
		st->Add(s);
	}

	// VERTEXES x VERTEXES ELEMENTS NEIGHBOURS
	st->Add("VELEMENTS");
	intToStr(m->vertexes->Count(),s);
    st->Add(s);

	TList<TObject*> *lneigh = new TList<TObject*>();
    
	for (i=0 ; i<m->vertexes->Count();i++)
	{
		v = m->vertexes->elementAt(i);
        intToStr(v->id,s);
		strcat(s, " ");

		v->getElemNeighbours(lneigh);    
				
		for (j = 0 ; j<lneigh->Count() ; j++)
		{
			TVertex* v3 = (TVertex*)(lneigh->elementAt(j));
			if (!v3) continue;
			intToStr(v3->id,s2);
	        strcat(s, s2);
			strcat(s, " ");
		}
		st->Add(s);
	}


	st->saveToFile(aMeshName);
	return true;
 //st.SaveToFile(fdir+fn+'.NEIGH');
 // .........................
	}


	bool TGIDLoad::save(char* aMeshName, TMesh* aMesh )
	{
		return false;
	}
	
	TMesh* TGIDLoad::load(char* aMeshName)
	{
		 FILE* fMesh;
		TVolumeMesh* m = new TVolumeMesh();
   
		TStringList* st = new TStringList();
		st->loadFromFile(aMeshName);
        
		// Read coords
		int i=0;

		while ( st->strings[i] != "Begin Nodes") i++;
		
		while ( true)
		{
			float x,y,z;
			std::string s = st->strings[i];

			fscanf_s(fMesh, "%f", &x);
			fscanf_s(fMesh, "%f", &x);
			fscanf_s(fMesh, "%f", &y);
			fscanf_s(fMesh, "%f", &z);			

			TVertex* v = new TVertex(x,y,z);
			v->id = i;
			m->addVertex(v);
		}

		// Find Elements position
		/*while ( strcmp( line , "Begin Elements" )== 0) fread (line,14,1,fMesh); 
		
		
		while ( strcmp( line , "End Elements" )== 0)
		{
			int iv0,iv1,iv2,iv3;
			// id element
			fscanf_s(fMesh, "%d", &iv0);
			// other
			fscanf_s(fMesh, "%d", &iv0);
			// elements
			fscanf_s(fMesh, "%d", &iv0);
			fscanf_s(fMesh, "%d", &iv1);
			fscanf_s(fMesh, "%d", &iv2);			
			fscanf_s(fMesh, "%d", &iv3);	

			TVertex* v0 = m->vertexes->elementAt(iv0-1);
			TVertex* v1 = m->vertexes->elementAt(iv1-1);
			TVertex* v2 = m->vertexes->elementAt(iv2-1);
			TVertex* v3 = m->vertexes->elementAt(iv3-1);

			TTetra* t = new TTetra(NULL, v0,v1,v2,v3);
			m->elements->Add(t);
		}
		*/
		fclose(fMesh);

		return m; 
	}

	