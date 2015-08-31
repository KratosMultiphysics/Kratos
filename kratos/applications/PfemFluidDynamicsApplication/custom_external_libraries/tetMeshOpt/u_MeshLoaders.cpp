
#include <iostream>
#include <string>
#include <cstdio>
#include "u_MeshLoaders.h"
#include "u_TetraFunctions.h"


bool TSurLoader:: save(const char* filename, TMesh* aMesh, int flags  )
{ 
	int  i,nG;
	TTriangle *tr;

	TStringList* st = new TStringList();
	st->Add((char*)("*ELEMENT GROUPS"));

	char* s = new char[120];

	if (aMesh->patchList== NULL)
		nG = 0 ;
	else
		nG = aMesh->patchList->Count();
	if (nG==0)
	{
		st->Add((char*)("1"));
		sprintf(s, "%s%d%s","1 ", aMesh->fFaces->Count(), " Tri3");   
		st->Add(s);
	}
	else
	{
		sprintf(s, "%d", aMesh->patchList->Count());   
		st->Add(s);	
	}

	st->Add((char*)("*SHELL GROUPS"));
	st->Add((char*)("1"));
	st->Add((char*)("*INCIDENCE"));
	for (i=0 ; i<aMesh->vertexes->Count() ; i++)
		aMesh->vertexes->elementAt(i)->setID( i+1 );
	for (i = 0 ; i< aMesh->fFaces->Count() ; i++)
	{
		tr = (TTriangle*)(  aMesh->fFaces->elementAt(i));
		std::string s2 = intToString(tr->vertexes[0]->getID())+ " "+intToString(tr->vertexes[1]->getID())+ " "+
			             intToString(tr->vertexes[2]->getID());

		st->Add(s2);
	}

	st->Add((char*)("*COORDINATES"));	
	sprintf(s, "%d", aMesh->vertexes->Count());
	st->Add(s);
	for (i=0; i<aMesh->vertexes->Count() ;i++)
	{
		TVertex* v = aMesh->vertexes->elementAt(i);
		sprintf(s, "%d%s%f%s%f%s%f",v->getID()," ", v->fPos.x," ",v->fPos.y," ",v->fPos.z);
		st->Add(s);
	}

	st->saveToFile(filename);
	delete st;
	return false;  
}



bool TVMWLoader::save(const char* aMeshName, TMesh* aMesh, int flags )
{ 
	//----------------------------
	TStringList* F;
	//S, S2 : String;
	int i;
	TVolumeMesh* m ;
	m = (TVolumeMesh*)(aMesh);
	F = new TStringList();
	F->Add((char*)("*COORDINATES"));
	char* s = new char[120];
	sprintf(s, "%d", m->vertexes->Count());
	F->Add(s);
	for (i=0; i<m->vertexes->Count() ;i++)
	{
		TVertex* v = m->vertexes->elementAt(i);
		if (flags == 1) v->setID( i);
		sprintf(s, "%d%s%f%s%f%s%f%s%f",v->getID()," ", v->fPos.x," ",v->fPos.y," ",v->fPos.z," ",v->expectedSize);
		F->Add(s);
	}

	F->Add((char*)("*ELEMENT GROUPS"));
	F->Add((char*)("1"));		
	sprintf(s, "%s%d%s","1 ", m->elements->Count()," Tetra4");
	F->Add( s);
	F->Add( (char*)("*INCIDENCE"));

	for (i=0; i<m->elements->Count();i++)
	{
		TTetra* t= (TTetra*)(m->elements->elementAt(i));
		sprintf(s, "%d%s%d%s%d%s%d",  t->vertexes[0]->getID()," ",t->vertexes[1]->getID()," ",
			t->vertexes[2]->getID()," ",t->vertexes[3]->getID());
		F->Add(s);
	}
	F->Add( (char*)("*END"));
	F->saveToFile(aMeshName);    

	delete F;

	return false;  
}
TMesh* TVMWLoader::load(const char* aMeshName)
{ 
	FILE* fMesh;
	TVolumeMesh* m = new TVolumeMesh();
	char *line = new char[200];		
	char *line2 = new char[200];	
	int nCoords;		
	size_t ctrlFlag , lSize;

	fMesh = fopen( aMeshName,"rb"); //xx = rb, wb, read and write binary, more 		
	if (fMesh==NULL) 
	{ 
		std::cout << "Invalid filename " << aMeshName;
		return NULL;
	}
	// obtain file size:
	fseek (fMesh , 0 , SEEK_END);
	lSize = ftell (fMesh);
	rewind (fMesh);


	// Read first Line
	//ctrlFlag = fread (,14,1,fMesh); 
	ctrlFlag = fscanf(fMesh, "%s", line);	
	// Read coords

	ctrlFlag = fscanf(fMesh, "%d", &nCoords);	
	if (ctrlFlag == 0) { std::cout << "Invalid read access " << aMeshName; 	return NULL; }

	for (int i=0; i<nCoords;i++)
	{
		int id;
		float x,y,z;
		ctrlFlag =fscanf(fMesh, "%d", &id);
		if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	break;			}

		ctrlFlag =fscanf(fMesh, "%f", &x);
		if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	break;			}

		ctrlFlag =fscanf(fMesh, "%f", &y);
		if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	break;			}

		ctrlFlag =fscanf(fMesh, "%f", &z);			
		if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	break;			}

		TVertex* v = new TVertex(x,y,z);
		v->setID( id );
		m->addVertex(v);
	}

	// Read *ELEMENT GROUPS
	int elGroups, numElements;

	//ctrlFlag =fread (line,17,1,fMesh); 
	ctrlFlag = fscanf(fMesh, "%s%s", line,line2);	
	if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

	ctrlFlag = fscanf(fMesh, "%s", line);	
	//if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

	ctrlFlag =fscanf(fMesh, "%d%d%s", &elGroups,&numElements,line);	
	if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

	ctrlFlag = fscanf(fMesh, "%s", line);	
	if (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

	for (int i=0; i<numElements;i++)
	{
		int iv0,iv1,iv2,iv3;
		ctrlFlag =fscanf(fMesh, "%d", &iv0);
		if  (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

		ctrlFlag =fscanf(fMesh, "%d", &iv1);
		if  (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

		ctrlFlag =fscanf(fMesh, "%d", &iv2);			
		if  (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

		ctrlFlag =fscanf(fMesh, "%d", &iv3);	
		if  (ctrlFlag == 0 )  { 		std::cout << "Invalid read access " << aMeshName; 	return NULL;			}

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


bool TElementTetraLoader::save(const char* aMeshName, TMesh* aMesh, int flags  )
{
	TStringList* st;
	int i,j;
	TVertex *v;
	TTriangle* tr;
	TTetra* t;
	TVolumeMesh *m;
	std::string s;
	std::string s2;
	// s ,fdir, fn : string;

	// ....................
	st = new TStringList();
	m = (TVolumeMesh*)(aMesh);
	st->Add((char*)("NEIGH"));
	s = intToStr(m->vertexes->Count());
	st->Add(s);
	/*
	for (i=0; i<m->elements->Count();i++)
	{
	((TTetra*)(m->elements->elementAt(i)))->id = i;
	}
	*/
	for (i=0; i<m->vertexes->Count();i++)
	{
		((TVertex*)(m->vertexes->elementAt(i)))->setID( i );
	}

	for (i=0; i<m->fFaces->Count();i++)
	{
		((TTriangle*)(m->fFaces->elementAt(i)))->setID( i );
	}


	// Tetra Vertex Neighbours



	for (i = 0 ; i<m->vertexes->Count() ; i++)
	{
		v = m->vertexes->elementAt(i);
		if (v->elementsList == NULL) continue;
		s = intToStr(v->getID());
		s = s + " ";

		for (j = 0 ; j<v->elementsList->Count() ; j++)
		{
			t = (TTetra*)(v->elementsList->elementAt(j));
			if (!t) continue;
			s2= intToStr(t->getID());
			s = s + s2 + " ";			
		}
		st->Add(s);
	}

	//delete st;
	st->Add((char*)("FACES"));
	s = intToString(m->fFaces->Count());
	st->Add(s);

	for (i=0 ; i<m->fFaces->Count();i++)
	{
		tr = (TTriangle*)(m->fFaces->elementAt(i));
		s = intToString(tr->getID());
		s = s + " ";

		for (j = 0 ; j<3 ; j++)
		{
			v = tr->vertexes[j];
			if (!v) continue;
			s2 = intToStr(v->getID());
			s = s + s2 + " ";
		}
		st->Add(s);
	}
	// VERTEXES
	st->Add((char*)("VERTEXES"));
	s = intToString(m->vertexes->Count());
	st->Add(s);

	for (i=0 ; i<m->vertexes->Count();i++)
	{
		v = (TVertex*)(m->vertexes->elementAt(i));
		s = intToString(v->getID())+ " " + floatToString(v->fPos.x) + " "  
			+ floatToString(v->fPos.y)+ " "  +floatToString(v->fPos.z);

		st->Add(s);
	}
	// ELEMENTS
	st->Add((char*)("ELEMENTS"));
	s = intToString(m->elements->Count());
	st->Add(s);

	for (i=0 ; i<m->elements->Count();i++)
	{
		t = (TTetra*)(m->elements->elementAt(i));
		s = intToString(t->getID()) + " ";


		for (j = 0 ; j<4 ; j++)
		{
			v = t->vertexes[j];
			if (!v) continue;
			s = s + intToString(v->getID()) +" ";			
		}
		st->Add(s);
	}

	// VERTEXES x VERTEXES ELEMENTS NEIGHBOURS
	st->Add((char*)("VELEMENTS"));
	s = intToString(m->vertexes->Count());
	st->Add(s);

	TList<TObject*> *lneigh = new TList<TObject*>();

	for (i=0 ; i<m->vertexes->Count();i++)
	{
		v = m->vertexes->elementAt(i);
		s = intToString(v->getID()) + " ";

		//v->getVertexNeighboursByElem(lneigh);    

		for (j = 0 ; j<lneigh->Count() ; j++)
		{
			TVertex* v3 = (TVertex*)(lneigh->elementAt(j));
			if (!v3) continue;
			s2= intToString(v3->getID());
			s = s + s2 + " ";	        
		}
		st->Add(s);
	}


	st->saveToFile(aMeshName);
	return true;
	//st.SaveToFile(fdir+fn+'.NEIGH');
	// .........................
}


bool TGIDLoad::save(const char* aMeshName, TMesh* aMesh, int flags  )
{
	return false;
}

TMesh* TGIDLoad::load(const char* aMeshName)
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
		float x = 0,y = 0,z = 0;
		std::string s = st->strings[i];

		TVertex* v = new TVertex(x,y,z);
		v->setID( i );
		m->addVertex(v);
	}

	// Find Elements position
	/*while ( strcmp( line , "Begin Elements" )== 0) fread (line,14,1,fMesh); 


	while ( strcmp( line , "End Elements" )== 0)
	{
	int iv0,iv1,iv2,iv3;
	// id element
	fscanf(fMesh, "%d", &iv0);
	// other
	fscanf(fMesh, "%d", &iv0);
	// elements
	fscanf(fMesh, "%d", &iv0);
	fscanf(fMesh, "%d", &iv1);
	fscanf(fMesh, "%d", &iv2);			
	fscanf(fMesh, "%d", &iv3);	

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

