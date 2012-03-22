// PrSpain.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "u_Types.h"
#include "u_Render.h"
#include "u_MeshLoaders.h"
#include "u_elementCluster.h"
#include "u_FormRender.h"
#include "u_TetraFunctions.h"

#include <gl\GLU.h>

#include <vector>

int main (int argc, char* argv[])
{
	//Create a std::vector of the six first primes
/*	TList<TTetra*>* tl = new TList<TTetra*>();
	for (int i=0;i<50;i++)
		if (i%2 == 0)
			tl->Add(new TTetra(NULL));
		else
			 tl->Add(NULL);
	tl->Pack();
 */
   TMeshLoader* ml = new TVMWLoader();
   TMesh* m = ml->load("D:/posDOC/cube_sphere2.vwm");
   m->normalizeMesh(1000);
  
   // init GLUT and create Window
    glutInit(&argc, argv);
	initRender((TVolumeMesh*)(m));
	addMeshToRender( (TObject*)(m));
	initForm();
	return 0;
} 
