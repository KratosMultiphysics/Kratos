#include "stdafx.h"
#include "apiPrimitives.h"
#include <gl\GL.h>
#include <gl\GLU.h>
#include <gl\glut.h>
#include "Math3D.h"

    double gbmodelMatrix[16];
	double gbprojMatrix[16];
    int gbviewport[4];

	

void ieStoreViewMatrix()
{
	glGetIntegerv(GL_VIEWPORT, gbviewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, gbmodelMatrix); // Load the matricies and viewport.
    glGetDoublev(GL_PROJECTION_MATRIX, gbprojMatrix);

}


bool iePointOnScreen(float4 pos, double*sx, double *sy, double *sz, bool useStoredMatrixes)
{
	float depth;
    int scrW;
    double *modelMatrix;
	double *projMatrix;
    int *viewport;
	bool result;
	if (useStoredMatrixes)
	{
	  viewport = gbviewport;
	  gbviewport[1] = 0;
	  modelMatrix = gbmodelMatrix;
	  projMatrix = gbprojMatrix;
	}
	else
	{
		viewport = new int[4];
		modelMatrix = new double[16];
		projMatrix = new double[16];
	glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix); // Load the matricies and viewport.
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	}
	scrW = viewport[2];

    gluProject(pos.x,pos.y,pos.z,modelMatrix,projMatrix,viewport,sx,sy,sz); // Find out where the light is on the screen.
    result = true;
	
    if ((*sx < viewport[2]) && (sx >= 0) &&
        (*sy < (viewport[3]+scrW/18)) && (*sy >= 20))
	{
   // El problema con el glReadPixels es que hace esperar al CPU
       glReadPixels(*sx,*sy,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&depth);
       if(depth < *sz)  // but it is behind something.
         result = false; // The light can't be seen.
	}
    else
	{
		// If the point isn't on the screen
      result = false; // The point can't be seen.
	}
	return result;
}
