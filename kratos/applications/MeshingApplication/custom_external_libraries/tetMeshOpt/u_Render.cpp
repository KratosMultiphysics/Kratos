#include "stdafx.h"
#include "u_Render.h"
#include "u_TetraFunctions.h"

// angle of rotation for the camera direction
float angle=0.0;
// actual vector representing the camera's direction
float lx=0.0f,lz=-1.0f;
// XZ position of the camera
float x=0.0f,z=5000.0f;

float deltaAngle = 0.0f;
float deltaAngleY = 0.0f;
int xOrigin = -1;
int yOrigin = -1;
float deltaMoveX = 0;
float deltaMoveY = 0;
float deltaStep = 100.0f;


int renderMode ;
TColorPalette* colorPaletteInstance;
TVolumeMesh* meshToRender;

void colorizeMesh()
{
	//TColorPalette->instance->colorizeMesh(meshToRender , relaxQuality);
	if (colorPaletteInstance == NULL)
		  colorPaletteInstance = new TColorPalette();
	colorPaletteInstance->colorizeMesh(meshToRender->elements , relaxQuality);
}

void processMenuEvents(int option) 
{
//Render_Properties* frm ;
if (option<=2)
	renderMode = option;
else
 if (option == 3)
  colorizeMesh();
 else
	 if (option == 5)
	 {
		 TetQuality* tq = new TetQuality(meshToRender);
		 tq->refresh();
		 tq->print();
	 }
 //else
	// frm = new Render_Properties();
}

void processMainMenu(int option) {

	// nothing to do in here
	// all actions are for submenus
}

void processFillMenu(int option) {

	switch (option) {

		case 0: glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); break;
		case 1: glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); break;
		
	}
}

void createGLUTMenus() {

	int menu;

	// create the menu and
	// tell glut that "processMenuEvents" will
	// handle the events
	int fillMenu = glutCreateMenu(processFillMenu);
	glutAddMenuEntry("Fill",0);
	glutAddMenuEntry("Line",1);

	menu = glutCreateMenu(processMenuEvents);

	//add entries to our menu
	glutAddMenuEntry("Surface",SURFACE);
	glutAddMenuEntry("Volume",VOLUME);	
	
	glutAddMenuEntry("Colorize",COLORIZE);
	glutAddMenuEntry("Properties",4);
	glutAddMenuEntry("Metrics",5);

	glutAddSubMenu("Polygon Mode", fillMenu);

	// attach the menu to the right button
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}



//- RENDERING

void ieLine(float4 fPoint,float4 lPoint, int width)
{
       glLineWidth( (float)width);
       glBegin(GL_LINES);
          glVertex3f(fPoint.x,fPoint.y,fPoint.z);
          glVertex3f(lPoint.x,lPoint.y,lPoint.z);
       glEnd();
       glLineWidth(1);
}

void glAxis(float4 center, int width)
{
   glColor3f(1.0,0.0,0.0);
   ieLine(center,(Float4)(3500.0,0.0,0.0),width);
   glColor3f(0.0,1.0,0.0);
   ieLine(center,(Float4)(0.0,3500.0,0.0),width);
   glColor3f(0,0,1);
   ieLine(center,(Float4)(0.0,0.0,3500.0),width);
}

void drawTriangle(TTriangle* tr, float4 color, long fillMode,float translateSc)
{
	int i;
	float4 tOff;
	
	  tOff = tr->normal*translateSc;
      glColor3f(color.x,color.y,color.z);

	  float4 vn =  Normal(tr->vertexes[0]->fPos,tr->vertexes[1]->fPos,tr->vertexes[2]->fPos);
        
      glBegin(GL_TRIANGLES);
       glNormal3f(vn.x,vn.y,vn.z);
	   for (i=0; i<3;i++)       
	      glVertex3f(tr->vertexes[i]->fPos.x+tOff.x,
                     tr->vertexes[i]->fPos.y+tOff.y,
                     tr->vertexes[i]->fPos.z+tOff.z);
      glEnd();
      glColor3f(1,0,0);

}

void drawTriangle(float4 v1,float4 v2,float4 v3,float4 color,long fillMode, float translateSc, bool wNormal )
{
 float4 vn;

  glColor3f(color.x,color.y,color.z);


   glPushMatrix();
      if (wNormal)
	  {
		  vn = Normal(v1,v2,v3);
          glNormal3f(vn.x,vn.y,vn.z);
	  }
      glTranslatef(vn.x*translateSc,vn.y*translateSc,vn.z*translateSc);


      glBegin(GL_TRIANGLES);

      glVertex3d(v1.x,v1.y,v1.z);
      glVertex3d(v2.x,v2.y,v2.z);
      glVertex3d(v3.x,v3.y,v3.z);
      glEnd();
     glPopMatrix();
  glColor3f(1,0,0);

}

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;

	float ratio =  w * 1.0f / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45.0f, ratio, 0.1f, 10000.0f);		

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

void processNormalKeys(unsigned char key, int xx, int yy) { 	

        if (key == 27)
              exit(0);
} 

void pressKey(int key, int xx, int yy) {

       switch (key) {
		   case GLUT_KEY_LEFT : deltaMoveX = -deltaStep; break;
             case GLUT_KEY_RIGHT : deltaMoveX = -deltaStep; break;
             case GLUT_KEY_UP : deltaMoveY = deltaStep; break;
             case GLUT_KEY_DOWN : deltaMoveY = -deltaStep; break;
       }
} 

void releaseKey(int key, int x, int y) { 	

        switch (key) {
			  case GLUT_KEY_LEFT :
             case GLUT_KEY_RIGHT : deltaMoveX = 0;break;

             case GLUT_KEY_UP :
             case GLUT_KEY_DOWN : deltaMoveY = 0;break;
        }
} 

void mouseMove(int x, int y) { 	

         // this will only be true when the left button is down
         if (xOrigin >= 0) {

		// update deltaAngle
		deltaAngle += (x - xOrigin) * 0.005f;
		deltaAngleY += (y - yOrigin) * 0.005f;


		// update camera's direction
		//lx = sin(angle + deltaAngle);
		//lz = -cos(angle + deltaAngle);
	}
}

void mouseButton(int button, int state, int x, int y) {

	// only start motion if the left button is pressed
	if (button == GLUT_LEFT_BUTTON) {

		// when the button is released
		if (state == GLUT_UP) {
			angle += deltaAngle;
			xOrigin = -1;
			yOrigin = -1;
		}
		else  {// state = GLUT_DOWN
			xOrigin = x;
			yOrigin = y;
		}
	}
}

void computePos(float dX,float dY) {

	x += dX * lx * 0.1f;
	z += dY * lz * 0.1f;
}


void drawTetraElement(TElement* e, float elemDisplacement,float4 meshCenter, float4 color,float elemScale  )
{
 float4 v0,v1,v2,v3, objCn ,vn,objDs;
  //------------    
      v0= e->vertexes[0]->fPos;    v1= e->vertexes[1]->fPos;
      v2= e->vertexes[2]->fPos;    v3= e->vertexes[3]->fPos;

	  if (!e->initialized )
	  {       
        e->normals[0] = Normal(v0,v1,v2)   ;
        e->normals[1] = Normal(v0,v2,v3)   ;
        e->normals[2] = Normal(v0,v3,v1)   ;
        e->normals[3] = Normal(v1,v3,v2)   ;
		e->initialized = true;
	  }
      
      glColor3f(color.x,color.y,color.z);
       if (elemDisplacement>0)
	   {
        objCn =(v0+v1+v2+v3)*0.25;
        objDs = objCn - meshCenter ;
		glPushMatrix();
        glTranslatef(objDs.x *elemDisplacement ,objDs.y *elemDisplacement,objDs.z *elemDisplacement);
	   }

       if (elemScale!=1.0f)
	   {
        objCn =(v0+v1+v2+v3)*0.25;
        v0 =  (v0-objCn)*elemScale+objCn;
        v1 =  (v1-objCn)*elemScale+objCn;
        v2 =  (v2-objCn)*elemScale+objCn;
        v3 =  (v3-objCn)*elemScale+objCn;
	   }
        glBegin(GL_TRIANGLES);
          vn =e->normals[0];
          glNormal3d(vn.x,vn.y,vn.z);
          glVertex3f(v0.x,v0.y,v0.z) ;
          glVertex3f(v2.x,v2.y,v2.z) ;
          glVertex3f(v1.x,v1.y,v1.z) ;

          vn =e->normals[1] ;
          glNormal3d(vn.x,vn.y,vn.z);
          glVertex3f(v0.x,v0.y,v0.z) ;
          glVertex3f(v3.x,v3.y,v3.z) ;
          glVertex3f(v2.x,v2.y,v2.z) ;

          vn =e->normals[2]  ;
          glNormal3d(vn.x,vn.y,vn.z);
          glVertex3f(v0.x,v0.y,v0.z) ;
          glVertex3f(v1.x,v1.y,v1.z) ;
          glVertex3f(v3.x,v3.y,v3.z) ;

          vn =e->normals[3]  ;
          glNormal3d(vn.x,vn.y,vn.z);
          glVertex3f(v1.x,v1.y,v1.z) ;
          glVertex3f(v2.x,v2.y,v2.z) ;
          glVertex3f(v3.x,v3.y,v3.z) ;
       glEnd();
     if (elemDisplacement>0)
       glPopMatrix();
}


void initRender(TMesh *m)
{
	meshToRender =(TVolumeMesh*)(m);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
//	glEnable(GL_CULL_FACE);
	// Somewhere in the initialization part of your program
glEnable(GL_LIGHTING);
glEnable(GL_LIGHT0);

// Create light components
float ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
float diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
float specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
float position[] = { -1.5f, 1.0f, -4.0f, 1.0f };

// Assign created components to GL_LIGHT0
glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
glLightfv(GL_LIGHT0, GL_POSITION, position);
}

void renderMesh(TMesh *mv, int renderMode, bool wireframe, float elementSize)
{
	TVolumeMesh *m = (TVolumeMesh*)(mv);
	if (wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	if (m)
    {
	  if (renderMode == SURFACE)
	  {
		 for (int i=0; i<m->fFaces->Count();i++)
	   {
		 TTriangle* tr =(TTriangle*)( m->fFaces->elementAt(i));
		 float4 cn =Float4(1.0f,0.0f,1.0f);
		 drawTriangle(tr,cn,GL_LINE,0.0f);
		 
	   }
	  }
	  else      
	  {
		for (int i=0; i<m->elements->Count();i++)
	   {
		 TTetra* t = (TTetra*)m->elements->elementAt(i);
		 if (t == NULL ) continue;
		 float4 cn =Float4(1.0f);
		 float4 elColor = t->rColor;
		 
		 drawTetraElement(t,0.0f,cn,elColor,elementSize);
	   }
	  }
	}
  
   

}

void display() 
{

	if (deltaMoveX || deltaMoveY)
		computePos(deltaMoveX,deltaMoveY);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);   // Clear the color buffer
	glLoadIdentity();
	// Set the camera
	// Set the camera
	gluLookAt(	x, 500.0f, z,
			x+lx, 500.0f,  z+lz,
			0.0f, 1.0f,  0.0f);
	
	glAxis(Float4(0.0f),2);
	
	glRotatef(deltaAngleY, 1,0,0);
	glRotatef(deltaAngle, 0,1,0);

	renderMesh(meshToRender,  renderMode,false,1);
	
    
   glutSwapBuffers();  // Render now
  
}
