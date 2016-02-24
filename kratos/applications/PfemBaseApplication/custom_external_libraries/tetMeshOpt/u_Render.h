
#include "stdafx.h"
#include "freeglut\include\GL\glut.h"
#include <gl\GL.h>

#include "u_colorScale.h"
#include "u_qualityMetrics.h"
#include "u_ShowMetrics.h"
#include "Math3D.h"
#include "apiPrimitives.h"

void colorizeMesh();

void processMenuEvents(int option) ;
void processMainMenu(int option) ;

void createGLUTMenus() ;


//- RENDERING

void ieLine(float4 fPoint,float4 lPoint, int width);

void glAxis(float4 center, int width);

void drawTriangle(TTriangle* tr, float4 color, long fillMode,float translateSc);

void drawTriangle(float4 v1,float4 v2,float4 v3,float4 color,long fillMode = GL_FILL, float translateSc=0.0f, bool wNormal= false );
void changeSize(int w, int h) ;

void processNormalKeys(unsigned char key, int xx, int yy) ;

void pressKey(int key, int xx, int yy) ;

void releaseKey(int key, int x, int y) ;

void mouseMove(int x, int y) ;
void mouseButton(int button, int state, int x, int y) ;

void computePos(float dX,float dY) ;

void drawTetraElement(TElement* e, float elemDisplacement,float4 meshCenter, float4 color,float elemScale  );
void initRender(TMesh *m);

void display() ;
void renderMesh(TMesh *mv, int renderMode, bool wireframe, float elementSize);
void glAxis(float4 center, int width);
