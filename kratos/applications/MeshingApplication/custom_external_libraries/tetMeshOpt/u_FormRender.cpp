/****************************************************************************
****************************************************************************/
#include <string.h>

#include "u_FormRender.h"
#include "stdafx.h"
#include <GL/glui.h>
#include "u_Render.h"
#include "u_TetraFunctions.h"

float xy_aspect;
int   last_x, last_y;
float rotationX = 0.0, rotationY = 0.0;

/** These are the live variables passed into GLUI ***/
int   wireframe = 0;
int   obj_type = 1;
int   segments = 8;
int   segments2 = 8;
int   light0_enabled = 1;
int   light1_enabled = 1;
float light0_intensity = 1.0f;
float light1_intensity = .4f;
int   main_window;
float scale = 1.0;
float elemenScale = 1.0f;
int   show_sphere=1;
int   show_torus=1;
int   show_axes = 1;
int   show_text = 1;
float sphere_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float torus_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float view_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float obj_pos[] = { 0.0, 0.0, 3000.0 };
char *string_list[] = { "Hello World!", "Foo", "Testing...", "Bounding box: on" };
int   curr_RenderMode = 0;


// Pointers to the windows and some of the controls we'll create 
GLUI *glui, *glui2;
GLUI_Spinner    *light0_spinner, *light1_spinner;
GLUI_RadioGroup *radio;
GLUI_Panel      *obj_panel;


//********** User IDs for callbacks ********
#define LIGHT0_ENABLED_ID    200
#define LIGHT1_ENABLED_ID    201
#define LIGHT0_INTENSITY_ID  250
#define LIGHT1_INTENSITY_ID  260
#define ENABLE_ID            300
#define DISABLE_ID           301
#define SHOW_ID              302
#define HIDE_ID              303


//********** Miscellaneous global variables **********
TList<TObject*>* meshesToRender ;

GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
GLfloat light0_position[] = {.5f, .5f, 1.0f, 0.0f};

GLfloat light1_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
GLfloat light1_diffuse[] =  {.9f, .6f, 0.0f, 1.0f};
GLfloat light1_position[] = {-1.0f, -1.0f, 1.0f, 0.0f};

GLfloat lights_rotation[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

void addMeshToRender(TObject* o)
{
	if (meshesToRender ==NULL) meshesToRender = new TList<TObject*>();
	meshesToRender->Add(o);
}
//**************************************** control_cb() *******************/
// GLUI control callback                                                 */

  void _myMotionFunc(int x, int y){}  
  void _myPassiveFunc(int x, int y){}  
  void _myEntryFunc(int x){}  
  void _myVisFunc(int x){}
  void _my_special_func(int key, int x, int y){}
  void _myGlutIdle( void ){}
  void _myGlutDisplay() {}
  void _myGlutReshape(int x, int y){}
  void _myGlutMotion(int x, int y){}
  void _myGlutKeyboard(unsigned char c , int x,int y){}
  void _myGlutMouse(int x, int y, int z, int w){}
  void _exitMenu(int x){}

  void innerInit()
  {
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition( 50, 50 );
  glutInitWindowSize( 300, 300 );
 
  main_window = glutCreateWindow( "GLUI Example 3" );
glutDisplayFunc( _myGlutDisplay );
  glutReshapeFunc( _myGlutReshape );  
  glutKeyboardFunc( _myGlutKeyboard );
  glutMotionFunc( _myGlutMotion );
  glutMouseFunc( _myGlutMouse );
   

	  // BUG!!
  
  glutSpecialFunc(_my_special_func);
  //glutReshapeWindow(300,300);

    glutPassiveMotionFunc( _myPassiveFunc );
  glutEntryFunc( _myEntryFunc );
  glutVisibilityFunc( _myVisFunc );
  glutBitmapWidth(NULL, 0);
  glutIdleFunc(_myGlutIdle);
  glutSetCursor(0);
  glutGet(0);
  glutBitmapCharacter(NULL,0);
  glutGetModifiers();
  glutSetWindow(0);
  glutHideWindow();
  glutShowWindow();
  glutReshapeWindow(0,0);
  glutCreateMenu(_exitMenu);
  glutCreateSubWindow(main_window,0,0,0,0);
  glutAddMenuEntry("as",0);
  
  glutAttachMenu(0);
  glutDetachMenu(0);
  glutDestroyMenu(1);
  glutPositionWindow(0,0);
  glutDestroyWindow(main_window);
  }



void control_cb( int control )
{
  if ( control == LIGHT0_ENABLED_ID ) {
    if ( light0_enabled ) {
      glEnable( GL_LIGHT0 );
      light0_spinner->enable();
    }
    else {
      glDisable( GL_LIGHT0 ); 
      light0_spinner->disable();
    }
  }
  else if ( control == LIGHT1_ENABLED_ID ) {
    if ( light1_enabled ) {
      glEnable( GL_LIGHT1 );
      light1_spinner->enable();
    }
    else {
      glDisable( GL_LIGHT1 ); 
      light1_spinner->disable();
    }
  }
  else if ( control == LIGHT0_INTENSITY_ID ) {
    float v[] = { 
      light0_diffuse[0],  light0_diffuse[1],
      light0_diffuse[2],  light0_diffuse[3] };
    
    v[0] *= light0_intensity;
    v[1] *= light0_intensity;
    v[2] *= light0_intensity;

    glLightfv(GL_LIGHT0, GL_DIFFUSE, v );
  }
  else if ( control == LIGHT1_INTENSITY_ID ) {
    float v[] = { 
      light1_diffuse[0],  light1_diffuse[1],
      light1_diffuse[2],  light1_diffuse[3] };
    
    v[0] *= light1_intensity;
    v[1] *= light1_intensity;
    v[2] *= light1_intensity;

    glLightfv(GL_LIGHT1, GL_DIFFUSE, v );
  }
  else if ( control == ENABLE_ID )
  {
    glui2->enable();
  }
  else if ( control == DISABLE_ID )
  {
    glui2->disable();
  }
  else if ( control == SHOW_ID )
  {
    glui2->show();
  }
  else if ( control == HIDE_ID )
  {
    glui2->hide();
  }
}

//**************************************** myGlutKeyboard() **********/

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
  case 27: 
  case 'q':
    exit(0);
    break;
  };
  
  glutPostRedisplay();
}


//***************************************** myGlutMenu() ***********/

void myGlutMenu( int value )
{
  myGlutKeyboard( value, 0, 0 );
}


//***************************************** myGlutIdle() ***********/

void myGlutIdle( void )
{
  /* According to the GLUT specification, the current window is 
     undefined during an idle callback.  So we need to explicitly change
     it if necessary */
  if ( glutGetWindow() != main_window ) 
    glutSetWindow(main_window);  

  /*  GLUI_Master.sync_live_all();  -- not needed - nothing to sync in this
                                       application  */

  glutPostRedisplay();
}

//***************************************** myGlutMouse() **********/

void myGlutMouse(int button, int button_state, int x, int y )
{
}


//***************************************** myGlutMotion() **********/

void myGlutMotion(int x, int y )
{
  glutPostRedisplay(); 
}

//**************************************** myGlutReshape() *************/

void myGlutReshape( int x, int y )
{
  int tx, ty, tw, th;
  GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
  glViewport( tx, ty, tw, th );

  xy_aspect = (float)tw / (float)th;

  glutPostRedisplay();
}


//************************************************** draw_axes() **********/
//* Disables lighting, then draws RGB axes                                */

void draw_axes( float scale )
{
  glDisable( GL_LIGHTING );

  glPushMatrix();
  glScalef( scale, scale, scale );

  glBegin( GL_LINES );
 
  glColor3f( 1.0, 0.0, 0.0 );
  glVertex3f( .8f, 0.05f, 0.0 );  glVertex3f( 1.0, 0.25f, 0.0 ); /* Letter X */
  glVertex3f( 0.8f, .25f, 0.0 );  glVertex3f( 1.0, 0.05f, 0.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 1.0, 0.0, 0.0 ); /* X axis      */

  glColor3f( 0.0, 1.0, 0.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 0.0, 1.0, 0.0 ); /* Y axis      */

  glColor3f( 0.0, 0.0, 1.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 0.0, 0.0, 1.0 ); /* Z axis    */
  glEnd();

  glPopMatrix();

  glEnable( GL_LIGHTING );
}


//***************************************** myGlutDisplay() *****************/

void myGlutDisplay( void )
{
  
  glClearColor( .9f, .9f, .9f, 1.0f );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum( -xy_aspect*.04, xy_aspect*.04, -.04, .04, .1, 15000.0 );

  glMatrixMode( GL_MODELVIEW );

  glLoadIdentity();
  glMultMatrixf( lights_rotation );
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  
  glLoadIdentity();
  glTranslatef( 0.0, 0.0, -2.6f );
  glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] ); 
  glMultMatrixf( view_rotate );

  glScalef( scale, scale, scale );
  
  ieStoreViewMatrix();
  
  /*** Now we render object, using the variables 'obj_type', 'segments', and
    'wireframe'.  These are _live_ variables, which are transparently 
    updated by GLUI ***/
  for (int i=0; i<meshesToRender->Count();i++)
  {
	  TVolumeMesh* vm = (TVolumeMesh*)(meshesToRender->elementAt(i));
	  renderMesh(vm,curr_RenderMode+1, wireframe == 1,elemenScale);
	  
  }
  
  glDisable( GL_LIGHTING );
  glAxis(Float4(0.0f),2);

  if ( show_text ) 
  {
    
	int tx, ty, tw, th;
    GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
	
	glDisable( GL_LIGHTING );  // Disable lighting while we render text 
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluOrtho2D( tx,tw,0,th );
    glMatrixMode( GL_MODELVIEW );
    
	glLoadIdentity();
    glColor3ub( 0, 0, 0 );
    

    //  printf( "text: %s\n", text );              

    //*** Render the live character array 'text' 
    int i;
	glRasterPos2i( 10,10);
	for ( i=0; i<meshesToRender->Count();i++)
    {
	  TVolumeMesh* vm = (TVolumeMesh*)(meshesToRender->elementAt(i));
	  char ids[50];

	  glDisable(GL_DEPTH_TEST);
	  
	  for (int j=0; j<vm->vertexes->Count();j++)
	  {
		  if (j>50) break;
		  TVertex* v = vm->vertexes->elementAt(j);
		  float4 vpos = v->fPos;
		  double sx,sy,sz;
	      //if (!iePointOnScreen(vpos,&sx,&sy,&sz,true)) continue;
		  if (!iePointOnScreen(vpos,&sx,&sy,&sz,true)) continue;
		  
		  sprintf(ids,"%d",v->id);
		  glRasterPos2d( sx,sy );
		  for( int k=0; k<(int)strlen( ids ); k++ )
		  {			 
		     glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, ids[k] );	
		  }


	  }
	  glEnable(GL_DEPTH_TEST);
    }
  }
  
  glEnable( GL_LIGHTING );
  glutSwapBuffers(); 
}


//**************************************** main() ********************/

void initForm()
{
  /****************************************/
  /*   Initialize GLUT and create window  */
  /****************************************/
  
  innerInit();

  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition( 50, 50 );
  glutInitWindowSize( 800, 600 );
 
  main_window = glutCreateWindow( "GLUI Example 5" );
  glutDisplayFunc( myGlutDisplay );
  
  GLUI_Master.set_glutReshapeFunc( myGlutReshape );  
  GLUI_Master.set_glutKeyboardFunc( myGlutKeyboard );
  GLUI_Master.set_glutSpecialFunc( NULL );
  GLUI_Master.set_glutMouseFunc( myGlutMouse );
  glutMotionFunc( myGlutMotion );

  /****************************************/
  /*       Set up OpenGL lights           */
  /****************************************/

  glEnable(GL_LIGHTING);
  glEnable( GL_NORMALIZE );

  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

  glEnable(GL_LIGHT1);
  glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

  /****************************************/
  /*          Enable z-buferring          */
  /****************************************/

  glEnable(GL_DEPTH_TEST);

  /****************************************/
  /*         Here's the GLUI code         */
  /****************************************/

  printf( "GLUI version: %3.2f\n", GLUI_Master.get_version() );

  /*** Create the side subwindow ***/
  glui = GLUI_Master.create_glui_subwindow( main_window, 
					    GLUI_SUBWINDOW_RIGHT );

  obj_panel = new GLUI_Rollout(glui, "Properties", false );

  /***** Control for object params *****/

  new GLUI_Checkbox( obj_panel, "Wireframe", &wireframe, 1, control_cb );

  GLUI_Spinner *scale_spinner = 
    new GLUI_Spinner( obj_panel, "Scale:", &scale);
  scale_spinner->set_float_limits( .2f, 4.0 );
  scale_spinner->set_alignment( GLUI_ALIGN_RIGHT );

  GLUI_Spinner *scaleElement_spinner = 
	  new GLUI_Spinner( obj_panel, "Elem Scale:", &elemenScale);
  scale_spinner->set_float_limits( .1f, 2.0 );
  scale_spinner->set_alignment( GLUI_ALIGN_RIGHT );

  /**** Add listbox ****/

  GLUI_Listbox *listModes = new GLUI_Listbox( obj_panel, "Modes:", &curr_RenderMode );    
  listModes->add_item( 0, "SURFACE" );
  listModes->add_item( 1, "VOLUME" );
  listModes->add_item( 2, "COLORIZE" );


  /******** Add some controls for lights ********/
  GLUI_Rollout *roll_lights = new GLUI_Rollout(glui, "Lights", false );

  GLUI_Panel *light0 = new GLUI_Panel( roll_lights, "Light 1" );
  GLUI_Panel *light1 = new GLUI_Panel( roll_lights, "Light 2" );

  new GLUI_Checkbox( light0, "Enabled", &light0_enabled,
                     LIGHT0_ENABLED_ID, control_cb );
  light0_spinner = 
    new GLUI_Spinner( light0, "Intensity:", 
                      &light0_intensity, LIGHT0_INTENSITY_ID,
                      control_cb );
  light0_spinner->set_float_limits( 0.0, 1.0 );
  GLUI_Scrollbar *sb;
  sb = new GLUI_Scrollbar( light0, "Red",GLUI_SCROLL_HORIZONTAL,
                           &light0_diffuse[0],LIGHT0_INTENSITY_ID,control_cb);
  sb->set_float_limits(0,1);
  sb = new GLUI_Scrollbar( light0, "Green",GLUI_SCROLL_HORIZONTAL,
                           &light0_diffuse[1],LIGHT0_INTENSITY_ID,control_cb);
  sb->set_float_limits(0,1);
  sb = new GLUI_Scrollbar( light0, "Blue",GLUI_SCROLL_HORIZONTAL,
                           &light0_diffuse[2],LIGHT0_INTENSITY_ID,control_cb);
  sb->set_float_limits(0,1);
  new GLUI_Checkbox( light1, "Enabled", &light1_enabled,
                     LIGHT1_ENABLED_ID, control_cb );
  light1_spinner = 
    new GLUI_Spinner( light1, "Intensity:",
                      &light1_intensity, LIGHT1_INTENSITY_ID,
                      control_cb );
  light1_spinner->set_float_limits( 0.0, 1.0 );
  sb = new GLUI_Scrollbar( light1, "Red",GLUI_SCROLL_HORIZONTAL,
                           &light1_diffuse[0],LIGHT1_INTENSITY_ID,control_cb);
  sb->set_float_limits(0,1);
  sb = new GLUI_Scrollbar( light1, "Green",GLUI_SCROLL_HORIZONTAL,
                           &light1_diffuse[1],LIGHT1_INTENSITY_ID,control_cb);
  sb->set_float_limits(0,1);
  sb = new GLUI_Scrollbar( light1, "Blue",GLUI_SCROLL_HORIZONTAL,
                           &light1_diffuse[2],LIGHT1_INTENSITY_ID,control_cb);
  sb->set_float_limits(0,1);


  /*** Add another rollout ***/
  GLUI_Rollout *options = new GLUI_Rollout(glui, "Options", true );
  new GLUI_Checkbox( options, "Draw Object", &show_sphere );
  new GLUI_Checkbox( options, "Draw text", &show_text );

  /**** Add listbox ****/
  new GLUI_StaticText( glui, "" );
  new GLUI_StaticText( glui, "" );


  /*** Disable/Enable buttons ***/
  new GLUI_Button( glui, "Disable movement", DISABLE_ID, control_cb );
  new GLUI_Button( glui, "Enable movement", ENABLE_ID, control_cb );
  new GLUI_Button( glui, "Hide", HIDE_ID, control_cb );
  new GLUI_Button( glui, "Show", SHOW_ID, control_cb );

  new GLUI_StaticText( glui, "" );

  /****** A 'quit' button *****/
  new GLUI_Button( glui, "Quit", 0,(GLUI_Update_CB)exit );


  /**** Link windows to GLUI, and register idle callback ******/
  
  glui->set_main_gfx_window( main_window );


  /*** Create the bottom subwindow ***/
  glui2 = GLUI_Master.create_glui_subwindow( main_window, 
                                             GLUI_SUBWINDOW_BOTTOM );
  glui2->set_main_gfx_window( main_window );

  GLUI_Rotation *view_rot = new GLUI_Rotation(glui2, "Objects", view_rotate );
  view_rot->set_spin( 1.0 );

  new GLUI_Column( glui2, false );
  GLUI_Rotation *lights_rot = new GLUI_Rotation(glui2, "Blue Light", lights_rotation );
  lights_rot->set_spin( .82f );
  new GLUI_Column( glui2, false );
  GLUI_Translation *trans_xy = 
    new GLUI_Translation(glui2, "Objects XY", GLUI_TRANSLATION_XY, obj_pos );
  trans_xy->set_speed( .5 );
  new GLUI_Column( glui2, false );
  GLUI_Translation *trans_x = 
    new GLUI_Translation(glui2, "Objects X", GLUI_TRANSLATION_X, obj_pos );
  trans_x->set_speed( .5 );
  new GLUI_Column( glui2, false );
  GLUI_Translation *trans_y = 
    new GLUI_Translation( glui2, "Objects Y", GLUI_TRANSLATION_Y, &obj_pos[1] );
  trans_y->set_speed( .5 );
  new GLUI_Column( glui2, false );
  GLUI_Translation *trans_z = 
    new GLUI_Translation( glui2, "Objects Z", GLUI_TRANSLATION_Z, &obj_pos[2] );
  trans_z->set_speed( .5 );

#if 0
  /**** We register the idle callback with GLUI, *not* with GLUT ****/
  GLUI_Master.set_glutIdleFunc( myGlutIdle );
#endif

  /**** Regular GLUT main loop ****/
  
  glutMainLoop();

  return ;
}

