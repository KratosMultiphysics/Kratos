/****************************************************************************

  A simple GLUT program using the GLUI User Interface Library

  This program sets up a checkbox and a spinner, both with live variables.
  No callbacks are used.

  -----------------------------------------------------------------------
  9/17/04 John Kew - Natural Selection, Inc. (jkew@natural-selection.com)   
  9/9/98 Paul Rademacher (rademach@cs.unc.edu)

****************************************************************************/
#include "freeglut\include\GL\glut.h"

#include <gl\GL.h>

#include <GL\glui.h>

/** These are the live variables passed into GLUI ***/
int main_window;
//int num_display  = 0;
//int num_format  = 0;
//int enable_textbox = 0;
GLUI_StaticText  *text;
GLUI             *tree;
GLUI_EditText    *bedit;



void initMenues()
{

  
  
}