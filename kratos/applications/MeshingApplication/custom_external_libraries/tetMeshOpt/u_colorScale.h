 #pragma once
#include "Math3D.h"
#include "u_delphiClasses.h"
#include "u_Types.h"

float4 TColorToRealColor(int ct);
int TVector4fToTColor(float4 aV );


class TColorValue 
{
public :
     float4 realColor;
     float pos;
	 int color;
     TColorValue( float4 rc, float p)
	 {
		 this->realColor = rc;
		 this->pos = p;
		 this->color = TVector4fToTColor(rc);
	 }
};

class TColorPalette
{
public :
   
   TList<TColorValue*>* colors ;
   TColorPalette()
   {
	    colors = new TList<TColorValue*>();
		colors->Add( new TColorValue( Float4(0.807843f,0.807843f,0.807843f),0.0f) );
		colors->Add( new TColorValue( Float4(0.75294f,0.75294f,0.75294f),0.1924f) );
		colors->Add( new TColorValue( Float4(0.33725f,0.33725f,0.25490f),0.446735f) );
		colors->Add( new TColorValue( Float4(0.50196f,0.50196f,0.0f),0.646446735f) );
		colors->Add( new TColorValue( Float4(0.95686274767f,0.97254f,0.42350f),0.81786f) );
		colors->Add( new TColorValue( Float4(0.75294f,0.862745f,0.7529411f),1.0f) );
   }

   float4 getColorValues(float p  )
   {
 
     int i;
     TColorValue* befValue;
	 TColorValue* nValue;
	 p = (float) Max(0,Min(1,p));
     float prop,distance;
 
      befValue= NULL;
	  nValue = colors->elementAt(0);
	  for (i = 1 ; i<colors->Count(); i++)
	  {
          TColorValue* tc = colors->elementAt(i);
		  
		  if (p<tc->pos)
		  {
            befValue = colors->elementAt(i-1);
            nValue =colors->elementAt(i);
            break;
		  }

	  }
     if (befValue == NULL)     
       return nValue->realColor  ;

     distance = nValue->pos - befValue->pos;
     prop = (p - befValue->pos) / distance;
     return befValue->realColor* (1-prop) + nValue->realColor*(prop);
     
   }

   void colorizeMesh(TList<TObject*>* elList , fncEvaluation fn  )
   {
	   double minV = 50000000;
	   double maxV = -500000;
	   for (int i=0; i<elList->Count() ; i++)
	   {
		   double value = fn(elList->elementAt(i));
		   minV = Min(value, minV);
		   maxV = Max(value, maxV);		   
	   }

	   for (int i=0; i<elList->Count() ; i++)
	   {
		   TTetra* el = (TTetra*)elList->elementAt(i);
		   double value = fn(el);
		   float4 c =this->getColorValues( (float)((value - minV)/(maxV - minV)));
		   el->rColor = c;
		   
	   }
   }
};

//TColorPalette* instance = NULL;   
