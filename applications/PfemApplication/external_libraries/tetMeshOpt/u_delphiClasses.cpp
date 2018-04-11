#include "u_delphiClasses.h"
#include <algorithm>
#include <vector>
#include <cstdio>

using namespace std;

#define CEILING_POS(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))
#define CEILING_NEG(X) ((X-(int)(X)) < 0 ? (int)(X-1) : (int)(X))
#define CEILING(X) ( ((X) > 0) ? CEILING_POS(X) : CEILING_NEG(X) )

TObject::TObject()
{
	// Nothing
}
void TObject::free()
{
	delete(this);
}

TObject::~TObject()
{
}

std::string intToString(int number)
{
	char* c = intToStr(number);

	return std::string(c);
}

std::string floatToString(float number)
{
	char* c = floatToStr(number);

	return std::string(c);
}

char* intToStr(int i)
{
	char* result = new char[100];       
	sprintf( result, "%d", i );
	return result ;
}

char* intToStr(int i, char* result)
{	
	sprintf( result, "%d", i );
	return result ;
}

char* floatToStr(float f)
{
	char* result = new char[100];       
	sprintf( result, "%f", f );
	return result ;
}

int dround(double d)
{
	return CEILING(d);
}


char* floatToStr(float f, char* result)
{	
	sprintf( result, "%f", f );
	return result ;
}


void freeAndNil(TList<TObject*>* l)
{	 
	for (int i=0 ; i<l->Count() ; i++)
	{
		TObject *o = l->structure[i];
		delete o;
	}
	delete(l);
}

