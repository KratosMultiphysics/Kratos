#include "u_delphiClasses.h"
#include <algorithm>
#include <vector>
using namespace std;



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
	return ceil(d);
}


char* floatToStr(float f, char* result)
{	
    sprintf( result, "%f", f );
	return result ;
}


void freeAndNil(TList<TObject*>* l)
{	 
	delete(l);
}
