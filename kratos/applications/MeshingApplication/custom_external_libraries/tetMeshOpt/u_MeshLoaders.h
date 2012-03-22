

#include <iostream>
#include <string>
#include "u_Types.h"


class TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TMeshLoader(){ return; }
	~TMeshLoader(){return; }

	virtual bool save(char* aMeshName, TMesh* aMesh ){ return false;  }

	virtual  TMesh* load(char* aMeshName) { return nil; }
};


class TSurLoader : public TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TSurLoader(){ return; }
	~TSurLoader(){return; }
	virtual bool save(char* filename, TMesh* aMesh );

	virtual TMesh* load(char* aMeshName)
	{ return nil; }
};

class TVMWLoader : public TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TVMWLoader(){ return; }
	~TVMWLoader(){return; }

	virtual bool save(char* aMeshName, TMesh* aMesh );

	virtual TMesh* load(char* aMeshName);

};

class TElementTetraLoader : public TMeshLoader
{
public :
	virtual TMesh* load(char* aMeshName){ return NULL;};
	virtual bool save(char* aMeshName, TMesh* aMesh );

};

class TGIDLoad : public TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TGIDLoad(){ return; }
	~TGIDLoad(){return; }

	virtual bool save(char* aMeshName, TMesh* aMesh );

	virtual TMesh* load(char* aMeshName);

};
