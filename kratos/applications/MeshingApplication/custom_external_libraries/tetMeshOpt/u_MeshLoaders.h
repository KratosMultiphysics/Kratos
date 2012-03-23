

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

	virtual bool save(const char* aMeshName, TMesh* aMesh )
	{ 
	    std :: cout << "Abstract save method "<< "\n";
		return false;  
	}
	virtual bool save(std::string aMeshFileName, TMesh* aMesh ) 
	{ 
		return save(aMeshFileName.data() , aMesh); 
	}

	virtual  TMesh* load(const char* aMeshName) 
	{ 
		std :: cout << "Abstract load method "<< "\n";
		return nil; 
	}
	virtual  TMesh* load(std::string aMeshFileName)
	{ 
		return load(aMeshFileName.data()); 
	}
};


class TSurLoader : public TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TSurLoader(){ return; }
	~TSurLoader(){return; }
	virtual bool save(const char* filename, TMesh* aMesh );

	virtual TMesh* load(const char* aMeshName)
	{ 
		std :: cout << "Sur format Loader not implemented "<< "\n";
		return nil; 
	}
};

class TVMWLoader : public TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TVMWLoader(){ return; }
	~TVMWLoader(){return; }

	virtual bool save(const char* aMeshName, TMesh* aMesh );

	virtual TMesh* load(const char* aMeshName);

};

class TElementTetraLoader : public TMeshLoader
{
public :
	virtual TMesh* load(const char* aMeshName){ return NULL;};
	virtual bool save(const char* aMeshName, TMesh* aMesh );

};

class TGIDLoad : public TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TGIDLoad(){ return; }
	~TGIDLoad(){return; }

	virtual bool save(const char* aMeshName, TMesh* aMesh );
	

	virtual TMesh* load(const char* aMeshName);

};
