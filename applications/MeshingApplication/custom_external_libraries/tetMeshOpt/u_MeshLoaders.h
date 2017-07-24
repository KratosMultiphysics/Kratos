

#include <iostream>
#include <string>
#include "u_Types.h"


class TMeshLoader
{
private :
	//char* fLastLoaded; //commented out to avoid warning (MA Celigueta 20-5-2016)
public :
	TMeshLoader(){ return; }
	virtual ~TMeshLoader(){return; }

	virtual bool save(const char* aMeshName, TMesh* aMesh, int flags = 0 )
	{ 
	    std :: cout << "Abstract save method "<< "\n";
		return false;  
	}
	virtual bool save(std::string aMeshFileName, TMesh* aMesh, int flags = 0 ) 
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
	~TSurLoader() override{return; }
	bool save(const char* filename, TMesh* aMesh , int flags = 0) override;

	TMesh* load(const char* aMeshName) override
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
	~TVMWLoader() override{return; }

	bool save(const char* aMeshName, TMesh* aMesh , int flags = 0) override;

	TMesh* load(const char* aMeshName) override;

};

class TElementTetraLoader : public TMeshLoader
{
public :
	TMesh* load(const char* aMeshName) override{ return nullptr;};
	bool save(const char* aMeshName, TMesh* aMesh , int flags = 0) override;

};

class TGIDLoad : public TMeshLoader
{
private :
	char* fLastLoaded;
public :
	TGIDLoad(){ return; }
	~TGIDLoad() override{return; }

	bool save(const char* aMeshName, TMesh* aMesh , int flags = 0) override;
	

	TMesh* load(const char* aMeshName) override;

};
