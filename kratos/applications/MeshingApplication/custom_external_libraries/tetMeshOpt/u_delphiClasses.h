#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdio>

#define nil  NULL; 

#define laAppend 10000
#define laOr   20000

char* intToStr(int i);
char* intToStr(int i, char* result);
char* floatToStr(float f);
char* floatToStr(float f, char* result);
int dround(double d);

std::string intToString(int number);
std::string floatToString(float number);

class TObject
{
public :	
	TObject();
	void free();
	~TObject();
};

typedef bool (*SortFunction) (TObject *a, TObject *b);





template <class T>
class TList 
{
public :
	bool ommitNULLSwhenAdd ;
	bool keepOrder ;
	std::vector<T> structure;
	T operator[](int v1) const
	{ 
		return structure[v1]; 
	}

	TList(bool ommitNULLSwhenAdd = true)  
	{  
		this->ommitNULLSwhenAdd = ommitNULLSwhenAdd;
		this->keepOrder = false;
	}
	~TList()
	{
		
		this->Clear();
		
	}

	void PackKeepingOrder()
	{
		int sSize = structure.size();
		int lIndex =sSize-1;
		
		for (int i=0; i<=lIndex ; i++)
		{
			if (structure[i] != NULL ) continue;
			int nextNotNull = -1;
			for (int j=i+1 ; j<=lIndex ; j++)
				if (structure[j] != NULL)
			    {
					nextNotNull = j;
					break;
				}
			if (nextNotNull >0)
			{
				structure[i] = structure[nextNotNull] ;
				structure[nextNotNull]  = NULL;				
			}
		}
		int lastIndex = 0;
		for ( lastIndex=0; lastIndex<=lIndex ; lastIndex++)
		{
			if (structure[lastIndex] == NULL) 
				break;
		}
		structure.resize(lastIndex);	  
	}

	void Pack() 
	{ 
		if (keepOrder)
		  return PackKeepingOrder();
		int sSize = structure.size();
		int lastIndex =sSize-1;
		while ((lastIndex>0) &&(structure[lastIndex] == NULL)  ) lastIndex--;
		int lIndex = lastIndex;

		for (int i=0; i<=lIndex ; i++)
		{
			if (i>=lastIndex) break;
			if (structure[i] == NULL)
			{			  
				structure[i] = structure[lastIndex] ;
				structure[lastIndex]  = NULL;
				while ((lastIndex>=i) &&(structure[lastIndex] == NULL)  ) lastIndex--;
			}

		}

		structure.resize(lastIndex+1);	  
	}
	int Count() { return structure.size(); }
	void Add(T elem)
	{ 
		if ((ommitNULLSwhenAdd) && (elem == NULL))
			return ; 
		structure.push_back(elem); 
	}
	void Remove(T elem) 
	{ 
		structure.erase( std::find( structure.begin(), structure.end(), elem ) );
	}

	void Clear()
	{ 
		structure.clear();
	}

	void Extract(T elem) 
	{ 
		structure.erase( std::find( structure.begin(), structure.end(), elem ) );
	}

	void setElementAt(int index, T elem)
	{
		structure[index] = elem;
	}

	T elementAt(int index)
	{
		return structure[index];
	}

	void Sort(SortFunction sortMethod)
	{
		// using function as comp
		std::sort (structure.begin(), structure.end(), sortMethod); 	  
	}

	bool contains(T elem)
	{
		for (int i=0 ; i<structure.size() ; i++)
			if (	  structure[i] == elem)
				return true;
		return false;
	}

	int indexOf(T elem)
	{
		unsigned int pos = std::find(structure.begin(), structure.end(), elem) - structure.begin();
		if( pos < structure.size() )
			return pos;
		else 
			return -1;	  // Not found  

	}

	void Assign(TList* f)
	{
		structure.clear();
		structure.assign(f->structure.begin(),f->structure.end());
	}

	void Assign(TList* f, int assignMode  )
	{ 
		if (assignMode == laAppend) 
			structure.assign(f->structure.begin(),f->structure.end());
		else
		{

			//Agregar no repetidos
			if (assignMode == laOr)
			{
				for (unsigned int i=0; i<f->structure.size();i++)
				{ 
					if (this->indexOf(f->structure[i])>=0) continue;
					this->structure.push_back(f->structure[i]);

				}
			}

		}
	}

};



class TStringList 
{
public :
	std::vector<std::string> strings;
	std::vector<TObject*> objects;
	TStringList() { }
	void Add(std::string str) 
	{ 
		strings.push_back(str); 
	}

	void splitString(std::string s, std::vector<std::string> sout)
	{
		// char * pch ;
		// const char *Data = s.data();

	}

	void Add(char* cc) 
	{ 
		std::string str(cc);
		strings.push_back(str); 
	}

	void Add(char* cc, char* cc2) 
	{ 		
		std::string str(cc);
		std::string str2(cc2);

		strings.push_back(str + " " + str2); 
	}

	void Add(char* cc, char* cc2, char* cc3) 
	{ 		  
		std::string str(cc);
		std::string str2(cc2);
		std::string str3(cc3);

		strings.push_back(str + " " + str2+ " " + str3); 

	}

	void saveToFile(const char* filename) 
	{
		FILE* fMesh;
		fMesh = fopen(filename,"wb"); //xx = rb, wb, read and write binary, more 		
		for (unsigned int i=0; i<strings.size(); i++)
		{ 
			std::string st = strings[i];
			fprintf( fMesh, "%s%s", st.data(),"\n");			 
		}
		fclose( fMesh );
	}

	void loadFromFile(const char *filename)
	{

		std::ifstream myfile (filename);

		while ( myfile.good() )
		{
			std::string line;
			getline (myfile,line);
			this->Add(line);
		}
		myfile.close();

	}

} ;

void freeAndNil(TList<TObject*>* l) ;

