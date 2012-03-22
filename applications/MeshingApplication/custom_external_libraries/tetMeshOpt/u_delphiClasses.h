#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#define nil  NULL; 

enum listOperations { laAppend, laOr   };

char* intToStr(int i);
char* intToStr(int i, char* result);
char* floatToStr(float f);
char* floatToStr(float f, char* result);
int dround(double d);

class TObject
{
public :	
	TObject(){}
	void free()
	{
		delete(this);
	}
};

typedef bool (*SortFunction) (TObject *a, TObject *b);





class u_delphiClasses
{
public:
	u_delphiClasses(void){};
	~u_delphiClasses(void){};
};



template <class T>
class TList 
{
public :
  std::vector<T> structure;
  	T operator[](int v1) const
	{ 
		 return structure[v1]; 
	}

  TList()  {  }
	~TList()
	{
		this->Clear();
	}

  void Pack() 
  { 
	  	  
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
  void Add(T elem){ structure.push_back(elem); }
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
      sort (structure.begin(), structure.end(), sortMethod); 	  
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
	  unsigned int pos = find(structure.begin(), structure.end(), elem) - structure.begin();
	  if( pos < structure.size() )
		   return pos;
	  else 
		  return -1;	  // Not found  

  }

  void Assign(TList* f)
  {
	  structure.assign(f->structure.begin(),f->structure.end());
  }

  void Assign(TList* f, listOperations assignMode  )
  { 
	 if (assignMode == laAppend) 
	  structure.assign(f->structure.begin(),f->structure.end());
	 else
	 {
		 
		 //Agregar no repetidos
		 if (assignMode == laOr)
		 {
			 for (int i=0; i<f->structure.size();i++)
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
		  strcat(cc," "); 
		  strcat(cc,cc2);
		  std::string str(cc);
		  
		  strings.push_back(str); 
	 }

	 void Add(char* cc, char* cc2, char* cc3) 
	 { 
		  strcat(cc," "); 
		  strcat(cc,cc2);
		  strcat(cc," "); 
		  strcat(cc,cc3);
		  std::string str(cc);
		  
		  strings.push_back(str); 
	 }

	 void saveToFile(char* filename) 
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

	  void loadFromFile(char *filename)
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


/*
template <class T>
class TList 
{
private:

	class Node
	{
	public:
		T info;
		Node *prev, *next;
	public:
		Node(T info): info(info) { prev = next = NULL; }
	} *head, *tail;

	
	int lastIndex ;
	Node* lastNode ;

private:
	void Prepend(Node *p)
	{
		if (this->head == NULL)
		{
			this->tail = p;
		}
		else
		{
			p->next = this->head;
			this->head->prev = p;
		}
		this->head = p;
		count++;
	}

	
	void Append(Node *p)
	{
		if (this->tail == NULL)
		{
			this->head = p;
		}
		else
		{
			p->prev = tail;
			this->tail->next = p;
		}
		this->tail = p;
		count++;
	}

	void InsertBefore(Node *q, Node *p)
	{
		p->prev = q->prev;
		p->next = q;
		q->prev = p;
		count++;
	}

	void InsertAfter(Node *q, Node *p)
	{
		p->next = q->next;
		p->prev = q;
		q->next = p;
		count++;
	}

	void Remove(Node *p)
	{
		Node *next = p->next;
		Node *prev = p->prev;

		if (p == this->head) this->head = next;
		if (p == this->tail) this->tail = prev;
		delete p;

		if (next != NULL) next->prev = prev;
		if (prev != NULL) prev->next = next;
		count--;
	}

public:
	int count;

	TList() : count(0) { head = tail = NULL; }
	~TList()
	{
		this->Clear();
	}

public:
	bool IsEmpty() { return this->head == NULL; }
	void Clear()
	{
		lastIndex = 0;
		lastNode = NULL;

		Node *p = this->head;
		while (p != NULL)
		{
			this->head = p->next;
			delete p;
			p = this->head;
		}
		this->head = NULL;
		this->tail = NULL;
		count = 0;
	}

public:
        void OrderInsert(T item)
          {
          ///////////// Insert code here
           ///////////// to insert data into the List sequentially
          }
	void PushBegin(T item) { Prepend(new Node(item)); }
	
	void PushEnd(T item) { Append(new Node(item)); }

	void Add(T item){ Append(new Node(item)); }
	
	int Count(){ return count; }	

	T elementAt(int v1) 
	{
		if ((lastIndex == v1-1)&& (lastNode!=NULL))
		{
             lastIndex = v1;
			 lastNode = lastNode->next;
			 return lastNode->info;
		}
	
	    Node* h = this->head; 
		int index =0;
		while (h!= NULL)
		{
			T item2 = h->info;
			
			lastNode = h;
			lastIndex = index;
			if (index == v1) return item2;
			h = h->next;
			index++;
		}
		lastIndex = -1;
		lastNode = NULL;
		return NULL;
	} 
	
	T operator[](int v1) const
	{ 
		 return elementAt(v1); 
	}

	T& operator[](int v1) 
	{ 
		 return elementAt(v1); 
	}
	
	T PopBegin()
	{
		T info = this->head->info;
		Remove(head);
		return info;
	}

	T PopEnd()
	{
		T info = this->tail->info;
		Remove(tail);
		return info;
	}

	void Assign( TList<T>* l)
	{
		// Hago una copia de esta lista
		this->Clear();
		Node* n = l->head;
		while (n!= l->tail)
		{
			this->Append(n->info);
			n = n->next;
		}
	}
	int indexOf(T item)
	{
		Node* h = this->head; 
		int index =0;
		while (h!= NULL)
		{
			T item2 = h->info;
			if (item2 == item) return index;
			h = h->next;
			index++;
		}
		return -1;
	};

	void InsertAt(T info, int index)
	{
		Node *p = new Node(info);

		if (index == 0 || q == NULL)
			Prepend(p);
		else
		{
			if (index >= this->count)
				Append(p);
			else
			{
				Node *q = this->head;
				while (index > 0)
				{
					index--;
					q = q->next;
				}
				InsertBefore(q, p);
			}
		}
	}

	void Extract(T info)
	{
		  int index =indexOf(info);
		  if (index>=0)
		    RemoveAt(index);
	}
	void RemoveAt(int index)
	{
		if (index == 0)
			Remove(head);
		else
		{
			if (index >= this->count)
				Remove(tail);
			else
			{
				Node *q = this->head;
				while (index > 0)
				{
					index--;
					q = q->next;
				}
				Remove(q);
			}
		}
	}
};
*/


