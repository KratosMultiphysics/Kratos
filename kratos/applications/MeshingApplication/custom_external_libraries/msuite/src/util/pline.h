// pline: array unidimensional de enteros
// cpline: cerrado (len->0)
// cpliner no importa la direccion (==)
// ordlist: se mantiene ordenado y sin repetir
// listnorep: ordenado, si se quiere agregar uno que estaba se borran ambos

#ifndef _PLINE_H
#define _PLINE_H

#include "utiles.h"

class cpliner; // necesto la definicion por inters

// Array de enteros.
class pline
{
public:
  int len;       // cantidad de indices ocupados
  int* vertex;   // el array de ints
  int  size;     // tamanio

  pline(){vertex=0;len=size=0;}   // default
  pline(int initsz);              // primer constructor
  pline(const pline &);           // segundo constructor
  pline(int, const int*);         // copia un array
  pline(int initsz, int initval); // inicializado (y con len)
  ~pline();                       // destructor

  inline int& operator[] (int n)    // indice
    const {return vertex[n];}
  int index(int elm,int start=0) const;    // indice de elem desde start (o len)
  int& last() const {return vertex[len-1];} // el ultimo (ojo si len==0)
  bool have(int elm) const {return index(elm)<len;} // esta en la pline ?
  pline& operator=(const pline &);         // deep copy
  pline& copia(int l, const int *v);       // copia una lista
  pline& roba(pline& p);                   // roba lista
  pline& roba(int l, int* v);              // roba lista
  void swap(pline &p);                     // intercambia listas
  friend void swap(pline &a, pline &b)     // intercambia listas
    {a.swap(b);}
  bool replace1(int viejo,int nuevo);      // cambia la primera aparicion (false=no estaba)
  int replaceall(int viejo,int nuevo);     // cambia en todas y devuelve la cant
  pline& remove(int qty, int place);       // elimina un pedazo
  pline& remove1(int element)              // elimina la primera aparicion 
    {return remove(1,index(element));}
  pline& removeall(int element);           // elimina en todas
  pline& remove1(const pline &intern);     // elimina cada elemento una vez
  pline& removeall(const pline &intern);   // toda aparicion de cada elemento
  pline& insert(const pline &intern, int place); // inserta un pedazo
  pline& insert(int element, int place);   // inserta un elemento
  pline& operator+=(int element)           // inserta un elemento al final
    {insert(element,len); return *this;}
  pline& operator+=(const pline &intern)   // inserta un pedazo al final
    {insert(intern,len); return *this;}
  pline& first(int ix);                    // ix (indice) al ppio
  bool operator==(const pline&) const;     // mismos elementos y orden
  bool operator!=(const pline& p)          // distintos idem
    const {return !(*this==p);}
  cpliner inters() const;                  // interseccion con sigo misma
  cpliner inters(const pline&) const;      // interseccion con otra pline
  pline& reverse();                        // invierte el orden
  int distance(int elm1,int elm2) const;   // distancia de elm1 a elm2 (>0)
  pline& fit(int newsize=-1);              // asigna espacio exacto (-1=len)
  pline& resize(int);                      // asigna espacio suficiente
  pline& clean() {len=0; return *this;}     // borra el contenido
  pline& ini();                            // la hace percha
  pline& set(int i=1, int setlen=0);       // asigna valor y setea len (0=> len=size)
  pline& unset(int setlen=0) {return set(0,setlen);}
  pline& natural(int qty);                 // numeros naturales (OJO construir con size 0)
  pline& random_index(int qty);            // numeros del 0 al n-1 en orden aleatorio
  operator int*()                          // para averiguar si existe
    const {return (len)? vertex : 0;}
};

// Array circular
class cpline : public pline
{
public:
  cpline(){vertex=0;len=size=0;}    // default
  cpline(int initsz)                // primer constructor
    : pline(initsz) {}
  cpline(const pline &p)            // segundo constructor
    : pline(p) {}
  cpline(int initsz, const int* ar) // copia un array
    : pline(initsz,ar) {}
  cpline(int initsz, int initval)    // inicializado
    : pline(initsz,initval) {}
  inline int& operator[] (int n)
    const {return vertex[ciclo(n,len)];} // indice circular (malditos negativos!)
//  int index(int elm,int start=0)     // indice de elem desde start (o len)
//    const {return (len) ? pline::index(elm,ciclo(start,len)) : 0;} NOOOO! (len->0)
  cpline& first(int new0)                         // new0 (indice) al ppio
    {new0=ciclo(new0,len); pline::first(new0); return *this;}
  cpline& insert(const pline &intern, int place);// inserta un pedazo
  cpline& insert(int element, int place);        // inserta un elemento
  cpline& remove(int qty, int place);            // elimina un pedazo
  bool operator==(const pline&) const;            // elementos y orden (circular)
  bool operator!=(const pline& c)                 // elementos u orden (circular)
    const {return !(*this==c);}
  int distance(int e1,int e2) const; // distancia del elemento 1 al 2
};

// Array circular pero no calienta la direccion (reversible)
class cpliner : public cpline
{
public:
  cpliner(){vertex=0;len=size=0;}     // default
  cpliner(int initsz)                 // primer constructor
    : cpline(initsz) {}
  cpliner(const pline &p)             // segundo constructor
    : cpline(p) {}
  cpliner(int initsz, const int* ar)  // copia un array
    : cpline(initsz,ar) {}
  cpliner(int initsz, int initval)     // inicializado
    : cpline(initsz,initval) {}

  bool operator==(const pline&) const;       // elementos y orden (circular reversible)
  bool operator!=(const pline& c)            // elementos u orden (circular reversible)
    const {return !(*this==c);}
};

// Array ordenado de ints distintos.
class ordlist : public pline
{
public:
  ordlist(){vertex=0;len=size=0;}  // default
  ordlist(int initsz)               // primer constructor
    : pline(initsz){}
  ordlist(const ordlist &o)        // segundo constructor
    : pline(o) {}
  ordlist(const pline &p);         // segundo constructor bis
  ordlist(int, const int*);        // copia un array

  inline int operator[](const int n) // copia del elemento, no lvalue (sin &)
    const {return vertex[n];}
  int last() const {return vertex[len-1];} // el ultimo (ojo si len==0)

  ordlist& operator=(const ordlist &);// deep copy
  ordlist& operator=(const pline &);// deep copy
  ordlist& copia(int l, const int *v); // copia una lista
  int index(int) const;  // indice si esta, -(donde debe ir)-1 si no esta
  bool have(int elm) const;       // esta en la ordlist?
  int replace(int viejo,int nuevo);  // cambia el primero por el segundo
  int operator +=(int element); // devuelve el indice
  ordlist& operator +=(const pline &intern);  // insercion de una pline
  ordlist& operator-=(int elm) {pline::remove1(elm); return *this;}
  ordlist& remove(int elm) {pline::remove1(elm); return *this;}
  ordlist& remove(const pline &intern) {pline::remove1(intern); return *this;}
};

// Array ordenado de ints distintos, pero al agregar uno repetido se elimina.
class listnorep : public ordlist
{
public:
  listnorep(){vertex=0;len=size=0;}  // default
  listnorep(int initsz)              // primer constructor
    : ordlist(initsz){}
  listnorep(const listnorep &o)      // segundo constructor
    : ordlist(o) {}
  listnorep(const ordlist &o)        // tercer constructor
    : ordlist(o) {}
  listnorep(const pline &p);         // cuarto constructor
  listnorep(int, const int*);        // copia un array

//  inline int operator[](const int n) // copia del elemento, no lvalue (sin &)
//    const {return vertex[n];}

  listnorep& operator =(const pline &);// deep copy
  bool operator +=(int element); // insercion de un elemento (si puso=>t)
  listnorep& operator +=(const pline &intern);  // insercion de una pline
};

#endif
