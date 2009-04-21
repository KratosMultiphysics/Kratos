#include <cstdlib> // alloc random qsort
#include <cstring> // memcpy, memset, memmove
#include "pline.h"

// primer constructor
pline::pline(int initsize)
{
  size=initsize;
  len=0;
  vertex=(size) ? (int *)malloc(size*SZI) : 0;
}

// segundo constructor (necesario para los return con plines)
pline::pline(const pline &a)
{
  size=a.len; len=a.len;
  if (!size) {vertex=0; return;}
  vertex=(int *)malloc(size*SZI);
  memcpy(vertex,a.vertex,len*SZI);
}
ordlist::ordlist(const pline &a)
{
  size=a.len; len=0;
  if (!size) {vertex=0; return;}
  vertex=(int *)malloc(size*SZI);
  for (int i=0;i<a.len;i++) *this+=a[i];
}
listnorep::listnorep(const pline &a)
{
  size=a.len; len=0;
  if (!size) {vertex=0; return;}
  vertex=(int *)malloc(size*SZI);
  for (int i=0;i<a.len;i++) *this+=a[i];
}

// copia de un array
pline::pline(int nv, const int* v)
{
  size=nv; len=nv;
  if (!size) {vertex=0; return;}
  vertex=(int *)malloc(size*SZI);
  memcpy(vertex,v,len*SZI);
}
ordlist::ordlist(int nv, const int* v)
{
  size=nv; len=0;
  if (!size) {vertex=0; return;}
  vertex=(int *)malloc(size*SZI);
  for (int i=0;i<nv;i++) *this+=v[i];
}
listnorep::listnorep(int nv, const int* v)
{
  size=nv; len=0;
  if (!size) {vertex=0; return;}
  vertex=(int *)malloc(size*SZI);
  for (int i=0;i<nv;i++) *this+=v[i];
}

pline::pline(int initsz, int n) // inicializado
{
  len=size=initsz;
  vertex=(int *)malloc(size*SZI);
  // pedazo de set()
  if (!len) return;
  if (!(n&&~n))
    memset(vertex,n,size*SZI); // funciona solo con chars!! => 0 o -1
  else
    for (int i=0;i<size;i++) vertex[i]=n;
}

// destructor
pline::~pline() {free(vertex);}

// roba lista
pline& pline::roba(pline& p)
{
  free(vertex);
  vertex=p.vertex; size=p.size; len=p.len;
  p.vertex=0; p.size=p.len=0; return *this;
}
pline& pline::roba(int l, int* v)
  {free(vertex); vertex=v; size=l; len=l; v=0; return *this;}

// swap de dos plines
void pline::swap(pline &p)
{
  int t, *T;
  T=vertex; vertex=p.vertex; p.vertex=T;
  t=len; len=p.len; p.len=t;
  t=size; size=p.size; p.size=t;
}


// asignacion deep copy
#define _EQP_ \
{\
  len=a.len;\
  if (len>size) {\
    free(vertex);\
    size=len;\
    vertex=(int *)malloc(len*SZI);\
  }\
  if (len) memcpy(vertex,a.vertex,len*SZI);\
  return *this;\
}

// no puede hacer a=a por el +=
#define _EQO_ \
{\
  if (&a==this) return *this;\
  if (a.len>size) {\
    free(vertex);\
    size=a.len;\
    vertex=(int *)malloc(a.len*SZI);\
  }\
  len=0; for (int i=0;i<a.len;i++) *this+=a[i];\
  return *this;\
}

pline& pline::operator=(const pline &a) _EQP_
ordlist& ordlist::operator=(const ordlist &a) _EQP_
ordlist& ordlist::operator=(const pline &a) _EQO_
listnorep& listnorep::operator=(const pline &a) _EQO_
#undef _EQO_
#undef _EQP_

pline& pline::copia(int l, const int *v){ // copia una lista
  len=l;
  if (len>size) {
    free(vertex);
    size=len;
    vertex=(int *)malloc(len*SZI);
  }
  if (len) memcpy(vertex,v,len*SZI);
  return *this;
}

ordlist& ordlist::copia(int l, const int *v){ // copia una lista
  if (l>size) {
    free(vertex);
    size=l;
    vertex=(int *)malloc(l*SZI);
  }
  len=0; for (int i=0;i<l;i++) *this+=v[i];
  return *this;
}

// asigna espacio suficiente
pline& pline::resize(int newsize)
{
  if (size>=newsize) return *this;
  return fit(newsize);
}

// asigna espacio exacto (len)
pline& pline::fit(int newsize)
{
  if (newsize==-1) newsize=len;
  else if (len>newsize) len=newsize;
  if (newsize==size) return *this;
  size=newsize;
  if (size) vertex=(int *)realloc(vertex,size*SZI);
  else {free(vertex); vertex=0;}
  return *this;
}

// asigna una constante
pline& pline::set(int n, int setlen)
{
  if (!size&&!setlen) return *this;
  if (!setlen) setlen=size;
  if (size<setlen) resize(setlen);
  if (!(n&&~n)) // memset funciona solo con chars!! => 0 o -1
    memset(vertex,n,setlen*SZI);
  else
    for (int i=0;i<setlen;i++) vertex[i]=n;
  len=setlen;
  return *this;
}

// desasigna espacio
pline& pline::ini()
{free(vertex); vertex=0; len=size=0; return *this;}

// indice de elm si no esta devuelve len
int pline::index(int elm, int start) const
{
  if (start>=len||start<0) return len;
  for (int i=start; i<len; i++) if (vertex[i]==elm) return i;
  return len;
}
// si no esta devuelve -(donde deberia ir)-1
int ordlist::index(int elm) const
{
  if (len==0) return -1;
  if (elm>=vertex[len-1]) return (elm==vertex[len-1] ? len-1 : -len-1);
  if (elm<=vertex[0]) return (elm==vertex[0] ? 0 : -1);
  // 0 < posicion de elm < len-1
  int min=0,max=len-1,i;
  while (max-min>1)
  {
    i=(max+min)/2;
    (elm<=vertex[i] ? max : min) = i;
  }
  return (vertex[max]==elm ? max : -max-1); // indice o -(donde_va)-1 (0=>-1)
}

// esta en la ordlist?
bool ordlist::have(int elm) const
{int i=index(elm); return (i>=0&&i<len);}

// cambia el viejo por el nuevo
bool pline::replace1(int viejo,int nuevo){
  int i=index(viejo);
  if (i==len) return false;
  vertex[i]=nuevo;
  return true;
}
int pline::replaceall(int viejo,int nuevo){
  int count=0,i=index(viejo);
  while (i<len){
    count++;
    vertex[i]=nuevo;
    i=index(viejo,i+1);
  }
  return count;
}

int ordlist::replace(int viejo,int nuevo){
  int i=index(viejo);
  if (i<0||i==len) return len;
  pline::remove(1,i);
  return operator+=(nuevo);
}

// inserta un vertice
pline& pline::insert(int element, int place)
{
  if (len==size) resize(len+(len>>1)+1);
  for(int i=len;i>place;i--) vertex[i]=vertex[i-1];
  vertex[place]=element;
  len++;
  return *this;
}
cpline& cpline::insert(int element, int place)
{
  if (!len||place==len) // para que len no vaya a 0 ni divida por 0
    pline::insert(element,len);
  else pline::insert(element,ciclo(place,len));
  return *this;
}

int ordlist::operator+=(int element)
{
  if (len==0) {pline::insert(element,0); return 0;}
  if (element>=vertex[len-1]) {
    if (element!=vertex[len-1]) pline::operator+=(element);
    return len-1;
  }
  if (element<=vertex[0]) {
    if (element!=vertex[0]) pline::insert(element,0);
    return 0;
  }
  // 0 < posicion de element < len-1
  int min=0,max=len-1,i;
  while (max-min>1)
  {
    i=(max+min)/2;
    (element<=vertex[i] ? max : min) = i;
  }
  if (element!=vertex[max]) pline::insert(element,max);
  return max;
}
bool listnorep::operator +=(int element){ //(si puso=>t)
  if (len==0) {pline::insert(element,0); return true;}
  if (element>=vertex[len-1]) {
    if (element!=vertex[len-1]) {
      pline::operator+=(element); return true;
    }
    len--; return false;
  }
  if (element<=vertex[0]) {
    if (element!=vertex[0]) {
      pline::insert(element,0); return true;
    }
    pline::remove(1,0); return false;
  }
  // 0 < posicion de element < len-1
  int min=0,max=len-1,i;
  while (max-min>1)
  {
    i=(max+min)/2;
    (element<=vertex[i] ? max : min) = i;
  }
  if (element!=vertex[max]) {pline::insert(element,max); return true;}
  pline::remove(1,max); return false;
}

// inserta una pline en otra
pline& pline::insert(const pline &intern, int place)
{
  int il=intern.len,nl=len+il;
  if (!il) return *this;
  if (nl>size) resize(nl);
  if (place<len) memmove(&vertex[place+il],&vertex[place],(len-place)*SZI);
  memcpy(&vertex[place],intern.vertex,il*SZI);
//  for(i=nl-1;i>=place+il;i--) vertex[i]=vertex[i-il];
//  for(i=0;i<il;i++) vertex[place+i]=intern[i];
  len=nl;
  return *this;
}
cpline& cpline::insert(const pline &intern, int place)
{
  if (!len||place==len) // para que len no vaya a 0 ni divida por 0
    pline::insert(intern,len);
  else pline::insert(intern,ciclo(place,len));
  return *this;
}
ordlist& ordlist::operator+=(const pline &intern)
  {for (int i=0;i<intern.len;i++) *this+=intern[i]; return *this;}
listnorep& listnorep::operator+=(const pline &intern)
  {for (int i=0;i<intern.len;i++) *this+=intern[i]; return *this;}

// elimina
// elimina un cacho de pline
pline& pline::remove(int qty, int place)
{
  if (qty<=0||place>len-qty||place<0) return *this;
  if (place!=len-qty)
    memmove(&vertex[place],&vertex[place+qty],(len-qty-place)*SZI);
  len-=qty;
  return *this;
}
cpline& cpline::remove(int qty, int place)
{
  if (len==0||qty<=0) return *this;
  place=ciclo(place,len); // positivo < len
  if ((place+qty)>len) // sobrepasa el final
  {
    qty-=(len-place); //los que hay que eliminar del ppio
    len=place; // elimina los ultimos
    place=0; // a partir del ppio
  }
  pline::remove(qty,place);
  return *this;
}

pline& pline::removeall(int element){
  int i=index(element);
  while (i<len){
    remove(1,i);
    i=index(element,i);
  }
  return *this;
}

pline& pline::remove1(const pline& intern)
  {for (int i=0;i<intern.len;i++) remove1(intern[i]); return *this;}
pline& pline::removeall(const pline& intern)
  {for (int i=0;i<intern.len;i++) removeall(intern[i]); return *this;}


// comparacion
bool pline::operator==(const pline& pl) const
{
  if (len != pl.len) return 0;
  int i;
  for (i=0;i<len;i++) if (vertex[i] != pl.vertex[i]) break;
  return (i==len);
}
bool cpline::operator==(const pline& pl) const
{
  if (len != pl.len) return 0;
  int i;
  for (i=0;i<len;i++) if (pl[0]==vertex[i]) break;
  if (i==len) return 0;
  int j=i;
  for (i=0;i<len;i++) if (vertex[ciclo(j+i,len)] != pl.vertex[i]) break;
  return (i==len);
}

bool cpliner::operator==(const pline& pl) const
{
  if (len != pl.len) return 0;
  int i;
  for (i=0;i<len;i++) if (pl[0]==vertex[i]) break;
  if (i==len) return 0; // el 0 no esta
  int j=i; //el 0 es = al j
  for (i=1;i<len;i++)
    if (vertex[j+i] != pl[i]) break;
  if (i==len) return true;
  for (i=1;i<len;i++)
    if (vertex[j-i] != pl[i]) break;
  return (i==len);
}

// interseccion
cpliner pline::inters() const // con sigo misma
{
  cpliner inter;
  for (int i=0;i<len-1;i++) if (index(vertex[i],i+1)<len) inter+=vertex[i];
  return inter;
}
cpliner pline::inters(const pline& pl)  const // con otra pline o con sigo misma
{
  if (&pl==this) return inters();
  cpliner inter;
  for (int i=0;i<len;i++) if (pl.have(vertex[i])) inter+=vertex[i];
  return inter;
}

// reordena la pline con el elemento n-esimo primero (indice n --> 0)
pline& pline::first(int new0)
{
  if (!(len&&new0)) return *this;
  int* aux=(int *)malloc(new0*SZI);
  memcpy(aux,vertex,new0*SZI);
  memcpy(vertex,vertex+new0,(len-new0)*SZI);
  memcpy(vertex+(len-new0),aux,new0*SZI);
  free(aux);
  return *this;
}

// invierte el orden de los elementos
pline& pline::reverse()
{
  for (int t,i=0;i<(len>>1);i++)
   {t=vertex[i];vertex[i]=vertex[len-1-i];vertex[len-1-i]=t;}
  return *this;
}

// distancia entre elm1 y elm2
int pline::distance(int elm1,int elm2) const
{
  int i1=index(elm1),i2=index(elm2);
  if (i1==len||i2==len) return -1;
  return (i2>i1)? (i2-i1) : (i1-i2);
}
int cpline::distance(int elm1,int elm2) const
{
  int i1=index(elm1),i2=index(elm2);
  if (i1==len||i2==len) return len;
  return (i2>i1)? (i2-i1) : (len+i2-i1);
}


// numeros naturales o a qty
// solo agrega pero no cambia los ya hechos (=> inicializar)
pline& pline::natural(int qty)
{
  if (size<qty) {
    int oldsz=size; fit(qty);
    for (int i=oldsz;i<size;i++) vertex[i]=i;
  }
  len=qty;
  return *this;
}

// randomiza el orden
// Primero hago una lista ordenada, luego intercambio cada uno con uno random
// RAND_MAX puede ser cualquier cosa desde short en adelante
// El problema es que no va a dar "todos" los nros entre 0 y len-1
//   eso puede provocar problemas de "resonancia" si hubiese un orden muy definido de subgrupos
// En sintesis: esto no es exacto pero sirve muy bien para todo lo que aqui se hace
pline& pline::random_index(int qty){
  int i,ix,t;
  if (size<qty) fit(qty);
  len=qty;
  for (i=0;i<len;i++) vertex[i]=i;
  const long double divisor=((long double)(len-1))/RAND_MAX; //[0,len)
  for (i=0;i<len;i++){
    ix=int(divisor*rand());
    t=vertex[i]; vertex[i]=vertex[ix]; vertex[ix]=t;
  }
  return *this;
}
             /*
// salida a archivo o window
pline& pline::print(ofstream* out)
{
  out << "# " << len << ":";
  for (int i=0;i<len;i++) out << vertex[i] << ",";
  out << "\t";
}
            */
