//                   VECTOR DE ACCESO ALEATORIO
//  Se pueden agregar o quitar elementos solo al final de la lista
//
// PARA LA CLASE DE ENTRADA DEBE HABER
//     CONSTRUCTOR DEFAULT, DEEP COPY Y COMPARACION (==)

// la lista se hace con malloc/realloc/free por la facilidad de usar realloc
// y los elementos se hacen con new/delete por el constructor/destructor

// array1h mantiene una lista de elementos borrados intermedios y agrega en los huecos

// modificado para aceptar two-stage name lookup (gcc 4)

#ifndef _ARRAY1_H
#define _ARRAY1_H

#include <cstdlib> // qsort, malloc/free
#include "pline.h" // plines para array1h
#include <cstring>

//template <class c> class array1h;

template <class c> class array1{
 public:
  int  len;  // cantidad de elementos
  c**  list; // array de elementos
  int  size; // tamanio maximo antes de realocar elementos

  array1(int initsz=0);      // constructor default
  array1(const array1<c>&);  // constructor por copia
//  array1(const array1h<c>&); // constructor por copia
  array1(const c*, int);     // constructor por copia
  array1(int initsz, const c& value);// inicializado y con len
  ~array1();  // destructor
  void ini(); // inicializa
  array1<c>& set(const c& value, int qty=0); // todos = a algo (0=size)
  array1<c>& operator=(const array1<c>&); // deep copy (copia)
//  array1<c>& operator=(const array1h<c>&); // deep copy (copia)
  array1<c>& roba(array1<c>&); // roba la lista
  inline c& operator[](int i) const {return *list[i];} // elemento
  c& last() const {return *list[len-1];}
  void clean() {len=0;} // borra el contenido
  int index(const c &element) const; // indice del elemento o len si no existe
  bool have(const c &element) const;  // si el elemento existe
  inline int next_index() const {return len;} // indice que corresponderia al proximo elemento a agregar
  int operator+=(const c &element);  // agrega un elemento y devuelve el indice
  array1<c>& operator+=(const array1<c>&); // agrega un array y devuelve este
  int add(const c &element, bool norepetir=0); // idem pero no repite si norepetir=1
  bool resize(int newsz);  // asigna espacio suficiente
  bool fit(int newsize=-1); // asigna espacio exacto (-1=len)
  operator bool() const {return (len)? true : false;} // para averiguar si tiene algo
  array1<c>& reverse(); // invierte el orden
  void swap(int i, int j){c* t=list[i]; list[i]=list[j]; list[j]=t;}
//  void sort(int compare(const void*,const void*)); // quick sort (hay que definir compare)
};

// la cantidad de elementos es len-borrado.len
// el ultimo elemento no debe ser borrado, se baja el len y listo
template <class c> class array1h : public array1<c>
{
 public:
  bool *es_borrado;
  ordlist borrado; // puede ser pline, pero asi es mas seguro y mas facil

  array1h(int initsz=0);        // constructor default
  array1h(const array1<c> &o);  // constructor por copia
  array1h(const array1h<c> &o); // constructor por copia
  array1h(const c* o, int sz);  // constructor por copia
  array1h(int initsz, const c& value); // inicializado y con len
  ~array1h();  // destructor
  void ini(); // inicializa
  array1h<c>& set(const c& value, int qty=0);   // todos = a algo (0=size)
  array1h<c>& operator=(const array1<c>&); // deep copy (copia)
  array1h<c>& operator=(const array1h<c>&); // deep copy (copia)
  array1h<c>& roba(array1<c>&);  // roba la lista (shallow copy)
  array1h<c>& roba(array1h<c>&); // roba la lista
  inline c& operator[](int i) const {return *(array1<c>::list)[i];} // elemento
  c& last() const {return *(array1<c>::list)[array1<c>::len-1];}
  void clean() {array1<c>::len=0; borrado.len=0;} // borra el contenido
  int index(const c &element) const; // indice del elemento o len si no existe o es borrado
  bool have(const c &element) const; // si el elemento existe (y no es borrado)
  inline int next_index() const; // indice que corresponderia al proximo elemento a agregar
  int operator+=(const c &element);  // agrega un elemento y devuelve el indice
  array1h<c>& operator+=(const array1<c>&);  // agrega un array1 y devuelve este
  array1h<c>& operator+=(const array1h<c>&);  // agrega un array1h y devuelve este
  int add(const c &element, bool flag=0); // idem pero no repite si flag=1
  bool resize(int newsz);  // asigna espacio suficiente
  bool fit(int newsize=-1); // asigna espacio exacto (-1=len+borrado.len)
  void swap(int i, int j);
  operator bool() const {return (array1<c>::len-borrado.len)? true : false;} // tiene algo?

  // propias
  array1h<c>& squeeze(); // elimina huecos
  array1h<c>& squeeze(pline &map); // y hace un mapa de indices antiguos
  array1h<c>& pack(); // squeeze+fit
  array1h<c>& pack(pline &map); // y hace un mapa de indices antiguos
  array1h<c>& operator-=(const c&e);  // elimina un elemento
  array1h<c>& remove(int ix);  // elimina un elemento
};


//===========================================================
//                  IMPLEMENTACION
//===========================================================
#define _INI_ES_BORRADO \
  if (!array1<c>::size) {es_borrado=0; return;}\
  es_borrado=(bool *)malloc(SZB*array1<c>::size);\
  memset(es_borrado,false,SZB*array1<c>::size);

// primer constructor
template <class c> array1<c>::array1(int initsz){
  len=0;
  size=initsz;
  if (size) {
    list=(c**)malloc(sizeof(c*)*size);
    for (int i=0;i<size;i++) list[i]=new c;
  }
  else list=0;
}
template <class c> array1h<c>::array1h(int initsz)
: array1<c>(initsz){
  _INI_ES_BORRADO
}

// constructor por copia
template <class c> array1<c>::array1(const array1<c>& r)
  {list=0; size=0; len=0; *this=r;}
//template <class c> array1<c>::array1(const array1h<c>& r)
//  {list=0; size=0; len=0; *this=r;}
template <class c> array1h<c>::array1h(const array1<c>& r)
  {array1<c>::list=0; array1<c>::size=0; array1<c>::len=0; es_borrado=0; *this=r;}
template <class c> array1h<c>::array1h(const array1h<c>& r)
  {array1<c>::list=0; array1<c>::size=0; array1<c>::len=0; es_borrado=0; *this=r;}

// constructor por copia de un c*
template <class c> array1<c>::array1(const c *r, int n){
  len=n;
  size=len;
  if (size){
    list=(c**)malloc(sizeof(c*)*size);
    for (int i=0;i<size;i++)
      {list[i]=new c; *list[i]=r[i];} // deep copy
  }
  else list=0;
}
template <class c> array1h<c>::array1h(const c *r, int n)
: array1<c>(r,n){
  _INI_ES_BORRADO
}

// constructor inicializado
template <class c> array1<c>::array1(int initsz, const c& value){
  len=initsz;
  size=initsz;
  if (size){
    list=(c**)malloc(sizeof(c*)*size);
    for (int i=0;i<size;i++)
      {list[i]=new c; *list[i]=value;} // deep copy
  }
  else list=0;
}
template <class c> array1h<c>::array1h(int initsz, const c& value)
: array1<c>(initsz,value){
  _INI_ES_BORRADO
}

#undef _INI_ES_BORRADO

// destructor (no hace falta poner en 0)
template <class c> array1<c>::~array1(){
  for (int i=0;i<size;i++) {delete list[i]; list[i]=0;}
  free(list);
}
template <class c> array1h<c>::~array1h(){
  free(es_borrado); // ojo: el borrador de array1 se llama solo
}

// inicializacion
template <class c> void array1<c>::ini(){
  for (int i=0;i<size;i++) {delete list[i]; list[i]=0;}
  free(list);list=0;
  len=0;
  size=0;
}
template <class c> void array1h<c>::ini(){
  array1<c>::ini();
  free(es_borrado); es_borrado=0;
  borrado.ini();
}

// todos igual a algo
template <class c> array1<c>& array1<c>::set(const c& value, int qty)
{
  int i;
  if (qty>size){
    list=(c**)realloc(list,sizeof(c*)*qty);
    for (i=0;i<qty;i++) list[size+i]=new c;
    size=qty;
  }
  if (!qty) qty=size;
  len=qty;
  for (i=0;i<qty;i++) *list[i]=value;
  return *this;
}
template <class c> array1h<c>& array1h<c>::set(const c& value, int qty){
  borrado.clean();
  if (qty>array1<c>::size) es_borrado=(bool *)realloc(es_borrado,SZB*qty);
  memset(es_borrado,false,SZB*array1<c>::size);
  array1<c>::set(value,qty);
  return *this;
}

// deep copy (solo len, pero incluyendo borrados)
template <class c> array1<c>& array1<c>::operator=(const array1<c>& r){
  int i;
  if (&r==this) return *this; // no puede hacer r=r por el delete
  for (i=r.len;i<size;i++) {delete list[i]; list[i]=0;}
  list=(c**)realloc(list,sizeof(c*)*r.len);
  for (i=size;i<r.len;i++) list[i]=new c;
  len=r.len;
  size=len;
  if (size) {
    for (int i=0;i<len;i++) *list[i]=*(r.list[i]);
  }
  else list=0;
  return *this;
}
//template <class c> array1<c>& array1<c>::operator=(const array1h<c>& r){
//  r.squeeze(); operator=((array1<c>)r);
//}
template <class c> array1h<c>& array1h<c>::operator=(const array1<c>& r){
  array1<c>::operator=(r);
  borrado.clean();
  es_borrado=(bool *)realloc(es_borrado,SZB*array1<c>::size);
  memset(es_borrado,false,SZB*array1<c>::size);
  return *this;
}
template <class c> array1h<c>& array1h<c>::operator=(const array1h<c>& r){
  array1<c>::operator=((array1<c>)r);
  borrado=r.borrado;
  es_borrado=(bool *)realloc(es_borrado,SZB*array1<c>::size);
  memcpy(es_borrado,r.es_borrado,SZB*array1<c>::size);
  return *this;
}

// roba
template <class c> array1<c>& array1<c>::roba(array1<c>& r){
  for (int i=0;i<size;i++) {delete list[i]; list[i]=0;}
  free (list);
  list=r.list; r.list=0;
  len=r.len; r.len=0;
  size=r.size; r.size=0;
  return *this;
}
template <class c> array1h<c>& array1h<c>::roba(array1<c>& r){
  array1<c>::roba(r);
  es_borrado=(bool *)realloc(es_borrado,SZB*array1<c>::size);
  memset(es_borrado,false,SZB*array1<c>::size);
  borrado.clean();
  return *this;
}
template <class c> array1h<c>& array1h<c>::roba(array1h<c>& r){
  array1<c>::roba(r); borrado.roba(r.borrado);
  free(es_borrado); es_borrado=r.es_borrado; r.es_borrado=0;
  return *this;
}

// asignacion de espacio suficiente
template <class c> bool array1<c>::resize(int newsize){
  if (newsize<=size) return true;
  return fit(newsize);
}
template <class c> bool array1h<c>::resize(int newsize){
  if (newsize<=array1<c>::size) return true;
  return fit(newsize);
}

// asignacion de espacio exacto
template <class c> bool array1<c>::fit(int newsize){
  int i;
  if (newsize==-1) newsize=len;
  if (newsize==size) return true;
  if (len>newsize) len=newsize;
  for (i=newsize;i<size;i++) {delete list[i]; list[i]=0;}
  list=(c**)realloc(list,sizeof(c*)*newsize);
  if (!list) return false;
  for (i=size;i<newsize;i++) list[i]=new c;
  size=newsize;
  return true;
}
template <class c> bool array1h<c>::fit(int newsize){
  if (newsize==-1) newsize=array1<c>::len+borrado.len;
  else {if (newsize<array1<c>::len+borrado.len) return false;}
  borrado.fit();
  if (!array1<c>::fit(newsize)) return false;
  es_borrado=(bool *)realloc(es_borrado,SZB*array1<c>::size);
  if (!es_borrado) return false;
  return true;
}

// indice de un elemento en la lista
template <class c> int array1<c>::index(const c& elm) const{
  int i;
  for (i=0;i<len&&*list[i]!=elm;i++);
  return i;
}
template <class c> int array1h<c>::index(const c& elm) const{
  int i;
  for (i=0;i<array1<c>::len;i++){
    if (!es_borrado[i]&&*(array1<c>::list)[i]==elm) break;
  }
  return i;
}

// si el elemento existe
template <class c> bool array1<c>::have(const c &element) const
{return (index(element)!=len);}
template <class c> bool array1h<c>::have(const c &element) const
{return (index(element)!=array1<c>::len);}

template <class c> int array1h<c>::next_index() const {
  return (borrado.len)? borrado.last() : array1<c>::len;
}

// agrega un elemento (debe haber deep copy para la clase c)
template <class c> int array1<c>::operator+=(const c& elm){
  if (len==size) fit(len+1+(len>>1));
  *list[len++]=elm;
  return len-1;
}
template <class c> int array1h<c>::operator+=(const c& elm){
  int i;
  if (borrado){
    i=borrado.last(); borrado.len--;
  }
  else{
    if (array1<c>::len==array1<c>::size) fit(array1<c>::len+1+(array1<c>::len>>1));
    i=array1<c>::len++;
  }
  *(array1<c>::list)[i]=elm; es_borrado[i]=false;
  return i;
}

// agrega un array
template <class c> array1<c>&
    array1<c>::operator+=(const array1<c>& o){
  int newlen=len+o.len;
  if (newlen>size) fit(newlen+1+(newlen>>1));
  for (int i=0;i<o.len;i++) *list[len+i]=o[i];
  len=newlen;
  return *this;
}
template <class c> array1h<c>&
    array1h<c>::operator+=(const array1<c>& o){
  int i,io=0,vivos=array1<c>::len-borrado.len+o.len;
  if (vivos>array1<c>::size) fit(vivos+1+(vivos>>1));
  while (borrado.len&&io<o.len) {
    i=borrado.last();
    *(array1<c>::list)[i]=o[io++];
    borrado.len--; es_borrado[i]=false;
  }
  while(io<o.len) *(array1<c>::list)[array1<c>::len++]=o[io++];
  return *this;
}
template <class c> array1h<c>&
    array1h<c>::operator+=(const array1h<c>& o){
  int i,io=0,vivos=array1<c>::len-borrado.len+o.len-o.borrado.len;
  if (vivos>array1<c>::size) fit(vivos+1+(vivos>>1));
  while (borrado.len&&io<o.len) {
    i=borrado.last();
    while (o.es_borrado[io]) io++;
    *(array1<c>::list)[i]=o[io++];
    borrado.len--; es_borrado[i]=false;
  }
  while(io<o.len) {
    while (o.es_borrado[io]) io++; 
    *(array1<c>::list)[array1<c>::len++]=o[io++];
  }
  return *this;
}

// idem pero no repite si flag=1
template <class c> int array1<c>::add(const c &elm, bool flag){
  int i;
  if (flag && (i=index(elm))<len) return i;
  return *this+=elm;
}
template <class c> int array1h<c>::add(const c &elm, bool flag){
  int i;
  if (flag && (i=index(elm))<array1<c>::len) return i;
  return *this+=elm;
}

// invierte el orden de los elementos
template <class c> array1<c>& array1<c>::reverse(){
  c *t; int i;
  for (i=0;i<(len>>1);i++)
   {t=list[i];list[i]=list[len-1-i];list[len-1-i]=t;}
  return *this;
}

// intercambia dos elementos
template <class c> void array1h<c>::swap(int i, int j){
  c* t=(array1<c>::list)[i]; 
  (array1<c>::list)[i]=(array1<c>::list)[j]; 
  (array1<c>::list)[j]=t;
  bool ebi=es_borrado[i],ebj=es_borrado[j];
  if (ebi&&!ebj) {borrado-=i; borrado+=j;}
  if (ebj&&!ebi) {borrado-=j; borrado+=i;}
  es_borrado[i]=ebj; es_borrado[j]=ebi;
  if (ebi||ebj) {
    while (array1<c>::len&&es_borrado[array1<c>::len-1]) {
      borrado.len--; array1<c>::len--;
    }
  }
}

/*
// quick sort
template <class c> void array1<c>::sort(int compare(const void*,const void*))
{std::qsort(list,len,sizeof(c),compare);}
// ejemplo de compare:
//int compare(const void *p1, const void *p2)
//{return ((*(punto *)p1)[2]<(*(punto *)p2)[2]) ? -1 : 1;}
*/

// compacta
// hay que swapear por el delete de ~array1h
template <class c> array1h<c>& array1h<c>::squeeze(){
  int ixo,ixn;
  while (borrado){
    while (array1<c>::len&&es_borrado[array1<c>::len-1]) {
      array1<c>::len--; borrado.len--;
    } 
    ixo=array1<c>::len-1; ixn=borrado.last();
    es_borrado[ixn]=false;
    c* t=(array1<c>::list)[ixo];
    (array1<c>::list)[ixo]=(array1<c>::list)[ixn]; 
    (array1<c>::list)[ixn]=t;
    array1<c>::len--; borrado.len--;
  }
  return *this;
}
template <class c> array1h<c>& array1h<c>::squeeze(pline &map){
  map.ini(); map.natural(array1<c>::len);
  int ixo,ixn;
  while (borrado){
    while (array1<c>::len&&es_borrado[array1<c>::len-1]) {
      array1<c>::len--; borrado.len--;
    } 
    ixo=array1<c>::len-1; ixn=borrado.last();
    es_borrado[ixn]=false;
    map[ixn]=ixo; map[ixo]=ixn;
    c* t=(array1<c>::list)[ixo];
    (array1<c>::list)[ixo]=(array1<c>::list)[ixn];
    (array1<c>::list)[ixn]=t;
    array1<c>::len--; borrado.len--;
  }
  return *this;
}
template <class c> array1h<c>& array1h<c>::pack(){
  squeeze(); fit();
  return *this;
}
template <class c> array1h<c>& array1h<c>::pack(pline &map){
  squeeze(map); fit();
  return *this;
}

// elimina un elemento
template <class c> array1h<c>& array1h<c>::remove(int ix){
  if (ix<0||ix>=array1<c>::len||es_borrado[ix]) return *this;
  es_borrado[ix]=true;
  if (ix==array1<c>::len-1) {
    array1<c>::len--; 
    while (array1<c>::len&&es_borrado[array1<c>::len-1]) {
      array1<c>::len--; borrado.len--;
    }
  }
  else borrado+=ix;
  return *this;
}
template <class c> array1h<c>& array1h<c>::operator-=(const c& e){
  return remove(index(e));
}

#endif
