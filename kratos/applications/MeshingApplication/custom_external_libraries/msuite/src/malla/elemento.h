#ifndef _ELEMENTO_H
#define _ELEMENTO_H

#include <iostream> // overload de <<

#include "flags.h"

//=================================================================================
// elemento
//=================================================================================

// Los elementos clasicos son estructuras definidas, los poligonos son simplemente
// un array circular de nodos con caras en correspondencia, pero los poliedros son
// totalmente distintos.
// Para un elemento clasico basta definir los nodos en un orden preestablecido y 
// todos los datos se deducen del orden: para todos los 2d cada cara va del nodo 
// al siguiente y todas las orientaciones son hacia fuera del elemento o girando a 
// la izquierda en 2d. (algun dia hay que ponerlas hacia dentro en 3d)
// Los poligonos necesitan ademas el dato de la cantidad de nodos.
// Los poliedros necesitan varios datos ordenados correctamente: cantidad de nodos, 
// de caras, nodos de cada cara, y puede ser caras de cada nodo (porque es caro hacer
// la lista ordenada)
// Se puede hacer la clase elemento puramente virtual o conteniendo un par de cosas 
// comunes, pero asi habria que tener arrays de punteros a elementos. Es preferible 
// proveer la simplicidad de un array de elementos directamente. Eso trae 
// complicaciones dentro de la clase, pero para el que la usa es mas facil.
// Para no poner todas las necesidades del poliedro a todos los elementos, se guarda
// solo una direccion de una estructura que contiene los datos del poliedro.

class punto;

// datos de elementos que son poligonos o poliedros
class poly_dat{
public:
  int nv,nc,nt;
  int *c; // caras triangulares (indices de n)
  int **cn; // caras de cada nodo en orden (1er nro: cantidad de caras del nodo)
  int *t; // tetraedros (indices de n) (para integracion o para dividir)

  poly_dat():nv(0),nc(0),nt(0),c(0),cn(0),t(0){};
  poly_dat(const poly_dat &p);
  ~poly_dat();
  
  void ini();

  void swap(); // invierte el orden
  void mk_cn(); // hace y ordena la lista de caras de nodo
  void rm_cn(); // elimina la lista de caras de nodo
};


// Para definir el tipo de elemento hay un enum estatico e_tipo
// y pueden pedirse:
//    Un string
//    Un entero n (el de e_tipo) 0 para indefinidos y despues 1,2,3....
//    Un flag 2^n, 1 para indefinidos y despues 2,4,8...
// El string puede servir de entrada si tiene al menos los primeros 
//    caracteres que definen el tipo, en ingles o castellano


// para lectura y escritura
static const char* s_tipo[]=
  {"undefined",
   "point",
   "segment",
   "triangle","quadriateral","polygon",
   "tetrahedron","brick","wedge","polyhedron"};

// de uso general para saber que tipo de elemento es
enum e_tipo {
   e_indefinido=0,
   e_punto,
   e_segmento,
   e_triangulo,e_cuadrilatero,e_poligono,
   e_tetraedro,e_cubo,e_wedge,e_poliedro,
   e_ntipos // como va desde 0 hasta aqui, este cuenta la cantidad de tipos
};

enum { 
  ef_indefinido   =1<<e_indefinido,
  ef_punto        =1<<e_punto,
  ef_segmento     =1<<e_segmento,
  ef_triangulo    =1<<e_triangulo,
  ef_cuadrilatero =1<<e_cuadrilatero,
  ef_poligono     =1<<e_poligono,
  ef_tetraedro    =1<<e_tetraedro,
  ef_cubo         =1<<e_cubo,
  ef_wedge        =1<<e_wedge,
  ef_poliedro     =1<<e_poliedro
};

class elemento{
  friend class malla; // para que acceda a los protected
  friend class voronoi; // para que acceda a los protected

  // variables y funciones de la clase
public:

  // Las primeras letras del string de entrada deben ser:  
  // "te"                            tetraedro
  // "polie" "polyh"                 poliedro
  // "cub"   "b"     "l"     "h"     cubo brick ladrillo o hexaedro
  // "tr"                            triangulo
  // "cua"   "q"                     cuadrilatero
  // "polig" "polyg"                 poligono
  // "s"                             segmento
  // "w"     "pr"                    wedge o prisma
  // "poi"   "pu"                    punto
  //  otra cosa                      indefinido
  // En mayusculas o minusculas.
  // Lee esas letras y nada mas, ej: "teta" lo toma como tetraedro
  static e_tipo etipo(const char *);

  // Salida:
  static const char* stipo(const e_tipo t) {return s_tipo[t];}
  static int         itipo(const e_tipo t) {return t;}
  // solo para usar en >> y <<
  std::istream &lee(std::istream &a);
  std::ostream &graba(std::ostream &a) const;

  // variables y funciones de cada elemento
  int* n;  // nodos (indices de nodos de la malla)
  e_tipo Tipo;
  poly_dat *pd; // datos y funciones especiales para poligonos y poliedros
  flagtype f; // flag multiroposito (material...)

  elemento(): n(0),Tipo(e_indefinido),pd(0),f(0){};
  elemento(const elemento& e);
  elemento(e_tipo tipo);
  elemento(const char *stipo);
  elemento(e_tipo tipo, const int* nin); // no poligono ni poliedro
  ~elemento();

  elemento &ini(); // inicializa y pasa a indefinido

  elemento &copy_nodes(const int* nin); // nv fijo de antemano

  e_tipo tipo(e_tipo); // fija el tipo con e_tipo
  e_tipo tipo(const char *s) {return tipo(etipo(s));} // con string
  e_tipo tipo() const {return Tipo;} // consulta el tipo
  int itipo() const {return itipo(Tipo);}
  const char* stipo() const {return stipo(Tipo);}
  
  elemento &operator=(const elemento& e);
  elemento &set(e_tipo t, const int* nn); //ni poligono ni poliedro

  bool operator==(const elemento &e) const; // mismos nodos
  bool operator!=(const elemento &e) const {return (!operator==(e));}
  
  int nv() const; // consulta cant de nodos
  int &operator [](int i) {return n[i];} // nodo (no circular!!)
  int operator [](int i) const{return n[i];}
  int &nant(int i, int d=1) const {return n[ciclo(i-d,nv())];} // anterior circular
  int &npos(int i, int d=1) const {return n[ciclo(i+d,nv())];} // posterior circular
  int index(int i) const; // indice (-1 si no esta)
  bool have(int i) const {return index(i)>=0;} // tiene el nodo?
  bool have(const elemento &e, bool strict=false) const; // todos los nodos (no caras)
  int replace(int iviejo,int inuevo); // reemplaza iviejo por inuevo

  int nc() const; // consulta cant de caras
  elemento cara(int) const; // cara (indices de la malla)
  void cara(int i, int &nvc, int *ix) const; // pone nvc indices de n en ix
  int icara(const elemento &c) const; // (-1 si no esta) (no importa el orden)
  
  int opuesto(int i, int ix=0); // opuesto ix-esimo al i-esimo nodo 

  inline int dim() const
    {return ((Tipo>=e_tetraedro)? 3 : 
            ((Tipo>=e_triangulo)? 2 : 
            ((Tipo>=e_segmento )? 1 : 0)));}  

  int nt() const; // consulta cantidad de tetraedros
  int ntetra() const {return nt();} // consulta cantidad de tetraedros
  const int *tetra(int) const; // tetraedro (indices de la malla)

  // overload de >> y <<
  // el default es con flag (-1), 0 no lo graba y otro numero sirve de mascara binaria
  static void flag_mask(int i=-1); // lee/graba el flag?
  // base es 0 o 1, el default es 0
  static int io_base(int i=-1);
  // Para los elementos no-variables se lee/escribe la lista de nodos
  // Para poligonos va NV y la lista de nodos
  // Para poliedros va NV, los nodos, NC y los indices de c (ccw desde afuera)
  friend std::istream& operator >>(std::istream &a, elemento &e)
    {return e.lee(a);}
  friend std::ostream& operator <<(std::ostream &a,const elemento &e)
    {return e.graba(a);}
  void lee(const int* linea); // lee numeros de nodos,caras,tetras, ni tipo ni flag, 

  // solo poligonos/poliedros
  elemento &p2t(); // poliedro -> tetraedro (si tiene 4 nodos)
  elemento &c2p(); // cubo -> poliedro
};

#endif
