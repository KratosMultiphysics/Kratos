#ifndef _MALLA_H
#define _MALLA_H

#include "utiles.h"
#include "mensajes.h"
#include "array1.h"
#include "pline.h"
#include "octree.h"
#include "flags.h"
#include "elemento.h"
#include "filename.h"

#ifdef _CALCULO_
#include "nodedata.h"
#else
#include "nodo.h"
#endif

//==============================================================
// flags
//==============================================================

// los primeros 16 son para el usuario, el resto son reservados

#define _FREE_FLAGS 16

enum { // malla
  m_vol=           1<<29,
  m_sup=           1<<28,
  m_lin=           1<<27,
  m_nodos=         1<<26,
  m_planaxy=       1<<25,
  m_cerrada=       1<<24,
  m_orientada=     1<<23,
  m_modificada=    1<<22,
  m_slivers_marked=1<<21
};

enum { // nodo
  n_h=             1<<29, //=e
  n_frontera=      1<<28, //=e
  n_virtual=       1<<27, //=e
  n_borrado=       1<<26, //=e
  n_offset=        1<<25, //=e
  n_simetria=      1<<24, //=e
  n_permanente=    1<<23,
  n_edge=          1<<22
};

enum { // elemento
  e_h=             1<<29, //=n
  e_frontera=      1<<28, //=n
  e_virtual=       1<<27, //=n
  e_borrado=       1<<26, //=n
  e_offset=        1<<25, //=n
  e_simetria=      1<<24, //=n
  e_exterior=      1<<23,
  e_sliver=        1<<22
};

enum { // para todo servicio
  flag1=           1<<20,
  flag2=           1<<19,
  flag3=           1<<18
};

// los que se leen/graban
static const int fmask=(1<<_FREE_FLAGS)-1+n_h+n_offset+n_simetria;


//==============================================================
// cascara
//==============================================================
class malla;
// Una cascara esta formada por nodos de la malla pero tiene sus
//   propios elementos (en gral. una dimension menos que la malla)
// Se usa como frontera o aristas de la malla o elementos sueltos

class cascara {
 public:
  malla *parent; // a quien pertenece (no necesariamente const, errores, etc.)
  pline n; // indices de los nodos de la malla
  array1<elemento> e;   // lista de elementos (nodos e indices de la malla)
  bool hayfe; // hay flags de elementos  

  //faltan constructores
  cascara(): parent(0) {}; //default

  bool graba_dat(const char *archivo=0);
  bool graba_con(const char *archivo=0);
  bool agrega_con(const char *stipo,const char *archivo=0);

//  double area() const;    // area
//  void bbox(punto[2]) const; // bbox de la cascara
//  bool inter(const punto&) const; // si un punto es interior
};

//==============================================================
// malla
//==============================================================
typedef class malla{

 public:

// customizacion para quantech

#ifdef _QUANTECH
#ifdef _QT
#include "quantech.h"
#endif
#endif
   
  // nombre y descripcion
  char nombre[_max_file_len],ext[_max_ext_len];
  flagtype tipo; // flag binario
  char *m_error,*m_warning; // mensajes (usa new/delete)
  static bool INFO_CL; // para rutinas largas si hay ventana del OS (default=false)

  // nodos
  array1<nodo> n; // lista de nodos
  int nodosh; // cantidad de nodos que solo son para definir h
  array1<cpline> nn; // vecinos naturales
  punto pmin,pmax; // bounding box
  octree *o;       // octree (* porque no hay default)
  double epsilon;  // igualacion de nodos (default: _prec de point.h)
  bool hayh,hayv,hayfn; 
  double hmin; // minimo h de toda la malla
  double nvmin,nvmax; // nodevalues

  // elementos
  array1<elemento> e;  // lista de elementos
  array1<cpline> vecino; // vecinos por cara (negativo=>frontera)
  array1<double> ve; // volumen
  array1<double> re; array1<punto> ce; // centro y radio
  bool hayfe; // hay flags de elementos  
  // normales hacia afuera por nodo/elemento (solo puede tener modulo uno, o cero si hay error)
  array1<punto> dir; bool ndir,edir; 

  // frontera
  // como un conjunto de cascaras conexas y cerradas
  array1<cascara> frontera; 

  malla();  // constructor default
  malla(const malla&); // constructor por copia
  malla( // a partir de arrays sueltos (conectividades en base 1)
    int ndim, int nlen, double* nl, // dim cantidad y lista de puntos (x0y0[z0]x1y1[z1].....)
    double* hl, // h de cada nodo (puede ser un array nulo)
    int edim, int elen, int* el   // dim, cantidad y lista de simplices (pueden ser 0 y nulo)
  );
  malla( // a partir de arrays sueltos (zmax-zmin define si es 3D o 2D, triangulos o segmentos)
    int ndim, int nlen, double* nl, // dim cantidad y lista de puntos (x0y0[z0]x1y1[z1].....)
    int* el  // lista de simplices base 1 terminada en 0 (pede ser nulo)
  );
  malla(const cascara&); // cascara->malla
  // corte (cada step elementos)
  malla(malla &m, const punto &p0, const punto &normal, pline *mapelm=0, int step=1); // plano
  malla(malla &m, double valor, pline *mapelm=0, int step=1); // isosuperficie
  
  ~malla();

  void ini();// inicializa
  void e_ini();// inicializa elementos y sus datos asociados
  malla &operator +=(malla &b); // agrega una malla (no const por el bbox)
  malla &operator =(const malla &b);
  malla &roba(malla &b); // roba el contenido de b que queda vacia
  bool agrega_elementos(const array1<elemento> &el);
  bool agrega_elementos(int nel,const elemento *el);

  // mensajes
  void ini_error();
  void add_error(const char*);
  void add_error(double);
  void add_error(punto);
  void add_error(int);
  void ini_warning();
  void add_warning(const char*);
  void add_warning(double);
  void add_warning(punto);
  void add_warning(int);

  // io
  void copy_nombre(const malla&);
  void renombra(const char* inombre);

  bool lee(const char *archivo=0);
  bool agrega(const char *archivo=0);
  bool graba(const char *archivo=0);
  // al grabar la frontera hay muchisimas opciones
  // todas las piezas o archivos separados por pieza o separados por flag
  // numeracion de nodos propia o la de la malla original
  // flags o no, tipo de archivo, etc. No esta todo implementado, ver
//  bool graba_frontera(const char *ext=0, bool todojunto=true); //uno o varios archivos
  bool graba_frontera(); //varios archivos
  bool graba_dxf_nnod(const char *archivo=0);
  bool graba_calidad(const char *archivo=0);

  // formatos de i/o
  bool lee_xy (const char*, array1<nodo>&); // puntos plano
  bool lee_xyz(const char*, array1<nodo>&); // puntos 3d
  bool lee_dat(const char*, array1<nodo>&, array1<elemento>&); // propio
  bool lee_msh(const char*, array1<nodo>&, array1<elemento>&); // GID
  bool lee_mai(const char*, array1<nodo>&, array1<elemento>&); // SAMCEF
  bool lee_e  (const char*, array1<nodo>&, array1<elemento>&); // OOFELIE
  bool lee_key(const char*, array1<nodo>&, array1<elemento>&); // LS-DYNA
  bool lee_dxf(const char*, array1<nodo>&, array1<elemento>&); // AutoCAD
  bool lee_nas(const char*, array1<nodo>&, array1<elemento>&); // NASTRAN
  bool lee_ans(const char*, array1<nodo>&, array1<elemento>&); // ANSYS
  bool lee_stl(const char*, array1<nodo>&, array1<elemento>&); // STL

  bool agrega_con(const char* tipo_de_elemento, const char* archivo); // solo conectividades

  bool graba_xyz(const char *archivo=0);
  bool graba_dat(const char *archivo=0);
  bool graba_msh(const char *archivo=0);
  bool graba_mai(const char *archivo=0);
  bool graba_key(const char *archivo=0);
  bool graba_con(const char *archivo=0);
  bool graba_dxf(const char *archivo=0, double esize=1);
  bool graba_wrl(const char *archivo=0); // vrml solo superficies

  // malla (cascaras) frontera
  // para saber si un nodo o elemento es frontera alcanza con mk_vecinos()
  bool mk_frontera(bool remake=false);
  bool es_frontera(
    int nod,
    int &ilim,
    int &ix) const; // es de algun limite?
  inline bool es_frontera(int i) const
    {return n[i].f.es(n_frontera);}
  void rm_frontera() {frontera.ini();}
  // elimina elementos exteriores a la frontera previa
  bool frontera_previa(const malla &previa);

  // vecinos por cara de cada elemento (negativo: frontera)
  // (2 y 3D: si hay malla de frontera: -frontera-1; si no: -1)
  bool mk_vecino(bool remake=false);
  bool mk_vecino_1d(bool remake=false); // -1 frontera; -2 bifurcacion
  void rm_vecino() {vecino.ini();}

  // graba piezas conexas por separado
  bool piezas_conexas();

  // unifica orientacion en cada pieza
  bool orienta(bool remake=false,bool mkdir_elm=false);
  // normaliza abiertas
  bool orienta_abierta
    (bool remake=false, bool mkdir_elm=false,int primernodo=0);
  // normaliza cerradas (numf es pieza conexa de cada nodo)
  bool orienta_cerrada
    (bool dirext=true,bool remake=false,int* numf=0,bool mkdir_elm=false,bool reentry=false);
  void swap(int); // invierte elemento y vecinos
  void swap(); // invierte todos los elementos y vecinos

  // normales
  bool mk_dir_elm(bool remake=false); // si esta orientada, solo hace versores normales por elemento
  void rm_dir();  // elimina normles (elemento y nodo)
  // versor normal por nodo para los nodos deben marcados con n_offset (modulo 0 para el resto)
  // marca con e_offset a los elementos cuyos nodos son todos n_offset
  // mata las componentes normales a los elementos no marcados con e_offset
  // delta + hacia afuera
  bool mk_dir_nod(bool remake=false);
  bool offset(double delta); // ofset segun dir de nodos
  // para blayer los deltas son distancias entre capas
  // positivo hacia adentro (al reves que la normal)
  bool blayer(
    int capas, const double *delta, // delta entre capas
    bool archiva_normales=false, // graba normales como malla de segmentos 
    array1<punto>* normal=0, pline* mapnod=0); // array de normales en nodos indicados por map

  // proyecta al plano definido por la normal, 
  // si la normal es 0 por el plano de la normal promedio (sin ponderar)
  // luego lleva la malla al plano xy
  void proyecta_xy(punto normal=pzero);

  // distancia entre nodos
  bool mk_h_nn(bool remake=false){return mk_h_nn_min(remake);}
  bool mk_h_nn_min(bool remake=false);  // distancia al nn mas cercano
  bool mk_h_nn_med(bool remake=false);  // distancia media a los nn
  bool mk_h_nn_max(bool remake=false);  // distancia al nn mas lejano (desrefina)
  double h_nn_min(int i,bool mknn=false) const;
  double h_nn_med(int i,bool mknn=false) const;
  double h_nn_max(int i,bool mknn=false) const;
  bool mk_h_min_max_arista(bool remake=false);// minima(maxima arista) de elementos del nodo
  double h_min_max_arista(int i) const;
  bool filter_nodes(double dmax=0); // elimina nodos cercanos

  // agrega un nodo interpolado y develve el indice
  int nodefrom(int n1, int n2, double f01=.5); // entre otros dos
  int nodefrom(int ie, const double *ff){ // de un elemento
    return nodefrom(e[ie].nv(),e[ie].n,ff);
  }
  int nodefrom( // dependiendo de len nodos de nn, con pesos ff
      int len,      // cantidad de nodos parent
      const int *nn,      // indices de los parent
      const double *ff);  // incidencia de cada parent
  void combine(int vuela, int queda); // de dos nodos proximos queda uno

  // nodos vecinos de cada nodo (~ vecinos naturales)
  void nn1(int in, pline *poner_aqui=0, bool remake=false); // nodos de los elementos de un nodo
  bool mk_nn(bool remake=false);
  void rm_nn(){nn.ini();}
  // vecinos y vecinos de los vecinos
  // limitado por cantidad de layers o cantidad de nodos
  bool vecindad_de_nodo(array1<pline> &vecindad,int nlayers=0,int nnodos=0);

  bool n_arista(double grados=40); //marca nodos de aristas geometricas
  bool graba_aristas(double grados=40, const char* extension="DXF"); //graba aristas geometricas
  bool smooth_surface(); // empareja el h real en fronteras

  // triangulos 2D y 3D
  // metodo 1=delaunay  2=maximo angulo  3=maxima minima altura
  bool diagonal_swap(int metodo=1);
  bool diagonal_swap(int ie1, int ie2, int metodo=1); // metodo 0=incondicional
  
  // suavizado laplaciano (factor=1/pasos)
  bool laplace_smooth(double factor=1);
  bool laplace_smooth(int in,double factor=1);
  
  // empareja el jacobiano, clim=cociente aceptable min/med
  // ps[0]=p0 y ps[0]=n(unitario) del plano de simetria
  // si hay plano de simetria los nodos deben tener n_simetria y los inmoviles n_permanente
//  bool suavej(double clim=.25,const punto *ps=0);
  void j_min_med(int in,double &jmin,double &jmed) const; //de los triedros del nodo

  // slivers
  bool es_sliver(int ie, double h) const; // no marca el elemento como sliver ni reordena
  bool es_sliver(int ie) const{
    return es_sliver(ie,Min(n[e[ie][0]].h,Min(n[e[ie][1]].h,Min(n[e[ie][2]].h,n[e[ie][3]].h))));
  }
  bool engorda_slivers(); // engorda slivers
  void rm_f_slivers(double min_vol_h3r=.1); // elimina slivers con dos caras de frontera (vol/h^3)
  bool mark_slivers(bool remake=false); // los estrictamente slivers: circulares
  bool mark_sliver(int ie); // marca y reordena

  bool diam_aspect(double cota_superior=2); // pone en nodevalues max/min diametro de elementos

  //caracteristicas de elementos
  void bbox(int ne, punto p[2]) const;
  double volumen(int i) const; // volumen del elemento
  double volumen(const elemento&) const; // volumen del elemento
  punto gp(const elemento &) const; // c. de g. de nodos
  punto gp(int ne) const {return gp(e[ne]);}
  punto gv(const elemento &) const; // c. de g. de volumen
  punto gv(int ne) const {return gv(e[ne]);}
  void normal(int ie, int ic, punto &gc, punto &nc) const; // centro y normal de cara
  double jaco(int ie, int in) const; // jacobiano en nodo triedrico (solo 3D)
  // modo de angulos: 0=grados, 1=rad, 2=sen, 3=cos
  double angulo(int ie,int i, int modo=0) const; // angulo 2d en nodo i [0,360)
  // diedros en 3d en rads (no esta para wedges ni cubos)
  double angulo(int ie,int ic1, int ic2, int modo=0) const; // diedro c1 c2 [0,360)
  bool angulos(int ie, double &amin, double &amax); // angulos o diedros (no const)
  double maxd(int ie,double &amax, int &c1, int &c2); // caras con maximo diedro
  double maxd(int ie); // maximo diedro
  double aspect(int) const; // aspect de un elemento (min/max dist e/nodos, o -1 si max=0)
  bool mk_esferas(bool remake=false, bool mkdir=false); // esfera (y normal?) de cada elemento
  bool esfera_e(int i, punto &c, double &r, double &v, punto* normal=0);// return r!=0
  void rm_esferas(){ce.ini(); re.ini();} // no borra volumen ni normales
  void vuelae(int,bool mark_f=true);// elimina swappeando con el ultimo y marca frontera
  bool fforma(   // funcion de forma
    const elemento &ei,      // elemento
    const punto &p,    // punto de calculo
    double *f,   // funciones de forma en el punto
    punto *Df=0  // gradientes (0 => no calcula gradientes)
    ) const;
  bool fforma(   // funcion de forma
    int ie,      // indice del elemento
    const punto &p,    // punto de calculo
    double *f,   // funciones de forma en el punto
    punto *Df=0  // gradientes (0 => no calcula gradientes)
    ) const {return fforma(e[ie],p,f,Df);}
  void pi(int ie, punto *pi, double *vi) const;  // punto y volumen de integracion

  // caracteristicas de la malla
  double volumen() ; // volumen de la malla
  void bbox(bool remake=false); // bbox
  double epsilon_bb(); // hace epsilon a partir del tamanio del boundig box

  // elementos que tienen qty nodos cumpliendo con f
  // en ex pone los elementos y en ix los indices (alguno puede ser 0)
  // los ix y ex se corresponden con los de la malla en forma monotona creciente
  // si qty es -1 son los que tienen todos sus nodos con f
  // si almenos=true testea >= y si es false testea que sea exactamente =
  int filtra_elms(array1<elemento> *ex,pline *ix,flagtype f,int qty=-1, bool almenos=false) const;
  // elementos que tienen todos sus nodos con algun bit de f encendido
  // pero todos el/los mismo/s bit/s
  int mask_elms(array1<elemento> *ex,pline *ix,flagtype f) const;

  // superficie a partir de puntos
  //  baja la dimension de los elementos de la malla 
  //  (convierte elementos solidos en elementos 2d o 2d en segmentos)
  bool baja_dim();
  // Elimina elementos con aristas largas (superficie a partir de puntos)
  bool elimina_aristas_largas();
  // Elimina aristas muy chichas sin nodos de frontera y solo con simplices (blayer)
  bool colapsa_aristas_chicas();

  // si un punto es interior a un elemento o a la malla
  bool inter(const punto &p, int ie);
  bool inter(const punto &p, int ie, int &nodo_cercano);
  bool inter(const punto &p) {return (dequien(p)>=0);} // interior de la malla
  // a que elemento pertenece (<0 -> ninguno) 
  // ie es un elemento cercano conocido (de partida) y el resultado
  // si no se conocen deben ser -1
  void dequien(const punto &p,int &ie, int &nodo_cercano);
  int  dequien(const punto &p);

  void mk_octree(bool remake=false);
  bool weld(double factor); // rehace el octree eliminando cercanos factor*hmin
  bool elm_rep(int in); // testea y elimina elementos repetidos del nodo in

  // swap de nodos y elementos
  void swap_e(int,int);
  void swap_n(int,int);

  // generacion por delaunay
  bool mesh1d(bool parcial=false); // subdivide segmentos
  bool f_delaunay();
  bool delaunay(bool eliminar_slivers=false);
  bool delaunay_refinado(malla &vieja,double alpha,bool eliminar_slivers=false);// delaunay adaptado a h
  bool delaunay_refinado(double alpha, array1<double>& new_h, pline& nodes_from, array1<double>& shape_functions);// delaunay adaptado a h elimina slivers agregando nodos
  // genera los puntos interiores
  //h normalmente es interpolado de los nodos de la esfera
  //si h_radio ==> h=factor*radio de la esfera (ver en voronoi.cpp)
  bool mk_puntos(bool eliminar_slivers=false, bool h_radio=false);
  // hace la frontera por alpha shape como multiplo de h
  void alpha_shape(double alpha=1.5);  
  // solo saca elementos "desde fuera hacia adentro" es decir que no genera huecos
  void alpha_shape_carving(double alpha);
  // no remalla una parte fija (solo triangulos),
  // pone n y e fijos al ppio, devuelve nro e fijos
  int mk_mfija(flagtype,bool primera=false);

  // interpola
  bool interpola(malla &bk);

  // poliedriza
  bool edt(double delta=.01);
  void juntae(int,int); // junta elementos
  bool p2t(int i); // poliedro a tetraedros (solo poligonos/poliedros por ahora)
  void p2t_all(); // todos los poliedros a tetraedros

  // cuadrilateros a triangulos (muy similar a swap2d)
  bool q2t(int i, int k); // quad i a triangulos por el vertice k
  bool q2t(int metodo=1); // todos (1=delaunay, 2=parte el angulo mayor, 3=maxima menor altura de triangulos)
  bool q2t_d(int i); // delaunay
  bool q2t_a(int i); // parte el angulo mayor
  bool q2t_h(int i); // maxima menor altura de triangulos

  bool h5t(); // hexa -> 5 tetras

  bool extrude(const punto &v); // extrusion con vector
  bool extrude(const char* tabfile); // con malla de segmentos

  bool custom(); // customizacion ad-hoc (normalmente no hace nada)

  bool fitsurface(); // superficie dentro de una polilinea 3d cerrada

} mesh;

#define _INI_MALLA\
  m_error(0), m_warning(0),\
  nodosh(0), o(0), epsilon(ERRADM),\
  hayh(false), hayv(false), hayfn(false),\
  hmin(MAXREAL), nvmin(0), nvmax(0),\
  hayfe(false), ndir(false), edir(false)

#endif
