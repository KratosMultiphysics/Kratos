#ifndef _VORONOI_H
#define _VORONOI_H

#include "array1.h"
#include "pline.h"
#include "malla.h"
#include "esfera.h"

// arma una clase voronoi para no dejar las funciones sueltas
class voronoi{
 public:
  malla *m;

  punto *pori; int nori; // cantidad y posicion original de los nodos

  // frontera original para discriminar el interior (hydir.cpp, orienta.cpp)
  int qef,qnf;  // cantidad de nodos y elementos
  int* numf;    // numero de frontera de cada nodo
  elemento *ef; // elementos de frontera
  cpline *efn;  // elementos de cada nodo
  punto *dir;   // normales hacia afuera por elemento
  punto *gef;   // centriodes de elementos
  double *hef;  // h por elemnto
  int nvirt;    // cantidad de nodos virtuales
  
  bool slivers_marked; // esferas marcadas con e_sliver
  
  void prepara(); // perturba quita elementos y flag de frontera
  void restaura(); // repone la posicion de los nodos

  // busca todas las esferas que contengan al nodo
  // y las reemplaza por esferas nuevas con el
  bool cavidad(
    int in, // indice del nodo
    int is=nada, // esfera que contiene al nodo
    int ic=nada, // nodo cercano (si no hay una esfera que lo tenga seguro)
    bool ini_o_clean=false  // asigna o devuelve memoria y sale
  );

  // interiores/exteriores y h
  bool mk_hydir();// genera h y una normal hacia afuera en cada nodo
  void frontera(); // elimina exteriores
  void borra_exteriores(); // borra las esferas externas
  // a que tetraedro pertenece el punto; nc nodo cercano; vlejano: vertice mas alejado
  int de_que_tetra(const punto &p, int &nc, int &vlejano);
  // idem con esfera cercana y ff (agregar nodos)
  int de_que_tetra(const punto &p, int ec, double *ff);
 
  void s2e(); // esferas de voronoi -> elementos (tetra) de la malla
  
  // elimina virtuales y esferas borradas
  void rm_virtuales();
  void squeeze();

//-public de verdad--------------------------

  array1h <esfera> s;
  int NV; // 3 o 4 vertices (2D o 3D)

  // hace la triangulacion o la generacion de puntos segun el caso
  voronoi(malla *M):
    m(M),
    pori(0),nori(0),
    numf(0),ef(0),efn(0),dir(0),gef(0),hef(0),
    slivers_marked(false),
    NV(M->tipo.es(m_planaxy)? 3 :4) {qnf=0;qef=0;};

  voronoi():
    m(0),
    pori(0),nori(0),
    numf(0),ef(0),efn(0),dir(0),gef(0),hef(0),
    slivers_marked(false),
    NV(-1) {qnf=0;qef=0;};

  ~voronoi(){
    delete [] pori;
    delete [] numf;
    delete [] ef;
    delete [] efn; 
    delete [] dir; 
    delete [] gef; 
    delete [] hef;
  }

  // esferas huecas de un grupo de puntos
  // mueve los nodos para evitar coesfericos
  bool delaunay();

  // generador de puntos
  // h normalmente es interpolado de los nodos de la esfera, 
  //   pero si h_radio => h=f(radio de la esfera) (ver en voronoi.cpp)
  // siempre elimina slivers de frntare
  bool mk_puntos(bool h_radio=false); 

  // delaunay de la frontera sin agregar puntos
  bool f_delaunay();

  // hace las esferas sacando los puntos muy cercanos y
  // agrega puntos en las esferas medianas
  bool refina_esferas(malla &vieja,double alpha);

  // Hace las esferas sacando los puntos muy cercanos y
  // agrega puntos en las esferas medianas y en los slivers
  bool refina_esferas(double alpha, array1<double>& new_h, pline& nodes_from, array1<double>& shape_functions);

  // elimina slivers
  bool mark_sliver(int is); // marca y reordena
  bool mark_slivers(bool remake=false);
  void rm_slivers(){rm_f_slivers();rm_i_slivers();}//el orden de los factores es este
  bool rm_i_sliver(int i);
  void rm_i_slivers();
  void rm_f_slivers();// no debe haber exteriores
  bool swap3x2(int,int,int);// arista -> cara
  bool swap4x4(int); //arista->arista

  // elimina esferas grandes
  void alpha_shape(double alpha=1.5, bool no_eliminar_esferas=false);
  void alpha_shape_carving(double alpha=1.5, bool no_eliminar_esferas=false); //sin huecos

  // hace h de con la minima distancia a los nn
  bool mk_h_nn(bool remake=false);

  //io
  bool graba_dat(const char *arch=0);

  friend class malla;
};
#endif
