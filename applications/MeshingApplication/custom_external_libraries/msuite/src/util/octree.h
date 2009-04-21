// El octree es la estructura madre que contiene cajaes
// El octree conoce los datos generales e interactua con el usuario
// En cada caja hay un maximo de nodos, cuando se llega al
// maximo el caja se divide en 2^dim childs.

// El principal problema es que pasa cuando hay muchos nodos
// coincidentes (mas que el maximo), si no se trata provoca
// stack overflow. La solucion implementada es:
//   void add(int indice_del_nodo_a_agregar,
//            int &nodo_cercano,
//            bool &puesto,
//            double prec=ERRADM) // ERRADM de utiles.h
// si la distancia al nodo mas cercano es < prec no lo agrega.
// Hay dos variantes: add(...) que calcula "el nodo puesto mas
// cercano" y add_no_rep(...) que calcula el mas cercano solo
// del caja y devuelve ese o, si no hay, uno cualquiera del parent.

// Si es seguro que no se producira stack overflow se puede usar
// otra variante, insegura pero mas rapida: int add(inodo)
// que devuelve un nodo cualquiera del mismo caja (o del parent)

#ifndef _OCTREE_H
#define _OCTREE_H

//#define _OADEBUG // hace un script para autocad de un cuadtree (2D)

#include "punto.h"

// esta implementado para un array1<nodo>
// si se quiere usar otra cosa ej:punto* se pueden cambiar
// los includes y el define de aca abajo
// lo que se use debe devolver un const punto& al pedir (*plist)[ix]

#include "array1.h"
#ifdef _CALCULO_
#include "nodedata.h"
#else
#include "nodo.h"
#endif
#define _OCTREE_PTARRAY array1<nodo>

//#define _OCTREE_PTARRAY punto*
//#define _OCTREE_PTARRAY array1<punto>

class caja;

class octree{

private:
  caja *child; // caja mayor

public:

  const int dim,nchild,maxlen; // dimension, 1<<dim y puntos/caja
  punto pmin, pmax; // bounding caja
  const _OCTREE_PTARRAY *plist; // puntero al array de puntos

  octree(const _OCTREE_PTARRAY &p, // si ya se conoce el bbox
         const punto &ipmin,const punto &ipmax,
         int idim=3, int imaxlen=15);
  octree(const _OCTREE_PTARRAY &plist, int plen, // calcula el bbox
         int idim=3, int imaxlen=15);
 ~octree(); // el destruuctor de caja borra todos los subboxes

  bool check_interior(const punto& p) const; // interior del bbox

  // cercano al punto o, si esta fuera del bbox, cercano al
  // punto mas cercano a p del bbox
  // las distancias son max_delta(x,y,z) (esferas cuadradas)
  int close(punto p) const; // misma caja (o del parent si no hay)
  int closest(punto p, double &d) const; // el punto (puesto) mas cercano
  void add // agrega solo si no coincide dentro de e con cualquiera (O??)
    (int in, int &nc, bool &puesto, double e=ERRADM);
  void add_no_rep // si no coincide dentro de e en la caja (O(log(n)))
    (int in, int &nc, bool &puesto, double e=ERRADM);
  int add(int in);// agrega y devuelve un nodo cercano (O(log(n)))

  void remove(int in); // si no estaba no hace nada
  bool swap_n(int in1, int in2);// intercambia indices
};

#endif
