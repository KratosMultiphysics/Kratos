#ifndef _ESFERA
#define _ESFERA

#include "utiles.h" // ciclo Swap
#include "array1.h" // container de nodos
#include "flags.h"
#ifdef _CALCULO_
#include "nodedata.h"
#else
#include "nodo.h"
#endif

// Si tiene radio 0 es porque los puntos son coplanares (colineales en 2d)

class esfera
{
 public:
  static const array1<nodo> *nod; // los nodos
  static double s_epsilon;
  static int NV; // la cantidad de nodos (dim+1)

  double r; // radio, si la esfera se degenera r=0
  punto c;  // centro
  int n[4]; // puntos de definicion (indices de nod)
  double vt;  // volumen del tetraedro de definicion
  flagtype f; // flag multiproposito
  int vecino[4]; // esferas vecinas

  // setea los static
  static void setclass(int dim, const array1<nodo> &p);

  esfera() : r(0),vt(0),f(0) {}; // n y vecino con garbage
//  esfera(const int ix[4]); // f=0 y vecinos=-1
  esfera(const esfera &s) {*this=s;}

  ~esfera() {};

  esfera &operator=(const esfera&); // deep copy
  int &operator[](int i) {return n[ciclo(i,NV)];} // supongo n
  int operator[](int i) const {return n[ciclo(i,NV)];} // supongo n

  // calcula centro, radio y volumen, pero no toca f
  esfera &define(const int ix[4]);
  esfera &define();

  //intercambia dos nodos y cambia la orientacion y el volumen
  void swap() {Swap(n[NV-1],n[NV-2]); Swap(vecino[NV-1],vecino[NV-2]); vt=-vt;}
  void swap(int i, int j) {Swap(n[i],n[j]); Swap(vecino[i],vecino[j]); vt=-vt;}

  // si el punto esta dentro de la esfera
  bool have(const punto& p) const; //depende de la precision numerica
  bool have_i(const punto& p) const; //si esta en la sup no esta dentro
  bool have_e(const punto& p) const; //si esta en la sup esta dentro
  // si esta cerca de la superficie (en terminos absolutos),
  // mueve el punto una distancia epsilon
  // hacia afuera o hacia adentro (justo => adentro)
  // divide epsilon/2 para no modificar las inclusiones previas
  bool have_moving(punto &p);

  // si tiene al nodo ix entre los nodos de definicion
  bool have(int ix) const
    {return (n[0]==ix||n[1]==ix||n[2]==ix||n[3]==ix);}
  bool tiene_vecino(int ix) const
    {return (vecino[0]==ix||vecino[1]==ix||vecino[2]==ix||vecino[3]==ix);}

  // si el punto p esta en el tetraedro (justo o dentro), lejano es el vertice mas alejado
  bool have_t(const punto &p, int *lejano=0);

  // indice de un nodo en la esfera (o -1 si no esta)
  int index(int ix) const
    {return (n[0]==ix)? 0 :
            (n[1]==ix)? 1 :
            (n[2]==ix)? 2 :
            (n[3]==ix)? 3 : -1;}
  int index_vecino(int ix) const
    {return (vecino[0]==ix)? 0 :
            (vecino[1]==ix)? 1 :
            (vecino[2]==ix)? 2 :
            (vecino[3]==ix)? 3 : -1;}

  int replace(int iviejo,int inuevo); // reemplaza iviejo por inuevo
  int replace_vecino(int iviejo,int inuevo); // reemplaza iviejo por inuevo

  // si la esfera es no degenerada
  inline operator bool() const {return (r!=0);}

  bool fforma( // funcion de forma
    const punto &p,    // punto de calculo
    double *f,   // funciones de forma en el punto
    punto *Df=0  // gradientes (0 => no calcula gradientes)
    ) const;
};

#endif //_ESFERA
