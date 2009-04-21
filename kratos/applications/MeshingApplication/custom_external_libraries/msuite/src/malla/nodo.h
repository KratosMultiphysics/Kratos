#ifndef _nodo_H_
#define _nodo_H_

#include "punto.h"
#include "flags.h"
#include "pline.h"

class nodo : public punto{
 public:
  flagtype f; // flag
  cpline e;   // elementos
  double h; // h esperado
  double v; // un real cualquiera asociado al nodo

  // constructores
  nodo() : f(0),h(0),v(0){}
  nodo(double xx, double yy, double zz=0):punto(xx,yy,zz),f(0),h(0),v(0){}
  nodo(const punto& o):punto(o),f(0),h(0),v(0){}
  nodo(const double* xx):punto(xx),f(0),h(0),v(0){}
  nodo(const nodo& o):punto(o.x),f(o.f),e(o.e),h(o.h),v(o.v){}
  // interpolados (flags oreados)
  nodo(const nodo& n1, const nodo &n2, double f01=.5){
    x[0]=(1-f01)*n1.x[0]+f01*n2.x[0]; 
    x[1]=(1-f01)*n1.x[1]+f01*n2.x[1]; 
    x[2]=(1-f01)*n1.x[2]+f01*n2.x[2];
    h=(1-f01)*n1.h+f01*n2.h; v=(1-f01)*n1.v+f01*n2.v;
    f=n1.f|n2.f;
  };
  nodo(int qn, const nodo** n, const double *ff){ // ff debe sumar 1
    x[0]=ff[0]*n[0]->x[0]; x[1]=ff[0]*n[0]->x[1]; x[2]=ff[0]*n[0]->x[2];
    h=ff[0]*n[0]->h; v=ff[0]*n[0]->v; f=n[0]->f;
    while(--qn){
      x[0]+=ff[qn]*n[qn]->x[0]; x[1]+=ff[qn]*n[qn]->x[1]; x[2]+=ff[qn]*n[qn]->x[2];
      h+=ff[qn]*n[qn]->h; v+=ff[qn]*n[qn]->v; f|=n[qn]->f;
    }
  };

  nodo& operator =(const nodo &o)
    {punto::operator=(o); e=o.e; f=o.f; h=o.h; v=o.v; return *this;}
  nodo& operator =(const punto &o) // cambia solo posicion
    {punto::operator=(o); return *this;}
  nodo& setpos (const nodo &o) // copia solo la posicion
    {punto::operator=(o); return *this;}
  nodo& setpos (const punto &o) // copia solo la posicion
    {punto::operator=(o); return *this;}

  void ini(){zero(); f=0; h=0; v=0; e.clean();}
};

#endif
