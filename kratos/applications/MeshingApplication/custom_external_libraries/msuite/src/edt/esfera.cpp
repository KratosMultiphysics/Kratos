//IMPLEMENTACION DE ESFERAS

#include <cstring> // memcpy
#include <cmath> // fabs

#include "esfera.h"

using namespace std;

int esfera::NV=4;
static const int SZAN=4*SZI; // 4 enteros
double esfera::s_epsilon=2048*ERRADM; // el record que obtube fue 1024 (10 divisiones)
const array1<nodo>* esfera::nod=0;

void esfera::setclass(int dim, const array1<nodo> &p){
  NV=dim+1;
  nod=&p;
}

//esfera::esfera(const int ix[4])
//{f=0; memset(vecino,~0,SZAN); define(ix);} // memset solo funciona asi con 0 y -1

// asignacion deep copy
esfera& esfera::operator=(const esfera &s){
  r=s.r;
  c=s.c;
  memcpy(n,s.n,SZAN);
  vt=s.vt;
  f=s.f;
  memcpy(vecino,s.vecino,SZAN);
  return *this;
}


//=================================================================================
// Punto dentro de la esfera:
//=================================================================================
// El punto se mueve epsilon y se divide epsilon por 2, sum(e/2^n)<e (si n>=1)
// Hay varias alternativas cada una con su historia
// El punto se mueve a e de la superficie o a e de la posicion original?
// Ambos casos son "aceptables" y son irreproducibles excepto en el mismo orden,
//   si se retestea en distinto orden solo se garantiza inclusion estricta (sin e)
// De todos modos se eligio mover a e de la superficie.
// Para minimizar slivers puede favorecerse la inclusion (justo => adentro).
//=================================================================================
/*
Esto es muy complicado:
  Si el punto esta justito por dentro de epsilon no se mueve, puede pasar
que los posteriores movimientos (e+e/2...) lo saquen (o viceversa).
  Por eso, si esta entre epsilon y 2e divido epsilon sin mover,
posteriores movimientos (e/2+e/4...) no lo pueden alterar.
  Si esta entre 0 y e, lo muevo a e y divido por 2, posteriores
movimientos (e/2+e/4...) no lo alteran
*/

bool esfera::have_moving(punto &p){
  punto v(p-c); double delta(v.mod()-r);
  // lejos
  if (delta>= 2*s_epsilon) return false; // fuera
  if (delta<=-2*s_epsilon) return true; // dentro

  // entre e y y 2e
  if (delta>= s_epsilon) {s_epsilon/=2; return false;} // fuera
  if (delta<=-s_epsilon) {s_epsilon/=2; return true; } // dentro

  // casi en la superficie
  if (delta>0) {  // fuera => fuera
    p=c+v*((r+s_epsilon)/(r+delta));
    s_epsilon/=2;
    return false;
  }
  else { // dentro o justo => dentro
    p=c+v*((r-s_epsilon)/(r+delta));
    s_epsilon/=2;
    return true;
  }
}

// si el punto esta dentro de la esfera
bool esfera::have(const punto& p) const{
  return (p.eq_b(c,r));
}
bool esfera::have_i(const punto& p) const{ //si esta en la sup no esta dentro
  return (p.eq_b(c,r-s_epsilon));
}
bool esfera::have_e(const punto& p) const{ //si esta en la sup esta dentro
  return (p.eq_b(c,r+s_epsilon));
}


// verifica si un punto es interior o limite del tetraedro
// jmin es el nodo mas alejado
bool esfera::have_t(const punto &pt, int *jmin){
  double f[4]; fforma(pt,f);
  int j=1; double fmin=f[0];
  if (jmin){
    for (*jmin=0;j<NV;j++) if (f[j]<fmin) {fmin=f[j]; *jmin=j;}
  }
  else {
    for (;j<NV;j++) if (f[j]<fmin) fmin=f[j];
  }
  return (fmin>=0);
}

//=================================================================================

esfera& esfera::define(const int *ix){
  n[0]=ix[0]; n[1]=ix[1]; n[2]=ix[2];
  if (NV==4) n[3]=ix[3];
  return define();
}

esfera& esfera::define(){
  if (NV==4){
    s4((*nod)[n[0]],(*nod)[n[1]],(*nod)[n[2]],(*nod)[n[3]],
       c,r,vt);
  }
  else {
    n[3]=-1; punto a;
    c3((*nod)[n[0]],(*nod)[n[1]],(*nod)[n[2]],c,r,a);
    vt=a[2];
  }
  return *this;
}

// reemplaza iviejo por inuevo
int esfera::replace(int iviejo,int inuevo){
  if (n[0]==iviejo) {n[0]=inuevo; return 0;}
  if (n[1]==iviejo) {n[1]=inuevo; return 1;}
  if (n[2]==iviejo) {n[2]=inuevo; return 2;}
  if (NV==4&&n[3]==iviejo) {n[3]=inuevo; return 3;}
  return -1;
}

int esfera::replace_vecino(int iviejo,int inuevo){
  if (vecino[0]==iviejo) {vecino[0]=inuevo; return 0;}
  if (vecino[1]==iviejo) {vecino[1]=inuevo; return 1;}
  if (vecino[2]==iviejo) {vecino[2]=inuevo; return 2;}
  if (NV==4&&vecino[3]==iviejo) {vecino[3]=inuevo; return 3;}
  return -1;
}
