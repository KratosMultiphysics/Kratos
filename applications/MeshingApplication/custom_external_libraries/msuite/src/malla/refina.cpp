// divide elementos hasta dar con h

#include <cstdlib> // alloc
#include "malla.h"
#include "voronoi.h"


using namespace std;

// elimina nodos a menos de d
// se hace con el octree pero se puede hacer via delaunay parcial
bool malla::filter_nodes(double d){
  if (e||(d<ERRADM&&!hayh)) return false;
  bool puesto;
  int i,nc;
  if (o) delete o;
  o=new octree(n,pmin,pmax,tipo.es(m_planaxy) ? 2 : 3);
  array1<nodo> nuevo(n.len);
  for (i=0;i<n.len;i++) {
    if (d<ERRADM) d=n[i].h;
    o->add_no_rep(i,nc,puesto,d);
    if (puesto) nuevo+=n[i];
    else combine(i,nc);
  }
  n.roba(nuevo);
  return true;
}

// ds/dx=1/h con h lineal => s=ln(h)/h' o bien h/h0=exp(sh')
// El problema es la asignacion de flags
//   Si los elementos tienen flags, cada elemento dividido adquiere el flag del padre
//   Si los nodos tienen flag hay un problema. Podria orear los flags, pero parece mas
//      razonable usar el flag de elemento si es que hay y asignarselo a los nodos.
//   En todo caso, el flag de nodo es "tocable".....
bool malla::mesh1d(bool parcial){
  if (!hayh) {add_error(No_h); return false;}
  ve.ini(); re.ini(); ce.ini();  nn.ini();

  int i,ie0,in0,in1,oldelen=e.len,pasos;
  double h0,h1,h,dh,s,fh,fn,l;
  nodo ni,nf; ni.e.resize(2); ni.e.len=2; nf.e.resize(2); nf.e.len=2; 
  elemento es(e_segmento);
  punto dn;
  bool hvar;

  for (ie0=0;ie0<oldelen;ie0++){
    elemento &e0=e[ie0];
    in0=e0[0];in1=e0[1];
    nodo &n0=n[in0],&n1=n[in1];
    dn=n1-n0;
    h0=n0.h; h1=n1.h; dh=h1-h0; hvar=(fabs(dh)>ERRADM);
    l=dn.mod(); if (hvar) s=l*log(h1/h0)/dh; else s=l/h0;
    pasos=(int)ceil(s); // s va de a uno, pero no da entero
    if (pasos<2) continue;
    es[1]=n.len; e0[1]=n.len; // el primero sigue siendo ei
    es.f=e0.f;
    if (hvar) fh=exp(s/pasos*dh/l); else fh=1;
    h=h0*fh; nf.h=h;
    if (hvar) fn=(h-h0)/dh; else fn=1.0/pasos;
    nf.setpos(n0+dn*fn);
    nf.e[0]=ie0; nf.e[1]=e.len;
    if (hayfe) nf.f=es.f; else nf.f=n0.f|n1.f;
    n+=nf;
    for(i=1;i<pasos-1;i++){
      es[0]=es[1]; es[1]=n.len; ni=nf;
      nf.e[0]=nf.e[1]; nf.e[1]=e.len+1;
      h*=fh; nf.h=h;
      if (hvar) fn=(h-h0)/dh; else fn=double(i+1)/pasos;
      nf.setpos(n0+dn*fn);
      n+=nf;
      e+=es;
    }
    // agrega el ultimo elemento sin el nodo
    es[0]=es[1]; es[1]=in1; n1.e.replace1(ie0,e.len);
    e+=es;
  }
  if (!mk_vecino_1d(true)) return false;

  if (!parcial) return true;
  // los elementos tienen a lo mas sqrt2 hmed de long
  // quita cercanos y genera elementos no menores que h/2
  //      d>(h0+h1)/2/2 = d^2>(h0+h1)^2/16
  ordlist nborrado(500);
  int ieant,iepos,inant,inpos;
  double dant,dpos;
  for (ie0=e.len-1;ie0>=0;ie0--){
    elemento &e0=e[ie0];
    in0=e0[0];in1=e0[1];
    nodo &n0=n[in0],&n1=n[in1];
    double hc=pown(n0.h+n1.h,2)/16;
    if (n0.distancia2(n1)>=hc) continue; // <.5 se quita
    // elimina un nodo
    // busca el vecino mas cercano
    ieant=vecino[ie0][0]; iepos=vecino[ie0][1];
    if (ieant>=0) {
      elemento &eant=e[ieant];
      inant=eant[1-eant.index(in0)]; dant=n0.distancia(n[inant]);
    }else {inant=-1; dant=MAXREAL;}
    if (iepos>=0) {
      elemento &epos=e[iepos];
      inpos=epos[1-epos.index(in1)]; dpos=n1.distancia(n[inpos]);
    }else {inpos=-1; dpos=MAXREAL;}
    if (inant<0&&inpos<0) continue; // ambos limite (podria volar el elm)
    if (dant<dpos) {
      // vuela n0
      e[ieant].replace(in0,in1);
      n1.e.replace1(ie0,ieant);
      n0.f.set(n_borrado); nborrado+=in0;
    }
    else {
      // vuela n1
      e[iepos].replace(in1,in0);
      n0.e.replace1(ie0,iepos);
      n1.f.set(n_borrado); nborrado+=in1;
    }
    // vuela ie0
    if (ieant>=0) vecino[ieant].replace1(ie0,iepos);
    if (iepos>=0) vecino[iepos].replace1(ie0,ieant);
    vuelae(ie0,false); //no marca frontera
  }

  // squeeze de los nodos
  ///////////////ojo: nodo borrado que apunta mal a un elemento cambiado
  if (nborrado&&o) {delete o; o=0;} // el octree quedaria mal
  while (nborrado.len){
    int ib=nborrado.last(),il=n.len-1;
    if (ib==il) {n.len--; nborrado.len--; continue;}
    nodo &nl=n[il]; // nborrado es ordlist (n.last() no es borrado)
    const cpline &enl=nl.e;
    for (int i=0;i<enl.len;i++) e[enl[i]].replace(il,ib);
    n[ib]=n[il]; nborrado.len--; n.len--;      
  }
//graba("borrar.dat");
  return true;
}

#ifndef _CALCULO_

// interpola en el medio de otros dos
// (la pregunta del millon es si orear o no los flags)
int malla::nodefrom(int i1, int i2, double f){
  nodo nm(n[i1],n[i2],f);
  nm.f.reset(n_permanente|n_frontera);
  return n+=nm;
}

int malla::nodefrom( // interpola a partir de otros
   int len,      // cantidad de nodos parent
   const int *v,       // indices de los parent
   const double *ff)  // incidencia de cada parent
{
  const nodo **nlist=(const nodo**)malloc(len*SZP);
  for (int i=0;i<len;i++) nlist[i]=&n[v[i]];
  nodo nm(len,nlist,ff); 
  free(nlist);
  nm.f.reset(n_permanente|n_frontera);
  return n+=nm;
}

// de dos nodos proximos queda uno
void malla::combine(int vuela, int queda){
  nodo &nv=n[vuela];
  nodo &nq=n[queda];
  if (nq.f.es(n_permanente)) return;
  nq.f|=nv.f;
  nq.setpos((nq+nv)/2);
  nq.h=Min(nq.h,nv.h);
  nq.v=(nq.v+nv.v)/2;
}

#endif // _CALCULO_

