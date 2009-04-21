//#define temporizar
#include "tiempo.h"

#include "voronoi.h"

using namespace std;

// elimina virtuales y hace squeeze a la vez;
void voronoi::rm_virtuales(){
  _initime;
  array1<nodo> &n=m->n;

  if (!nvirt||n.last().f.noes(n_virtual)) return;

  int j,is,isl,in,v,newlen=n.len-nvirt;
  static const int f=e_borrado|e_virtual;
  ordlist &sb=s.borrado;

  for (is=s.len-1;is>=0;is--){
    esfera &si=s[is]; 
    if (si.f.noes_ninguno(f)) continue;
    if (si.f.es(e_virtual)){
      const int *ns=si.n;
      for (j=0;j<NV;j++) {
        in=ns[j];
        if (in<newlen) {n[in].e.remove1(is); n[in].f.set(n_frontera);}
        v=si.vecino[j]; if (v<0) continue; //ojo si hay indices de esferas exteriores <-1
        esfera &sv=s[v]; if (sv.f.es_alguno(f)) continue;
        sv.f.set(e_frontera);
        sv.replace_vecino(is,-1);
      }
    }
    // cambia por el ultimo (si no es este)
    if (is!=s.len-1) {
      isl=s.len-1; esfera &sl=s[isl]; // la ultima no es borrada
      const int *nsl=sl.n;
      for (j=0;j<NV;j++){
        n[nsl[j]].e.replace1(isl,is);
        v=sl.vecino[j]; if (v>=0) s[v].replace_vecino(isl,is);
      }      
      if (si.f.es(e_borrado)) {// si es borrado es el ultimo borrado
        sb.len--; s.es_borrado[is]=false;
      }
      si=sl;
    }
    s.len--;
    while (sb.len&&sb.last()==s.len-1) {sb.len--; s.len--; is--;}
  }

  n.len=newlen;
  _savetime(rm_virtuales);
}

/*
void voronoi::rm_virtuales(){
  _initime;
  array1<nodo> &n=m->n;

  if (n.last().f.noes(n_virtual)) return;

  int i,j,is,in,v,snlen,newlen=n.len-NV;

  pline nf(n.last().e.len);

  // los NV ultimos nodos son virtuales
  for (i=newlen;i<n.len;i++){
    cpline &sn=n[i].e;
    while (sn.len){
      is=sn.last(); sn.len--; esfera &si=s[is]; 
      if (si.f.es(e_borrado)) continue;
      const int *ns=si.n;
      for (j=0;j<NV;j++) {
        in=ns[j];
        if (in<newlen) {
          if (n[in].f.noes(n_frontera)) {
            nf+=in;
            n[in].f.set(n_frontera);
          }
        }
        else {// solo si el nodo es virtual cambia vecino
          v=si.vecino[j]; if (v<0) continue; //ojo si hay indices de esferas exteriores <-1
          s[v].f.set(e_frontera);
          s[v].replace_vecino(is,-1);
        }
      }
      si.f.set(e_borrado); s.remove(is);
    }
    sn.ini();
  }

  cpline newsn(32);
  for(i=0;i<nf.len;i++){
    cpline &sn=n[i].e; snlen=sn.len;
    newsn.clean();
    for (j=0;j<snlen;j++)
      if(s[sn[j]].f.noes(e_borrado)) newsn+=sn[j];
    sn=newsn;
  }

  n.len=newlen;
  _savetime(rm_virtuales);
}

*/
// elimina esferas borradas
void voronoi::squeeze(){
  _initime;
  array1<nodo> &n=m->n;

  int k,is,isl,v;

  ordlist &sb=s.borrado;
  while (sb.len){
    is=sb.last(); isl=s.len-1; esfera &sl=s[isl]; // la ultima no es borrada
    const int *nsl=sl.n;
    for (k=0;k<NV;k++) {
      n[nsl[k]].e.replace1(isl,is);
      v=sl.vecino[k]; if (v>=0) s[v].replace_vecino(isl,is);
    }
    s[is]=sl;
    s.len--; sb.len--; s.es_borrado[is]=false;
    while (sb.len&&sb.last()==s.len-1) {sb.len--; s.len--;}
  }
  _savetime(s_squeeze);
}

