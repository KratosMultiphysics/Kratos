// frontera de la malla o de la nube

//#define temporizar
#include "tiempo.h"

#include "voronoi.h"

using namespace std;

// esferas externas y nodos frontera por alpha shapes
// para h variable, alpha es un factor
void voronoi::alpha_shape(double alpha, bool nodelete){
  _initime;
  array1<nodo> &n=m->n;
  mk_h_nn(); // si no hay h lo hace

  int i,j,v;
  double hi=1,r,rc;

//  for (i=0;i<n.len;i++) n[i].f.reset(n_frontera); !!no, ya hay frontera

  // vuela las grandes
  for (i=0;i<s.len;i++) {
    esfera &ei=s[i];
    if (ei.f.es(e_borrado)) continue;
    const int *nei=ei.n;
    for (hi=n[nei[0]].h,j=1;j<NV;j++) hi+=n[nei[j]].h; hi/=NV; // h medio
//    for (hi=n[nei[0]].h,j=1;j<NV;j++) set_min(hi,n[nei[j]].h);// h min
    r=ei.r; if (r==0) rc=MAXREAL; else rc=r/hi; // (0=>colineales o coplanares => vuela)
    if (rc<=alpha) continue;
    // la esfera es grande
    ei.f.reset(e_frontera);
    if (nodelete) {
      for (j=0;j<NV;j++) n[nei[j]].f.set(n_frontera);
      ei.f.set(e_exterior);
      for (j=0;j<NV;j++) {
        v=ei.vecino[j];
        if (v==-1) continue;
        if (v<0) v=-v-2;
        else s[v].f.set(e_frontera);
        s[v].replace_vecino(i,-i-2);
      }
      continue;
    }
    // vuela
    for (j=0;j<NV;j++) {
      nodo &nj=n[nei[j]];
      nj.e.remove1(i);
      nj.f.set(n_frontera);
      v=ei.vecino[j]; if (v<0) continue;
      s[v].replace_vecino(i,-1);
      s[v].f.set(e_frontera);
    }
    // swap con el ultimo (la ultima no puede ser borrada)
    int last=s.len-1;
    if (i==last){s.remove(last); continue;}
    const esfera &el=s[last];
    const int *nel=el.n;
    for (j=0;j<NV;j++) {
      n[nel[j]].e.replace1(last,i);
      v=el.vecino[j];
        if (v==-1) continue;
        if (v<0) v=-v-2;
        s[v].replace_vecino(last,i);
    }
    s[i]=s[last]; s.remove(last);
    i--;
  }
  _savetime(alpha);
}

// solo saca elementos "desde fuera hacia adentro" es decir que no genera huecos
void voronoi::alpha_shape_carving(double alpha, bool nodelete){
  _initime;
  array1<nodo> &n=m->n;
  mk_h_nn(); // si no hay h lo hace

  int i,j,v;
  double hi=1,r,rc;

  // listas de mantenimiento, flag1 rehacer, flag2 procesado
  pline rehacer(s.len),procesado(s.len); 
  // presupongo frontera marcada en elementos => e_frontera -> rehacer
  for (i=0;i<s.len;i++) if (s[i].f.es(e_frontera)) {rehacer+=i; s[i].f.set(flag1);}

  // vuela las grandes
  
  while (rehacer){
    i=rehacer.last(); rehacer.len--;     
    esfera &ei=s[i]; 
    if (ei.f.es(e_borrado)) continue;
    ei.f.reset(flag1);
    const int *nei=ei.n;
    for (hi=n[nei[0]].h,j=1;j<NV;j++) hi+=n[nei[j]].h; hi/=NV; // h medio
//    for (hi=n[nei[0]].h,j=1;j<NV;j++) set_min(hi,n[nei[j]].h);// h min
    r=ei.r; if (r==0) rc=MAXREAL; else rc=r/hi; // (0=>colineales o coplanares => vuela)
    if (rc<=alpha) continue;
    // la esfera es grande
    ei.f.reset(e_frontera);
    if (nodelete) {
      for (j=0;j<NV;j++) n[nei[j]].f.set(n_frontera);
      ei.f.set(e_exterior);
      for (j=0;j<NV;j++) {
        v=ei.vecino[j];
        if (v==-1) continue; // frontera sin nada
        if (v<0) v=-v-2; // tetraedro exterior
        else {// tetraedro interior contra la frontera
          s[v].f.set(e_frontera);
          if (s[v].f.noes_ninguno(flag1|flag2)) {rehacer+=v; s[v].f.set(flag1);}
        }
        s[v].replace_vecino(i,-i-2);
      }
      continue;
    }
    // vuela
    for (j=0;j<NV;j++) {
      nodo &nj=n[nei[j]];
      nj.e.remove1(i);
      nj.f.set(n_frontera);
      v=ei.vecino[j]; if (v<0) continue;
      s[v].replace_vecino(i,-1);
      if (s[v].f.noes_ninguno(flag1|flag2)) {rehacer+=v; s[v].f.set(flag1);}
      s[v].f.set(e_frontera);
    }
    // swap con el ultimo (la ultima no puede ser borrada)
    int last=s.len-1;
    if (i==last){s.remove(last); continue;}
    const esfera &el=s[last];    
    if (el.f.es(flag1)) rehacer.replace1(last,i); // si el esta en rehacer reemplaza
    const int *nel=el.n;
    for (j=0;j<NV;j++) {
      n[nel[j]].e.replace1(last,i);
      v=el.vecino[j];
        if (v==-1) continue;
        if (v<0) v=-v-2;
        s[v].replace_vecino(last,i);
    }
    s[i]=s[last]; s.remove(last);
    i--;
  }
  _savetime(alpha);
}

/////////////////////////////////////////////////////////////////////
// idem pero en malla
void malla::alpha_shape(double alpha){
  _initime;
  mk_h_nn(); // si no hay h lo hace de los nn
  if (!re) mk_esferas();

  int i,j,nv;
  double hi,r,rc;
  bool alguno=false;

//  for (i=0;i<n.len;i++) n[i].f.reset(n_frontera); !!no, ya hay frontera

  // vuela las grandes
  for (i=e.len-1;i>=0;i--) { // al reves porque vuelae swappea
    elemento &ei=e[i];
    const int *nei=ei.n; nv=ei.nv();
    for (hi=n[nei[0]].h,j=1;j<nv;j++) hi+=n[nei[j]].h; hi/=nv; // h medio
//    for (hi=h_nn[nei[0]],j=1;j<nv;j++) set_min(hi,h_nn[nei[j]]);// h min
    r=re[i]; if (r==0) rc=MAXREAL; else rc=r/hi; // (0=>colineales o coplanares => vuela)
    if (rc<alpha) continue;
    vuelae(i,true); // marca frontera
    alguno=true;
    // i++; //ojo si va de 0 a len hay que retestear
  }
  if (alguno) {frontera.ini();nn.ini();}
  _savetime(alpha);
}

// solo saca elementos "desde fuera hacia adentro" es decir que no genera huecos
void malla::alpha_shape_carving(double alpha){
  _initime;
  mk_h_nn(); // si no hay h lo hace de los nn
  if (!re) mk_esferas();

  int i,j,nv,nc;
  double hi,r,rc;
  bool alguno=false;
  // listas de mantenimiento, flag1 rehacer, flag2 procesado
  pline rehacer(e.len),procesado(e.len); 
  // presupongo frontera marcada en elementos => e_frontera -> rehacer
  for (i=0;i<e.len;i++) if (e[i].f.es(e_frontera)) {rehacer+=i; e[i].f.set(flag1);}

  // vuela las grandes
  while (rehacer){
    i=rehacer.last(); rehacer.len--; elemento &ei=e[i]; ei.f.reset(flag1);
    const int *nei=ei.n; nv=ei.nv();
    hi=n[nei[0]].h; for (j=1;j<nv;j++) hi+=n[nei[j]].h; hi/=nv; // h medio
//    for (hi=h_nn[nei[0]],j=1;j<nv;j++) set_min(hi,h_nn[nei[j]]);// h min
    r=re[i]; if (r==0) rc=MAXREAL; else rc=r/hi; // (0=>colineales o coplanares => vuela)
    if (rc<alpha) {procesado+=i; ei.f.set(flag2); continue;}
    const cpline &v=vecino[i]; nc=ei.nc();
    for(j=0;j<nc;j++){
      if (v[j]<0) continue; // frontera
      if (e[v[j]].f.es(flag1|flag2)) continue; // ya estaba
      rehacer+=v[j]; e[v[j]].f.set(flag1);
    }
    // vuelae va a swappear este con el ultimo
    // si el ultimo esta en rehacer hay que cambiar el numero por el actual.
    if (i!=e.len-1&&e[e.len-1].f.es(flag1))
      rehacer.replace1(e.len-1,i);
    vuelae(i,true); // marca frontera
    alguno=true;
  }
  // restaura flag2  
  for (i=0;i<procesado.len;i++) e[procesado[i]].f.reset(flag2);

  if (alguno) {frontera.ini();nn.ini();}
  _savetime(alpha);
}

