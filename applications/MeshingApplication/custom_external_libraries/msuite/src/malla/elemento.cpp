//=================================================================================
// elemento
//=================================================================================
#include <cstdlib> // alloc
#include <cstring> // memcmp memcp .....
#include "utiles.h" // ciclo swap
#include "case.h" //minusc
#include "elemento.h"

using namespace std;
/*e_tipo
e_indefinido,e_punto,e_segmento,e_triangulo,e_cuadrilatero,e_poligono,
e_tetraedro,e_cubo,e_wedge,e_poliedro,e_ntipos 
*/
static const int
  tnv[]={0,1,2,3,4,0,4,8,6,0}; // cantidad de nodos
static const int
  tnc[]={0,0,2,3,4,0,4,6,5,0}; // cantidad de caras
static const int
  tns[]={0,0,0,1,2,0,1,6,3,0}; // cantidad de simplices

e_tipo elemento::etipo(const char* is){
  char s[6]; memcpy(s,is,5); s[5]=0;
  minusc(s); // ordenado por probabilidad y dificultad
  if      (s[0]=='t'){
    if      (s[1]=='e') return e_tetraedro; // tetraedro o tetra o tetrahedral
    else if (s[1]=='r') return e_triangulo; // triangulo o tri o triangle
  }
  else if (s[0]=='p'){
    if      (s[1]=='r') return e_wedge;     // prisma o prism
    else if (s[4]=='e') return e_poliedro;  // poliedro
    else if (s[4]=='h') return e_poliedro;  // polyhedron
    else if (s[4]=='g') return e_poligono;  // polygon o poligono
    else if (s[1]=='u') return e_punto;     // punto
    else if (s[2]=='i') return e_punto;     // point
  }
  else if (s[0]=='c'){
    if      (s[2]=='b') return e_cubo;         // cubo o cube
    else if (s[2]=='a') return e_cuadrilatero; // cuadrilatero o cuadrado o cuadrangular
  }
  else if (s[0]=='h')   return e_cubo;         // hexaedro o hexahedral
  else if (s[0]=='q')   return e_cuadrilatero; // quad o quadrilateral
  else if (s[0]=='s')   return e_segmento;     // segmento o segment
  else if (s[0]=='w')   return e_wedge;        // wedge
  else if (s[0]=='h')   return e_cubo;         // hexaedro o hexahedral
  else if (s[0]=='q')   return e_cuadrilatero; // quad o quadrilateral
  else if (s[0]=='b')   return e_cubo;         // brick o box
  else if (s[0]=='l')   return e_cubo;         // ladrillo
                        return e_indefinido;
}

elemento::elemento(const elemento& e){
  Tipo=e.Tipo;
  int NV=e.nv();

  if (NV>0){
    n=new int[NV]; memcpy(n,e.n,NV*SZI);
  } else n=0;
  f=e.f;
  if (e.pd){
    pd=new poly_dat(*(e.pd));
  } else pd=0;
}

poly_dat::poly_dat(const poly_dat &p){
  nv=p.nv; nc=p.nc; nt=p.nt;
  int NV=3;
  if (p.c){
    c=new int[3*nc]; memcpy(c,p.c,3*nc*SZI); 
    NV=4; // los poligonos no tienen (c=0 es indicador de 2d/3d)
  } else c=0;
  if (p.cn){
    cn=new int*[nv];
    int i,nci;
    for (i=0;i<nv;i++){
      nci=p.cn[i][0]+1;
      cn[i]=new int[nci]; memcpy(cn[i],p.cn[i],nci*SZI);
    }
  } else cn=0;
  
  t=new int[NV*nt]; memcpy(t,p.t,NV*nt*SZI);
}

elemento::elemento(e_tipo t) {
  Tipo=t;
  int NV=tnv[Tipo];
  if (NV) n=new int[NV]; else n=0;
  f=0; pd=0;
}
elemento::elemento(const char* s) {
  Tipo=etipo(s);
  int NV=tnv[Tipo];
  if (NV) n=new int[NV]; else n=0;
  f=0; pd=0;
}

elemento::elemento(e_tipo t, const int* nn) {
  Tipo=t; int NV=tnv[Tipo];
  if (NV) {
    n=new int[NV]; 
    memcpy(n,nn,NV*SZI);
  } else n=0;
  f=0; pd=0;
}

elemento::~elemento(){
  if (n) delete[] n;
  if (pd) delete pd;
}

elemento &elemento::ini(){
  if (n) delete[] n; n=0;
  if (pd) delete pd; pd=0;
  Tipo=e_indefinido;
  f=0;
  return *this;
}

elemento& elemento::copy_nodes(const int* nin){ // nv fijo de antemano  
    memcpy(n,nin,tnv[Tipo]*SZI);
    return *this;
}

e_tipo elemento::tipo(e_tipo t){
  if (Tipo==t&&t!=e_poligono&&t!=e_poliedro) return t;
  // distinto de poly
  if (n) delete[] n; n=0;
  if (pd) delete pd; pd=0;
  Tipo=t; int NV=tnv[t];
  if (NV>0) n=new int[NV];
  return t;
}

elemento& elemento::operator=(const elemento& e){
  if (&e==this) return *this;
  int NV=e.nv();
  if (nv()!=NV){
    if (n) delete[] n; n=0;
    if (NV) n=new int[NV];
  }
  Tipo=e.Tipo;
  if (n) memcpy(n,e.n,NV*SZI);
  if (pd) delete pd; pd=0;
  if (e.pd){
    pd=new poly_dat(*(e.pd));
  }
  f=e.f;
  return *this;
}

// no poli
elemento& elemento::set(e_tipo t, const int* nn){
  int NV=tnv[t];
  if (nv()!=NV){
    if (n) delete[] n; n=0;
    if (NV) n=new int[NV];
  }
  Tipo=t;
  if (n) memcpy(n,nn,NV*SZI);
  if (pd) delete pd; pd=0;
  f=0;
  return *this;
}

bool elemento::operator==(const elemento &e) const{
  if (Tipo!=e.Tipo) return false;
  int NV=nv();
  if (NV!=e.nv()) return false;
  for (int i=0;i<NV;i++) if (!e.have(n[i])) return false;
  return true;
}

int elemento::nv() const {
  return (pd)? pd->nv : tnv[Tipo];
}

int elemento::index(int i) const{
  int j=nv();
  while (--j>=0) if (n[j]==i) break;
  return j;
}

bool elemento::have(const elemento& e, bool strict) const{
  int i=nv(),j=e.nv(); if (i<j||(strict&&i==j)) return false;
  while (--j>=0) if (index(e.n[j])<0) return false;
  return true;
}

// reemplaza iviejo por inuevo
int elemento::replace(int iviejo,int inuevo){
  int j=nv();
  while (--j>=0) if (n[j]==iviejo) {n[j]=inuevo; break;}
  return j;
}

int elemento::nc() const {
  return (pd)? pd->nc : tnc[Tipo];
}

elemento elemento::cara(int i) const{
  static elemento 
    ept(e_punto), es(e_segmento), 
    et(e_triangulo), ec(e_cuadrilatero);
  // para todos los 2d va del nodo al siguiente
  // orientaciones 3D hacia fuera del elemento
  if (Tipo==e_segmento) { // los nodos
    ept.n[0]=n[ciclo(i,2)];
    ept.f=f; return ept;
  }
  if (dim()==2) {
    int NV=nv(); i=ciclo(i,NV);    
    es.n[0]=n[i];
    es.n[1]=n[(i+1)%NV];
    es.f=f; return es;
  }
  else if (Tipo==e_tetraedro) { // cara i opuesta al vertice i
    i=ciclo(i,4);
    et.n[0]=n[(i+1)%4];
    if (i&1){ // impar
      et.n[1]=n[(i+3)%4];
      et.n[2]=n[(i+2)%4];
    }
    else {
      et.n[1]=n[(i+2)%4];
      et.n[2]=n[(i+3)%4];
    }
    et.f=f; return et;
  }
  else if (Tipo==e_cubo) {   
    i=ciclo(i,6);
         if (i==0){ec.n[0]=0; ec.n[1]=3; ec.n[2]=2; ec.n[3]=1;} // 0: plano xy z=0
    else if (i==1){ec.n[0]=0; ec.n[1]=1; ec.n[2]=5; ec.n[3]=4;} // 1: plano xz y=0
    else if (i==2){ec.n[0]=0; ec.n[1]=4; ec.n[2]=7; ec.n[3]=3;} // 2: plano yz x=0
    else if (i==3){ec.n[0]=6; ec.n[1]=7; ec.n[2]=4; ec.n[3]=5;} // 3: plano xy z=1
    else if (i==4){ec.n[0]=6; ec.n[1]=2; ec.n[2]=3; ec.n[3]=7;} // 4: plano xz y=1
    else if (i==5){ec.n[0]=6; ec.n[1]=5; ec.n[2]=1; ec.n[3]=2;} // 5: plano yz x=1
    ec.n[0]=n[ec.n[0]]; ec.n[1]=n[ec.n[1]]; ec.n[2]=n[ec.n[2]]; ec.n[3]=n[ec.n[3]];
    ec.f=f; return ec;
  }
  else if (Tipo==e_wedge) {
    i=ciclo(i,5);
         if (i==0){ec.n[0]=1; ec.n[1]=2; ec.n[2]=5; ec.n[3]=4;} // 0: cuadrilatero opuesto al 0
    else if (i==1){ec.n[0]=2; ec.n[1]=0; ec.n[2]=3; ec.n[3]=5;} // 1: cuadrilatero opuesto al 1
    else if (i==2){ec.n[0]=0; ec.n[1]=1; ec.n[2]=4; ec.n[3]=3;} // 2: cuadrilatero opuesto al 2
    else if (i==3){et.n[0]=0; et.n[1]=2; et.n[2]=1;}            // 3: triangulo inferior
    else if (i==4){et.n[0]=3; et.n[1]=4; et.n[2]=5;}            // 4: triangulo superior
    if (i<3){
      ec.n[0]=n[ec.n[0]]; ec.n[1]=n[ec.n[1]]; ec.n[2]=n[ec.n[2]]; ec.n[3]=n[ec.n[3]];
      ec.f=f; return ec;
    }
    else {
      et.n[0]=n[et.n[0]]; et.n[1]=n[et.n[1]]; et.n[2]=n[et.n[2]];
      et.f=f; return et;
    }
  }
  else if (Tipo==e_poliedro&&pd->nc){
    elemento e(e_triangulo); 
    i=3*ciclo(i,pd->nc);
    e.n[0]=n[pd->c[i]];e.n[1]=n[pd->c[i+1]];e.n[2]=n[pd->c[i+2]];
    e.f=f; return e;
  }
  else return elemento();
}

// indice de una cara dada por los nodos de la malla
int elemento::icara(const elemento &c) const{
  int nvc=c.nv(),i=nc(),j;
  while (--i>=0){
    elemento ci=cara(i);
    if (ci.nv()!=nvc) continue;
    j=nvc; while (--j>=0) if (!ci.have(c[j])) break;
    if (j<0) break;
  }
  return i;
}

void elemento::cara(int i, int &nvc, int *ix) const{
  // para todos los 2d va del nodo al siguiente
  // orientaciones 3D hacia fuera del elemento
  if (Tipo==e_segmento) { // los nodos
    nvc=1; ix[0]=ciclo(i,2);
    return;
  }
  if (dim()==2) {
    int NV=nv(); i=ciclo(i,NV);
    nvc=2; ix[0]=i; ix[1]=(i+1)%NV;
    return;
  }
  else if (Tipo==e_tetraedro) { // cara i opuesta al vertice i
    nvc=3; i=ciclo(i,4);
    ix[0]=(i+1)%4;
    if (i&1){ // impar
      ix[1]=(i+3)%4;
      ix[2]=(i+2)%4;
    }
    else {
      ix[1]=(i+2)%4;
      ix[2]=(i+3)%4;
    }
    return;
  }
  else if (Tipo==e_cubo) {
    nvc=4; i=ciclo(i,6);
         if (i==0){ix[0]=0; ix[1]=3; ix[2]=2; ix[3]=1;} // 0: plano xy z=0
    else if (i==1){ix[0]=0; ix[1]=1; ix[2]=5; ix[3]=4;} // 1: plano xz y=0
    else if (i==2){ix[0]=0; ix[1]=4; ix[2]=7; ix[3]=3;} // 2: plano yz x=0
    else if (i==3){ix[0]=6; ix[1]=7; ix[2]=4; ix[3]=5;} // 3: plano xy z=1
    else if (i==4){ix[0]=6; ix[1]=2; ix[2]=3; ix[3]=7;} // 4: plano xz y=1
    else if (i==5){ix[0]=6; ix[1]=5; ix[2]=1; ix[3]=2;} // 5: plano yz x=1
    return;
  }
  else if (Tipo==e_wedge) {
    i=ciclo(i,5);
         if (i==0){ix[0]=1; ix[1]=2; ix[2]=5; ix[3]=4;} // 0: cuadrilatero opuesto al 0
    else if (i==1){ix[0]=2; ix[1]=0; ix[2]=3; ix[3]=5;} // 1: cuadrilatero opuesto al 0
    else if (i==2){ix[0]=0; ix[1]=1; ix[2]=4; ix[3]=3;} // 2: cuadrilatero opuesto al 0
    else if (i==3){ix[0]=0; ix[1]=2; ix[2]=1;}          // 3: triangulo inferior
    else if (i==4){ix[0]=3; ix[1]=4; ix[2]=5;}          // 4: triangulo superior
    if (i<3) nvc=4; else nvc=3;
    return;
  }
  else if (Tipo==e_poliedro&&pd->nc){
    i=3*ciclo(i,pd->nc); nvc=3;
    ix[0]=pd->c[i]; ix[1]=pd->c[i+1]; ix[2]=pd->c[i+2];
    return;
  }
  else nvc=0;
}

void poly_dat::swap(){ // no cambia los nros de nodo (n)
  int i,j,nci;
  for (i=0;i<nc;i++) Swap(c[3*i+1],c[3*i+2]);
  for (i=0;i<nv;i++) { // reversa
    nci=cn[i][0];
    for (j=1;j<1+(nci>>1);j++) Swap(cn[i][j],cn[i][nci-j]);
  }
}

int elemento::opuesto(int i, int k){
  // para todos los 2d va del nodo al siguiente
  // orientaciones hacia fuera del elemento o izquierda en 2d
  if (Tipo==e_segmento)  // los nodos
    return (k)? -1 : ciclo(i,2);
  else if (dim()==2)
    return (k<nv()-2)? ciclo(i+k+1,nv()) : -1;
  else if (Tipo==e_tetraedro)  // cara i opuesta al vertice i
    return (k)? -1 : i;
  else if (Tipo==e_cubo) {
         if (i==0){
           if (k==0) return 0;
      else if (k==1) return 1;
      else if (k==2) return 2;
      else return -1;
    }
    else if (i==1){
           if (k==0) return 0;
      else if (k==1) return 5;
      else if (k==2) return 1;
      else return -1;
    }
    else if (i==2){
           if (k==0) return 0;
      else if (k==1) return 4;
      else if (k==2) return 5;
      else return -1;
    }
    else if (i==3){
           if (k==0) return 0;
      else if (k==1) return 2;
      else if (k==2) return 4;
      else return -1;
    }
    else if (i==4){
           if (k==0) return 1;
      else if (k==1) return 3;
      else if (k==2) return 2;
      else return -1;
    }
    else if (i==5){
           if (k==0) return 1;
      else if (k==1) return 5;
      else if (k==2) return 3;
      else return -1;
    }
    else if (i==6){
           if (k==0) return 3;
      else if (k==1) return 5;
      else if (k==2) return 4;
      else return -1;
    }
    else if (i==7){
           if (k==0) return 2;
      else if (k==1) return 3;
      else if (k==2) return 4;
      else return -1;
    }
    else return -1;
  }
  else if (Tipo==e_wedge) {
         if (i==0){
           if (k==0) return 3;
      else if (k==1) return 2;
      else if (k==2) return 1;
      else return -1;
    }
    else if (i==1){
           if (k==0) return 3;
      else if (k==1) return 0;
      else if (k==2) return 2;
      else return -1;
    }
    else if (i==2){
           if (k==0) return 3;
      else if (k==1) return 1;
      else if (k==2) return 0;
      else return -1;
    }
    else if (i==3){
           if (k==0) return 4;
      else if (k==1) return 1;
      else if (k==2) return 2;
      else return -1;
    }
    else if (i==4){
           if (k==0) return 4;
      else if (k==1) return 2;
      else if (k==2) return 0;
      else return -1;
    }
    else if (i==5){
           if (k==0) return 4;
      else if (k==1) return 0;
      else if (k==2) return 1;
      else return -1;
    }
    else return -1;
  }
  else if (Tipo==e_poliedro){
    if (!pd->cn) pd->mk_cn();
    int NV=nv(),j,ixj,l,nc,ic;
    for (j=0; j<nv()-1; j++){
      ixj=ciclo(i+j,NV); nc=pd->cn[ixj][0];
      // solo caras en que ixj sea el menor
      for (l=0;l<nc;l++){
        ic=pd->cn[ixj][l+1];
        if(pd->c[3*ic]>=ixj&&pd->c[3*ic+1]>=ixj&&pd->c[3*ic+2]>=ixj) k--;
        if (k==-1) return ic;
      }
    }
    return -1;
  }
  else return -1;
}

// simplices
int elemento::nt() const {
  return (pd)? pd->nt : tns[Tipo];
}


const int* elemento::tetra(int i) const{
  static int t[4];
  if (Tipo==e_triangulo||Tipo==e_tetraedro) return n;
  if (Tipo==e_wedge){
    if (i==0) return n;
    if (i==1) {t[0]=n[1];t[1]=n[4];t[2]=n[5];t[3]=n[3];return t;}
    if (i==2) {t[0]=n[2];t[1]=n[1];t[2]=n[5];t[3]=n[3];return t;}
    return tetra(ciclo(i,3));
  }
  if (Tipo==e_cubo){ 
    // division en 6, z=0,2,4,6
    switch (i) {
      case 0: {t[0]=n[0];t[1]=n[2];t[2]=n[4];t[3]=n[1];return t;}
      case 1: {t[0]=n[2];t[1]=n[6];t[2]=n[4];t[3]=n[1];return t;}
      case 2: {t[0]=n[1];t[1]=n[6];t[2]=n[4];t[3]=n[5];return t;}
    // simetrico
//      case 3: {t[0]=n[0];t[1]=n[4];t[2]=n[2];t[3]=n[3];return t;}
//      case 4: {t[0]=n[4];t[1]=n[6];t[2]=n[2];t[3]=n[3];return t;}
//      case 5: {t[0]=n[3];t[1]=n[4];t[2]=n[6];t[3]=n[7];return t;}
    // antisimetrico
      case 3: {t[0]=n[0];t[1]=n[4];t[2]=n[2];t[3]=n[7];return t;}
      case 4: {t[0]=n[4];t[1]=n[6];t[2]=n[2];t[3]=n[7];return t;}
      case 5: {t[0]=n[7];t[1]=n[2];t[2]=n[0];t[3]=n[3];return t;}
      default: {return tetra(ciclo(i,6));}
    }
  }
  if (Tipo==e_cuadrilatero){
    if (i==0) return n;
    if (i==1) {t[0]=n[0];t[1]=n[2];t[2]=n[3];return t;}
    return tetra(ciclo(i,2));
  }
  if (pd) {
    int* nt=&(pd->t[4*ciclo(i,pd->nt)]);
    {t[0]=n[nt[0]];t[1]=n[nt[1]];t[2]=n[nt[2]];t[3]=n[nt[3]];return t;}
  }
  return 0;
}


//========================================
//io
static int _flag_mask=0;
static int iobase=0;

void elemento::flag_mask(int i){_flag_mask=i;}

int elemento::io_base(int i){
  if (i==0||i==1) iobase=i;
  return iobase;
}

istream& elemento::lee(istream &a){
  int i=0,NV,NC,NT;
  if (Tipo==e_indefinido) return a;
  else if (Tipo==e_poligono){
    if (pd) delete pd; pd=new poly_dat; if (n) delete [] n;
    a >> NV; pd->nv=pd->nc=NV; n=new int[NV]; while (i<NV) a >> n[i++];
    pd->nt=NT=NV-2; NT*=3; pd->t=new int[NT];
    i=0; while (i<NT) a >> pd->t[i++];
  }
  else if (Tipo==e_poliedro){
    if (pd) delete pd; pd=new poly_dat; if (n) delete [] n;
    a >> NV; pd->nv=NV; n=new int[NV]; while (i<NV) a >> n[i++];
    a >> NC; pd->nc=NC; NC*=3; pd->c=new int[NC];
    i=0; while (i<NC) a >> pd->c[i++];
    a >> NT; pd->nt=NT; NT*=4; pd->t=new int[NT];
    i=0; while (i<NT) a >> pd->t[i++];
  }
  else
    {NV=tnv[Tipo]; while (i<NV) a >> n[i++];}

  if (_flag_mask) {int flag=0; a >> flag; f=flag&_flag_mask;}
  if (iobase) {i=0; while (i<NV) n[i++]--;}
  return a;
}

void elemento::lee(const int *linea){
  int NV,NC,NT;
  if (Tipo==e_indefinido) return;
  else if (Tipo==e_poligono){
    if (pd) delete pd; pd=new poly_dat; if (n) delete [] n;
    NV=linea[0]; pd->nv=pd->nc=NV; n=new int[NV]; memcpy(n,&linea[1],NV*SZI);
    pd->nt=NT=NV-2; NT*=3; pd->t=new int[NT]; memcpy(pd->t,&linea[1+NV],NT*SZI);
  }
  else if (Tipo==e_poliedro){
    if (pd) delete pd; pd=new poly_dat; if (n) delete [] n;
    NV=linea[0]; pd->nv=NV; n=new int[NV]; memcpy(n,&linea[1],NV*SZI);
    NC=linea[1+NV]; pd->nc=NC; NC*=3; pd->c=new int[NC]; memcpy(pd->c,&linea[2+NV],NC*SZI);
    NT=linea[2+NV+NC]; pd->nt=NT; NT*=4; pd->t=new int[NT]; memcpy(pd->t,&linea[3+NV+NC],NT*SZI);
  }
  else
    {NV=tnv[Tipo]; memcpy(n,linea,NV*SZI);}
  if (iobase) {int i=0; while (i<NV) n[i++]--;}
}

ostream& elemento::graba(ostream &a) const{
  int i=0,NV=nv(),NC,NT;
  if (Tipo==e_indefinido) return a;
  else if (Tipo==e_poligono){
    a << NV; while (i<NV) a << "\t" << n[i++]+iobase;
    NT=NV-2; i=0; NT*=3; a << "\t"; while (i<NT) a << " " << pd->t[i++];
  }
  else if (Tipo==e_poliedro){
    a << NV; while (i<NV) a << "\t" << n[i++]+iobase;
    NC=nc(); a << "\t" << NC; i=0; NC*=3; while (i<NC) a << " " << pd->c[i++];
    NT=nt(); a << "\t" << NT; i=0; NT*=4; while (i<NT) a << " " << pd->t[i++];
  }
  else
    {a << n[0]+iobase; i=1; while (i<NV) a << "\t" << n[i++]+iobase;}

  if (_flag_mask) a << "\t" << (f&_flag_mask);
  return a;
}


//========================================================
// solo poligonos y poliedros

void poly_dat::rm_cn() {
  if (cn) {for (int i=0;i<nv;i++) delete [] cn[i]; delete[] cn; cn=0;}
}

poly_dat::~poly_dat(){
  rm_cn();
  if (c) delete[] c;
  delete[] t;
}

void poly_dat::ini(){
  if (c) delete[]c; c=0;
  delete[]t; t=0;
  rm_cn();
  nv=0; nc=0; nt=0;
}


// convierte en tetra un poliedro de cuatro nodos
elemento& elemento::p2t(){
  if (nv()!=4||Tipo!=e_poliedro) return *this;
  // orientacion de la cara 0 (invertida!)
  int *newn=new int[4];
  newn[0]=n[pd->c[0]];newn[1]=n[pd->c[2]];newn[2]=n[pd->c[1]];
  int i=4; while(i--)
    {if (n[i]!=newn[0]&&n[i]!=newn[1]&&n[i]!=newn[2]) break;}
  newn[3]=n[i];
  delete pd; delete []n; n=newn;
  Tipo=e_tetraedro;
  return *this;
}

//convierte un cubo en un poliedro
elemento& elemento::c2p(){
  if (Tipo!=e_cubo) return *this;
  pd=new poly_dat;
  pd->nv=8;pd->nc=12;pd->nt=6;
  // division en 5, z=0,2,4,6, comentada la simetrica
  static const int 
    c[]={
      0,2,1, 0,3,2,
      0,1,4, 5,4,1,
      /*0,4,3, 7,3,4,*/ 0,4,7, 3,0,7,
      6,4,5, 6,7,4,
      /*6,2,3, 7,6,3,*/ 6,2,7, 3,7,2,
      6,1,2, 5,1,6
    },
    t[]={
      0,2,4,1,
      2,6,4,1,
      1,6,4,5,
//      0,4,2,3, // simetrica
//      4,6,2,3,
//      3,4,6,7,
      0,4,2,7, // antisimetrica
      4,6,2,7,
      7,2,0,3
    };
  pd->c=new int[36]; memcpy(pd->c,c,36*SZI);
  pd->t=new int[24]; memcpy(pd->t,t,24*SZI);
  Tipo=e_poliedro;
  return *this;
}

// caras de cada vertice (solo para funciones de forma en poliedros)
void poly_dat::mk_cn(){
  int i,j;
  // cuenta caras/nodo
  int *cpn=(int*)calloc(nv,SZI);
  for (i=0;i<3*nc;i++) cpn[c[i]]++;
  // construye cn
  cn=new int* [nv]; for (i=0;i<nv;i++) {cn[i]=new int[cpn[i]+1]; cn[i][0]=cpn[i];}
  free(cpn);
  // arma vecinos de cara (horrible, pero como sea son 6 comparaciones)
  int *vecino=new int [3*nc];
  for (i=0;i<3*nc-3;i+=3){
    for (j=i+3;j<3*nc;j+=3){
      if (c[i]==c[j  ]){ // 0=0
        if (c[i+1]==c[j+2]) {vecino[i  ]=j/3; vecino[j+2]=i/3;} //1=2
        if (c[i+2]==c[j+1]) {vecino[i+2]=j/3; vecino[j  ]=i/3;} //2=1
        continue;
      }
      if (c[i]==c[j+1]){ // 0=1
        if (c[i+1]==c[j  ]) {vecino[i  ]=j/3; vecino[j  ]=i/3;} //1=0
        if (c[i+2]==c[j+2]) {vecino[i+2]=j/3; vecino[j+1]=i/3;} //2=2
        continue;
      }
      if (c[i]==c[j+2]){ // 0=2
        if (c[i+1]==c[j+1]) {vecino[i  ]=j/3; vecino[j+1]=i/3;} //1=1
        if (c[i+2]==c[j  ]) {vecino[i+2]=j/3; vecino[j+2]=i/3;} //2=0
        continue;
      }
      if (c[i+1]==c[j  ]){ // 1=0
        if (c[i+2]==c[j+2]) {vecino[i+1]=j/3; vecino[j+2]=i/3;} //2=2
        continue;
      }
      if (c[i+1]==c[j+1]){ // 1=1
        if (c[i+2]==c[j  ]) {vecino[i+1]=j/3; vecino[j  ]=i/3;} //2=0
        continue;
      }
      if (c[i+1]==c[j+2]){ // 1=2
        if (c[i+2]==c[j+1]) {vecino[i+1]=j/3; vecino[j+1]=i/3;} //2=1
        continue;
      }
    }
  }
  // arma cn
  // asigna a cada nodo una cara que lo contenga (para usar como primera cara)
  for (i=0;i<3*nc;i++) cn[c[i]][1]=i/3;
  // arma el resto por vecindad
  int qc,ic,iv; 
  for (i=0;i<nv;i++){
    qc=cn[i][0];
    for (j=1;j<qc;j++){
      ic=cn[i][j];
      if (i==c[3*ic]) iv=2; else if (i==c[3*ic+1]) iv=0; else iv=1; // el anterior
      cn[i][j+1]=vecino[3*cn[i][j]+iv];
    }
  }
  delete [] vecino;
}
