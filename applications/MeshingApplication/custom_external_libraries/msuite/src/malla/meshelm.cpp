// elementos de malla
#include <cstdlib> // qsort
//#define temporizar
#include "tiempo.h"
#include "malla.h"
#include "esfera.h" // sliverlim

using namespace std;

static const int t4[]={0,1,2,3}; // para pd trucho de tetraedros

//caracteristicas del elemento

//    _      ^      _
//   |\    ni|      /|
//      \nd  |    /
//        \  |  /
//          \|/
//   <-------+ a (hacia dentro de la pantalla)

// diedro dadas las normales hacia dentro y la arista (todos unitarios!!)
// izquierda y derecha desde el interior en la dir de arista
static double diedro_(const punto &ni, const punto &nd, const punto &a){
  double
    s=triple(nd,ni,a), // sin(180-d)
    c=-ni*nd,          // cos(180-d)
    ang=r2g(atan2(s,c)); // -180 a 180
  if (ang<0) ang+=360;
  return ang;
}

// angulo 2d o diedro 3d max y min (no esta para cubos)
// OJO LAS ORIENTACIONES DE LAS CARAS
bool malla::angulos(int ie,double &amin,double &amax){
  bool retval=true;
  amin=MAXREAL; amax=-MAXREAL;
  elemento &ei=e[ie];// no cacheo por cn
  int i,nv=ei.nv(); punto *ni=new punto[nv]; for (i=0;i<nv;i++) ni[i]=n[ei[i]];
  e_tipo tipo=ei.tipo();
  if (tipo==e_tetraedro){
    // aristas de menor a mayor y normales hacia dentro
    punto
      a[6]={(ni[3]-ni[0]).dir(),(ni[3]-ni[1]).dir(),(ni[3]-ni[2]).dir(),
            (ni[1]-ni[0]).dir(),(ni[2]-ni[1]).dir(),(ni[0]-ni[2]).dir()},
      x[4]={(a[1]%a[4]).dir(),(a[2]%a[5]).dir(),(a[0]%a[3]).dir(),(a[3]%a[4]).dir()};
    set_min_max(amin,amax,diedro_(x[2],x[1],a[0]));
    set_min_max(amin,amax,diedro_(x[0],x[2],a[1]));
    set_min_max(amin,amax,diedro_(x[1],x[0],a[2]));
    set_min_max(amin,amax,diedro_(x[3],x[2],a[3]));
    set_min_max(amin,amax,diedro_(x[3],x[0],a[4]));
    set_min_max(amin,amax,diedro_(x[3],x[1],a[5]));
  }
  else if (tipo==e_poliedro){
    if (!ei.pd->cn) ei.pd->mk_cn();
    int j;
    // caras y caras del nodo i
    int **cn=ei.pd->cn,ncni,*cni,nvc,nca[3],ncp[3],
      ino,ina,ica,inp,icp,ix;
    punto a,xa,xp;
    for (i=0;i<nv;i++){
      ncni=cn[i][0]; cni=&(cn[i][1]);
      for (j=0;j<ncni;j++){
        // cara anterior (derecha desde dentro)
        ica=cni[j]; ei.cara(ica,nvc,nca);
        if (i==nca[0]) ix=0; else if (i==nca[1]) ix=1; else ix=2;
        ino=nca[(ix+2)%3]; if (ino<i) continue; // de menor a mayor
        a=(ni[i]-ni[ino]).dir();
        ina=nca[(ix+1)%3]; xa=((ni[ina]-ni[ino])%a).dir(); // hacia dentro
        // cara posterior (izquierda)
        icp=cni[(j+1)%ncni]; ei.cara(icp,nvc,ncp);
        if (i==ncp[0]) ix=0; else if (i==ncp[1]) ix=1; else ix=2;
        inp=ncp[(ix+2)%3]; xp=(a%(ni[inp]-ni[ino])).dir();
        set_min_max(amin,amax,diedro_(xp,xa,a));
      }
    }
  }
  else if (tipo==e_triangulo){
    punto a[3]={(ni[1]-ni[0]).dir(),(ni[2]-ni[1]).dir(),(ni[0]-ni[2]).dir()};
    double x[3];
    x[0]=r2g(acos(-a[0]*a[1])); set_min_max(amin,amax,x[0]);
    x[1]=r2g(acos(-a[1]*a[2])); set_min_max(amin,amax,x[1]);
    x[2]=180-x[0]-x[1]; set_min_max(amin,amax,x[2]);
  }
  else if (tipo==e_cuadrilatero){
    punto a[4]={(ni[1]-ni[0]).dir(),(ni[2]-ni[1]).dir(),(ni[3]-ni[2]).dir(),(ni[0]-ni[3]).dir()};
    double x[4]={
      r2g(acos(-a[0]*a[1])),r2g(acos(-a[1]*a[2])),
      r2g(acos(-a[2]*a[3])),r2g(acos(-a[3]*a[0]))};
    set_min_max(amin,amax,x[0]); set_min_max(amin,amax,x[1]);
    set_min_max(amin,amax,x[2]); set_min_max(amin,amax,x[3]);
  }
  else if (tipo==e_poligono){
    punto ap=(ni[0]-ni[nv-1]).dir(),aa;
    double x;
    for (i=0;i<nv;i++){
      aa=ap;
      ap=(ni[(i+1)%nv]-ni[i]).dir();
      x=r2g(acos(-aa*ap)); set_min_max(amin,amax,x);
    }
  }
  else if (tipo==e_wedge){
    // normales hacia dentro y aristas
    punto
      a[9]={(ni[3]-ni[0]).dir(),(ni[4]-ni[1]).dir(),(ni[5]-ni[2]).dir(),
            (ni[1]-ni[0]).dir(),(ni[2]-ni[1]).dir(),(ni[0]-ni[2]).dir(),
            (ni[4]-ni[3]).dir(),(ni[5]-ni[4]).dir(),(ni[3]-ni[5]).dir()},
      x[6]={((ni[4]-ni[2])%(ni[5]-ni[1])).dir(),
            ((ni[5]-ni[0])%(ni[3]-ni[2])).dir(),
            ((ni[3]-ni[1])%(ni[4]-ni[0])).dir(),
            (a[5]%a[3]).dir(),(a[6]%a[8]).dir()};
    set_min_max(amin,amax,diedro_(x[2],x[1],a[0]));
    set_min_max(amin,amax,diedro_(x[0],x[2],a[1]));
    set_min_max(amin,amax,diedro_(x[1],x[0],a[2]));
    set_min_max(amin,amax,diedro_(x[3],x[2],a[3]));
    set_min_max(amin,amax,diedro_(x[3],x[0],a[4]));
    set_min_max(amin,amax,diedro_(x[3],x[1],a[5]));
    set_min_max(amin,amax,diedro_(x[2],x[4],a[6]));
    set_min_max(amin,amax,diedro_(x[0],x[4],a[7]));
    set_min_max(amin,amax,diedro_(x[1],x[4],a[8]));
  }
  else retval=false;
  delete [] ni;
  return retval;
}

// maximo diedro en 3d (no esta para wedges ni cubos)
double malla::maxd(int ie){double amax; int ic1,ic2; return maxd(ie,amax,ic1,ic2);}
double malla::maxd(int ie,double &amax, int &ic1, int &ic2){
  amax=-MAXREAL;
  elemento &ei=e[ie];// no cacheo por cn
  int i,nv=ei.nv(); punto *ni=new punto[nv]; for (i=0;i<nv;i++) ni[i]=n[ei[i]];
  e_tipo tipo=ei.tipo();
  if (tipo==e_tetraedro){
    // aristas de menor a mayor y normales hacia dentro
    punto
      a[6]={(ni[3]-ni[0]).dir(),(ni[3]-ni[1]).dir(),(ni[3]-ni[2]).dir(),
            (ni[1]-ni[0]).dir(),(ni[2]-ni[1]).dir(),(ni[2]-ni[0]).dir()},
      x[4]={(a[1]%a[4]).dir(),(a[5]%a[0]).dir(),(a[0]%a[3]).dir(),(a[3]%a[5]).dir()};
    if (set_max(amax,diedro_(x[2],x[1],a[0]))) {ic1=1; ic2=2;}
    if (set_max(amax,diedro_(x[0],x[2],a[1]))) {ic1=2; ic2=0;}
    if (set_max(amax,diedro_(x[1],x[0],a[2]))) {ic1=0; ic2=1;}
    if (set_max(amax,diedro_(x[3],x[2],a[3]))) {ic1=2; ic2=3;}
    if (set_max(amax,diedro_(x[3],x[0],a[4]))) {ic1=0; ic2=3;}
    if (set_max(amax,diedro_(x[1],x[3],a[5]))) {ic1=3; ic2=1;}
  }
  else if (tipo==e_poliedro){
    if (!ei.pd->cn) ei.pd->mk_cn();
    int j;
    // caras y caras del nodo i
    int **cn=ei.pd->cn,ncni,*cni,nvc,nca[3],ncp[3],
      ino,ina,ica,inp,icp,ix;
    punto a,xa,xp;
    for (i=0;i<nv;i++){
      ncni=cn[i][0]; cni=&(cn[i][1]);
      for (j=0;j<ncni;j++){
        // cara anterior (derecha desde dentro)
        ica=cni[j]; ei.cara(ica,nvc,nca);
        if (i==nca[0]) ix=0; else if (i==nca[1]) ix=1; else ix=2;
        ino=nca[(ix+2)%3]; if (ino<i) continue; // de menor a mayor
        a=(ni[i]-ni[ino]).dir();
        ina=nca[(ix+1)%3]; xa=((ni[ina]-ni[ino])%a).dir(); // hacia dentro
        // cara posterior (izquierda)
        icp=cni[(j+1)%ncni]; ei.cara(icp,nvc,ncp);
        if (i==ncp[0]) ix=0; else if (i==ncp[1]) ix=1; else ix=2;
        inp=ncp[(ix+2)%3]; xp=(a%(ni[inp]-ni[ino])).dir(); // hacia dentro
        if (set_max(amax,diedro_(xp,xa,a))) {ic1=ica; ic2=icp;}
      }
    }
  }
  delete [] ni;
  return amax;
}

// modo de angulos: 0=grados, 1=rad, 2=sen, 3=cos

// diedro 3d entre caras c1 y c2
//
//       j+1   i-1
//           |
//           |
//           |
// j-1 \  c2 | c1  / i+1
//      \    |    /
//        \  |  /
//          \|/
//         j   i
//
double malla::angulo(int ie,int ic1, int ic2, int modo) const{
  elemento &ei=e[ie];
  if (ei.dim()!=3) return -MINREAL;
  int i,nc1[4],nvc1,j,nc2[4],nvc2;
  ei.cara(ic1,nvc1,nc1); ei.cara(ic2,nvc2,nc2);
  for (i=0;i<nvc1;i++){
    for (j=0;j<nvc2;j++) {
      if ((nc1[i]!=nc2[j])||(nc1[(i+nvc1-1)%nvc1]!=nc2[(j+1)%nvc2])) continue;
      punto
        n0=n[ei[nc2[j]]],
        a=(n[ei[nc2[(j+1)%nvc2]]]-n0).dir(),
        xc1=((n[ei[nc1[(i+1)%nvc1]]]-n0)%a).dir(),
        xc2=(a%(n[ei[nc2[(j+nvc2-1)%nvc2]]]-n0)).dir();
      if (modo<=1){// angulo
        double
          s=triple(xc1,xc2,a), // sin
          c=-xc1*xc2,            // cos  (cos(180-d)=-cos(d))
          ang=atan2(s,c);      // -180 a 180
        if (ang<0) ang+=DOSPI;
        if (modo==0) return ang; // en radianes
        return r2g(ang); // en grados
      }
      else if (modo==2) return triple(xc1,xc2,a);  // sen
      else if (modo==3) return -(xc1*xc2);  // cos
    }
  }
  return -MINREAL;
}

// jacobiano de un elemento en un nodo triedrico
// bah.... modulo del triple producto en 3d 
// (solo 3d, en 2d hay que ver que se hace con la orientacion)
double malla::jaco(int ie, int in) const{
  const elemento &ei=e[ie];
  int nv=ei.nv(),dim=ei.dim(); if (dim!=3) return MAXREAL;
  int ix=ei.index(in); if (ix==nv) return MAXREAL;
  e_tipo tipo=ei.tipo();
  int a[3]; punto ni=n[in];
  if (tipo==e_tetraedro) {
    if (ix<3){a[0]=(ix+1)%3;   a[1]=(ix+2)%3;   a[2]=3;   } // abajo
    else     {a[0]=2;          a[1]=1;          a[2]=0;   } // arriba
  }
  else if (tipo==e_wedge) {
    if (ix<3){a[0]=(ix+1)%3;   a[1]=(ix+2)%3;   a[2]=ix+3;} // abajo
    else     {a[0]=(ix+2)%3+3; a[1]=(ix+1)%3+3; a[2]=ix-3;} // arriba
  }
  else if (tipo==e_cubo) {
    if (ix<4){a[0]=(ix+1)%4;   a[1]=(ix+3)%4;   a[2]=ix+4;} // abajo
    else     {a[0]=(ix+3)%4+4; a[1]=(ix+1)%4+4; a[2]=ix-4;} // arriba
  }
  else return MAXREAL; // si hace falta, implementar para poliedros
  return triple(n[ei[a[0]]]-ni,n[ei[a[1]]]-ni,n[ei[a[2]]]-ni);
}

// angulo de un elemento 2d en el nodo i
double malla::angulo(int ie,int i, int modo) const {
  elemento &ei=e[ie];
  if (ei.dim()!=2) return -MINREAL;
  int nv=ei.nv();
  punto
    ni=n[ei[i]],
    aa=(n[ei[(i+nv-1)%nv]]-ni).dir(),
    ap=(n[ei[(i+1)%nv]]-ni).dir();
  if (modo<=1){// angulo
    double
      s=(ap%aa).mod(), // sin
      c=ap*aa, // cos
      ang=atan2(s,c);    // -180 a 180
    if (ang<0) ang+=DOSPI; // <0 -> >180
    if (modo==0) return ang; // en radianes
    return r2g(ang); // en grados
  }
  else if (modo==2) return ap.pv2d(aa);  // sen
  else if (modo==3) return ap*aa;  // cos
  return -MINREAL;
}

// convierte el quad i en triangulos con una diagonal por el vertice k
// (muy parecido a diagonal_swap en smooth.cpp)
bool malla::q2t(int i, int k){
  if (e[i].tipo()!=e_cuadrilatero) return false;
  elemento &ei=e[i];
  int ie, ine[4], vec[4];
  elemento e3(e_triangulo);
  memcpy(ine,ei.n,4*SZI);
  e3[0]=ine[(k+2)%4]; e3[1]=ine[(k+3)%4]; e3[2]=ine[k]; ie=e+=e3;
  e3[0]=ine[k]; e3[1]=ine[(k+1)%4]; e3[2]=ine[(k+2)%4]; ei=e3;
  if (vecino){
    memcpy(vec,vecino[i].vertex,4*SZI);
    cpline vec3(3); vec3.len=3;
    vec3[0]=vec[(k+2)%4]; vec3[1]=vec[(k+3)%4]; vec3[2]=i; vecino+=vec3;
    if (vec3[0]>=0) vecino[vec3[0]].replace1(i,ie);
    if (vec3[1]>=0) vecino[vec3[1]].replace1(i,ie);
    vec3[0]=vec[k]; vec3[1]=vec[(k+1)%4]; vec3[2]=ie; vecino[i]=vec3;
  }
  if (ce||re||ve||(edir&&dir)){// esfera y dir
    if (!ce&&!re&&!(edir&&dir)){ // solo volumen
      ve+=volumen(ie);
      ve[i]=volumen(ei);
    }
    else{
      double r,v; punto c,pdir;
      if (edir&&dir){ // esfera y normal
        esfera_e(ie,c,r,v,&pdir);
        re+=r; ce+=c; ve+=v; dir+=pdir;
        esfera_e(i,ce[i],re[i],ve[i],&dir[i]);
      }
      else { //solo esfera
        esfera_e(ie,c,r,v);
        re+=r; ce+=c; ve+=v;
        esfera_e(i,ce[i],re[i],ve[i]);
      }
    }
  }
  // elementos de nodo
  n[ine[(k+3)%4]].e.replace1(i,ie);
  n[ine[k]].e+=ie; n[ine[(k+2)%4]].e+=ie;
  // vecinos naturales
  if (nn){
    nn[ine[(k+1)%4]].remove1(ine[(k+3)%4]);
    nn[ine[(k+3)%4]].remove1(ine[(k+1)%4]);
  }
  tipo.set(m_modificada);
  return true;
}

// convierte cuad en tri por delaunay
bool malla::q2t_d(int i){
  if (e[i].tipo()!=e_cuadrilatero) return false;
  if (tipo.noes(m_orientada)) orienta(); // orientada
  punto c,a; double r;
  c3(n[e[i][3]],n[e[i][0]],n[e[i][1]],c,r,a); int idiag=0; double a2=a.mod2();
  if (a2<ERRADM){
    c3(n[e[i][0]],n[e[i][1]],n[e[i][2]],c,r,a); idiag=1; a2=a.mod2();
    if (a2<ERRADM) return q2t(i,idiag);
  }
  if (tipo.es(m_planaxy)) { // plana
    if (n[e[i][idiag+2]].distancia(c)>r) return q2t(i,idiag); return q2t(i,idiag+1);
  }
  // 3D
  punto d=n[e[i][idiag+2]]-c;
  d-=a*((d*a)/a.mod2()); //proyectado en el plano del 1er triangulo
  if (d.mod()>r) return q2t(i,idiag); return q2t(i,idiag+1);
}

// convierte cuad en tri por el ang mas obtuso
bool malla::q2t_a(int i){
  if (e[i].tipo()!=e_cuadrilatero) return false;
  orienta(); // orientada
  punto ai[4],ntmp; int j,jmax,ine[4];
  double a,amax=0;
  elemento &ei=e[i];
  // versores arista
  memcpy(ine,ei.n,4*SZI); // cache
  ntmp=ai[3]=n[ine[0]];
  ai[3]-=(ai[2]=n[ine[3]]); ai[3].dir(); // 3->0
  ai[2]-=(ai[1]=n[ine[2]]); ai[2].dir(); // 2->3
  ai[1]-=(ai[0]=n[ine[1]]); ai[1].dir(); // 1->2
  ai[0]-=ntmp;              ai[0].dir(); // 0->1
  // busca el mayor angulo
  // uso angulo y no coseno por si hay invertidos
  jmax=-1;
  for (j=0;j<4;j++){
    // la segunda arista debe ir opuesta ==> tangente no cambia
    a=atan2((ai[j]%ai[(j+3)%4]).mod(),ai[j]*ai[(j+3)%4]); // -180 a 180
    if (a<0) a+=DOSPI; // <0 -> >180
    if (set_max(amax,a)) jmax=j;
  }
  return q2t(i,jmax);
}

//convierte quad en tri generando la mayor minima altura e triangulos
bool malla::q2t_h(int i){
  if (e[i].tipo()!=e_cuadrilatero) return false;
  if (tipo.noes(m_orientada)) orienta(); // orientada

  elemento &ei=e[i];
  punto n0=n[ei[0]],n1=n[ei[1]],n2=n[ei[2]],n3=n[ei[3]];
  
  double 
    d01=(n0-n1).mod2(),
    d02=(n0-n2).mod2(),
    d03=(n0-n3).mod2(),
    d13=(n1-n3).mod2(),
    d12=(n1-n2).mod2(),
    d23=(n2-n3).mod2(),
  
  // areas antes y despues del swap
    aa1=((n0-n2)%(n1-n2)).mod2(),
    aa2=((n0-n3)%(n1-n3)).mod2(),
    ad1=((n2-n0)%(n3-n0)).mod2(),
    ad2=((n2-n1)%(n3-n1)).mod2(),
  
    mha1 = aa1/Max3(d01,d02,d12),
    mha2 = aa2/Max3(d03,d01,d13),
    mhd1 = ad1/Max3(d02,d03,d23),
    mhd2 = ad2/Max3(d12,d13,d23);
  
  if (Min(mha1,mhd1)>Min(mha2,mhd2)) return q2t(i,0); return q2t(i,1);
}

bool malla::q2t(int metodo){
  int i,lenori=e.len;
  bool cambio_alguno=false;
       if (metodo==1) {for (i=0;i<lenori;i++) if (q2t_d(i)) cambio_alguno=true;}
  else if (metodo==2) {for (i=0;i<lenori;i++) if (q2t_a(i)) cambio_alguno=true;}
  else if (metodo==3) {for (i=0;i<lenori;i++) if (q2t_h(i)) cambio_alguno=true;}
  return (cambio_alguno);
}

// hexaedros en 5 tetraedros
// trabaja por nodos, de modo que lo hace desordenado y solo para mallas conectadas por cara o al menos arista
// ordena para que el tetraedro n venga del hexa int(n/5)
bool malla::h5t(){
  int in,ie,i,ix;
  for (i=0;i<e.len;i++) if (e[i].tipo()!=e_cubo) return false;
  if (!mk_vecino()) return false;

  // marco un bit a los nodos con flag1
  n[0].f.set(flag1);
  ordlist procesar(n.len); procesar+=0;
  elemento t(e_tetraedro);
  array1<elemento> et(5*e.len); et.len=5*e.len;
  while (procesar){
    in=procesar.last(); procesar.len--;
    if (n[in].f.es(flag2)) continue; // ya fue procesado
    n[in].f.set(flag2);
    const cpline& en=n[in].e;
    for (i=0;i<en.len;i++){
      ie=en[i]; elemento &ei=e[ie];
      if (ei.f.es(flag2)) continue;  // ya fue procesado
      t.f=ei.f&fmask; ei.f.set(flag2);
      ix=ei.index(in);
      switch (ix){
        case 0: case 2: case 5: case 7:
          if (n[ei[0]].f.noes(flag1)) {n[ei[0]].f.set(flag1); procesar+=ei[0];}
          if (n[ei[2]].f.noes(flag1)) {n[ei[2]].f.set(flag1); procesar+=ei[2];}
          if (n[ei[5]].f.noes(flag1)) {n[ei[5]].f.set(flag1); procesar+=ei[5];}
          if (n[ei[7]].f.noes(flag1)) {n[ei[7]].f.set(flag1); procesar+=ei[7];}
          t[0]=ei[0]; t[1]=ei[4]; t[2]=ei[1]; t[3]=ei[3]; et[5*ie  ]=t;
          t[0]=ei[2]; t[1]=ei[6]; t[2]=ei[3]; t[3]=ei[1]; et[5*ie+1]=t;
          t[0]=ei[5]; t[1]=ei[1]; t[2]=ei[4]; t[3]=ei[6]; et[5*ie+2]=t;
          t[0]=ei[7]; t[1]=ei[4]; t[2]=ei[3]; t[3]=ei[6]; et[5*ie+3]=t;
          t[0]=ei[1]; t[1]=ei[3]; t[2]=ei[4]; t[3]=ei[6]; et[5*ie+4]=t;
          break;
        default:
          if (n[ei[1]].f.noes(flag1)) {n[ei[1]].f.set(flag1); procesar+=ei[1];}
          if (n[ei[3]].f.noes(flag1)) {n[ei[3]].f.set(flag1); procesar+=ei[3];}
          if (n[ei[4]].f.noes(flag1)) {n[ei[4]].f.set(flag1); procesar+=ei[4];}
          if (n[ei[6]].f.noes(flag1)) {n[ei[6]].f.set(flag1); procesar+=ei[6];}
          t[0]=ei[1]; t[1]=ei[2]; t[2]=ei[0]; t[3]=ei[5]; et[5*ie  ]=t;
          t[0]=ei[3]; t[1]=ei[0]; t[2]=ei[2]; t[3]=ei[7]; et[5*ie+1]=t;
          t[0]=ei[4]; t[1]=ei[7]; t[2]=ei[5]; t[3]=ei[0]; et[5*ie+2]=t;
          t[0]=ei[6]; t[1]=ei[5]; t[2]=ei[7]; t[3]=ei[2]; et[5*ie+3]=t;
          t[0]=ei[0]; t[1]=ei[5]; t[2]=ei[2]; t[3]=ei[7]; et[5*ie+4]=t;
      }
    }
  }
  
  rm_nn(); rm_vecino(); rm_dir(); rm_frontera(); ve.ini(); re.ini(); ce.ini();

  e.roba(et); // reemplaza los elementos

  // arregla elementos y flags de nodos  
  int f=(flag1|flag2);
  for (in=0;in<n.len;in++) {
    n[in].f.reset(f);
    n[in].e.ini();
  }
  for (i=0;i<e.len;i++){
    const elemento &ei=e[i];
    for (ix=0;ix<4;ix++) n[ei[ix]].e+=i;
  }

  return true;
}


// aspect como min/max distancia e/nodos
double malla::aspect(int ix) const{
  elemento &ei=e[ix];
  double min=MAXREAL,max=0;
  int i,j,nv=ei.nv();
  punto p; // cache
  for (i=0;i<nv-1;i++){
    p=n[ei[i]];
    for (j=i+1;j<nv;j++) set_min_max(min,max,p.distancia(n[ei[j]]));
  }
  return (max<ERRADM)? -1 : min/max;
}

//bbox de un elemento
void malla::bbox(int ie, punto p[2]) const {
  const elemento &el=e[ie];
  p[0]=p[1]=n[el[0]];
  int i,nv=el.nv();
  for (i=1;i<nv;i++) n[el[i]].set_min_max(p[0],p[1]);
}

//Centro de gravedad
punto malla::gp(const elemento &el) const {
  int nv=el.nv();
  punto g(0.,0.,0.);
  for (int i=0;i<nv;i++) g+=n[el[i]];
  g/=nv;
  return g;
}

// g de simplice pesado con el area o voumen orientados
///////////////////////////////////////////////////////////////hay algo raro
static void gvol(const punto *p, int i1, int i2, int i3, int i4, punto &g, double &vt){
  double v=triple(p[i2]-p[i1],p[i3]-p[i1],p[i4]-p[i1]); vt+=v;
  g+=(p[i1]+p[i2]+p[i3]+p[i4])*v;
}
static void garea(const punto *p, int i1, int i2, int i3,const punto &vert, punto &g, double &at){
  double a=triple(p[i2]-p[i1],p[i3]-p[i1],vert); at+=a;
  g+=(p[i1]+p[i2]+p[i3])*a;
}

punto malla::gv(const elemento &ei) const {
  punto g(0.,0.,0.),vert;
  e_tipo t=ei.tipo();
  if(t==e_segmento||t==e_triangulo||t==e_tetraedro) return gp(ei);

  int i,nv=ei.nv();
  double vt=0;
  punto *ni=new punto[nv]; for(i=0;i<nv;i++) ni[i]=n[ei[i]];
  if (t==e_cuadrilatero){ // supuesto convexo!!!!!!!!!!!!!!!!!!!!!!!!!!!
    vert=((ni[2]-ni[0])%(ni[3]-ni[1])).dir();
    ///////////////////////////////////////////////////////////////hay algo raro
    garea(ni,0,1,2,vert,g,vt);
    garea(ni,0,2,3,vert,g,vt);
    g/=(3*vt);
  }
  else if (t==e_cubo){     // supuesto convexo!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gvol(ni,0,1,3,5,g,vt);
    gvol(ni,0,5,3,4,g,vt);
    gvol(ni,3,4,5,7,g,vt);
    gvol(ni,2,3,1,5,g,vt);
    gvol(ni,2,3,5,7,g,vt);
    gvol(ni,5,2,7,6,g,vt);
    g/=(4*vt);
  }
  else if (t==e_wedge){    // supuesto convexo!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gvol(ni,0,1,2,4,g,vt);
    gvol(ni,0,4,2,3,g,vt);
    gvol(ni,2,3,4,5,g,vt);
    g/=(4*vt);
  }
  else if (t==e_poligono){
    int *t=ei.pd->t;
    vert=((ni[t[1]]-ni[t[0]])%(ni[t[2]]-ni[t[0]])).dir();
    for (i=0;i<nv-2;i++) garea(ni,t[3*i],t[3*i+1],t[3*i+2],vert,g,vt);
    g/=(3*vt);
  }
  else if (t==e_poliedro){
    int nt=ei.pd->nt,*t=ei.pd->t;
    for (i=0;i<nt;i++) gvol(ni,t[4*i],t[4*i+1],t[4*i+2],t[4*i+3],g,vt);
    g/=(4*vt);
  }
  delete [] ni;
  return g;
}

// centroide y normal hacia afuera de la cara ic del elemento ie
void malla::normal(int ie, int ic, punto &gc, punto &nc) const {
  elemento c=e[ie].cara(ic);
  e_tipo t=c.tipo();
  int i,nv=c.nv();
  double vt=0;
  punto *ni=new punto[nv]; for(i=0;i<nv;i++) ni[i]=n[c[i]];

  if (t==e_segmento){
    gc=(ni[1]+ni[0])/2;
    nc=(ni[1]-ni[0]).dir().giro90();
  }
  else if (t==e_triangulo) {
    gc=(ni[2]+ni[1]+ni[0])/3;
    nc=((ni[1]-ni[0])%(ni[2]-ni[0])).dir();
  }
  else if (t==e_cuadrilatero){
    nc=((ni[2]-ni[0])%(ni[3]-ni[1])).dir();
    garea(ni,0,1,2,nc,gc,vt);
    garea(ni,0,2,3,nc,gc,vt);
    gc/=(3*vt);
  }
  else if (t==e_poligono){
    int *t=c.pd->t;
    nc=((ni[2]-ni[0])%(ni[nv-1]-ni[1])).dir();
    for (i=0;i<nv-2;i++) garea(ni,t[3*i],t[3*i+1],t[3*i+2],nc,gc,vt);
    gc/=(3*vt);
  }
  delete [] ni;
}

// verifica si un punto es interior o limite de un elemento
bool malla::inter(const punto &p, int i){
  int nc;
  return inter(p,i,nc);
}
// en nc devuelve el nodo mas cercano
bool malla::inter(const punto &p, int ie, int &nc){
  const elemento &ei=e[ie];
  int j,k,nv=ei.nv(),dim=ei.dim();
  double f[4];

  if (nv==dim+1){ //simplice
    fforma(ie,p,f);
    double fj,fmin=f[0],fmax=fmin; nc=ei[0];
    for (j=1;j<nv;j++) {
      fj=f[j];
      if (fj<fmin) fmin=fj;
      else if (fj>fmax) {fmax=fj; nc=ei[j];}
    }
    if (fmin>0&&dim==1){ // admite solo cierta tolerancia fuera
      punto pf=n[ei[0]]*f[0]+n[ei[1]]*f[1];
      double l=(n[ei[0]]-n[ei[1]]).mod2();
      if (p.distancia2(pf)>.04*l) return false;
    }
    return (fmin>0);
  }

  elemento et((dim==3)? e_tetraedro : e_triangulo);
  int nf=ei.nc(),in0=ei[0];
  bool dentro=false;
  double dmin;
  // considero el elemento convexo y hago simplices 
  // entre el nodo 0 y cada cara triangular o cada triangulo de cara cuadrilatera
  for(j=0;j<nf;j++){ // caras sin el nodo 0
    elemento cj=ei.cara(j); // hacia fuera del elemento en 3d
    if (cj.have(in0)) continue;
    // triangulo o primer triangulo de un quad
    et[0]=in0; for (k=0;k<dim;k++) et[k+1]=cj[k]; 
    fforma(et,p,f);
    for (k=0;k<nv;k++) if (f[k]<0) break;
    if (k==nv) {dentro=true; break;}
    else if (k==0) break; // esta fuera del elemento (detras de la cara original)
    if (cj.tipo()==e_cuadrilatero){ // segundo triangulo de la cara
      et[2]=et[3]; et[3]=cj[3];
      fforma(et,p,f);
      for (k=0;k<nv;k++) if (f[k]<0) break;
      if (k==nv) {dentro=true; break;}
      else if (k==0) break; // esta fuera del elemento (detras de la cara original)
    }
  }
  // dentro o fuera, busco el nodo mas cercano
  nc=ei[0]; dmin=p.distancia2(n[nc]);
  for (j=1;j<nv;j++) if (set_min(dmin,p.distancia2(n[ei[j]]))) nc=ei[j];

  return dentro;
}

// a que elemento pertenece un punto (nada = a ninguno)
// esta rutina solo es correcta para conjuntos convexos
// (el octree puede dar un nodo cercano pero no garantias)
// (hay una implementacion similar para esferas)
// requiere conexion recta entre un nodo cercano (octree) y el punto
// es decir que anda mal para mallas de lineas!!
int malla::dequien(const punto &p){
  int nc=-1,ie=-1;
  if (o) nc=o->close(p);
  dequien(p,ie,nc); return ie;
}
void malla::dequien(const punto &p, int &ie, int &nc){
  int j,inmin;
  double dmin,hi=0;
  static pline nni(32);

  if (ie>=0) { // verifico que realmente sea cercano (5h)
    if (inter(p,ie,nc)) return; // ya estaba dentro (probable)
    dmin=p.distancia2(n[nc]);
    if (hayh)
    hi=n[nc].h;
    if (hi==0) hi=n[e[ie][0]].distancia2(n[e[ie][e[ie].nv()/2]]);
    if (dmin>25*hi*hi) ie=-1; //d>5h
  }
  else { // si ie no viene dado, usa nc o busca un nodo cercano
    if (nc>=0&&n[nc].e) { // verifico que realmente sea cercano (5h)
      dmin=p.distancia2(n[nc]);
      hi=n[nc].h;
      if (hi==0) {
        const elemento &ei=e[n[nc].e[0]];
        hi=n[ei[0]].distancia2(n[ei[ei.nv()/2]]);
      }
      if (dmin>25*hi*hi) ie=-1; //d>5h
    }
    else {
      if (!o) mk_octree();
      nc=o->close(p); // nodo cercano
    }
    if (nc==-1||!n[nc].e) return; // sin elementos!
  }

  // busco el nodo mas cercano
  dmin=p.distancia2(n[nc]);
  while (nc>=0){
    // busca entre los nodos de los elementos si hay uno mas cercano
    nn1(nc,&nni); inmin=nc;
    for (j=0;j<nni.len;j++) if (set_min(dmin,p.distancia2(n[nni[j]]))) inmin=nni[j];
    if (inmin==nc) break;
    nc=inmin;
  }
  // nc es el mas cercano (puede que no, si no habia condicion convexa)
  // verifica si esta en alguno de sus elementos
  const cpline enc=n[nc].e;
  for (j=0;j<enc.len;j++) {ie=enc[j]; if (inter(p,ie)) return;}
  ie=-1; // ninguno!!  
}

//=====================================================================================
// nodos vecinos 

// nodos de los elementos de cada nodo
// sin orden ojo: tambien contiene diagonales
bool malla::mk_nn(bool remake){
  if (nn&&!remake) return true;
  if (!e) return false;
  _initime;

  int i,j,k,nv,in,nlen=n.len;
  nn.clean(); nn.resize(nlen);
  ordlist nni;
  for (i=0;i<nlen;i++){
    nni.clean();
    const nodo &ni=n[i];
    const cpline &eni=ni.e;
    for (j=0;j<eni.len;j++){
      const elemento &ej=e[eni[j]];
      nv=ej.nv();
      for (k=0;k<nv;k++){
        in=ej[k];
        if (in!=i) nni+=in;
      }
    }
    nn+=nni;
  }
  _savetime(nn);
  return true;
}

// vecinos y vecinos de los vecinos
// limitado por cantidad de layers o cantidad de nodos
///////////////////
// para ordenar
struct dupla{int n;double d;};
static int compare(const void *d1, const void *d2){
  double d=(*((dupla**)d1))->d-(*((dupla**)d2))->d;
  return (d<0)? -1: ((d>0)? 1 :0);
}
///////////////////
// flag1: puesto
// flag2: puestos los nn
bool malla::vecindad_de_nodo(array1<pline> &vecindad,int nlayers,int nnodos){
  if (!nlayers&!nnodos) return false;
  if (!mk_nn()) return false;
  vecindad.resize(n.len);
  pline v1(1000);
  array1<dupla> nd(1000);
  dupla d1; 
  int i,j,k,lastlen,newlen,layers,nlen=n.len,inj,ink;
  static const int flag12=flag1|flag2;
  for (i=0;i<nlen;i++){ // para cada nodo
    v1.clean(); // limpia la lista
    // agrega la primera capa 
    n[i].f.set(flag1|flag2); 
    const cpline nni=nn[i]; v1+=nni;
    for(j=0;j<v1.len;j++) n[nni[j]].f.set(flag1);
    newlen=0; layers=1; 
    // agrega mas capas
    while ((nlayers&&layers<nlayers)||(nnodos&&v1.len<nnodos)){
      layers++; lastlen=newlen; newlen=v1.len;
      for(j=lastlen;j<newlen;j++){
        inj=v1[j]; if (n[inj].f.es(flag2)) continue; // ya se pusieron los vecinos
        n[inj].f.set(flag2);
        // agrega los vecinos del nodo inj
        const cpline nnj=nn[inj];
        for(k=0;k<nnj.len;k++) {
          ink=nnj[k];
          if (n[ink].f.es(flag1)) continue; // ya se puso
          n[ink].f.set(flag1);
          v1+=ink;
        }
      }
    }    
    for(j=0;j<v1.len;j++) n[v1[j]].f.reset(flag12); // resetea flags
    n[i].f.reset(flag12);
    if (nlayers) {
      vecindad+=v1;
      continue; // ya esta
    }
    // limitado por cantidad => se tiran los mas lejanos (bonito quilombo!!)
    punto ni=n[i]; // para abaratar el calculo de distancias
    nd.clean(); nd.resize(v1.len);// vector de nodos y distancias
    for(j=0;j<v1.len;j++){ 
      d1.n=inj=v1[j]; d1.d=ni.distancia2(n[inj]); nd+=d1;
    }
    qsort((void*)nd.list,(size_t)nd.len,SZP,compare);
    // rearma la lista
    v1.clean();for(j=0;j<nnodos;j++) v1+=nd[j].n;
    vecindad+=v1;
  }
  return true;
}

// sin orden
void malla::nn1(int i, pline *poneraca, bool remake){
  if (!remake&&nn&&poneraca) {*poneraca=nn[i]; return;}

  int j,k,nv,in;
  ordlist nni; nni.clean();
  const nodo &ni=n[i]; const cpline &eni=ni.e;
  for (j=0;j<eni.len;j++){
    const elemento &ej=e[eni[j]]; nv=ej.nv();
    for (k=0;k<nv;k++) {in=ej[k]; if (in!=i) nni+=in;}
  }
  if (!poneraca) nn[i]=nni; else (*poneraca)=nni;
}

// genera una lista de h a partir del tamanio de elementos
// modo: 0 min    1 med    2 max
static bool _mk_h_nn(malla &m, bool remake, int modo){
  array1<nodo> &n=m.n; array1<cpline> &nn=m.nn; array1<elemento> &e=m.e;
  if (!e) return false;
  if (m.hayh&&!remake) return true;
  _initime;
  int i,j,k,nv,in,nlen=m.n.len-m.nodosh;
  double hi,hmin=m.hmin=MAXREAL;
  ordlist nni(64);

  nn.clean(); nn.resize(nlen); nn.len=nlen;
  for (i=0;i<nlen;i++) {
    const nodo &p=n[i]; const cpline &eni=p.e; nni.clean();
    for (j=0;j<eni.len;j++){
      const elemento &ej=e[eni[j]];
      nv=ej.nv();
      for (k=0;k<nv;k++){
        in=ej[k];
        if (in!=i) nni+=in;
      }
    }
    nn[i]=nni;

    if (modo==1){ // medio
      hi=0;
      for (j=0;j<nni.len;j++) hi+=p.distancia(n[nni[j]]);
      hi/=nni.len;
    }
    else if (modo==0){ // minimo
      hi=MAXREAL;
      for (j=0;j<nni.len;j++) p.distancia2_menor(n[nni[j]],hi);
      hi=sqrt(hi);
    }
    else if (modo==2){ // maximo
      hi=0;
      for (j=0;j<nni.len;j++) set_max(hi,p.distancia2(n[nni[j]]));
      hi=sqrt(hi);
    }

    n[i].h=hi; set_min(hmin,hi);
  }
  // nodos de h
  nlen=m.n.len; for(;i<nlen;i++) set_min(hmin,n[i].h);
  m.hayh=true; m.hmin=hmin;
  _savetime(mk_h_nn);
  return true;
}

bool malla::mk_h_nn_min(bool remake){
  return _mk_h_nn(*this, remake, 0);
}
bool malla::mk_h_nn_med(bool remake){
  return _mk_h_nn(*this, remake, 1);
}
bool malla::mk_h_nn_max(bool remake){
  return _mk_h_nn(*this, remake, 2);
}

double malla::h_nn_med(int i,bool mknn) const{
  const nodo &p=n[i]; const pline &eni=p.e;
  int j,k,nv,in;
  ordlist nni(64);
  for (j=0;j<eni.len;j++){
    const elemento &ej=e[eni[j]];
    nv=ej.nv();
    for (k=0;k<nv;k++){
      in=ej[k];
      if (in!=i) nni+=in;
    }
  }
  if (mknn) nn[i]=nni;

  double hmed=0;
  for (j=0;j<nni.len;j++) {
    hmed+=p.distancia(n[nni[j]]);
  }hmed/=nni.len;

  return hmed;
}

double malla::h_nn_min(int i,bool mknn) const{
  const nodo &p=n[i]; const pline &eni=p.e;
  int j,k,nv,in;
  ordlist nni(64);
  for (j=0;j<eni.len;j++){
    const elemento &ej=e[eni[j]];
    nv=ej.nv();
    for (k=0;k<nv;k++){
      in=ej[k];
      if (in!=i) nni+=in;
    }
  }
  if (mknn) nn[i]=nni;

  double hmin_i=MAXREAL;
  for (j=0;j<nni.len;j++) p.distancia2_menor(n[nni[j]],hmin_i);
  return sqrt(hmin_i);
}

double malla::h_nn_max(int i,bool mknn) const{
  const nodo &p=n[i]; const pline &eni=p.e;
  int j,k,nv,in;
  ordlist nni(64);
  for (j=0;j<eni.len;j++){
    const elemento &ej=e[eni[j]];
    nv=ej.nv();
    for (k=0;k<nv;k++){
      in=ej[k];
      if (in!=i) nni+=in;
    }
  }
  if (mknn) nn[i]=nni;

  double hmax=0;
  for (j=0;j<nni.len;j++) set_max(hmax,p.distancia2(n[nni[j]]));
  return sqrt(hmax);
}

// La minima (maxima longitud de arista) de entre sus elementos
double malla::h_min_max_arista(int in) const {
  const cpline &eni=n[in].e;
  int i,j,k,nv,nc;
  double min=MAXREAL,max;
  if (tipo.es(m_sup)){
    for (i=0;i<eni.len;i++){
      const elemento &c=e[eni[i]]; nv=c.nv();
      for(max=0,j=0;j<nv;j++){
        set_max(max,n[c[j]].distancia2(n[c[(j+1)%nv]]));
      }
      set_min(min,max);
    }
  }
  else if (tipo.es(m_vol)){
    for (i=0;i<eni.len;i++){
      const elemento &ei=e[eni[i]]; nc=ei.nc();
      for(max=0,j=0;j<nc;j++){
        const elemento c=ei.cara(j); nv=c.nv();
        for(k=0;k<nv;k++){
          set_max(max,n[c[k]].distancia2(n[c[(k+1)%nv]]));
        }
      }
      set_min(min,max);
    }
  }
  else if (tipo.es(m_lin)){
    for (i=0;i<eni.len;i++){
      const elemento &s=e[eni[i]];
      n[s[0]].distancia2_menor(n[s[1]],min);
    }
  }
  return sqrt(min);
}
bool malla::mk_h_min_max_arista(bool remake){
  if (!e) return false;
  if (hayh&&!remake) return true;
  for(int i=0;i<n.len;i++) {
    n[i].h=h_min_max_arista(i);
    set_min(hmin,n[i].h);
  }
  hayh=true;
  return true;
}


//========================================
// esfera del elemento
// si pdir<>0 (y es 1D o 2D) calcula la normal y la pone ahi
// ojo elementos planos en el espacio (area vector)
bool malla::esfera_e(int ie, punto &c, double &r, double &v, punto *pdir){
  elemento &ei=e[ie];
  e_tipo etipo=ei.tipo();
  if (etipo==e_tetraedro){
    s4(n[ei[0]],n[ei[1]],n[ei[2]],n[ei[3]],c,r,v);
    if (r<ERRADM||v<ERRADM){
      punto a;
      c3(n[ei[0]],n[ei[1]],n[ei[2]],c,r,a);
    }
    return (r!=0);
  }
  else if (etipo==e_triangulo){
    punto a;
    c3(n[ei[0]],n[ei[1]],n[ei[2]],c,r,a);
    if (tipo.es(m_planaxy)) v=a[2];
    else v=a.mod();
    if (pdir) *pdir=a/v;
    return (r!=0);
  }
  else if (etipo==e_cuadrilatero){
    punto a1,a2,c1,c2; double v1,v2,r1,r2;

    c3(n[ei[0]],n[ei[1]],n[ei[2]],c1,r1,a1);
    if (tipo.es(m_planaxy)) v1=a1[2];
    else v1=a1.mod();

    c3(n[ei[2]],n[ei[3]],n[ei[0]],c2,r2,a2);
    if (tipo.es(m_planaxy)) v2=a2[2];
    else v2=a2.mod();

    c=(c1+c2)/2; r=r1+r2/2; v=v1+v2;
    if (pdir) *pdir=(a1+a2).dir();;
    return (r!=0);
  }
  else if (etipo==e_punto){c=n[ei[0]]; r=0; v=0; return false;}
  else if (etipo==e_segmento){
    c=(n[ei[0]]+n[ei[1]])/2;
    r=n[ei[0]].distancia(c); v=2*r;    
    if (pdir&&tipo.es(m_planaxy)) *pdir=(((c-n[ei[0]]).giro90())*v).dir();
    return (r!=0);
  }
  else if (etipo==e_indefinido) {r=0; return false;}
  else{// poligono/poliedro .... no me gusta    
    // esfera por varios puntos:
    // origen en P0 (r^2=C^2)
    // sumi[(Pi-C)^2-C^2]^2=sumi[Pi^2-2*Pi*C]^2=min
    // Dj{sumi[Pi^2-2*Pi*C]^2}=-4*sumi(Pi^2-2Pi*C)*Pi(j)=0
    // sumi(Pi^2*Pi(j))-2*sumk[sumi(Pi(j)*Pi(k))*C(k)]=0
    // sumk[sumi(Pi(j)*Pi(k))*C(k)]=sumi(Pi^2*Pi(j))/2
    bool es2d=ei.dim()==2;
    int i,j,k,nv=ei.nv()-1; if (es2d) nv++;
    punto P0=n[ei[0]];
    static const int maxv=32;
#ifndef _DEBUG
    static 
#endif
      punto P[maxv],sumiPijPik[3],sumiPi2Pij;
    double P2[maxv];

    for (i=0;i<nv-1;i++){// punto
      P[i]=n[ei[i+1]]-P0; P2[i]=P[i].mod2();
    }
    if (es2d) {P[i]=(P[0]%P[1]).dir(); P2[i]=1;}
    else {P[i]=n[ei[i+1]]-P0; P2[i]=P[i].mod2();}

    for (j=0;j<3;j++){ // ecuacion
      sumiPi2Pij[j]=0;
      sumiPijPik[j].zero();
      for (i=0;i<nv;i++) {// punto
        sumiPi2Pij[j]+=P2[i]*P[i][j];
        for (k=0;k<3;k++) { // componente
          sumiPijPik[j][k]+=P[i][j]*P[i][k];
        }
      }
    }
    sumiPi2Pij/=2;
    c.solve(sumiPijPik,sumiPi2Pij,v);
    if (es2d) c-=P[nv-1]*(c*P[nv-1]);
    if (v==0) r=0;
    else {r=c.mod(); c+=P0; v=0;} // no es volumen
    return (r!=0);
  }
}

bool malla::mk_esferas(bool remake, bool mkdir){
  _initime;
  if (ce&&re&&ve&&(mkdir&&edir&&dir)&&!remake) return true;
  int elen=e.len; punto c; double r,v; e_tipo t;
  ce.ini(); ce.resize(elen);
  re.ini(); re.resize(elen);
  ve.ini(); ve.resize(elen);
  if (mkdir){
    dir.ini(); 
    dir.resize(elen); edir=true; ndir=false;
  }
  punto *pdir=(mkdir&&tipo.es_alguno(m_sup|m_planaxy))? new punto : 0;
  for (int i=0;i<elen;i++){
    esfera_e(i,c,r,v,pdir);
    ce+=c; re+=r;
    t=e[i].tipo();
    if (t!=e_segmento&&t!=e_triangulo&&
      t!=e_tetraedro&&t!=e_cuadrilatero) v=volumen(e[i]);
    ve+=v;
    if(mkdir) dir+=(*pdir);
  }
  delete pdir;
  _savetime(esferas);
  return true;
}

// deben estar orientados y haber vecinos
// es de malla y no de elemento por los elementos de nodo y por los vecinos
void malla::juntae(int ie1, int ie2){
  nn.ini();
  elemento &e1=e[ie1], &e2=e[ie2];
  if (e1.pd) e1.pd->rm_cn();
  bool es2d=tipo.es(m_planaxy);

  int nv1=e1.nv(), nv2=e2.nv();
  int j,k,nt,nt1,nt2;
  const int *t1,*t2;
  cpline &v1=vecino[ie1], &v2=vecino[ie2];

  if (es2d){
    // los nodos forman un loop
    // las caras tienen el mismo loop que los nodos
    int i1=0,i2,in,iv;
    cpline n1(nv1,e1.n),n2(nv2,e2.n),
      newn(nv1+nv2-2),newv(nv1+nv2-2);
    while (i1<nv1){
      in=n1[i1];
      if (v1[i1]==ie2){  // cara comun
        i2=n2.index(in);
        while ((iv=v2[i2])!=ie1) {
          newn+=in=n2[i2]; newv+=v2[i2]; i2++;
          // arregla vecino del vecino
          if (iv>=0) { // si no es frontera
            cpline &vv=vecino[iv]; vv.replaceall(ie2,ie1);
          }
          // arregla eldenod
          cpline &ne=n[in].e;
          if (!n1.have(in))
            ne[ne.index(ie2)]=ie1;
          else
            ne.remove1(ie2);
        }
        i1=n1.index(in=n2[i2]);
        n[in].e.remove1(ie2); // arregla eldenod
        if (!i1) break;
      }
      newn+=in; newv+=v1[i1]; i1++;
    }
    // dos triangulos dan un poligono (por los puntos de integracion)
    e1.Tipo=e_poligono;
    if (!e1.pd) {e1.pd=new poly_dat; nt1=e1.pd->nt=1; t1=t4;}
    else {nt1=e1.pd->nt; t1=e1.pd->t;}
    e1.pd->nv=e1.pd->nc=int(newn.len);
    vecino[ie1]=newv;

    // tetraedros del elemento
    if (!e2.pd) {nt2=1; t2=t4;} else {nt2=e2.pd->nt; t2=e2.pd->t;}
    nt=nt1+nt2; int *t=new int[3*nt];
    for (k=0;k<nt1;k++) {
      for (j=0;j<3;j++) t[3*k+j]=newn.index(n1[t1[3*k+j]]);
    }
    for (k=0;k<nt2;k++) {
      for (j=0;j<3;j++) t[3*(nt1+k)+j]=newn.index(n2[t2[3*k+j]]);
    }
    if (e1.pd->t) delete[] e1.pd->t; e1.pd->t=t; e1.pd->nt=nt;

    // nodos
    delete [] e1.n; e1.n=new int[newn.len]; memcpy(e1.n,newn.vertex,newn.len*SZI);
  }
  else{
    // poliedro, no importa el orden de los nodos
    // no suponer una sola cara comun puede haber mas de una o ninguna!
    // newv (vecinos) tiene la cantidad de caras
    // nodos
    int i,olen;
    ordlist newn(nv1+nv2-3);
    for (i=0;i<nv1;i++) newn+=e1.n[i];
    for (i=0;i<nv2;i++) {
      olen=newn.len;
      newn+=e2.n[i];
      // arregla eldenod
      cpline &ne=n[e2.n[i]].e;
      if (newn.len>olen) // nodo de e2 que no estaba en e1
        ne[ne.index(ie2)]=ie1;
      else ne.remove1(ie2); // estaba en ambos
    }
    // caras
    cpline &v1=vecino[ie1], &v2=vecino[ie2];
    int nc1=v1.len, nc2=v2.len;
    int *newc=new int[3*(nc1+nc2)];
    int iv,ixc=0;
    cpline newv(nc1+nc2);
    for (i=0;i<nc1;i++){
      if (v1[i]==ie2) continue; // cara comun
      newv+=v1[i];
      elemento ec(e1.cara(i));
      newc[ixc++]=newn.index(ec.n[0]);
      newc[ixc++]=newn.index(ec.n[1]);
      newc[ixc++]=newn.index(ec.n[2]);
    }
    for (i=0;i<nc2;i++){
      if ((iv=v2[i])==ie1) continue; // cara comun
      newv+=iv;
      if (iv>=0) vecino[iv].replaceall(ie2,ie1);// si no es frontera
      elemento ec(e2.cara(i));
      newc[ixc++]=newn.index(ec.n[0]);
      newc[ixc++]=newn.index(ec.n[1]);
      newc[ixc++]=newn.index(ec.n[2]);
    }
    e1.Tipo=e_poliedro; 
    if (!e1.pd) {e1.pd=new poly_dat; nt1=e1.pd->nt=1; t1=t4;}
    else {nt1=e1.pd->nt; t1=e1.pd->t;}
    e1.pd->nv=int(newn.len); e1.pd->nc=newv.len;
    if (e1.pd->c) delete[] e1.pd->c; e1.pd->c=newc;
    vecino[ie1].roba(newv);

    // tetraedros del elemento
    if (!e2.pd) {nt2=1; t2=t4;} else {nt2=e2.pd->nt; t2=e2.pd->t;}
    nt=nt1+nt2; int *t=new int[4*nt];
    for (k=0;k<nt1;k++) for (j=0;j<4;j++) t[4*k+j]=newn.index(e1.n[t1[4*k+j]]);
    for (k=0;k<nt2;k++) for (j=0;j<4;j++) t[4*(nt1+k)+j]=newn.index(e2.n[t2[4*k+j]]);
    if (e1.pd->t) delete[] e1.pd->t; e1.pd->t=t; e1.pd->nt=nt;

    // nodos
    delete [] e1.n; e1.n=new int[newn.len]; memcpy(e1.n,newn.vertex,newn.len*SZI);
  }

  // arregla esferas
  double vol1=1,vol2=1,vol=2;
  if (ve) {vol1=ve[ie1]; vol2=ve[ie2]; vol=vol1+vol2;}
  if (ce) ce[ie1]=(ce[ie1]*vol1+ce[ie2]*vol2)/vol;
  if (re) re[ie1]=(re[ie1]*vol1+re[ie2]*vol2)/vol;
  e1.f.set(e2.f); // orea flags
  vuelae(ie2,false); // borra sin marcar frontera
}

// desarma poligonos o poliedros
bool malla::p2t(int ie){
  bool desarmo=false,hayv=vecino,hays=(ce||ve||re);
  // copia porque el i vuela
  cpline vi; if (hayv) vi=vecino[ie];
  elemento ei=e[ie];
  e_tipo etipo=ei.tipo();
  if (etipo==e_poliedro||etipo==e_poligono){
    int i,j,k,ic,icv,nt=ei.pd->nt,*t=ei.pd->t,nv;
    int ix,iv;
    double r,v; punto c;
    pline newe(nt); newe.len=nt;
    cpline newv(4); newv.len=4;
    elemento et;
    if (etipo==e_poliedro) {nv=4; et.tipo(e_tetraedro);}
    else {nv=3; et.tipo(e_triangulo);}

    vuelae(ie,false); // de movida nomas

    for (i=0;i<nt;i++,t+=nv){
      j=nv; while (--j>=0) et[j]=ei[t[j]];
      newe[i]=ix=e+=et; // agrega el tetraedro
      j=nv; while (--j>=0) n[et[j]].e+=ix;
      // esferas
      if (hays) {
        esfera_e(ix,c,r,v);
        if (ce) ce+=c; if (re) re+=r; if (ve) ve+=v;
      }
      // vecinos (por cara porque puede haber vecinos por dos caras)
      if (hayv){
        j=nv; while (--j>=0){ // cada cara del tetra
          elemento ct=et.cara(j);
          ic=ei.icara(ct); // la ubica en ei
          icv=iv=-1;
          if (ic<0)  //interior
            for (k=0;k<i;k++) // busca en los nuevos tetras
              if ((icv=e[iv=newe[k]].icara(ct))>=0) break;
          else {
            iv=vi[ic];
            if (iv>=0) icv=e[iv].icara(ct);
          }
          newv[j]=iv; if (iv>=0) vecino[iv][(int)icv]=ix;
        }
        vecino+=newv;
      }
    }    
    desarmo=true;
  }
  if (desarmo) nn.ini();
  return desarmo;
}

// volumen
double malla::volumen(int ie) const{return volumen(e[ie]);}

// areas o voumen orientados de simplices
static double svol(const punto *p, int i1, int i2, int i3, int i4, double &vt){
  double v=triple(p[i2]-p[i1],p[i3]-p[i1],p[i4]-p[i1]); vt+=v;
  return v;
}
static double sarea(const punto *p, int i1, int i2, int i3,const punto &vert, double &at){
  double a=triple(p[i2]-p[i1],p[i3]-p[i1],vert); at+=a;
  return a;
}

// volumen(elemento)
double malla::volumen(const elemento &ei) const{
  e_tipo t=ei.tipo();
  if (t==e_tetraedro){
    punto p0=n[ei[0]];
    return triple(n[ei[1]]-p0,n[ei[2]]-p0,n[ei[3]]-p0)/6;
  }
  else if (t==e_triangulo){ // area2d solo sirve en 2d
    punto p0=n[ei[0]],x=(n[ei[1]]-p0)%(n[ei[2]]-p0);
    double m=fabs(x.mod2());
    if (m<ERRADM) return 0;
    return sqrt(m)/2;
  }
  else if (t==e_segmento){
    return n[ei[0]].distancia(n[ei[1]]);
  }
  else { // no simplice
    int i,nv=ei.nv();
    double vt=0;
    punto vert,*ni=new punto[nv]; for(i=0;i<nv;i++) ni[i]=n[ei[i]];

    if (t==e_cuadrilatero){
      vert=((ni[1]-ni[0])%(ni[2]-ni[0])); vt=vert.mod(); vert/=vt;
      sarea(ni,0,2,3,vert,vt);
      vt/=2;
    }
    else if (t==e_cubo){
      svol(ni,0,1,3,5,vt);
      svol(ni,0,5,3,4,vt);
      svol(ni,3,4,5,7,vt);
      svol(ni,2,3,1,5,vt);
      svol(ni,2,3,5,7,vt);
      svol(ni,5,2,7,6,vt);
      vt/=6;
    }
    else if (t==e_wedge){
      svol(ni,0,1,2,4,vt);
      svol(ni,0,4,2,3,vt);
      svol(ni,2,3,4,5,vt);
      vt/=6;
    }
    else if (t==e_poligono){
      vert=((ni[1]-ni[0])%(ni[2]-ni[0])); vt=vert.mod(); vert/=vt;
      for (i=2;i<nv-1;i++) sarea(ni,0,i,i+1,vert,vt);
      vt/=2;
    }
    else if (t==e_poliedro){
      int i1,i2,i3,nc=3*ei.nc();
      for (i=0;i<nc;i+=3){ // caras
        i1=ei.pd->c[i  ]; if (!i1) continue;
        i2=ei.pd->c[i+1]; if (!i2) continue;
        i3=ei.pd->c[i+2]; if (!i3) continue;
        svol(ni,0,i1,i2,i3,vt);
      }
      vt/=6;
    }
    delete [] ni;
    return vt;
  }
}

double malla::volumen() {
  double vtot=0;
  int i;
  if (ve.len!=e.len){
    double v;
    ve.resize(e.len);
    for (i=0;i<e.len;i++) {vtot+=v=volumen(e[i]); ve+=v;}
  }
  else
    for (i=0;i<e.len;i++)
      vtot+=ve[i];

  return vtot;
}

void malla::swap(){for (int i=0;i<e.len;i++)  swap(i);}


void malla::swap(int i){
  elemento & ei=e[i];
  int *ne=ei.n;
  const e_tipo Tipo=ei.tipo();
  if (ei.dim()<3) { // en 2 d invierte el orden
    int nv=ei.nv();
    for (int j=0;j<(nv>>1);j++) {
      Swap(ne[j],ne[nv-1-j]);
      if (vecino) Swap(vecino[i][j],vecino[i][nv-1-j]);
    }
    if (vecino) vecino[i].first(1);
  }
  else if (Tipo==e_tetraedro) {
    Swap(ne[1],ne[2]); if (vecino) Swap(vecino[i][1],vecino[i][2]);
  }
  else if (Tipo==e_cubo)  {
    Swap(ne[1],ne[3]);Swap(ne[5],ne[7]);
    if (vecino) {Swap(vecino[i][1],vecino[i][2]);Swap(vecino[i][5],vecino[i][4]);}
  }
  else if (Tipo==e_wedge) {
    Swap(ne[1],ne[2]);Swap(ne[4],ne[5]);
    if (vecino) {Swap(vecino[i][1],vecino[i][2]);Swap(vecino[i][4],vecino[i][5]);}
  }
  else if (ei.pd) ei.pd->swap();
  if (ve) ve[i]=-ve[i];
}

// punto y volumen de integracion
void malla::pi(int ie, punto *pi, double *vi) const{
  const poly_dat *pd=e[ie].pd;
  const int *in=e[ie].n;
  const int *ti,*t=(pd)? pd->t : t4;
  int nt=(pd) ? pd->nt : 1;
  int i;
  punto n0,n1,n2,n3;
  if (tipo.es(m_planaxy)){
    for (i=0;i<nt;i++) {
      ti=&t[3*i];
      n0=n[in[ti[0]]]; n1=n[in[ti[1]]]; n2=n[in[ti[2]]];
      pi[i]=(n0+n1+n2)/3;
      vi[i]=(ve&&nt==1)? ve[ie] : area2d(n1-n0,n2-n0);
    }
  }
  else {
    for (i=0;i<nt;i++) {
      ti=&t[4*i];
      n0=n[in[ti[0]]]; n1=n[in[ti[1]]]; n2=n[in[ti[2]]]; n3=n[in[ti[3]]];
      pi[i]=(n0+n1+n2+n3)/4;
      vi[i]=(ve&&nt==1)? ve[ie] : triple(n1-n0,n2-n0,n3-n0)/6;
    }
  }
}

// intercambia dos elementos
void malla::swap_e(int ie1,int ie2){
  if (ie1==ie2) return;
  elemento &e1=e[ie1]; elemento &e2=e[ie2];
  int j,v;
  bool hayv=vecino;

  // elementos de los nodos
  // creo que no importa si un nodo esta en los dos, salvo por el orden
  for (j=0;j<e1.nv();j++) { // nodos de e1
    n[e1[j]].e.replace1(ie1,ie2);
  }
  for (j=0;j<e2.nv();j++) { // nodos de e2
    n[e2[j]].e.replace1(ie2,ie1);
  }

  // vecinos
  if(hayv){
    cpline &v1=vecino[ie1];
    for (j=0;j<v1.len;j++) {
      v=v1[j];
      if(v>=0)
        vecino[v].replace1(ie1,ie2);
    }
    cpline &v2=vecino[ie2];
    for (j=0;j<v2.len;j++) {
      v=v2[j];
      if(v>=0)
        vecino[v].replace1(ie2,ie1);
    }
    vecino.swap(ie1,ie2);
  }

  // esferas y volumen
  if (re) re.swap(ie1,ie2);
  if (ce) ce.swap(ie1,ie2);
  if (ve) ve.swap(ie1,ie2);

  if (edir&&dir) dir.swap(ie1,ie2);

  // los elementos
  e.swap(ie1,ie2);
}

// elimina un elemento swappeando con el ultimo
// (hay rutinas que esperan el swap (en vez de la eliminacion o copia)
//  son change alpha y delta)
// los replace son replace1 y no replaceall (si hay mas de uno pasara mas de una vez)
void malla::vuelae(int i,bool mark_f){
//  nn.ini(); OJO!!! no lo borra pero queda invalido
  elemento &ei=e[i]; if (mark_f) ei.f.reset(e_frontera);
  int j,v,nv=ei.nv(),last=e.len-1;
  bool hayv=vecino;
  // elimina i de los nodos del elemento y los marca como frontera
  for (j=0;j<nv;j++) {
    nodo &nj=n[ei[j]];
    nj.e.remove1(i);
    if (mark_f) nj.f.set(n_frontera);
    if (hayv){
      v=vecino[i][j];
      if(v<0) continue;
      vecino[v].replace1(i,-i-2);
      if (mark_f) e[v].f.set(e_frontera);
    }
  }

  // si es el ultimo, simplemente lo elimina
  if (i==last){
    if (re) re.len--; if (ce) ce.len--; if (ve) ve.len--;
    if (hayv) vecino.len--;
    if (edir&&dir&&dir.len==e.len) dir.len--;
    e.len--;
    return;
  }
  // el i no es el ultimo  
  // swap con el ultimo
  const elemento &el=e[last]; nv=el.nv();
  for (j=0;j<nv;j++) n[el[j]].e.replace1(last,i);
  if(hayv){
    cpline &vl=vecino[last];
    for (j=0;j<vl.len;j++) {v=vl[j]; if(v>=0) vecino[v].replace1(last,i);}
    vecino.swap(i,last); vecino.len--;
  }
  if (re) {re.swap(i,last); re.len--;}
  if (ce) {ce.swap(i,last); ce.len--;}
  if (ve) {ve.swap(i,last); ve.len--;}
  if (edir&&dir){
    if (dir.len==e.len) {dir.swap(i,last); dir.len--;}
    else rm_dir();
  }
  e.swap(i,last); e.len--;
}

// extrae elementos que tienen todos sus nodos con un mismo bit de f encendido
// en ex pone los elementos y en ix los indices (alguno puede ser 0)
// si qty es -1 son los que tienen todos sus nodos con f
// si almenos=true testea >= y si es false testea que sea exactamente =
// los ix y ex se corresponden con los de la malla en forma monotona creciente
int malla::filtra_elms(array1<elemento> *ex,pline *ix,flagtype f,int qty, bool almenos) const{
  if(!qty||!f||(!ix&&!ex)) return 0;
  int i,j,qi,nv,elen=e.len;
  bool todos=(qty==-1);
  if (ex) ex->clean(); if (ix) ix->clean();
  for(i=0;i<elen;i++){ // creciente!
    const elemento& ei=e[i]; nv=ei.nv();
    if (todos) qty=nv;
    qi=j=0;
    do if(n[ei[j]].f.es_todos(f)) qi++; while(++j<nv&&qi<qty);
    if(j!=nv&&!almenos) do if(n[ei[j]].f.es_todos(f)) qi++; while(++j<nv);
    if(qi!=qty) continue;
    if (ex) *ex+=ei;
    if (ix) *ix+=i;
  }
  if (ex) return ex->len;
  return ix->len;
}

// elementos que tienen todos sus nodos con algun bit de f encendido
// pero todos el/los mismo/s bit/s
int malla::mask_elms(array1<elemento> *ex,pline *ix,flagtype f) const{
  if(!f||(!ix&&!ex)) return 0;
  int i,j,nv,elen=e.len;
  flagtype mask;
  if (ex) ex->clean(); if (ix) ix->clean();
  for(i=0;i<elen;i++){ // creciente!
    const elemento& ei=e[i]; nv=ei.nv();
    for (mask=f,j=0;j<nv&&mask;j++) mask&=n[ei[j]].f;
    if(!mask) continue;
    if (ex) *ex+=ei;
    if (ix) *ix+=i;
  }
  if (ex) return ex->len;
  return ix->len;
}

// convierte elementos 3d en 2d o 2d en 1d
bool malla::baja_dim(){
  if (!mk_vecino()) return false;
  array1<elemento> newe(e.len);
  array1<cpline> newen(n.len); newen.len=n.len;
  int ie,ic,iv,nc,nv,ix,in,fdim=0;  
  for (ie=0;ie<e.len;ie++){
    const elemento &ei=e[ie];
    nc=ei.nc(); if (nc<3) continue; // ni segmentos ni nodos
    for (ic=0;ic<nc;ic++){
      if (vecino[ie][ic]>ie) continue; // una sola vez, incluyendo frontera
      elemento c=ei.cara(ic); nv=c.nv();
      ix=newe+=c;
      fdim|=(1<<c.dim());
      for (iv=0;iv<nv;iv++) newen[c[iv]]+=ix;
    }
  }
  // reemplaza
  e=newe;  
  for (in=0;in<n.len;in++) {
    n[in].f.reset(n_frontera);
    n[in].e=newen[in];
  }

  // mantenimiento
  rm_vecino(); rm_frontera(); 
  ve.ini(); re.ini(); ce.ini(); rm_dir();
  nn.ini();

  tipo.reset(m_vol|m_sup|m_lin|m_nodos|m_cerrada|m_orientada|m_slivers_marked);
  if (fdim&2) tipo.set(m_lin);
  if (fdim&4) tipo.set(m_sup);
  return true;
}


// elimina aristas largas o triangulos de aristas largas
bool malla::elimina_aristas_largas(){
  if (tipo.es_alguno(m_nodos|m_vol)) return true;
  _initime;
//  mk_h_min_max_arista(true); ???????????????
  mk_h_nn_min();
  rm_vecino(); rm_nn();

  int i,j,k,nv,kmax=-1;
  double h1,h2,d,dmax=3;
  bool volar;

  // primer paso vuela las grandes
  for (i=e.len-1;i>=0;i--) { // al reves porque vuelae swappea
    elemento &ei=e[i]; volar=false;
    const int *nei=ei.n; nv=ei.nv();
    for (j=0;j<nv-1;j++){
      nodo &n1=n[nei[j]]; h1=n1.h;
      for (k=j+1;k<nv;k++){// simplices (si no de uno al siguiente)
        nodo &n2=n[nei[k]]; h2=n2.h;
        d=n1.distancia(n2);
        if ((2*d)/(h1+h2)<dmax) continue;
        volar=true; break;
      }
      if (volar) break;
    }
    if (!volar) continue;
    vuelae(i,false); // no marca frontera
    // i++; //ojo si va de 0 a len hay que retestear
  };

  //segundo paso elimina hasta quedar con dos vecinos
  if (tipo.es(m_sup)){ // malla de superficie
    int l,nnlen;
    bool unovolado;
    cpline nni;
    for (i=0;i<n.len-1;i++){
      nn1(i,&nni,true); nnlen=nni.len;
      cpline &nei=n[i].e;
      for (j=0;j<nnlen;j++){
        if (nni[j]<i) continue;
        // busca elementos con los dos nodos
        cpliner ea=nei.inters(n[nni[j]].e);
        unovolado=false;
        while (ea.len>2) {
          // busca el elemento con la suma de aristas mas larga y lo elimina
          dmax=0; kmax=-1;
          for (k=0;k<ea.len;k++){
            elemento &ek=e[ea[k]]; nv=ek.nv(); d=0;
            for (l=0;l<nv;l++){
              d+=n[ek[l]].distancia(n[ek[(l+1)%nv]]);
            }
            if (set_max(dmax,d)) kmax=k;
          }
          vuelae(ea[kmax],false);
          ea.remove(1,kmax);
          unovolado=true;
        }
        if (!unovolado) continue;
        nn1(i,&nni,true); nnlen=nni.len; j=0; //testea de nuevo
      }
    }
  }
  else{ // malla de lineas
    for (i=0;i<n.len;i++){
      while (n[i].e.len>2){
        // busca el elemento mas largo y lo elimina
        dmax=0;
        cpline &en=n[i].e;
        for (k=0;k<en.len;k++){
          elemento &ek=e[en[k]]; d=n[ek[0]].distancia(n[ek[1]]);
          if (set_max(dmax,d)) kmax=k;
        }
        vuelae(en[kmax],false);   
      }
    }
  }

  return true;
  _savetime(aristas);
}

