////////////////////////////////////////////////////////////
// Calculo de funciones de forma laplacianas en poliedros //
////////////////////////////////////////////////////////////

#include "esfera.h"
#include "malla.h"

using namespace std;

// Funciones de forma con origen en el punto de calculo

// segmento
static void fforma2(const punto *n,double *f,punto *Df){
  punto l=n[1]-n[0]; double lt=l.mod();
  f[0]=(n[1]*l)/lt; f[1]=1-f[0];
  if (Df) {Df[0][0]=-1/lt; Df[1][0]=1/lt;}
}

// triangulo
static void fforma3(const punto *n,double *f,punto *Df){
  int i;
  double at=(n[1]-n[0]).pv2d(n[2]-n[0]),a[3]; // areas
  for (i=0;i<3;i++) a[i]=n[(i+1)%3].pv2d(n[(i+2)%3]);

  //Funcion de Forma y derivada
  //f[p]=a[p]/at
  //Da[p]/Dx[i]=e[i][j](n[p+1][j]-n[p+2][j])
  for (i=0;i<3;i++) {
    f[i]=a[i]/at;
    if (Df) Df[i]=(n[(i+2)%3]-n[(i+1)%3]).giro90()/at;
  }
}

// tetraedro
static void fforma4(const punto *n,double *f,punto *Df){
  int i,j;
  static punto xv[4][4]; // productos vectoriales
  for (i=0;i<3;i++) for (j=i+1;j<4;j++) {
    xv[i][j]=n[i]%n[j]; xv[j][i]=-xv[i][j];
  }

  // ojo: en tetraedros las caras no ciclan bien, depende de la paridad
  double vt=triple(n[1]-n[0],n[2]-n[0],n[3]-n[0]),v[4]; // volumenes
  for (i=0;i<4;i++) {
    v[i]=n[(i+1)%4]*xv[(i+2)%4][(i+3)%4];
    if (i%2) v[i]=-v[i]; // impar
  }

  //Funcion de Forma y derivada
  //f[p]=v[p]/vt
  //Dv[p]/Dx[i]=-(xv[p+1][p+2]+xv[p+2][p+3]+xv[p+3][p+1])[i]
  punto Dfi;
  for (i=0;i<4;i++) {
    f[i]=v[i]/vt;
    if (Df) {
      Dfi=(xv[(i+2)%4][(i+1)%4]+
           xv[(i+3)%4][(i+2)%4]+
           xv[(i+1)%4][(i+3)%4])/vt;
      if (i%2) Dfi=-Dfi;
      Df[i]=Dfi;
    }
  }
}

// centro y derivada del centro de una esfera
// pi*c=pi^2/2
// D[i][j]= dc[j]/dx[i] = D[i](c[j])
// pj*D[i]=c[i]
// l22 son los mod2()/2 de los vectores
static bool define_s(
  const punto *p, const punto l22,
  punto &c, punto *D
  ){
  static punto xv[3];
  double vt;
  c.solve(p,l22,xv,vt);
  if (vt==0)
    return false;
  /*
  punto c_p[3]={c-p[0],c-p[1],c-p[2]};
  transpose(c_p);
  D[0][0]=(xv[0]*c_p[0])/vt;
  D[0][1]=(xv[1]*c_p[0])/vt;
  D[0][2]=(xv[2]*c_p[0])/vt;

  D[1][0]=(xv[0]*c_p[1])/vt;
  D[1][1]=(xv[1]*c_p[1])/vt;
  D[1][2]=(xv[2]*c_p[1])/vt;

  D[2][0]=(xv[0]*c_p[2])/vt;
  D[2][1]=(xv[1]*c_p[2])/vt;
  D[2][2]=(xv[2]*c_p[2])/vt;
  */

  punto xvt(xv[0].trace(),xv[1].trace(),xv[2].trace());
  D[0]=xvt*(c[0]/vt); D[0][0]-=1;
  D[1]=xvt*(c[1]/vt); D[1][1]-=1;
  D[2]=xvt*(c[2]/vt); D[2][2]-=1;
  return true;
}

// poligono
static bool fformapg(const punto *n,int nv,double *f,punto *Df){
  int i;
  double *l2=new double[nv];// modulos al cuadrado
  for (i=0;i<nv;i++) l2[i]=n[i].mod2();
  // esferas de las caras
  punto ps[3],l22; ps[2]=_ez; l22[2]=.5;
  punto *Ds=new punto [nv*3], *c=new punto[nv];
  for (i=0;i<nv;i++) {
    ps[0]=n[i]; ps[1]=n[(i+1)%nv];
    l22[0]=l2[i]/2; l22[1]=l2[(i+1)%nv]/2;
    if (!define_s(ps,l22,c[i],&(Ds[3*i]))){
      delete[] l2; delete[] Ds; delete[]c;
      return false;
    }
  }

  double fsuma=0,fsuma2; // suma de funciones de forma individuales
  punto Dsuma(0,0,0); // derivada de la suma
  punto v,Dfi; //temporarios
  double ai;
  int iant;
  for (i=0;i<nv;i++){
    iant=(i-1+nv)%nv;
    // calculo de la funcion de forma individual
    ai=c[iant].pv2d(c[i]);
    fsuma+=f[i]=ai/l2[i];
    if (Df){
      Dfi.zero();
      // componentes k del gradiente
      Dfi[0]+=Ds[3*iant].pv2d(c[i])+c[iant].pv2d(Ds[3*i]);
      Dfi[1]+=Ds[3*iant+1].pv2d(c[i])+c[iant].pv2d(Ds[3*i+1]);
      Dsuma+=Df[i]=(Dfi+2*n[i]*f[i])/l2[i];
    }
  }
  fsuma2=fsuma*fsuma;
  for (i=0;i<nv;i++){
    if (Df) Df[i]=(Df[i]*fsuma-f[i]*Dsuma)/fsuma2;
    f[i]/=fsuma; // funcion de forma
  }
  delete[] l2; delete[] Ds; delete[]c;
  return true;
}

// poliedro
static bool fformapl(const punto *n,int nv,
                     const int *cara,int nc,const int **cn,
                     double *f,punto *Df){
  int i,j,k,qc,ic1,ic2;
  bool retval=true;

  double *l2=new double[nv];// modulos al cuadrado
  for (i=0;i<nv;i++) l2[i]=n[i].mod2();
  // esferas de las caras
  punto ps[3],l22;
  punto *Ds=new punto [nc*3], *c=new punto[nc];
  for (i=0;i<3*nc;i+=3) {
    ps[0]=n[cara[i]]; ps[1]=n[cara[i+1]]; ps[2]=n[cara[i+2]];
    l22[0]=l2[cara[i]]/2; l22[1]=l2[cara[i+1]]/2; l22[2]=l2[cara[i+2]]/2;
    if (!define_s(ps,l22,c[i/3],&(Ds[i]))){
      delete[] l2; delete[] Ds; delete[]c;
      return false;
    }
  }

  // caras del poliedro de Voronoi (centros de esferas por cada nodo y el punto)
  // el pie de la perpendicular al plano es el medio entre el punto y el nodo
  double l2i,l4i,fij,fsuma=0,fsuma2; // suma de funciones de forma individuales
  punto Dsuma; Dsuma.zero(); // derivada de la suma
  punto v,Dfi; //temporarios
  for (i=0;i<nv;i++){ // nodo i
    // calculo de la funcion de forma individual
    f[i]=0; Dfi.zero(); qc=cn[i][0]; l2i=l2[i]; l4i=l2i*l2i;
    for (j=0;j<qc;j++){ // arista j
      ic1=cn[i][1+j]; ic2=cn[i][1+(j+1)%qc];
      v=c[ic1]%c[ic2];
      f[i]+=fij=v*n[i];
/*
#ifdef _DEBUG
  // da 0 en los slivers (dos triangulos contiguos mismo plano)
  // puede dar negativo cuando la cara del poliedro de voronoi esta a un lado
  // del segmento del punto al nodo
  assert(fij>-ERRADM);
#endif
*/
      if (Df) for (k=0;k<3;k++)  // componente k del gradiente
        Dfi[k]+=((Ds[3*ic1+k]%c[ic2]+c[ic1]%Ds[3*ic2+k])*n[i]-v*_e[k])/l2i
            +2*fij*n[i][k]/l4i;
    }
    fsuma+=f[i]/=l2i;
    if (Df) Dsuma+=Df[i]=Dfi;
  }
  fsuma2=fsuma*fsuma;
  for (i=0;i<nv;i++){
    if (Df) Df[i]=(Df[i]*fsuma-f[i]*Dsuma)/fsuma2;
    f[i]/=fsuma; // funcion de forma
  }
  delete[] l2; delete[] Ds; delete[]c;
  return retval;
}

// funciones de forma de un punto en un elemento
bool malla::fforma(
  const elemento &ei,    // elemento
  const punto &p,  // punto de calculo (modificable si hay error en concavos)
  double *f, // funciones de forma en el punto
  punto *Df  // gradientes
  )const {

  const int nv=ei.nv();
  if (nv<2) return false; // indefinido o nodo suelto
  const e_tipo tipo=ei.tipo();

  // origen en p
  int i;
  punto *v=new punto [nv];
  for (i=0;i<nv;i++) v[i]=n[ei[i]]-p;

  bool retval=true;

       if(tipo==e_segmento)  fforma2(v,f,Df);
  else if(tipo==e_triangulo) fforma3(v,f,Df);
  else if(tipo==e_tetraedro) fforma4(v,f,Df);
  else if (tipo==e_cuadrilatero||tipo==e_poligono) fformapg(v,nv,f,Df);
  else if (tipo==e_poliedro){
    if (!ei.pd->cn) ei.pd->mk_cn();
    retval=fformapl(v,nv,ei.pd->c,ei.pd->nc,(const int**)(ei.pd->cn),f,Df);
  // verificar en linux!!!!!!!!!!!!!!!!!
    // puedo modificar pd (usando mk_cn) sin que el const de elemento zapatee (es puntero),
    // pero no me deja pasar cn de int** a const int** !!??
  }
  else retval=false;
  delete [] v;
  if (retval) return true;
  // poligonos o poliedros fallados
  for (i=0;i<nv;i++){
    if (Df) Df[i][0]=Df[i][1]=Df[i][2]=MAXREAL;
    f[i]=0;
  }

  return false;
}

bool esfera::fforma(
  const punto &p,  // punto de calculo
  double *f, // funciones de forma en el punto
  punto *Df  // gradientes
  )const {

  if (vt==0) return false;

  // origen en p
  int i;
  punto v[4];
  for (i=0;i<NV;i++) v[i]=(*nod)[n[i]]-p;

       if(NV==2) fforma2(v,f,Df);
  else if(NV==3) fforma3(v,f,Df);
  else if(NV==4) fforma4(v,f,Df);
  return true;
}

