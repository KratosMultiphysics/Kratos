// superficie interior de una polilinea

#include "malla.h"

using namespace std;


#define _add1 {\
    ni.e[0]=nlen-1; ni.e[1]=nlen;\
    n+=ni;\
    s[0]=nlen-1; s[1]=nlen; e+=s;\
}

bool malla::fitsurface(){

  if ((tipo&(m_vol|m_sup|m_lin|m_nodos))!=m_lin)
    return false;
  if (!mk_vecino()) return false;
  int i;
  if (tipo.noes(m_cerrada)) {
    add_error(Open_Mesh);
    for (i=0;i<n.len;i++) if (n[i].f.es(n_frontera)) break;
    add_error(i);
    return false;
  }
  bbox();
  if (!hayh&&!mk_h_nn_min()) return false;
  for (i=0;i<n.len;i++) n[i].h=hmin;
  mesh1d();

  int ndiv,j,len=n.len,&nlen=n.len;

  array1<punto> nbase; // nodos ordenados
  array1<double> vbase; // doubles asignados (para interpolar)
  int in=e[0][0],ie=0;
  // array circular ordenado de nodos
  for(i=0;i<len;i++){
    nbase+=n[in]; vbase+=n[in].v;
    const elemento &ei=e[ie];
    j=(ei.index(in)+1)%2;
    in=ei[j];
    ie=vecino[ie][j];
  }

  // inicializa la malla conservando el nombre y el bbox
  malla copia(*this); // solo por si hay algun return false despues de inicializar
  ini();
  copy_nombre(copia);
  pmin[0]=pmin[1]=-1;pmax[0]=pmax[1]=1;

  double t,factor=DOSPI/len;
  nodo ni,nj;
  ni.h=.03;
  double edgelength=2*sin(factor/2);
  if (edgelength<ni.h) {ndiv=1; ni.h=edgelength;}
  else ndiv=int(edgelength/ni.h+.5);
  ni.e.resize(2); ni.e.len=2;
  nj.e.resize(1); nj.e.len=1; nj.e[0]=0;
  punto p0,p1=_ex;
  elemento s(e_segmento), b(e_poligono);
  { // define el poligono
    int *pd=new int[1+len+3*(len-2)];
    pd[0]=len; for (i=0;i<len;i++) pd[1+i]=i;
    for (i=0;i<len-2;i++) {pd[1+len+3*i]=0;pd[1+len+3*i+1]=i+1;pd[1+len+3*i+2]=i+2;}
    b.lee(pd);
    delete [] pd;
  }

  pline en(2);

  n.resize(len*ndiv);
  e.resize(len*ndiv);
  tipo.set(m_lin|m_planaxy|m_cerrada|m_orientada);
  
  malla mb; mb.e.resize(1); mb.e+=b;
  mb.tipo.set(m_sup|m_planaxy|m_orientada);

  for(i=0;i<len;i++){
    ni=p0=p1; _add1; mb.n+=nj=p0; b[i]=i;
    t=factor*(i+1); p1[0]=cos(t); p1[1]=sin(t);
    for (j=1;j<ndiv;j++) {ni=((ndiv-j)*p0+j*p1)/ndiv; _add1;}
  }
  e.last()[1]=0; n[0].e[0]=nlen-1;  // (mas barato)
  ni.zero(); ni.e.ini(); ni.h*=5; ni.f.set(n_h); n+=ni; nodosh=1;

  if (!mk_puntos()) {
    copia.add_error(m_error);
    *this=copia;
    add_error(copia.m_error);
    return false;
  }

  // posicion de los nodos de frontera
  for(i=0;i<len;i++){
    n[ndiv*i].setpos(nbase[i]);
    for (j=1;j<ndiv;j++) {
      n[ndiv*i+j].setpos(((ndiv-j)*nbase[i]+j*nbase[(i+1)%len])/ndiv);
      n[ndiv*i+j].v=((ndiv-j)*vbase[i]+j*vbase[(i+1)%len])/ndiv;
    }
  }

  // posicion de los nodos interiores
  double *ff=new double[len];
  for (i=ndiv*len;i<nlen;i++){
    mb.fforma(0,n[i],ff,0); n[i].zero(); n[i].v=0;
    for (j=0;j<len;j++) {
      n[i]+=nbase[j]*ff[j];
      n[i].v+=vbase[j]*ff[j];
    }
  }

  if (copia.tipo.es(m_planaxy)) tipo.set(m_planaxy); else tipo.reset(m_planaxy);
  epsilon=copia.epsilon;
  if (o) delete o; o=0;
  pmin=copia.pmin;pmax=copia.pmax;
  hayh=false; hmin=MAXREAL;
  rm_dir();
  return true;
}
