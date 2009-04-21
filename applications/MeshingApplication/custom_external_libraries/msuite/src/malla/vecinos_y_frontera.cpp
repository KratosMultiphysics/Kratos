// caras de poliedros no compartidas

//#define temporizar
#include "tiempo.h"
#include "malla.h"

using namespace std;

bool malla::mk_vecino_1d(bool remake){
  if (vecino&&!remake) return true;
  _initime;
  int in,j,ie,ie0,ie1,qn=n.len-nodosh;

  vecino.clean(); vecino.resize(e.len); vecino.len=e.len;

  cpline tmp(2); tmp.len=2;
  for (ie=0;ie<e.len;ie++) vecino[ie]=tmp;

  tipo.set(m_cerrada);

  for (in=0;in<qn;in++){
    const cpline &en=n[in].e;
    int elen=en.len;
    if (elen==0) continue; // nodo suelto
    if (elen==1) {
      ie=en[0];
      vecino[ie][int(e[ie].index(in))]=-1;
      n[in].f.set(n_frontera);
      tipo.reset(m_cerrada);
    }
    else if (elen>2) { // non-manifold
      add_warning(NonManifold);
      add_warning(" at node: ");
      add_warning(in); 
      for (j=0;j<elen;j++){
        ie=en[j];
        vecino[ie][int(e[ie].index(in))]=-2;
      }
      tipo.reset(m_cerrada);
    }
    else{
      ie0=en[0]; ie1=en[1];
      vecino[ie0][int(e[ie0].index(in))]=ie1;
      vecino[ie1][int(e[ie1].index(in))]=ie0;
    }
  }  
  _savetime(vecinos);
  return true;
}


// un nodo suelto es nodo frontera
bool malla::mk_vecino(bool remake){
  if (vecino&&!remake) return true;
  if (!e) return false;
  bool lin=tipo.es(m_lin),sup=tipo.es(m_sup),vol=tipo.es(m_vol);
  if ((lin&vol)||(lin&sup)||(sup&vol)) {//mezcla
    add_error(Mixed);
    return false; 
  }
  if (lin) return mk_vecino_1d(remake); // lineas

  _initime;

  int i,j,k,l,m,elen=e.len,nlen=n.len,nc,ie,enlen,nv,iv,iv3,ie_v;
  bool tieneali;
  elemento c;

  vecino.clean(); vecino.resize(elen); vecino.len=elen;

  for (i=0;i<elen;i++){
    nc=e[i].nc();
    vecino[i].resize(nc); vecino[i].len=nc;
    e[i].f.reset(e_frontera);
  }
  for (i=0;i<nlen;i++) n[i].f.reset(n_frontera);

  // Para cada nodo, junta las caras de cada elemento del nodo en
  // las cuales el nodo es el menor (una sola vez por cara)
  // Cuenta en la cantidad de elementos que comparten 
  // la cara (1 frontera, 2 interior, mas de dos non manifold)

  // Es mas rapido si no testea si es no manifold (mas de dos vecinos por cara) 
  // (en 0.5Mnodos hay una diferencia del 1% en tiempo => testea)

  array1<elemento> cn(32); // caras del nodo
  pline elc(32),ic(32),qv(32); // elemento e indice correspondientes cantidad de vecinos por esa cara
  bool cerrada=true;

  for (i=0;i<nlen;i++) { // cada nodo
    const cpline &en=n[i].e; enlen=en.len;
    if (!enlen) {n[i].f.set(n_frontera); continue;} // nodo aislado
    cn.clean(); elc.clean(); ic.clean(); qv.clean();
    for(j=0;j<enlen;j++){ // cada elemento del nodo
      ie=en[j]; const elemento &ej=e[ie]; nc=ej.nc();
      for(k=0;k<nc;k++){ // cada cara
        c=ej.cara(k);
        nv=c.nv();
        tieneali=false;
        for (l=0;l<nv;l++) {if (c[l]<i) break; if (c[l]==i) tieneali=true;}
        if (l<nv||!tieneali) continue; // otra cara
        // la cara tiene al nodo y es el menor de la cara
        for (l=0;l<cn.len;l++) {if (cn[l]==c) break;}
        if (l==cn.len) { // agrega la cara
          qv+=1; cn+=c; elc+=ie; ic+=k; // el flag es para contar
        }
        else {// la cara ya fue puesta
          qv[l]++;
          iv=elc[l]; ie_v=ic[l]; // indice del vecino y de la cara en iv
          if (qv[l]==2){//situacion normal
            vecino[ie][k]=iv;
            vecino[iv][ie_v]=ie;
            //cn[l].tipo(e_indefinido); //para no comparar (si para testear manifold)
          }
          else{// non manifold (abierta aunque sea par)
            add_warning(NonManifold);
            for(m=0;m<cn[l].nv();m++) n[cn[l][m]].f.set(n_frontera);
            if (qv[l]==3){// tres con la misma cara
              // ie, el guardado en iv=elc[l] y su vecino
              iv3=vecino[iv][ie_v];
              add_warning(" on elements: ");
              add_warning(ie);
              add_warning(", ");
              add_warning(iv);
              add_warning(", ");
              add_warning(iv3);
              vecino[iv3].replace1(iv,-2);
              vecino[ie][k]=vecino[iv][ie_v]=-2;
              e[ie].f.set(e_frontera);
              e[iv].f.set(e_frontera);
              e[iv3].f.set(e_frontera);
              cerrada=false;
            }
            else{//mas de tres con una cara
              add_warning(" on element: ");
              add_warning(ie);
              vecino[ie][k]=-2;
              e[ie].f.set(e_frontera);
            }
          }
        }
      }
    }
    // si queda alguna cara puesta solo una vez es frontera
    // es cerrada si no tiene caras de frontera (puede tener nodos sueltos)
    for (j=0;j<cn.len;j++){
      if (qv[j]!=1) continue;
      const elemento &cj=cn[j];
      e[elc[j]].f.set(e_frontera);
      vecino[elc[j]][ic[j]]=-1;
      nv=cj.nv();
      for (k=0;k<nv;k++) n[cj[k]].f.set(n_frontera);
      cerrada=false;
    }
  }

  if (cerrada) tipo.set(m_cerrada);  

  _savetime(vecinos);
  return true;
}


// piezas conexas de la frontera
bool malla::mk_frontera(bool remake){
  if ((!e||frontera||tipo.es(m_cerrada))&&!remake) return true;

  //no va para mezcla de tipos
  flagtype fdim=tipo&(m_vol|m_sup|m_lin);
  if (fdim!=m_vol&&fdim!=m_sup&&fdim!=m_lin) {
    add_error(Mixed);
    return false;
  }

  frontera.ini();tipo.reset(m_cerrada);

  if (!mk_vecino()) return false;
  if (tipo.es(m_cerrada)) return true;
  _initime;

  // flag1=puesto, flag2=hecho
  int in,i,j,k,l,nc,nv,inc,iej,lastfound;
  static const int flag12=flag1|flag2;
  int tmpsize=int(1.1*pow(double(e.len),.75));
  flagtype f;
  bool tienealnodo,menor;
  pline
    npuesto(tmpsize), // nodos puestos
    nf(tmpsize);      // nodos puestos en total

  // busca un nodo de frontera
  for (i=0;i<n.len;i++) if (n[i].f.es(n_frontera)&&n[i].e) break;
  if (i==n.len) {tipo.set(m_cerrada); _savetime(frontera); return true;}
  npuesto+=i; nf+=i; n[i].f.set(flag1); lastfound=i;
  double maxx=n[i][0],lastmax=n[i][0]; int imax=i; // punto de maximo x
  // puede haber piezas sueltas, tangentes, elementos contra dos fronteras, etc. (alpha)
  while (npuesto){
    frontera+=cascara();
    cascara &fronti=frontera.last();
    fronti.parent=this;
    for (i=0;i<npuesto.len;i++){
      in=npuesto[i]; n[in].f.set(flag2); // hecho
      cpline &en=n[in].e;
      for (j=0;j<en.len;j++){
        iej=en[j]; 
        _revienta(iej>e.len);
        elemento &ej=e[iej];
        cpline vj=vecino[iej]; nc=vj.len;
        for (k=0;k<nc;k++) {
          if (vj[k]>=0) continue;
          elemento ck=ej.cara(k); 
          nv=ck.nv();
          for(tienealnodo=false,menor=true,l=0;l<nv;l++){
            inc=ck[l];
            if (inc==in) {tienealnodo=true; continue;}
            if (inc<in) menor=false;
            if(n[inc].f.noes(flag1)){
              npuesto+=inc; nf+=inc; n[inc].f.set(flag1);
              if (set_max(maxx,n[inc][0])) imax=inc;
            }
          }
          if (!tienealnodo||!menor) continue;
          ck.f.reset(e_frontera); 
          fronti.e+=ck; 
        }
      }
    }
    fronti.n=npuesto; npuesto.clean();

    // verifica si es la mas exterior
    if (maxx>lastmax&&frontera.len>1){
      frontera.swap(frontera.len-1,0);
    }
    lastmax=maxx;

    // busca nodos de frontera sin agregar
    for(i=lastfound+1;i<n.len;i++){
      f=n[i].f;
      if (f.es(n_frontera)&&f.noes(flag1)&&n[i].e)
        break;
    }
    if (i==n.len) break;
    nf+=i; npuesto+=i; n[i].f.set(flag1); lastfound=i;
  }

  // limpia el flag de nodos
  for(i=0;i<nf.len;i++) {n[nf[i]].f.reset(flag12);}

  _savetime(frontera);
  return true;
}


malla::malla(const cascara &c):_INI_MALLA{
  nombre[0]=ext[0]=0;

  const malla &m=(*c.parent);
  strcpy(nombre,m.nombre); strcat(nombre,"_c");
  if (!c.n) return;
  e=c.e; // copia los elementos
  if (m.tipo.es(m_orientada)) tipo.set(m_orientada);
  tipo.set(m_modificada);
  if (m.hayh) hayh=true;

  int i,j,in,nv,d=0;

  // copia los nodos
  n.resize(c.n.len);
  pline mapnod(m.n.len,nada);
  pmin=m.pmax; pmax=m.pmin;
  if (hayh) hmin=MAXREAL;
  for (i=0;i<c.n.len;i++) {
    in=c.n[i]; nodo ni=m.n[in];
    ni.f.reset(n_frontera);
    if (ni.f&fmask) hayfn=true;
    if (hayh) set_min(hmin,ni.h);
    ni.set_min_max(pmin,pmax);
    ni.e.ini();
    mapnod[in]=n+=ni;
  }
  if (pmin[2]==pmax[2]) tipo.set(m_planaxy);
  if (hayh) epsilon=hmin/1000; else epsilon_bb();

  // arregla los elementos
  for (i=0;i<c.e.len;i++) {
    elemento &ei=e[i]; nv=ei.nv(); d|=(1<<(ei.dim()));
    if (ei.f&fmask) hayfe=true;
    for (j=0;j<nv;j++) {
      in=mapnod[ei[j]];
      ei[j]=in;
      n[in].e+=i;
    }
  }
  if (d&2) tipo.set(m_lin);
  if (d&4) tipo.set(m_sup);
  if (d&8) tipo.set(m_vol);
}

// Separa componentes conexas
bool malla::piezas_conexas(){
  int qn=n.len-nodosh,qe=e.len,in,ie,lastin=0,pieza=0,qetotal=0,i,j,d;
  int *mapn=new int[qn],*mape=new int[qe];
  pline procesar(qn);
  while (qetotal<qe){
    // busca el 1er nodo (sin flag1 y con elementos)
    for (in=lastin;in<qn;in++) {if (n[in].f.noes(flag1)&&n[in].e) break;}
    if (in==qn) break; // puede habar sin elementos
    malla mt; d=0;
    lastin=in+1; procesar.clean(); procesar+=in;
    while (procesar){
      in=procesar.last(); procesar.len--;
      const cpline nie=n[in].e;
      for (i=0;i<nie.len;i++){
        ie=nie[i];
        if (e[ie].f.noes(flag1)) {
          e[ie].f.set(flag1);
          elemento ei=e[ie]; d|=(1<<(ei.dim()));
          for (j=0;j<ei.nv();j++) {
            in=ei[j];
            if (n[in].f.noes(flag1)){
              procesar+=in;
              n[in].f.set(flag1);
              mapn[in]=mt.n+=n[in];
            }
            ei[j]=mapn[in];
          }
          mape[ie]=mt.e+=ei;
        }        
        mt.n[mapn[in]].e[i]=mape[ie];
      }
    }
    if (d&2) mt.tipo.set(m_lin);
    if (d&4) mt.tipo.set(m_sup);
    if (d&8) mt.tipo.set(m_vol);
    sprintf(mt.nombre,"%s_part_%i",nombre,pieza);
    qetotal+=mt.e.len;
    mt.graba();
    pieza++;
  }
  delete [] mape; delete [] mapn; 
  for (i=0;i<qn;i++) n[i].f.reset(flag1);
  for (i=0;i<qe;i++) e[i].f.reset(flag1);
  return true;
}
