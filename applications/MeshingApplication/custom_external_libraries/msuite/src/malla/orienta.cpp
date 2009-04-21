// Las normales deben tener modulo unitario asi no hay que renormalizar cuando se usan
// Pueden ser por nodo o por elemento, cualquiera va en dir[] pero se indica
//   cual es mediante un bool: uno de ndir o edir debe estar en true y el otro en false
// Se asume que hay tantas normales como nodos o elementos
// Los elementos "bien orientados" son:
//    3d: Volumen positivo con los elementos numerados con el estandar de FEM
//    2d: Cerrado: En cada frontera giran a la izquierda vistos desde fuera
//        Malla abierta: todos como el primero que encuentra de cada pieza conexa
//    1d: Cerrado: Giro a la izquierda la exterior, a la derecha los huecos
//        Malla abierta: todos como el primero que encuentra de cada pieza conexa

// Reconoce como frontera exterior al nodo con maximo x. Si hay mas piezas conexas
//   las considera huecos interiores => NO ADMITE MAS DE UNA FRONTERA EXTERNA

// NO ADMITE FRONTERAS DISCONEXAS en mallas 2D o 1D (en 3D solo hace volumenes positivos)

//#define temporizar
#include "tiempo.h"
#include "cone.h" // dirnod
#include "malla.h"

using namespace std;

//===========
// defines
static const int flag12=flag1|flag2;
#define _repone {\
  for (in=0;in<qn;in++) n[in].f.reset(flag12);\
  if (tipo.noes(m_lin)) for (ie=0;ie<e.len;ie++) e[ie].f.reset(flag1);\
}

#define _NONMANIFOLD(i) {\
  add_error(NonManifold);\
  add_error(" at node: ");\
  add_error(i); add_error(": ");\
  add_error(n[i]);\
  _repone;\
  _savetime(error);\
  return false;\
}

#define _MOPEN(i) {\
  add_error(Open_Mesh);\
  add_error(i); add_error(": "); \
  add_error(n[i]);\
  _repone;\
  _savetime(error);\
  return false;\
}

#define _ENULL(i) {\
  add_error(Null_Element);\
  add_error(i);\
  hayerror=true;\
  break;\
}
// defines
//===========


void malla::rm_dir(){
  dir.ini(); 
  ndir=edir=false;
}


//=============================================================================================
// ELEMENTOS
//=============================================================================================

// orienta elementos y calcula normales por elemento

// implementacion especial para cerradas por el generador
bool malla::orienta(bool remake,bool mkdir){
  if (!e) return false; // no hay elementos
  // si ya esta orientada y 
  //   no se quiere dir o ya fue hecho por elemento y
  //     no se pide que se rehaga
  //       ==> sale
  if (mkdir&&tipo.es(m_vol)) mkdir=false;
  if (tipo.es(m_orientada)&&!remake){
    if (!mkdir||(dir&&edir)) return true;
    return mk_dir_elm(false);  // solo normales
  }
  if (!mk_vecino(remake)) return false;
  rm_frontera(); // cambia la orientacion de la frontera
  if (tipo.es(m_cerrada))
    return orienta_cerrada(true,remake,0,mkdir); // hacia fuera
  return orienta_abierta(remake,mkdir,0);
}

// orienta de a una pieza conexa abierta
// se supone hecha la lista de vecinos
// bits:
//   flag1 ordenado en nodos y elementos
//   flag2 nodo ordenable
// primera vez => marca nodos sueltos
// empieza a buscar vecinos por el nodo in

bool malla::orienta_abierta
    (bool remake, bool mkdir, int in){
  if (mkdir&&tipo.es(m_vol)) mkdir=false;
  if (tipo.es(m_orientada)&&(!mkdir||dir)&&!remake) return true;

  tipo.reset(m_orientada);

  //no orienta mezcla de tipos
  flagtype fdim=tipo&(m_vol|m_sup|m_lin);
  if (fdim!=m_vol&&fdim!=m_sup&&fdim!=m_lin) {
    add_error(Mixed);
    return false;
  }

  if (mkdir){
    if (dir.len!=e.len) {
      dir.clean(); dir.fit(e.len); dir.len=e.len; 
    }
    edir=true; ndir=false;
  }

  _initime;

  int ie,qn=n.len-nodosh;

  // primera vez ==> marca nodos sueltos
  if (in==0){for (;in<qn;in++) if (!n[in].e) n[in].f.set(flag1);in=0;};

  if (tipo.es(m_lin)){ // lineas
    int il,np;

    // ordena la arista posterior de cada nodo y eldenod
    // primer nodo: nodo con un solo elemento
    for (;in<qn;in++) if (n[in].f.noes(flag1)&&n[in].e.len==1) break;
    if (in==qn){
      // puede que sea un cacho cerrado (ademas del abierto)
      // porque si vino aca con in es que hay un nodo sin flag1
      bool retval=orienta_cerrada();
      tipo.reset(m_cerrada); // orienta_cerrada lo setea
      if (retval) tipo.set(m_orientada); else tipo.reset(m_orientada);
      _savetime(orienta); 
      return retval;
    }
    nodo &n0=n[in];  cpline &e0=n0.e;
    // ambas partiendo del nodo (puede revertir la orientacion original)
    elemento &lp=e[il=e0[0]]; np=lp[1];
    if (np==in) {swap(il); np=lp[1];}
    n0.f.set(flag1);
    if (mkdir)
      dir[il]=((n[np]-n0).giro90()).dir();
    in=np;
    while (n[in].e.len>1){
      cpline &ln=n[in].e;
      if (ln.len!=2) // 3 o mas
        _NONMANIFOLD(in);
      if (ln[1]==il) Swap(ln[0],ln[1]);
      il=ln[1]; elemento &lp=e[il]; np=lp[1];
      if (np==in) {swap(il);np=lp[1];}
      n[in].f.set(flag1);
      if (mkdir)
        dir[il]=((n[np]-n[in]).giro90()).dir();
      in=np;
    }
    n[in].f.set(flag1);
  }
  else if (tipo.es(m_sup)){ // superficie

    int i,v,ix,ipuesto,nv,in1;
    int j;
    pline nordenable(100);   

    for (;in<qn;in++)
      if (n[in].f.noes(flag1)) break;
    if (in==qn){_repone;_savetime(orienta);return true;}
    i=n[in].e[0]; elemento &e0=e[i]; e0.f.set(flag1); nv=e0.nv();
    // el primero queda como estaba (por si es un cacho de malla orientada)
    if (mkdir) dir[i]=((n[e0[1]]-n[e0[0]])%(n[e0[2]]-n[e0[0]])).dir();
    for (j=0;j<nv;j++) {n[e0[j]].f.set(flag2); nordenable+=e0[j];}
    while (nordenable.len){
      in=nordenable.last(); nordenable.len--; n[in].f.set(flag1);
      cpline &eni=n[in].e;
      if (eni.len==1) continue; // vertice
      // busca un elemento ordenado
      for (i=0;i<eni.len;i++) if (e[eni[i]].f.es(flag1)) break;
      if (i) Swap(eni[i],eni[0]); // lo pone al principio
      ipuesto=0; ie=eni[0]; elemento &e1=e[ie]; nv=e1.nv();
      // busca el nodo en el elemento y el vecino anterior
      ix=e1.index(in); v=vecino[ie][(ix+nv-1)%nv];
      // ordena el resto
      while(ipuesto<eni.len-1){
        ipuesto++;
        if (v<0) {
          for (i=ipuesto;i<eni.len;i++){
            v=eni[i]; elemento &ev=e[v];
            ix=ev.index(in); nv=ev.nv();
            if (vecino[v][ix]<0) break;
            if (vecino[v][(ix+nv-1)%nv]<0) {swap(v); ix=ev.index(in); break;}
          }
          if (i==eni.len)
            _NONMANIFOLD(in);
          if (i!=ipuesto) Swap(eni[i],eni[ipuesto]);
        }
        else {
          elemento &ev=e[v]; ix=ev.index(in); nv=ev.nv();
          if (vecino[v][ix]!=ie) {swap(v);ix=ev.index(in);}
          i=eni.index(v,ipuesto);
          if (i!=ipuesto) Swap(eni[i],eni[ipuesto]);
        }
        ie=v; elemento &ei=e[ie]; v=vecino[ie][(ix+nv-1)%nv];
        if (ei.f.es(flag1)) continue;
        ei.f.set(flag1);
        // agrega a nordenable
        for (j=0;j<nv;j++) {
          in1=ei[j]; if (n[in1].f.noes(flag2))
            {n[in1].f.set(flag2); nordenable+=in1;}
        }
        if (mkdir){
          // punto medio y puntos medios de las dos mitades (aprox)
          j=nv/2; // 3=01x02 4=02x13 5=02x13 6=03x14 .....
          dir[ie]=((n[ei[j]]-n[ei[0]])%(n[ei[(nv+j)/2]]-n[ei[j/2]])).dir();
        }
      }
    }
  }
  else if (tipo.es(m_vol)){
    int i;
    if (!ve) volumen();
    for (i=0;i<e.len;i++)
      if (ve[i]<0) swap(i);
  }

  if (tipo.noes(m_vol)){
    // verifica si restan nodos
    for (in=0;in<qn;in++) {if (n[in].f.noes(flag1)) break;}
    if (in<qn) if (!orienta_abierta(true,mkdir,in)) return false;
    // listo
    _repone;
  }

  tipo.set(m_orientada);
  _savetime(orienta);
  return true;
}

// normaliza
// ordena tambien los eldenod
// ojo con piezas disconexas (varias fronteras exteriores)
// habria que testear si las piezas disconexas son interiores
// bits:
// flag1 ordenado en nodos y elementos
// flag2 nodo ordenable
// dirext=true significa giro a la izq visto desde afuera/arriba
// se presupone que la malla es cerrada
// y no se supone que haya vecinos (si viene del generador)

bool malla::orienta_cerrada
    (bool dirext,bool remake,int* numf, bool mkdir, bool reentry){
  if (tipo.es(m_orientada)&&(!mkdir||dir)&&!remake) return true;

  tipo.reset(m_orientada);

  if (mkdir){
    if (dir.len!=e.len) {
      dir.clean(); dir.fit(e.len); dir.len=e.len;
    }
    edir=true; ndir=false;
  }
  _initime;

  int in,ie,qn=n.len-nodosh;
  static int numero; static bool dir_ori;
  if (!reentry) {numero=0; dir_ori=dirext;} else numero++;

  // busca max x
  int imax=-1; double xmax=-MAXREAL;
  for (in=0;in<qn;in++) // el primero
    if (n[in].f.noes(flag1)) {xmax=n[in][0]; imax=in; break;}
  for (in++;in<qn;in++) {
    const nodo &ni=n[in];
    if (ni.f.noes(flag1)&&ni[0]>xmax) {xmax=ni[0]; imax=in;}
  }

  if (tipo.es(m_lin)){ // lineas
    int il,na,np;

    // ordena la arista posterior de cada nodo y eldenod
    // primer nodo
    nodo &nx=n[imax];  cpline &ex=nx.e;  if (ex.len!=2) _MOPEN(imax);
    if (numf) numf[imax]=numero;
    // ambas partiendo del nodo
    elemento &la=e[ex[0]]; na=la[1]; if (na==imax) {swap(ex[0]); na=la[1];}
    elemento &lp=e[ex[1]]; np=lp[1]; if (np==imax) {swap(ex[1]); np=lp[1];}
    // posterior % anterior debe dar hacia adentro o afuera segun dir
    bool izq;
    if ((fabs(n[na][0]-xmax)/xmax<.01)&&(fabs(n[np][0]-xmax)/xmax<.01))
      izq=(n[na][1]<n[np][1]); // puntos colineales
    else
      izq=((n[np]-nx).pv2d(n[na]-nx)>0);

    if ((!dirext&&izq)||(dirext&&!izq)) Swap(ex[0],ex[1]);
    swap(ex[0]);
    n[imax].f.set(flag1);
    if (mkdir)
      dir[ex[1]]=(n[imax]-n[e[ex[1]][1]]).giro90().dir(); // hacia afuera

    il=ex[1]; in=e[il][1];
    while (n[in].f.noes(flag1)){
      if (numf) numf[in]=numero;
      cpline &ln=n[in].e;
      if (ln.len!=2)
        _MOPEN(in);
      if (ln[1]==il) Swap(ln[0],ln[1]); // elementos del nodo
      il=ln[1];
      elemento &l=e[il];
      if (l[1]==in) swap(il); // el elemento siguiente
      n[in].f.set(flag1);
      if (mkdir)
        dir[il]=(n[in]-n[l[1]]).giro90().dir(); // hacia afuera
      in=l[1];
    }
  }
  else { // superficie

    int i,j,ix,na,np,in1,ipuesto,nv;
    pline nordenable(100);

    // ordena el elemento mas perpendicular al eje x
    // puede que alguna cara este en un plano paralelo al eje x o mas o menos por perturbacion
    cpline &enx=n[imax].e; if (enx.len<2) _MOPEN(imax);
    int iemax=-1; xmax=-1; double vx,vxmax=-MAXREAL;
    for (i=0;i<enx.len;i++){
      ie=enx[i]; elemento& e0=e[ie]; ix=e0.index(imax);
      // direccion del producto vectorial de la arista anterior por la posterior
      na=e0.nant(ix); np=e0.npos(ix);
      vx=((n[np]-n[imax])%(n[na]-n[imax])).dir()[0];
      if (fabs(vx)>xmax) {xmax=fabs(vx); vxmax=vx; iemax=ie;}
    }
    e[iemax].f.set(flag1); elemento &emax=e[iemax]; nv=emax.nv();
    // vxmax debe ser <0 si debe apuntar para adentro
    if ((!dirext&&vxmax>0)||(dirext&&vxmax<0)) swap(iemax); // invierte el orden
    if (mkdir)
      dir[iemax]=((n[emax[1]]-n[emax[0]])%(n[emax[2]]-n[emax[0]])).dir();

    // ya tenemos un elemento ordenado => todos sus nodos son ordenables
    for (j=0;j<nv;j++) {n[emax[j]].f.set(flag2); nordenable+=emax[j];}

    while (nordenable.len){
      in=nordenable.last(); nordenable.len--; n[in].f.set(flag1);
      if (numf) numf[in]=numero;
      cpline &eni=n[in].e; if (eni.len<2)
        _MOPEN(in);
      // busca un elemento ordenado
      for (i=0;i<eni.len;i++) if (e[eni[i]].f.es(flag1)) break;
      if (i) Swap(eni[i],eni[0]); // lo pone al principio
      ipuesto=0; ie=eni[0]; elemento &e1=e[ie];
      // busca el nodo en el elemento
      ix=e1.index(in); na=e1.nant(ix);
      // ordena el resto
      while(ipuesto<eni.len-1){
        // busca un elemento con el nodo anterior
        for (i=ipuesto+1;i<eni.len;i++) { // ojo, cuadratico eni.len
          ix=e[eni[i]].index(na);
          if (ix>=0) break;
        }
        if (i==eni.len)
          _MOPEN(in);
        ie=eni[i]; elemento &e2=e[ie];
        // coloca este siguiente en la lista
        ipuesto++; if (i!=ipuesto) Swap(eni[i],eni[ipuesto]);
        // ordena y busca el nodo posterior
        if (e[ie].f.es(flag1)) {na=e2.nant(ix,2); continue;} // para el proximo
        e[ie].f.set(flag1); nv=e2.nv();
        // agrega a nordenable
        for (j=0;j<nv;j++) {
          in1=e2[j];
          if (n[in1].f.noes(flag2))
            {n[in1].f.set(flag2); nordenable+=in1;}
          if (in1==in) ix=j;
        }
        // orienta
        if (e2.nant(ix)==na) {
          // desordenado
          na=e2.npos(ix); // para el proximo
          swap(ie); // invierte el orden
        }
        else na=e2.nant(ix); // para el proximo
        if (mkdir){
          // punto medio y puntos medios de las dos mitades (aprox)
          j=nv/2; // 3=01x02 4=02x13 5=02x13 6=03x14 .....
          dir[ie]=((n[e2[j]]-n[e2[0]])%(n[e2[(nv+j)/2]]-n[e2[j/2]])).dir();
        }
      }
    }
  }

  // verifica si restan piezas conexas
  for (in=0;in<qn;in++) {if (n[in].f.noes(flag1)) break;}
  if (in<qn) // las piezas interiores al reves
    if (!orienta_cerrada(!dir_ori,true,numf,mkdir,true)) return false;
  // listo, limpia
  if (reentry) return true;
  _repone;
  tipo.set(m_orientada|m_cerrada);
  _savetime(orienta);
  return true;
}

bool malla::mk_dir_elm(bool remake){
  if (tipo.es(m_vol)) return false;
  if (dir&&edir&&!remake) return true;
  if (tipo.noes(m_orientada)) return orienta(true,true);
  int i,elen=e.len;
  if (dir.len!=elen) {
    dir.clean(); dir.fit(elen); dir.len=elen;
  }
  edir=true; ndir=false;
  if (tipo.es(m_lin)) {
    for (i=0;i<elen;i++) dir[i]=(n[e[i][1]]-n[e[i][0]]).giro90().dir();
    return true;
  }
  int j,nv;
  for (i=0;i<elen;i++){
    nv=e[i].nv(); // punto medio y puntos medios de las dos mitades (aprox)
    j=nv/2; // 3=01x02 4=02x13 5=02x13 6=03x14 .....
    dir[i]=((n[e[i][j]]-n[e[i][0]])%(n[e[i][(nv+j)/2]]-n[e[i][j/2]])).dir();
  }
  return true;
}

//=============================================================================================
// NODOS
//=============================================================================================

// normales por nodo
// mata las componentes normales a los elementos marcados sin e_offset
bool malla::mk_dir_nod(bool remake){  
  if (!remake && dir && ndir && dir.len==n.len-nodosh) return true;
  if (!orienta(remake)) return false;

  _initime;
  int i,nlen=n.len-nodosh; // los nodos de h no tienen dir
  dir.clean(); dir.fit(nlen); dir.len=nlen; ndir=true; edir=false;
  bool hayerror=false;

  if (e[0].dim()==1){// lineas
    punto d0,d1; double l;
    for (i=0;i<nlen;i++) {
      nodo &ni=n[i]; const cpline &eni=ni.e;
      if (eni.len==2) { // nodo interior
        if (e[eni[0]].f.es(e_offset)) {
          d0=n[e[eni[0]][0]]-ni; l=d0.mod(); if (l<ERRADM) _ENULL(eni[0]);
          d0/=l;
        } else d0.zero();
        if (e[eni[1]].f.es(e_offset)) {          
          d1=ni-n[e[eni[1]][1]]; l=d1.mod(); if (l<ERRADM) _ENULL(eni[1]);
          d1/=l;
        } else d1.zero();
        dir[i]=(d0+d1).giro90();
        l=dir[i].mod(); if (l>ERRADM) dir[i]/=l; else dir[i].zero();
      }
      else if (eni.len==1) { // nodo terminal
        if (e[eni[0]].f.noes(e_offset)) continue; // normal nula
        int n0=e[eni[0]][0], n1=e[eni[0]][1];
        d1=n[n0]-n[n1]; l=d1.mod(); if (l<ERRADM) _ENULL(eni[0]);
        if (n1!=i) l=-l; dir[i]=(d1/l).giro90(); // supuesto que i es el primer nodo
      }
      // else la normal queda nula
    }
  }
  else if (e[0].dim()==2){// superficie
    // minimiza el maximo angulo con las normales
    bool haysim;
    double m,cc;
    int j,ix,elen,maxelen=32,realelen;
    punto ai,ad,dir1,dirsim,mdir,vc,*adir=new punto[maxelen],*edir=new punto[maxelen],d1;
    for (i=0;i<nlen;i++){
      nodo &ni=n[i];
      const cpline &eni=ni.e; elen=eni.len;
      if (elen>maxelen) { // mantenimiento de arrays
        maxelen=int(1.5*elen); 
        delete [] edir; edir=new punto [maxelen];
        delete [] adir; adir=new punto [maxelen];
      }
      // aristas
      for (j=0;j<elen;j++) {
        const elemento &ei=e[eni[j]]; ix=ei.index(i);
        dir1=n[e[eni[j]].npos(ix)]-ni;
        m=dir1.mod(); if (m<ERRADM) _ENULL(eni[j]); adir[j]=dir1/=m;
      } 
      if (hayerror) break; //ENULL
      // normales por cara y media (NOponderada por angulos)
      realelen=0; mdir=pzero; haysim=false;
      for (j=0;j<elen;j++) {
        dir1=adir[j]%adir[(j+1)%elen];
        m=dir1.mod(); if (m<ERRADM) _ENULL(eni[j]); // seno
        if (e[eni[j]].f.noes(e_offset)){
          if (!haysim) {haysim=true; dirsim=dir1/=m;}
        }
        else {
          mdir+=edir[realelen++]=dir1/=m;
//          edir[realelen++]=dir1/=m;
//          mdir+=dir1*acos(adir[j]*adir[(j+1)%elen]);        
        }
      } 
      if (hayerror) break; //ENULL
      if (!realelen) continue;
      if (realelen==1) {vc=edir[0]; cc=1;}
      else {
        mdir.dir();
        if (!cone(edir, realelen, vc, cc)) vc=mdir;
      }
      if (haysim) vc=(vc-dirsim*(vc*dirsim)).dir(); // mata la componente normal
      if (cc<=0){ // vc no sirve
        // puede servir para offset pero no para blayer
        add_warning(Void_Kernel); add_warning(i);
      }
      dir[i]=vc;
    }
    delete [] edir;
    delete [] adir;
  }
  else hayerror=true;
  
  if (hayerror){
    rm_dir();
    _savetime(error);
    return false;
  }

  return true;
}

// offset hacia afuera (delta +)
bool malla::offset(double delta){
  if (fabs(delta)/10<ERRADM) return false;
  if (fabs(delta)/10<epsilon) epsilon=fabs(delta)/10;
  int i,nlen=n.len-nodosh;

  if (!ndir){
    // si no estan marcados marco todos
    bool deshacer=false;
    for(i=0;i<e.len;i++) {
      if (e[i].f.es(e_offset)) {deshacer=true; break;} // oops ya habia marcas!
      e[i].f.set(e_offset);
    }
    if (deshacer) for(i--;i>=0;i--) e[i].f.reset(e_offset);
    if (!mk_dir_nod()) return false; // calcula normales por nodo
  }
  if (o) delete o; o=0; ve.ini(); re.ini(); ce.ini();

  // offestea y recalcula el bbox;
  pmin=pmax=n[0]+=dir[0]*delta;
  for (i=0;i<nlen;i++) {
    n[i]+=dir[i]*delta;
    n[i].set_min_max(pmin,pmax);
  }

  return true;
}

// Hace una malla extrudando la superficie "hacia adentro" 
//   (delta + => adentro => opuesto al estandar 2D de elemento orientado hacia afuera)
// Los nodos que se extrudan vienen con n_offset
// Los delta son distancias entre capas
// Los nodos originales quedan donde estaban, 
// Los elementos atachados a estos se comprimen hacia dentro sin suavizado
//   y cambian de nodos, pasan a tener los nodos del limite interior de la blayer
// Los elementos extrudados contra el plano de simetria dan rectangulitos
//   si la deteccion es por nodos se generan falsos internos que pueden estar repetidos 
//   o no dependiendo de cuantos elementos del plano de simetria vengan dados ?????
bool malla::blayer(
  int capas, const double *dcapa, // delta acumulado entre capas
  bool archiva_normales, // graba normales como malla de segmentos 
  array1<punto>* normal, pline* mapnod // array de normales en cada nodo indicado por map  
  ){
  int i,j,k,in,io,nv;
  bool eoffset,noffset,haynoffset,hayrec,expande;

  // verifica deltas
  for(i=0,expande=false;i<capas;i++) {
    if (dcapa[i]<0) expande=true; //hay que recalcular bbox
    double delta=fabs(dcapa[i]-((i)? dcapa[i-1] : 0))/10;
    if (delta<ERRADM) return false;
    if (delta<epsilon) epsilon=delta;
  }
  
  // frontera
  orienta(); // la necesito orientada porque 2d puede invertir las lineas de frontera
  if (tipo.es(m_lin)||!mk_frontera()||!frontera||!frontera[0].e.len) return false;

  // Si no hay flag de offset => marca todos los nodos de frontera
  for (k=0,haynoffset=false;k<frontera.len;k++){
    // al buscar marco y si habia, desmarco los que marque aca
    const pline &nf=frontera[k].n;
    for (i=0;i<nf.len;i++) {
      if (n[nf[i]].f.es(n_offset)) {haynoffset=true;break;} 
      else n[nf[i]].f.set(n_offset);
    }
    if (haynoffset){ // desmarca los marcados aca
      for (i--;i>=0;i--) n[nf[i]].f.reset(n_offset); // esta frontera
      for (k--;k>=0;k--){ // las anteriores
        const pline &nfant=frontera[k].n;
        for (i=0;i<nfant.len;i++) n[nfant[i]].f.reset(n_offset);
      }
      break;
    }
  }

  // malla de elementos frontera que tienen al menos un nodo a offsetear
  // los que tienen todos, son elementos a offsetear, el resto define plano de simetria
  malla mo; 
  if (tipo.es(m_vol)) mo.tipo=m_sup;
  else if (tipo.es(m_sup)) mo.tipo=m_lin;
  pline map(frontera[0].n.len); pline pam(n.len,-1);
  for (k=0,hayrec=false;k<frontera.len;k++){
    array1<elemento> efk=frontera[k].e;
    for (i=0;i<efk.len;i++){
      elemento ei=efk[i]; nv=ei.nv(); // copia
      // eoffset=todos offset;  noffset=alguno offset
      eoffset=true; noffset=false;
      for (j=0;j<nv;j++) if (n[ei[j]].f.es(n_offset)) noffset=true; else eoffset=false;
      if (!noffset) continue; // ninguno
      for (j=0;j<nv;j++){
        in=ei[j]; io=pam[in]; 
        nodo nj=n[in]; // copia
          if (io==-1) { // primera vez que aparece
          map+=in;
          nj.e.clean(); io=pam[in]=mo.n+=nj;
        }
        ei[j]=io;
        mo.n[io].e+=mo.e.len;
      }
      if (eoffset) ei.f.set(e_offset); else hayrec=true;
      mo.e+=ei;
    }
  }
  if (!(mo.n&&mo.e)) return false;

  // si es 1d, orienta() puede revertir el orden
  int n00_ori=mo.e[0][0];
  if (!mo.orienta()) return false;
  if (mo.e[0][0]!=n00_ori) swap();

  if(!mo.mk_dir_nod()) return false;

  // guarda las normales (en la dir de offset) si se piden
  // las longitudes no son 1 sino el maximo offset
  // mo contiene nodos de simetria ademas de offsets
  if (archiva_normales) {
    malla mn; mn.n.resize(2*mo.n.len); mn.e.resize(mo.n.len);
    mn.hayh=false; mn.hayv=true; mn.hayfn=false; mn.hayfe=false;
    elemento en(e_segmento);
    for (i=0;i<mo.n.len;i++){
      if (mo.n[i].f.noes(n_offset)) continue;
      in=map[i]; // el original
      nodo newn=mo.n[i]; newn.f.reset(); newn.v=in; 
      newn.e.clean(); newn.e+=mn.e.len;
      en[0]=mn.n.len; en[1]=en[0]+1; mn.e+=en;
      mn.n+=newn; mn.n+=newn-=dcapa[capas-1]*mo.dir[i];
    }
    mn.copy_nombre(*this); strcat(mn.nombre,"_n");
    mn.graba_dat();
  }
  // arma el array de normales (si se pide)
  if (normal && mapnod){
    normal->resize(mo.n.len); mapnod->resize(mo.n.len);
    for (i=0;i<mo.n.len;i++){
      in=map[i]; // el original
      if (n[in].f.noes(n_offset)) continue;
      (*normal)+=mo.dir[i]; (*mapnod)+=in;
    }
  }

//  if (mo.m_warning&&strstr(mo.m_warning,Void_Kernel)) return false; // normal imposible

  // propaga los nodos 
  nodo newn; newn.f.set(n_offset);
  int oldnlen=n.len; // para rehacer el bbox
//  double r/*,minr=1*/; ///?????????????
  pline newmap(map.len,-1); // nodos nuevos
  for (i=0;i<mo.n.len;i++) {
    in=map[i]; nodo &ni=n[in]; // el original
    if (ni.f.noes(n_offset)) continue;
//    r=1; if (hayv&&ni.v!=0&&dcapa[capas-1]>ni.v/2) { ///?????????????
//      r=ni.v/deltatotal[capas-1]/2; r*=(2-r);
//    }
//newn.v=ni.v;
    newmap[i]=n.len; // primer nodo agregado de la serie
    if (ni.f.es(n_simetria)) newn.f.set(n_simetria); else newn.f.reset(n_simetria);
    for (j=0;j<capas;j++) n+=newn.setpos(ni-mo.dir[i]*dcapa[j]/**r*/); // offsets
    // los elementos que estaban con in ahora van con el ultimo
    for (j=0;j<ni.e.len;j++) e[ni.e[j]].replace(in,n.len-1);
    n.last().e.roba(ni.e);
//set_min(minr,r);
  }

  if (dcapa[capas-1]<0) {// hacia afuera => bbox
    for (i=oldnlen+capas;i<n.len;i+=capas) n[i].set_min_max(pmin,pmax);
  }
  delete o; o=0; // hay mas nodos
  rm_nn(); rm_vecino(); rm_esferas(); ve.ini();
  if (hayrec) rm_frontera(); // ya no es valida por los rectangulos
 
  // elementos
  elemento s(e_segmento),q(e_cuadrilatero), c(e_cubo), w(e_wedge);
  int elen=mo.e.len;
  for (i=0;i<elen;i++) {
    elemento &ei=mo.e[i]; nv=ei.nv();
    if (ei.f.noes(e_offset)) continue; // no es offset 
    // agrega prismas
    if (nv==2){ // segmento -> cuadrilatero (orden circular)
      // el exterior
      q[0]=map[ei[0]]; q[1]=map[ei[1]];
      q[2]=newmap[ei[1]]; q[3]=newmap[ei[0]];
      n[q[0]].e+=e.len; n[q[1]].e+=e.len; n[q[2]].e+=e.len; n[q[3]].e+=e.len; 
      e+=q;
      for (j=1; j<capas; j++){ // interiores
        q[0]=q[3]; q[1]=q[2]; q[2]++; q[3]++; 
        n[q[0]].e+=e.len; n[q[1]].e+=e.len; n[q[2]].e+=e.len; n[q[3]].e+=e.len; 
        e+=q;
      }
    }
    else if (nv==3){ // triangulo -> wedge
      // el exterior
      w[0]=map[ei[0]]; w[1]=map[ei[2]]; w[2]=map[ei[1]]; // triangulos hacia afuera
      w[3]=newmap[ei[0]]; w[4]=newmap[ei[2]]; w[5]=newmap[ei[1]]; 
      n[w[0]].e+=e.len; n[w[1]].e+=e.len; n[w[2]].e+=e.len;
      n[w[3]].e+=e.len; n[w[4]].e+=e.len; n[w[5]].e+=e.len;
      e+=w;
      for (j=1; j<capas; j++){ // interiores
        w[0]=w[3]; w[1]=w[4]; w[2]=w[5]; 
        w[3]++; w[4]++; w[5]++;
        n[w[0]].e+=e.len; n[w[1]].e+=e.len; n[w[2]].e+=e.len;
        n[w[3]].e+=e.len; n[w[4]].e+=e.len; n[w[5]].e+=e.len;
        e+=w;
      }
    }
    else if (nv==4){ // cuadrilatero -> cubo
      // exterior
      c[0]=map[ei[1]]; c[1]=map[ei[0]]; c[2]=map[ei[3]]; c[3]=map[ei[2]]; // quads hacia afuera
      c[4]=newmap[ei[1]]; c[5]=newmap[ei[0]]; c[6]=newmap[ei[3]]; c[7]=newmap[ei[2]]; 
      n[c[0]].e+=e.len; n[c[1]].e+=e.len; n[c[2]].e+=e.len; n[c[3]].e+=e.len;
      n[c[4]].e+=e.len; n[c[5]].e+=e.len; n[c[6]].e+=e.len; n[c[7]].e+=e.len;
      e+=c;      
      for (j=1; j<capas; j++){// interiores
        c[0]=c[4]; c[1]=c[5]; c[2]=c[6]; c[3]=c[7];
        c[4]++; c[5]++; c[6]++; c[7]++;
        n[c[0]].e+=e.len; n[c[1]].e+=e.len; n[c[2]].e+=e.len; n[c[3]].e+=e.len;
        n[c[4]].e+=e.len; n[c[5]].e+=e.len; n[c[6]].e+=e.len; n[c[7]].e+=e.len;
        e+=c;
      }
    }
    else{
_revienta(true);
    }
  }

  return true;
}

// proyecta sobre un plano normal a norm
// si norm es 0,0,0 proyecta sobre la normal media
// Puede quedar orientada al reves o aun pueden solaparse partes (depende de norm)
void malla::proyecta_xy(punto norm) {
	int i;
	if (norm==pzero) {
    orienta(false,true);
		for (i=0;i<e.len;i++) norm+=dir[i];
	}
  
	norm.dir(); // normaliza la normal
	punto vx=norm%_ez; 
  double m=vx.mod(); if (m>ERRADM) vx/=m; else vx=(_ey%norm).dir();
  punto vy=norm%vx;
  
	for (i=0;i<n.len;i++) {
		punto aux=n[i];
		n[i][0]=aux*vx;
		n[i][1]=aux*vy;
		n[i][2]=0;
	}
	
  tipo.set(m_planaxy);
  rm_dir(); rm_esferas(); ve.ini();
	if (o) { delete o; o=NULL; }
	bbox(true);
}

#undef _MOPEN
#undef _ENULL
#undef _repone
#undef _NONMANIFOLD
#undef sueltos
