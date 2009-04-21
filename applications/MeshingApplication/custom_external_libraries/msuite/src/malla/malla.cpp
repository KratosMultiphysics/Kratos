// operaciones generales de la clase malla

#include <cstring> // strcpy/cat...
#include <cstdio> // sprintf
//#define temporizar
#include "tiempo.h"
#include "malla.h"
#include "mensajes.h" // mensajes
#include "filename.h"
#include "case.h"

using namespace std;

bool malla::INFO_CL=false;

malla::malla():_INI_MALLA {nombre[0]=ext[0]=0;}

malla::malla(const malla &m1):_INI_MALLA {operator=(m1);}

// listas simples
// conectividades base 1
// ojo: ndim=3 no indica 3d sino que hay z (puede ser 0)
malla::malla(
  int ndim, int nlen, double* nl, // cantidad y lista de puntos (x0y0[z0]x1y1[z1].....)
  double* hl, // h de cada nodo
  int edim, int elen, int* el   // cantidad y lista de elementos
  ):_INI_MALLA
{
  nombre[0]=ext[0]=0;

  int i,j;
  nodo p;

  // nodos, bounding box y h minimo
  n.resize(nlen);
  p[0]=nl[0]; p[1]=nl[1]; if (ndim==3) p[2]=nl[2];
  n+=p;
  pmin=pmax=p;
  hayh=(hl!=0); if (hayh) hmin=hl[0];
  for (i=1,j=ndim;i<nlen;i++) {
    p[0]=nl[j++]; p[1]=nl[j++]; if (ndim==3) p[2]=nl[j++];
    p.set_min_max(pmin,pmax);
    if (hayh) {set_min(hmin,hl[i]); p.h=hl[i];}
    n+=p;
  }

  // hay z?
  bool es3D=pmax[2]-pmin[2]>ERRADM;
  if (!es3D) tipo.set(m_planaxy);

  // epsilon (para varios usos)
  if (hayh)
    epsilon=hmin/1000;
  else
    epsilon_bb();

  o=new octree(n,pmin,pmax,es3D ? 3 : 2);// octree
  for (i=0; i<nlen; i++) o->add(i);

  // elementos y elementos de cada nodo
  if (el) {
    elemento ei((edim==1)? e_segmento : ((edim==2)? e_triangulo : e_tetraedro));
    int nv=edim+1;
    e.resize(elen);
    for (i=0;i<elen;i++){
      for (j=0;j<nv;j++) {ei[j]=el[nv*i+j]-1; n[ei[j]].e+=i;}
      e+=ei; // presupongo que es e.len-1
    }
    tipo.set((edim==1)? m_lin : ((edim==2)? m_sup : m_vol));
  }
  else tipo.set(m_nodos);
}


// listas mas simples
// ojo:
// ndim=3 no indica 3d sino que hay z (puede ser 0)
// elementos base 1
malla::malla(
  int ndim, int nlen, double* nl, // cantidad y lista de puntos (x0y0[z0]x1y1[z1].....)
  int* el   // lista de elementos base 1 terminada en 0
  ):_INI_MALLA
{
  nombre[0]=ext[0]=0;

  int i,j;
  punto p;

  // nodos y bounding box
  n.resize(nlen);
  p[0]=nl[0]; p[1]=nl[1]; if (ndim==3) p[2]=nl[2];
  n+=p;
  pmin=pmax=p;
  hayh=false;
  for (i=1,j=ndim;i<nlen;i++) {
    p[0]=nl[j++]; p[1]=nl[j++]; if (ndim==3) p[2]=nl[j++];
    p.set_min_max(pmin,pmax);
    n+=p;
  }

  // hay z?
  bool es3D=pmax[2]-pmin[2]>ERRADM;
  if (!es3D) tipo.set(m_planaxy);

  // epsilon (para varios usos)
  epsilon_bb();

  o=new octree(n,pmin,pmax,es3D ? 3 : 2);// octree
  for (i=0; i<nlen; i++) o->add(i);

  // elementos y elementos de cada nodo
  if (el) {
    elemento ei((!es3D)? e_segmento : e_triangulo);
    int nv=ei.nv();
    e.resize((!es3D)? nlen : 7*nlen); // prealoca (no es necesario)
    for(i=0;el[nv*i];i++){
      for (j=0;j<nv;j++) {ei[j]=el[nv*i+j]-1; n[ei[j]].e+=i;}
      e+=ei; // presupongo que es e.len-1
    }
    tipo.set((!es3D)? m_lin : m_sup);
  }
  else tipo.set(m_nodos);
}

malla::~malla(){  // destructor
  _initime;
  delete [] m_error; m_error=0;
  delete [] m_warning; m_warning=0;
  delete o; o=0;
  _savetime(~malla);
}

#define _INI(array){\
  _initime;\
  array.ini();\
  _savetime(malla::ini array);\
}


void malla::ini(){// inicializacion
  _initime;

  nombre[0]=ext[0]=0;
  delete [] m_error; delete [] m_warning; m_error=m_warning=0;
  tipo=0;
  nodosh=0;
  pmin=pmax=pzero;

  _initime;
  delete o; o=0;
  _savetime(malla::ini octree);

  epsilon=ERRADM;
  hayh=false; hmin=MAXREAL;
  hayfn=false; hayfe=false;
  hayv=false; nvmin=MAXREAL; nvmax=-MAXREAL;
  _INI(frontera);
  _INI(dir); ndir=false; edir=false;
  _INI(e);
  _INI(ce);
  _INI(re);
  _INI(ve);
  _INI(vecino);
  _INI(n);
  _INI(nn);

  _savetime(malla::ini);
}

void malla::e_ini(){
  _initime;
  tipo.reset(m_vol|m_sup|m_lin|m_cerrada|m_orientada|m_slivers_marked);
  tipo.set(m_nodos);
  if (n[0].e||n[n.len/2].e||n.last().e) { // para no hacerlo al cuete
    for (int i=0;i<n.len;i++) n[i].e.ini();
  }
  hayfe=false;
  frontera.ini();
  rm_dir(); e.ini(); ce.ini(); re.ini(); ve.ini(); vecino.ini();
  _savetime(malla::e_ini);
}

void malla::copy_nombre(const malla &m){
  strcpy(nombre,m.nombre);
  strcpy(ext,m.ext);
}

malla& malla::operator=(const malla& m){
  if (&m==this) return *this;
  _initime;

  if (o==m.o) o=0; // sin borrar
  ini();

  copy_nombre(m);
  tipo=m.tipo;

  n=m.n; nn=m.nn; nodosh=m.nodosh;
  pmin=m.pmin; pmax=m.pmax;
  epsilon=m.epsilon;
  hayh=m.hayh; hmin=m.hmin;
  hayv=m.hayv; nvmin=m.nvmin; nvmax=m.nvmax;
  hayfn=m.hayfn; hayfe=m.hayfe;
  e=m.e; vecino=m.vecino; ve=m.ve; re=m.re; ce=m.ce;
  frontera=m.frontera;
  dir=m.dir; ndir=m.ndir; edir=m.edir;
  // no copia o ni b ni error ni warning
  delete [] m_error; delete [] m_warning; m_error=m_warning=0;

  _savetime(malla=);
  return *this;
}

malla& malla::roba(malla& m){
  if (&m==this) return *this;

  _initime;

  if (o==m.o) o=0; // sin borrar
  ini();

  copy_nombre(m);
  tipo=m.tipo;
  m_error=m.m_error; m.m_error=0;
  m_warning=m.m_warning; m.m_warning=0;
  hayh=m.hayh; hmin=m.hmin;
  hayv=m.hayv; nvmin=m.nvmin; nvmax=m.nvmax;
  hayfn=m.hayfn; hayfe=m.hayfe;

  n.roba(m.n); nn.roba(m.nn); nodosh=m.nodosh;
  o=m.o; o->plist=&n; m.o=0;
  pmin=m.pmin; pmax=m.pmax;
  epsilon=m.epsilon;

  e.roba(m.e); vecino.roba(m.vecino); ve.roba(m.ve); re.roba(m.re); ce.roba(m.ce);
  frontera.roba(m.frontera);
  dir.roba(m.dir); ndir=m.ndir; edir=m.edir;

  m.ini();

  _savetime(malla_roba);
  return *this;
}

#define _mk_msg(type)\
void malla::ini_##type() {\
  if (!m_##type) return;\
  delete[] m_##type; m_##type=0;\
}\
void malla::add_##type(const char *msg){\
  if (INFO_CL) cout << msg;\
  size_t len=((m_##type)? strlen(m_##type) : 0) + strlen(msg) + 1;\
  char *nuevo=new char[len]; *nuevo=0;\
  if (m_##type) {strcpy(nuevo,m_##type);delete[] m_##type;}\
  strcat(nuevo,msg);\
  m_##type=nuevo;\
}\
void malla::add_##type(double d){\
  char msg[32]; sprintf(msg,"%.4g",d);\
  add_##type(msg);\
}\
void malla::add_##type(int i){\
  char msg[32]; sprintf(msg,"%i",i);\
  add_##type(msg);\
}\
void malla::add_##type(punto p){\
  char msg[32]; sprintf(msg,"%.4g,%.4g,%.4g",p[0],p[1],p[2]);\
  add_##type(msg);\
}

_mk_msg(error)
_mk_msg(warning)
#undef _mk_msg

void malla::renombra(const char* inombre){
  // nombre y extension
  if (inombre[0]) strcpy(nombre,inombre);
  else return;
  char* extaddr=ext_begin(nombre);
  if (extaddr) {
    *(extaddr-1)=0;
    strcpy(ext,extaddr);
  }
  if (!ext[0]) strcpy(ext,"dat");
}

#define _graba(_ext) \
  else if (!strcmp_nocase(ext,#_ext)) {\
    if (!graba_##_ext(filename)) {_savetime(error);return false;}\
  }

bool malla::graba(const char* inombre) {
  _initime;

  // si la malla tiene nodos de h pero el resto de los nodos no tiene h es un problema
  if (nodosh){
    if (!hayh)
      mk_h_nn_min();
  }

  // nombre y extension
  char filename[_max_file_len];
  if (inombre&&inombre[0])
    strcpy(filename,inombre);
  else
    strcpy(filename,nombre);
  char *extaddr=ext_begin(filename);
  if (extaddr) {
    *(extaddr-1)=0;
    strcpy(ext,extaddr);
  }
  if (!ext[0]) strcpy(ext,"dat");
  if (!strcmp_nocase(ext,"xy")){
    if (!graba_xyz(filename)) {_savetime(error);return false;}
  }
//  else if (strlen(ext)==1&&ext[0]=='e') {if (!graba_e(filename)) {_savetime(error);return false;}}
  _graba(dat)
  _graba(dxf)
  _graba(msh)
  _graba(mai)
  _graba(key)
  _graba(con)
  _graba(wrl)
  else { // error pero de todos modos graba dat
    add_error(Unknown_Extension); add_error(ext);
    if (!graba_dat(filename)) {_savetime(error);return false;}
  }
  strcpy(nombre,filename);

  tipo.reset(m_modificada);
  _savetime(graba);
  return true;
}
#undef _graba

/*
// aca graba frontera haciendo una malla, pero se pierde la numeracion de nodos
bool malla::graba_frontera(const char *exti,bool todojunto) {
  _initime;
  if (!e||!mk_frontera()) _return_error
  if (!frontera) {_savetime(graba frontera);return true;}
  if (todojunto) {
    malla m(frontera[0]);
    // nombre y extension
    sprintf(m.nombre,"%s_bdry",nombre);
    if (exti&&exti[0]){
      if (exti[0]=='.') strcat(m.nombre,exti);
      else sprintf(m.nombre+strlen(m.nombre),".%s",exti);
    }
    else strcat(m.nombre,".dat");
    for (int i=1;i<frontera.len;i++) {
      malla mt(frontera[i]);
      m+=mt;
    }
    if(!(m.graba())) _return_error
  }
  else for (int i=0;i<frontera.len;i++){
    malla m(frontera[i]);
    // nombre y extension
    sprintf(m.nombre,"%s_%i",nombre,i);
    if (exti&&exti[0]){
      if (exti[0]=='.') strcat(m.nombre,exti);
      else sprintf(m.nombre+strlen(m.nombre),".%s",exti);
    }
    else strcat(m.nombre,".dat");
    if(!(m.graba())) _return_error
  }
  _savetime(graba frontera);

  return true;
}
*/
bool malla::graba_frontera() {
  _initime;
  if (!e||!mk_frontera()) {_savetime(error);return false;}
  if (!frontera) {_savetime(graba frontera);return true;}
  for (int i=0;i<frontera.len;i++){
    // nombre y extension
    char fnombre[_max_file_len];
    if (frontera.len==1) sprintf(fnombre,"%s_bdry.dat",nombre);
    else sprintf(fnombre,"%s_bdry%i.dat",nombre,i);
    if(!(frontera[i].graba_dat(fnombre))) {_savetime(error);return false;}
  }
  _savetime(graba frontera);
  return true;
}


double malla::epsilon_bb(){
  return epsilon=Max(pmin.distancia(pmax)*1e-8,ERRADM);
}


#define _lee(_ext) \
  else if (!strcmp_nocase(ext,#_ext)) {\
    if (!lee_##_ext(nombre,nl,el)) {_savetime(error);return false;}\
  }

bool malla::lee(const char* const_inombre){
  _initime;

  // nombre y extension
  char inombre[_max_file_len];
  if (const_inombre&&const_inombre[0]) strcpy(inombre,const_inombre);
  else { // usa el que ya tiene
    if (!nombre[0]){
      add_error(No_Mesh); _savetime(error); return false;
    }
    else {
      strcpy(inombre,nombre);
      strcat(inombre,".");
      if (ext[0]) strcat(inombre,ext); else strcat(inombre,"dat");
    }
  }

  ini();  renombra(inombre);

  // lee de acuerdo a filext
  if (INFO_CL) cout << "Reading" << endl;
  array1<nodo> nl; array1<elemento> el;
  if (!strcmp_nocase(ext,"xyz")) {if (!lee_xyz(nombre,nl)) {_savetime(error);return false;}}
  else if (!strcmp_nocase(ext,"xy")) {if (!lee_xy(nombre,nl)) {_savetime(error);return false;}}
  else if (strlen(ext)==1&&ext[0]=='e') {if (!lee_e(nombre,nl,el)) {_savetime(error);return false;}}
  _lee(dxf)
  _lee(dat)
  _lee(mai)
  _lee(msh)
  _lee(key)
  _lee(nas)
  _lee(ans)
  _lee(stl)
  else {add_error(Unknown_Extension); add_error(ext); _savetime(error);return false;}

  if (!nl) {add_error(No_Mesh); _savetime(error); return false;}

  if (INFO_CL) cout << "Analizing" << endl;

  int i,j,in,nc; flagtype f1;
  int &nlen=nl.len;
  n.resize(nlen);

  // bounding box y elimina borrados
  pmin=pmax=nl[0];
  for (i=1;i<nlen;i++) {
    if (nl[i].f.es(n_borrado)) continue;
    nl[i].set_min_max(pmin,pmax);
  }

  // hay z?
  bool hayz=pmax[2]-pmin[2]>ERRADM;
  if (!hayz) tipo.set(m_planaxy);
  else tipo.reset(m_planaxy);

  // octree (para no poner nodos repetidos)
  if (hayh) epsilon=hmin/1000; else epsilon_bb();
  o=new octree(n,pmin,pmax,hayz ? 3 : 2);
  bool puesto,hayrep=false;
  double d,dmax=0;
  pline map(nlen); map.len=nlen;
  // los nodos frontera van primero y los
  // de h (si hay) despues
  bool haynh=false;
  for (i=0; i<nlen; i++){
    f1=nl[i].f;
    if (f1.es(n_borrado)) continue;
    if (f1.es(n_h)) {haynh=true; continue;}
    if (f1) hayfn=true;
    in=map[i]=n+=nl[i]; // presupongo que es n.len-1
    o->add_no_rep(in,nc,puesto,epsilon);
    if (puesto) continue;
    if (nc<0||nc>=nlen)
      {add_error(Bad_Format); _savetime(error);return false;}
    hayrep=true;
    d=n[nc].distancia(nl[i]);
    if (d>dmax) dmax=d;
    n.len--; map[i]=nc;// si no es n.len-1 usar remove
    n[nc].f.set(nl[i].f); // oreo flags
    if (n[nc].h>ERRADM&&nl[i].h>ERRADM)
      set_min(n[nc].h,nl[i].h); // h es el menor
    else
      set_max(n[nc].h,nl[i].h); // h > 0
    n[nc].v=(n[nc].v+nl[i].v)/2;
  }
  if (hayrep&&strcmp_nocase(ext,"dxf")){
    add_warning(Repeated_Nodes); add_warning(dmax);
  }
  // los nodos de h
  if (haynh) for (i=0; i<nlen; i++){
    if (nl[i].f.es(n_borrado)) continue;
    if (nl[i].f.noes(n_h)) continue;
    in=map[i]=n+=nl[i]; // presupongo que es n.len-1
    o->add_no_rep(in,nc,puesto,epsilon);
    if (puesto) {
      nodosh++;
      if (n[in].h<ERRADM) n[in].h=MAXREAL;
      continue;
    }
    n.len--; map[i]=nc;// si no es n.len-1 usar remove
    // no oreo flags
    if (n[nc].h>ERRADM&&nl[i].h>ERRADM)
      set_min(n[nc].h,nl[i].h); // h es el menor
    else
      set_max(n[nc].h,nl[i].h); // h > 0
    if (n[nc].h<ERRADM) n[nc].h=MAXREAL;
    n[nc].v=(n[nc].v+nl[i].v)/2;
  }
  // minimo h
  if (hmin<ERRADM||hmin==MAXREAL) {hmin=MAXREAL; hayh=false;}
  else hayh=true;


  // modifica los elementos y hace eldenod
  if (el) {
    tipo.reset(m_nodos);
    int nv,ie,fdim=0;
    hayrep=false;
    e.resize(el);
    for (i=0;i<el.len;i++){
      elemento &ei=el[i]; nv=ei.nv();
      if (ei.f) hayfe=true;
      for (j=0;j<nv;j++) ei[j]=map[ei[j]];
      // verifica nodos h
      for (haynh=false,j=0;j<nv&&nodosh;j++) {
        if (n[ei[j]].f.es(n_h)) haynh=true;
      }
      // verifica si el elm tiene nodos repetidos
      for (j=0;j<nv;j++)
        if (ei[j]==ei.npos(j))
          break;
      if (j<nv) {
        if (ei.tipo()==e_cuadrilatero){
          elemento e3(e_triangulo);
          for (j=0,nv=0;j<4;j++) if (ei[j]!=ei.npos(j)) e3[nv++]=ei[j];
          ei=e3;
        }
        else if (ei.tipo()==e_cubo && ei[j+4]==ei.npos(j+4)){
          elemento e3(e_wedge);
          for (j=0,nv=0;j<4;j++) if (ei[j]!=ei.npos(j)) {e3[nv]=ei[j];e3[nv+3]=ei[j+4];nv++;}
          ei=e3;
        }
        else
          continue; // no lo agrega
      }
      // verifica si el elm esta repetido
      cpline &en=n[ei[0]].e;
      for (j=0;j<en.len;j++) if (ei==e[en[j]]) break;
      if (j<en.len) { // repetido
        hayrep=true;
        e[en[j]].f|=ei.f;
        continue;
      }
      // nodos h
      if (haynh){ // no agrega el elemento, solo saca el h de cada nodo
        for (j=0;j<nv;j++) {
          nodo &ni=n[ei[j]];
          if (ni.f.noes(n_h)) continue;
          if (ni.h!=MAXREAL) continue; // ya tiene h
          d=ni.distancia(n[ei.nant(j)]);
          if (d>ERRADM) ni.h=d;
          d=ni.distancia(n[ei.npos(j)]);
          if (d>ERRADM&&d<ni.h) ni.h=d;
        }
        continue;
      }
      // agrega el elemento
      ie=e+=ei;
      for (j=0;j<nv;j++) n[ei[j]].e+=ie;
      fdim|=(1<<ei.dim());
    }

    if (fdim&1) tipo.set(m_nodos);
    if (fdim&2) tipo.set(m_lin);
    if (fdim&4) tipo.set(m_sup);
    if (fdim&8) tipo.set(m_vol);

    if (hayrep)
      add_warning(Repeated_Elements);
  }
  else tipo.set(m_nodos);

  // paso de customizacion
  custom();

  _infotime(nombre);
  _infotime(endl);
  _savetime(lee);
  _infotime("nodos: ");
  _infotime(n.len);
  _infotime("; elementos: ");
  _infotime(e.len);
  _infotime(endl);

  return true;
}
#undef _lee

bool malla::agrega(const char* inombre){
// el octree se usa para soldar al hacer merge (+=)
  malla ml;
  if (!ml.lee(inombre))
    {add_error(ml.m_error); return false;}
  operator+=(ml);
  return true;
}

bool malla::agrega_con(const char* stipo, const char* archivo){ // solo conectividades
  if (!n) return false;
  cascara c; c.parent=this;
  if (!c.agrega_con(stipo,archivo)||!c.e) return false;
  rm_vecino(); rm_frontera(); rm_dir(); rm_nn(); rm_esferas(); ve.ini();
  tipo.reset(m_nodos|m_cerrada|m_orientada);
  int i,j,nv=c.e[0].nv(),d=c.e[0].dim();
  if      (d==1) tipo.set(m_lin);
  else if (d==2) tipo.set(m_sup);
  else if (d==3) tipo.set(m_vol);
  else return false;
  e.resize(e.len+c.e.len);
  for (i=0;i<c.e.len;i++){
    const elemento &ei=c.e[i];
    for(j=0;j<nv;j++) n[ei[j]].e+=e.len;
    e+=ei;
  }
  if (c.hayfe) hayfe=true; // si no no se sabe
  return true;
}

bool malla::agrega_elementos(const array1<elemento> &el){
  if (!n||!el) return false;
  int i,j,nv=el[0].nv(),d=0,nel=el.len;
  // a partir de aqui debe ser igual a la rutina que sigue
  rm_vecino(); rm_frontera(); rm_dir(); rm_nn(); rm_esferas(); ve.ini();
  tipo.reset(m_nodos|m_cerrada|m_orientada);
  e.resize(e.len+nel);
  for (i=0;i<nel;i++){
    const elemento &ei=el[i]; d|=(1<<ei.dim());
    for(j=0;j<nv;j++) n[ei[j]].e+=e.len;
    e+=ei;
  }
  if (d&2) tipo.set(m_lin); if (d&4) tipo.set(m_sup); if (d&8) tipo.set(m_vol);
  return true;
}
bool malla::agrega_elementos(int nel,const elemento *el){
  if (!n||!nel||!el) return false;
  int i,j,nv=el[0].nv(),d=0;
  // a partir de aqui debe ser igual a la rutina anterior
  rm_vecino(); rm_frontera(); rm_dir(); rm_nn(); rm_esferas(); ve.ini();
  tipo.reset(m_nodos|m_cerrada|m_orientada);
  e.resize(e.len+nel);
  for (i=0;i<nel;i++){
    const elemento &ei=el[i]; d|=(1<<ei.dim());
    for(j=0;j<nv;j++) n[ei[j]].e+=e.len;
    e+=ei;
  }
  if (d&2) tipo.set(m_lin); if (d&4) tipo.set(m_sup); if (d&8) tipo.set(m_vol);
  return true;
}

malla& malla::operator +=(malla &m){// suma de mallas
  _initime;
  int i,j;
  array1 <nodo> nodo_h(nodosh);

  // se puede sumar fronteras, pero...
  frontera.ini(); vecino.ini();
  nn.ini(); rm_dir();
  bool
    hayce=(!e||ce)&&(!m.e||m.ce),
    hayre=(!e||re)&&(!m.e||m.re),
    hayve=(!e||ve)&&(!m.e||m.ve);

  // bounding box
  punto oldpmin=pmin,oldpmax=pmax;
  m.bbox();
  m.pmin.set_min_max(pmin,pmax);
  m.pmax.set_min_max(pmin,pmax);

  // hay z?
  bool hayz=pmax[2]-pmin[2]>ERRADM;
  if (!hayz) tipo.set(m_planaxy);
  else tipo.reset(m_planaxy);

  // flags
  hayfn|=m.hayfn; hayfe|=m.hayfe;

  // octree
  epsilon_bb();
  if (!o||
      (pmin[0]<oldpmin[0]||pmax[0]>=oldpmax[0]) ||
      (pmin[1]<oldpmin[1]||pmax[1]>=oldpmax[1]) ||
      (pmin[2]<oldpmin[2]||pmax[2]>=oldpmax[2]) ||
      nodosh) {
    if (o) delete o;
    o=new octree(n,pmin,pmax,hayz ? 3 : 2);
    for (i=0;i<n.len-nodosh;i++) o->add(i); // no hay que usar epsilon
    for (;i<n.len;i++) nodo_h+=n[i]; n.len-=nodosh;
  }

  if ((!n||hayh)&&m.hayh) {hayh=true; set_min(hmin,m.hmin);} else hayh=false;
  if ((!n||hayv)&&m.hayv) {hayv=true; set_min(nvmin,m.nvmin); set_max(nvmax,m.nvmax);}
  else {hayv=false; nvmin=MAXREAL; nvmax=-MAXREAL;}

  // agrego los nodos nuevos (eldenod despues)
  bool puesto; int nc,oldnlen=n.len;
  pline map(m.n.len),vacia;
  n.resize(n.len+m.n.len);
  for (i=0; i<m.n.len-m.nodosh; i++){
    nodo &ni=m.n[i];
    if (ni.f.es(n_h)) continue; // los de h al final
    ni.e.clean(); map+=n+=ni; // presupongo que es n.len-1
    o->add_no_rep(n.len-1,nc,puesto,epsilon);
    if (puesto) continue;
    n.len--; map[i]=nc;// ojo si no es n.len-1
    n[nc].f.set(ni.f);
    if (n[nc].h>ERRADM&&ni.h>ERRADM)
      set_min(n[nc].h,ni.h); // h es el menor
    else
      set_max(n[nc].h,ni.h); // h > 0
    n[nc].v=(n[nc].v+ni.v)/2;
  }

  // nodos de h (no tienen elementos)
  nodosh=0;
  for (i=0;i<nodo_h.len;i++){ // de esta malla
    nodo &ni=nodo_h[i];
    map+=n+=ni; // presupongo que es n.len-1
    o->add_no_rep(n.len-1,nc,puesto,epsilon);
    if (puesto) {nodosh++; continue;}
    n.len--; map[i]=nc;// ojo si no es n.len-1
    if (n[nc].h>ERRADM&&ni.h>ERRADM)
      set_min(n[nc].h,ni.h); // h es el menor
    else
      set_max(n[nc].h,ni.h); // h > 0
    if (n[nc].h<ERRADM) n[nc].h=MAXREAL;
    n[nc].v=(n[nc].v+ni.v)/2;
  }
  for (i=m.n.len-m.nodosh;i<m.n.len;i++){ // de b
    nodo &ni=m.n[i];
    map+=n+=ni; // presupongo que es n.len-1
    o->add_no_rep(n.len-1,nc,puesto,epsilon);
    if (puesto) {nodosh++; continue;}
    n.len--; map[i]=nc;// ojo si no es n.len-1
    if (n[nc].h>ERRADM&&ni.h>ERRADM)
      set_min(n[nc].h,ni.h); // h es el menor
    else
      set_max(n[nc].h,ni.h); // h > 0
    if (n[nc].h<ERRADM) n[nc].h=MAXREAL;
    n[nc].v=(n[nc].v+ni.v)/2;
  }

  // agrega los elementos nuevos y arregla eldenod
  int k,in,nv,ix,fdim=0;
  e.resize(e.len+m.e.len);
  bool rep,hayrep=false;
  for (i=0;i<m.e.len;i++){ // b es const
    ix=e+=m.e[i]; elemento &ei=e[ix]; nv=ei.nv();
    for (j=0;j<nv;j++) ei[j]=map[ei[j]]; // indices nuevos
    // chequea si no esta repetido
    in=ei[0]; rep=false; if (in<oldnlen) {// si no es un nodo nuevo
      const cpline &eni=n[in].e; // este todavia no esta en eni
      for (k=0;k<eni.len;k++)
        if (ei==e[eni[k]]) {
          rep=true; break;}
    }
    if (rep) {
      hayrep=true;
      e.len--; //(supongo ix=len-1)
      continue;
    }
   for (j=0;j<nv;j++) n[ei[j]].e+=ix;
    fdim|=(1<<ei.dim());
    if (hayce) ce+=m.ce[i]; else ce.ini();
    if (hayre) re+=m.re[i]; else re.ini();
    if (hayve) ve+=m.ve[i]; else ve.ini();
  }
  if (hayrep)
    add_warning(Repeated_Elements);

  // tipo
  tipo.reset(m_orientada|m_cerrada);
  tipo.set(m_modificada);
  if (fdim&1) tipo.set(m_nodos);
  if (fdim&2) tipo.set(m_lin);
  if (fdim&4) tipo.set(m_sup);
  if (fdim&8) tipo.set(m_vol);

  if (tipo.noes(m_slivers_marked)||m.tipo.noes(m_slivers_marked))
    tipo.reset(m_slivers_marked);

  _savetime(malla+=);
  return *this;
}

void malla::mk_octree(bool remake){// octree
  if (o&&!remake) return;
  if (epsilon<=ERRADM) epsilon_bb();
  if (epsilon==0) return;
  _initime;
  if (o) delete o;
  o=new octree(n,pmin,pmax,tipo.es(m_planaxy) ? 2 : 3);
  for (int i=0;i<n.len;i++) o->add(i);
  _savetime(octree);
}

bool malla::weld(double weldfactor){ // suelda
  _initime;
  bbox();
  if (!hayh&&!mk_h_nn_min()) return false;
  if (o) delete o;
  o=new octree(n,pmin,pmax,tipo.es(m_planaxy) ? 2 : 3);
  int i,j,ivuela,nc,iqueda,ilast;
  double minh;
  bool puesto,retval=false,hayerep=false;
  for (i=0;i<n.len;i++) {
    // intenta con el h del nodo i
    o->add_no_rep(i,nc,puesto,n[i].h*weldfactor);
    if (puesto) continue;
    // puede que uno de los nodos tenga un h mucho mayor que el otro
    // intenta, por las dudads, con minh
    do{
      iqueda=nc;
      minh=Min(n[i].h,n[iqueda].h);
      o->add_no_rep(i,nc,puesto,minh*weldfactor);
      if (puesto) break;
    } while(iqueda!=nc);
    if (puesto) continue;

    ivuela=i;
    if (n[i].f.es(n_permanente)) Swap(ivuela,iqueda);
    // verifica si algun elemento tiene a ambos
    cpline &enq=n[iqueda].e;
    for (j=0;j<enq.len;j++)
      if (e[enq[j]].have(ivuela)) break;
    if (j<enq.len)
      continue;
    // vuela nomas
    combine(ivuela,iqueda); retval=true;
    if (iqueda==i) {o->remove(ivuela); o->add(iqueda);}
    if (ivuela<n.len-nodosh||iqueda<n.len-nodosh) n[iqueda].f.reset(n_h);
    nn.ini(); tipo.reset(m_orientada); frontera.ini(); vecino.ini();
    // reemplaza en los elementos
    const cpline &env=n[ivuela].e;
    for (j=0;j<env.len;j++) e[env[j]].replace(ivuela,iqueda);
    enq+=env;
    // cambia por el ultimo segun sea h o no-h (solo el que vuela)
    ilast=n.len-1; if (ivuela>=n.len-nodosh) {ilast-=nodosh; --nodosh;}
    if (ivuela!=ilast){
      const cpline &enl=n[ilast].e;
      for (j=0;j<enl.len;j++) e[enl[j]].replace(ilast,ivuela);
      n[ivuela]=n[ilast];
      if (ilast<n.len-1) // borro un nodo no-h pero hay nodosh
        n[ilast]=n.last(); // llena el hueco con el ultimo nodo de h
    }
    n.len--; i--;
    // verifica si hay elementos repetidos del nodo soldado
    hayerep|=elm_rep(iqueda);
  }
  if (hayerep)
    add_warning(Repeated_Elements);

  _savetime(weld);
  return retval;
}

bool malla::elm_rep(int in){
// verifica si hay elementos repetidos de un nodo
  cpline &en=n[in].e;
  int j,k,l,iqueda,ivuela,ilast;
  bool hayerep=false;

  for (j=0;j<en.len-1;j++) {
    for (k=j+1;k<en.len;k++) {
      if (e[en[j]]!=e[en[k]]) continue;
      // repetido => lo swappea con el ultimo
      ilast=e.len-1; ivuela=en[k]; iqueda=en[j];
      const elemento &eq=e[iqueda];
       // elimina ivuela de las listas de elementos de sus nodos
      for (l=0;l<eq.nv();l++) n[eq[l]].e.remove1(ivuela);
      // swapea con el ultimo
      if (ivuela!=ilast) {
        const elemento &elast=e[ilast];
        for (l=0;l<elast.nv();l++) n[elast[l]].e.replace1(ilast,ivuela);
        e[ivuela]=elast;
      }
      // elimina el elemento
      e.len--; hayerep=true;
      vecino.ini();ve.ini();ce.ini();re.ini(); if (edir) rm_dir();
      if (k<en.len-1) {en[k]=en.last(); k--;} en.len--;
    }
  }
  return hayerep;
}


void malla::bbox(bool remake){
// bounding box
  if (!remake&&(pmin.distanciac(pmax)>ERRADM)) return;
  pmin=pmax=n[0];
  for (int i=1;i<n.len;i++) n[i].set_min_max(pmin,pmax);
}

bool malla::extrude(const punto &v){
  tipo.reset(m_planaxy|m_cerrada|m_sup|m_lin|m_nodos);
  tipo.set(m_modificada);

  rm_nn(); // vecinos naturales
  rm_esferas(); ve.ini(); // esferas y volumen
  rm_dir(); // normales
  rm_frontera(); // frontera
  if (o) delete o;  o=0; // octree

  // bounding box
  (pmin+v).set_min_max(pmin,pmax);
  (pmax+v).set_min_max(pmin,pmax);

  //copia los nodos
  punto p;
  int i,nlen=n.len;
  for (i=0;i<nlen;i++) {
    p=n[i]+v;
    n+=n[i]; n.last().setpos(p); // mismo h, flag y valor
  }

  //extruda los elementos (faltan....)
  elemento s(e_segmento),q(e_cuadrilatero),w(e_wedge),c(e_cubo);
  e_tipo etipo; int fdim=0;
  for(i=0;i<e.len;i++){
    elemento &ei=e[i]; etipo=ei.tipo();
    if (etipo==e_segmento){
      q[0]=ei[0];q[1]=ei[1];
      q[2]=ei[0]+nlen;q[3]=ei[1]+nlen;
      ei=q;
    }
    else if (etipo==e_triangulo){
      w[0]=ei[0];w[1]=ei[1];w[2]=ei[2];
      w[3]=ei[0]+nlen;w[4]=ei[1]+nlen;w[5]=ei[2]+nlen;
      ei=w;
    }
    else if (etipo==e_cuadrilatero){
      c[0]=ei[0];c[1]=ei[1];c[2]=ei[2];c[3]=ei[3];
      c[4]=ei[0]+nlen;c[5]=ei[1]+nlen;c[6]=ei[2]+nlen;c[7]=ei[3]+nlen;
      ei=c;
    }
    fdim|=(1<<ei.dim());
  }
  if (fdim&1) tipo.set(m_nodos);
  if (fdim&2) tipo.set(m_lin);
  if (fdim&4) tipo.set(m_sup);
  if (fdim&8) tipo.set(m_vol);

  // segmentos con los nodos sueltos
  bool haysegs=false;
  for (i=0;i<nlen;i++) {
    if (n[i].e) continue;
    s[0]=i; s[1]=i+nlen;
    n[i].e+=e.len; n[i+nlen].e+=e.len;
    e+=s;
    haysegs=true;
  }
  if (haysegs) tipo.set(m_lin);

  return true;
}

bool malla::extrude(const char* tabfile){

  malla m; if (!m.lee(tabfile)) return false;

  punto p,v;
  elemento s(e_segmento),q(e_cuadrilatero),w(e_wedge),c(e_cubo);
  e_tipo etipo;

  tipo.reset(m_planaxy|m_cerrada|m_sup|m_lin|m_nodos);
  tipo.set(m_modificada);
  rm_nn(); // vecinos naturales
  rm_esferas(); ve.ini(); // esferas y volumen
  rm_dir(); // normales
  rm_frontera(); // frontera
  if (o) delete o;  o=0; // octree

  int i,j,k,start=0,nlen=n.len,elen=e.len,enlen;

  n.resize(n.len*m.n.len);
  array1<elemento> eori=e;
  e_ini(); e.resize(elen*m.e.len);

  for (j=0;j<m.e.len;j++,start+=nlen){
    v=m.n[m.e[j][1]]-m.n[m.e[j][0]];
    //copia los nodos
    for (i=start;i<start+nlen;i++) {
      p=n[i]+v;
      n+=n[i];
      nodo& ni=n.last(); ni.setpos(p);
      enlen=ni.e.len;
    }

    //extruda los elementos
    for(i=0;i<elen;i++){
      elemento &ei=eori[i]; etipo=ei.tipo();
      if (etipo==e_segmento){
        q[0]=ei[0]+start;q[1]=ei[1]+start;
        q[2]=q[0]+nlen;q[3]=q[1]+nlen;
        e+=q;
        tipo.set(m_sup);
      }
      else if (etipo==e_triangulo){
        w[0]=ei[0]+start;w[1]=ei[1]+start;w[2]=ei[2]+start;
        w[3]=w[0]+nlen;w[4]=w[1]+nlen;w[5]=w[2]+nlen;
        e+=w;
        tipo.set(m_vol);
      }
      else if (etipo==e_cuadrilatero){
        c[0]=ei[0]+start;c[1]=ei[1]+start;c[2]=ei[2]+start;c[3]=ei[3]+start;
        c[4]=c[0]+nlen;c[5]=c[1]+nlen;c[6]=c[2]+nlen;c[7]=c[3]+nlen;
        e+=c;
        tipo.set(m_vol);
      }
      const elemento &elast=e.last();
      for (k=0;k<elast.nv();k++) n[elast[k]].e+=e.len-1;
    }
  }

  // segmentos con los nodos sueltos
  for (i=0;i<nlen;i++) {
    if (n[i].e) continue;
    s[1]=i;
    for (j=0;j<m.e.len;j++){
      s[0]=s[1]; s[1]=s[0]+nlen;
      n[s[0]].e+=e.len; n[s[1]].e+=e.len;
      e+=s;
    }
    tipo.set(m_lin);
  }

  bbox(true);// bounding box
  return true;
}

bool malla::diam_aspect(double cota_superior){
  bool haye=e;
  if (!haye&&!delaunay(true)) return false;; // sin slivers
  if (!re) {mk_esferas(); if (!re) return false;}
  int i,j,nlen=n.len,enlen;
  double minr,maxr;
  hayv=true;
  for (i=0;i<nlen;i++){
    if (n[i].f.es(n_frontera)){
      n[i].v=cota_superior;
      continue;
    }
    const cpline& en=n[i].e; enlen=en.len;
    if (!enlen) {n[i].v=0; continue;}
    minr=maxr=re[en[0]];
    for (j=0;j<enlen;j++)
      set_min_max(minr,maxr,re[en[j]]);
    n[i].v=Min(maxr/minr,cota_superior);
  }
  for (i=0;i<nlen;i++) set_min_max(nvmin,nvmax,n[i].v);

//  if (!haye) e_ini();
  return true;
}
// intercambia dos nodos
void malla::swap_n(int in1,int in2){
  int j;

  // elementos de los nodos
  // si un elemento tiene los dos y no hago esto (n.len) puede quedar invertido
  const cpline &en1=n[in1].e; // elementos de n1
  for (j=0;j<en1.len;j++) {
    e[en1[j]].replace(in1,n.len);
  }
  const cpline &en2=n[in2].e; // elementos de n2
  for (j=0;j<en2.len;j++) {
    e[en2[j]].replace(in2,in1);
  }
  for (j=0;j<en1.len;j++) {// elementos de n1
    e[en1[j]].replace(n.len,in2);
  }

  // listas
  if (nn) nn.swap(in1,in2);
  if (o&&!o->swap_n(in1,in2)){
    delete o; o=0;
  }
  if (ndir&&dir) dir.swap(in1,in2);

  // los elementos
  n.swap(in1,in2);
}

bool malla::custom(){
///////////////
// aca se pueden agregar rutinas para cambiar o deformar la malla leida
// no olvidarse de rehacer el ocree, el bbox y lo que haga falta
///////////////

/*
#pragma message("<---------------- ***   OJO, CUSTOMIZADO!!!  ***\n")

  bool tomarselas=true;
  __asm int 3;
  if (tomarselas) return false;
*/

  /*
  //elimino nodos de uno o dos elementos (y los elementos)
  int i,j;
  rm_vecino();
  delete o; o=0;
  for (i=0;i<n.len;i++){
    if (n[i].e.len>1) continue;
    cpline &en=n[i].e;
    for (j=0;j<en.len;j++) vuelae(en[j],false);
    if (i==n.len-1) {n.len--; continue;}
    nodo &nl=n.last();
    cpline &enl=nl.e;
    for (j=0;j<enl.len;j++) e[enl[j]].replace(n.len-1,i);
    n[i]=nl; n.len--;
  }
  */

  /*
  // elimino nodos sin elementos
  pline map(1000); array1<nodo> newn(1000);
  for (int in=0;in<n.len;in++) {
    if (n[in].e) {map+=in; newn+=n[in];}
  }
  for (int ie=0;ie<e.len;ie++) {
    for (int iv=0;iv<e[ie].nv();iv++) {
      e[ie][iv]=map.index(e[ie][iv]);
    }
  }
  n=newn;
  */

/*
  mk_vecino();
  for (int i=0;i<e.len;i++){    
    if (e[i].tipo()==e_wedge) continue;
    if (vecino[i][0]>=0&&e[vecino[i][0]].tipo()==e_wedge) continue;
    if (vecino[i][1]>=0&&e[vecino[i][1]].tipo()==e_wedge) continue;
    if (vecino[i][2]>=0&&e[vecino[i][2]].tipo()==e_wedge) continue;
    if (vecino[i][3]>=0&&e[vecino[i][3]].tipo()==e_wedge) continue;
    vuelae(i,true);
  }
  graba_dat();
*/

/*
  // soldadura
  int i;
  rm_vecino(); rm_frontera();
  malla m(*this);
  for (i=0;i<59;i++){
    for (int i=0;i<m.n.len;i++) {
      m.n[i][2]+=.015;
    }
    m.bbox(true);
    *this+=m;
  }
  bbox(true);
  graba_dat("media_placa");
  m=*this;
  for (i=0;i<m.n.len;i++) {
    m.n[i][0]=-m.n[i][0];
  }
  m.bbox(true);
  *this+=m;
  orienta(true);
  renombra("placa");
  graba_dat();
  graba_frontera();
*/

/*
  // nodevalues
  if (hayv) return true; hayv=true;
  for (int i=0;i<n.len;i++){
    //n[i].v=n[i][0]*n[i][1]+n[i][1]*n[i][1];
    n[i].v=n[i].h;
    set_min_max(nvmin,nvmax,n[i].v);
  }
*/

/*
// test de delaunay refinado
   malla vieja;
   char v_nombre[255]; strcpy(v_nombre,nombre); int len=strlen(v_nombre);
   char* last=&(v_nombre[len-1]);
   (*last)--; while ((*last)=='/') {
     (*last)='9';
     last--; (*last)--;
   }

   vieja.lee(v_nombre);
   delaunay_refinado(vieja,true);
*/

/*
  //escala y mirror y -> -y
  int i;
  for (i=0;i<n.len;i++) n[i][1]=3*n[i][1];
  pmin[1]=3*pmin[1]; pmax[1]=3*pmax[1];
  delete o; o=0; mk_octree();
  malla m(*this);
  for (i=0;i<n.len;i++) m.n[i][1]=-m.n[i][1];
  m.swap();
  m.pmin[1]=-pmax[1]; pmax[1]=-pmin[1];
  delete m.o; m.o=0; m.mk_octree();
  *this+=m;
*/
/*
  for (int i=0;i<n.len;i++) {n[i].h=2+n[i][1]*3/100;}
  hayh=true; hmin=2;
*/
/*
  //cilindro -> plano
  delete o; o=0;
  double ang,r=44;
  for (int i=0;i<n.len;i++) {
    nodo &ni=n[i];
    ang=atan2(ni[1],ni[0]);
    ni[1]=ni[2];
    ni[0]=r*ang;
    ni[2]=0;
  }
  tipo.set(m_planaxy);
  bbox(true); mk_octree();
*/
/*
  //plano -> cilindro
  delete o; o=0;
  double ang,r=44;
  for (int i=0;i<n.len;i++) {
    nodo &ni=n[i];
    ang=ni[0]/r;
    ni[2]=ni[1];
    ni[0]=r*cos(ang);
    ni[1]=r*sin(ang);
  }
  tipo.reset(m_planaxy);
  bbox(true); mk_octree();
*/

/*
  // plano no-xy
  punto norm(-0.00018166,0.47128425,-0.88198136);
  punto xv=norm%_ez; double ang=asin(xv.mod());

  delete o; o=0;
  for (int i=0;i<n.len;i++) n[i].setpos(n[i].giro(xv,ang));
  bbox(true); mk_octree();
*/
/*
 // fija h
 for (i=0; i<n.len;i++){
//   double y=n[i][1];
   n[i].h=.75;
      //1+(y-20)/(150-20)*(10-1);
 }
 hayh=true; hmin=.75;
*/
/*
 // perturba
 delete o; o=0;
 int i; punto pr;
 for (i=0; i<n.len;i++){
    n[i]+=pr.randdir()*5;
    n[i].set_min_max(pmin,pmax);
 }
 mk_octree();
*/

//  bbox();
  return true;

}
