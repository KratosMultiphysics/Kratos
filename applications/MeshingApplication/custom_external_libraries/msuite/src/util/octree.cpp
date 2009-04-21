// Las comparaciones de puntos se hacen con <=
// un punto esta en un caja si esta dentro
// o en la frontera "inferior"

// los childs se numeran de acuerdo al bit:
//                  2 | 3
//                  --+--
//                  0 | 1

// La distancia entre puntos es max_delta(x,y,z) (punto::distanciac)

#include "utiles.h" // MAXREAL
#include "octree.h"

using namespace std;

static double epsilon=1e-5;

//---------------------------------------------------------------
//saca un script de autocad con un cuadtree (si, solo 2D)
#ifdef _OADEBUG
#include <fstream>
static ofstream oad;
static void _oaini(const punto &pmin, const punto &pmax){
  oad.open("octree.scr");
  punto::o_separator(",");
  oad << "osmode 0\nline "
      << pmin << " "
      << punto(pmax,pmin) << " "
      << pmax << " "
      << punto(pmin,pmax) << " c\n";
}
static void _oadivide(const punto &pmin, const punto &pmax){
  oad << "line "
      << punto(pmin[0],(pmin[1]+pmax[1])/2) << " "
      << punto(pmax[0],(pmin[1]+pmax[1])/2) << " \n"
      << "line "
      << punto((pmin[0]+pmax[0])/2,pmin[1]) << " "
      << punto((pmin[0]+pmax[0])/2,pmax[1]) << " \n";
}
static void _oaclose(){
  punto::io_reset();
  oad.close();
}
#endif
//---------------------------------------------------------------

// caja es la clase de las ramas del octree
class caja{
public:
  octree *o; // el octree
  caja *parent; // un parent o 0 si es el mayor
  caja *child;  // 2^dim childs o 0 si es terminal
  int len; // cantidad de puntos dentro (o dentro de los childs)
  int* ip; // los puntos que hay dentro (solo si es terminal)
  int prof;// profundidad, solo para calcular cercano :-(
  bool pasado; // idem
  punto pmin,pmax; // bbox

  caja() :
    o(0),parent(0),child(0),
    len(0),ip(0),prof(0),pasado(false){}

  ~caja(){
    if (child) delete [] child;
    if (ip)    delete [] ip;
  }

  bool check_interior(const punto& p); // verifica con el bbox
  // chlidp y cajap presuponen que p es interior
  int childp(const punto& p); // en que child va?
  caja* cajap(const punto&); // en que sub-caja va? (O(log(n)))
  int otro(int ix); // devuelve otro nodo del caja o el parent
  bool divide(); // hace los childs y reparte los puntos
  caja* add(int ix); // agrega
  void closest_en_caja // el mas cercano a p de este caja
    (const punto &p, int &nc, double &d, double e);
  int closest(int ix, double &d, double e);
  int closest(const punto &p, double &d);
};

// calcula el bbox
static void bbox(const _OCTREE_PTARRAY &p, int plen,
                 punto &pmin, punto &pmax){
  pmin=pmax=p[0];
  for (int i=1;i<plen;i++) p[i].set_min_max(pmin,pmax);
}

/*
// sube el techo de pmax
static void add_epsilon(const punto &pmin, punto &pmax, int dim){
  while (dim--) {Max(pmax[dim]+=Max(pmax[dim]-pmin[dim],fabs(pmax[dim]))*ERRADM,ERRADM);}
}
*/

// costructor de octree sabiendo de antemano el bbox
octree::octree(const _OCTREE_PTARRAY &p,
               const punto &ipmin,const punto &ipmax,
               int idim, int imaxlen)
: dim(idim), nchild(1<<dim), maxlen(imaxlen), plist(&p){
  pmin=ipmin; pmax=ipmax;
  double delta=Max(pmin.distancia(pmax)*epsilon,ERRADM);
  pmax+=punto(delta,delta,delta); /*add_epsilon(pmin,pmax,dim);*/
  pmin-=punto(ERRADM,ERRADM,ERRADM); 
  child=new caja;
  child->o=this;
  child->parent=0;
  child->pmin=pmin;
  child->pmax=pmax;
#ifdef _OADEBUG
  _oaini(pmin,pmax);
#endif
}

octree::octree(const _OCTREE_PTARRAY &p, int plen,
               int idim, int imaxlen)
:dim(idim), nchild(1<<dim), maxlen(imaxlen), plist(&p){
  bbox(p,plen,pmin,pmax);
  double delta=(pmax-pmin).mod()*epsilon;
  pmax+=punto(delta,delta,delta); /*add_epsilon(pmin,pmax,dim);*/
  pmin-=punto(ERRADM,ERRADM,ERRADM); 
  child=new caja;
  child->o=this;
  child->parent=0;
  child->pmin=pmin;
  child->pmax=pmax;
#ifdef _OADEBUG
  _oaini(pmin,pmax);
#endif
}

octree::~octree(){
  delete child;
#ifdef _OADEBUG
  _oaclose();
#endif
}

// verifica que este dentro del bbox
bool octree::check_interior(const punto& p) const{
  if (p[0]<pmin[0]||p[0]>=pmax[0]) return false;
  if (dim==1) return true;
  if (p[1]<pmin[1]||p[1]>=pmax[1]) return false;
  if (dim==2) return true;
  if (p[2]<pmin[2]||p[2]>=pmax[2]) return false;
  return true;
}

bool caja::check_interior(const punto& p) {
  if (p[0]<pmin[0]||p[0]>=pmax[0]) return false;
  if (o->dim==1) return true;
  if (p[1]<pmin[1]||p[1]>=pmax[1]) return false;
  if (o->dim==2) return true;
  if (p[2]<pmin[2]||p[2]>=pmax[2]) return false;
  return true;
}

// se supone que entro por un caja grande que incluye al punto
int caja::childp(const punto& p) {
  int nch=0;
  punto pcen=(pmin+pmax)/2;
  if (pcen[0]<=p[0]) nch+=1;
  if (o->dim==1) return nch;
  if (pcen[1]<=p[1]) nch+=2;
  if (o->dim==2) return nch;
  if (pcen[2]<=p[2]) nch+=4;
  return nch;
}
caja* caja::cajap(const punto& p){
  if (!child)  // va en este
/*    if (
      pmax[0]==pmin[0]||
      ((o->dim>1)&&(pmax[1]==pmin[1]))||
      ((o->dim>2)&&(pmax[2]==pmin[2]))|
      )
      return 0;*/
    return this;
  return (child[childp(p)].cajap(p));
}

// busca un nodo distinto de ix en el caja o sus hermanos
// no asumo que ix es el ultimo
int caja::otro(int ix){
  if (pasado) return -1;
  pasado=true;
  int nc=-1;
  if (ip){ // es terminal con nodos
    if (ip[0]!=ix) nc=ip[0];
    else if (len>1) nc=ip[1];
    else if (parent) nc=parent->otro(ix);
  }
  // es parent o terminal sin nodos
  else if (!len){ // no tiene nodos
    if (parent) nc=parent->otro(ix);
  }
  // es parent con len (y len debe ser > 1 sino no estaria dividido)
  else {
    for (int i=0;i<o->nchild;i++){
      nc=child[i].otro(ix);
      if (nc!=-1) break;
    }
  }
  pasado=false; return nc;
}

int octree::close(punto p) const{
  for (int i=0;i<dim;i++){
    if (p[i]<pmin[i])  p[i]=pmin[i];
    if (p[i]>=pmax[i]) p[i]=pmax[i]-(pmax[i]-pmin[i])*epsilon;
  }
  caja &t=*(child->cajap(p));
  return t.otro(-1);
}

// hace los childs y reparte los puntos
bool caja::divide(){
  int i,j,k,profc=prof+1;
  punto dp=(pmax-pmin)/2;
  if (
      fabs(dp[0])<ERRADM ||
      ((o->dim>1)&&(fabs(dp[1])<ERRADM)) ||
      ((o->dim>2)&&(fabs(dp[2])<ERRADM))
      ){
    return false;
  }
  // hace los  child
  child=new caja[o->nchild];
#ifdef _OADEBUG
  _oadivide(pmin,pmax);
#endif
  for (i=0;i<o->nchild; i++){
    child[i].o=o; child[i].parent=this; child[i].prof=profc;
    child[i].pmin=pmin;
    for (j=0;j<o->dim;j++) if (i&(1<<j)) child[i].pmin[j]+=dp[j];
    child[i].pmax=child[i].pmin+dp;
  }
  // reparte los nodos
  const _OCTREE_PTARRAY &plist=*(o->plist);
  for (k=0;k<len;k++){
    const punto &pk=plist[ip[k]];
    caja &ch=*cajap(pk);
    if (!ch.ip) ch.ip=new int [o->maxlen];
    ch.ip[ch.len++]=ip[k];}
  delete [] ip; ip=0;
  return true;
}

// agrega el nodo
caja* caja::add(int ix){
  const punto &p=(*(o->plist))[ix];
  caja *t=this;
  while ((t=t->parent)) t->len++; // aumenta len a todos los parent
  t=this;
  while (t->len==o->maxlen) { // subdivide
    // aca puede haber stack overflow
    if(!t->divide()) { // maxlen puntos iguales??
      while ((t=t->parent)) t->len--; // deshace
      return 0;
    }
    t->len++;
    t=t->cajap(p);
  }
  if (!t->ip) t->ip=new int [o->maxlen];
  t->ip[t->len++]=ix;
  return t;
}

// agrega un punto y devuelve un punto cercano
int octree::add(int ix){
  const punto &p=(*plist)[ix];
  if (!check_interior(p)) return -1;
  caja *t=child->cajap(p); /// en que box va
  t=t->add(ix);
  if (!t) return -1;
  return t->otro(ix); // busca otro nodo para devolver
}

// punto mas cercano
// p no esta necesariamente en este caja
void caja::closest_en_caja
    (const punto &p, int &nc, double &d, double e){
  int i;
  if (!ip&&!len) return; // aca no hay nada
  if (!ip) { // hay childs
    for (i=0;i<o->nchild;i++) {
      child[i].closest_en_caja(p,nc,d,e);
      if (d<e) return;
    }
    return;
  }
  // es terminal con nodos
  const _OCTREE_PTARRAY &plist=*(o->plist);
  for (i=0;i<len;i++){
    if (p.distanciac_menor(plist[ip[i]],d)) {
      nc=ip[i]; if (d<e) return;
    }
  }
}

// el while es porque puede dmin abarcar varios cajaes
static bool _test_mayor
(caja *vin, const punto &p, punto &pt,
 double &d, double e, int &nc, int i0){
  int i,dim=vin->o->dim;
  for (i=i0;i<dim;i++){
    caja *v=vin;
    while (v && p[i]+d>=v->pmax[i]){
      pt[i]=v->pmax[i];
      while ((v=v->parent) && !v->check_interior(pt)) {};//sube
      if (v) {
        v=v->cajap(pt); // baja
        while (v->prof>vin->prof) v=v->parent;// sube a la misma profundidad
        v->closest_en_caja(p,nc,d,e);
        if (d<=e) return true;
        if (_test_mayor(v,p,pt,d,e,nc,i+1)) return true;
      }
    }
    pt[i]=p[i];
  }
  return false;
}

static bool _test_menor
(caja *vin, const punto &p, punto &pt,
 double &d, double e, int &nc, int i0){
  int i,dim=vin->o->dim;
  for (i=i0;i<dim;i++){
    caja *v=vin;
    while (v && p[i]-d<v->pmin[i]){
      pt[i]=v->pmin[i]-(v->pmax[i]-v->pmin[i])*epsilon;
      while ((v=v->parent) && !v->check_interior(pt)) {}; // sube
      if (v) {
        v=v->cajap(pt); // baja
        while (v->prof>vin->prof) v=v->parent;// sube a la misma profundidad
        v->closest_en_caja(p,nc,d,e);
        if (d<=e) return true;
        if (_test_menor(v,p,pt,d,e,nc,i+1)) return true;
      }
    }
    pt[i]=p[i];
  }
  return false;
}

// p va en este caja
// Para buscar vecino,no se si es mas rapido subir hasta encontrar
// en que caja esta o empezar desde el mas grande para abajo
// sobretodo si hay mucha frontera
// Lo hago asi para abaratar lo mas caro (3d con interiores)
int caja::closest(int ix, double &d, double e){
  const _OCTREE_PTARRAY &plist=*(o->plist);
  const punto &p=plist[ix];
  int nc=otro(ix); if (nc==-1) return -1;
  d=p.distanciac(plist[nc]);
  closest_en_caja(p,nc,d,e);
  if (d<=e) return nc;
  // busca en los cajaes vecinos
  punto pt(p);
  if (_test_mayor(this,p,pt,d,e,nc,0)) return nc;
  pt=p;
  if (_test_menor(this,p,pt,d,e,nc,0)) return nc;
  return nc;
}

// p no necesariamente va en este caja
int caja::closest(const punto &p, double &d){
  int nc=otro(-1); if (nc==-1) return nc; // no hay otro
  d=p.distanciac((*(o->plist))[nc]);
  closest_en_caja(p,nc,d,0);
  // busca en los cajaes vecinos
  punto pt(p);
  if (_test_mayor(this,p,pt,d,0,nc,0)) return nc;
  pt=p;
  if (_test_menor(this,p,pt,d,0,nc,0)) return nc;
  return nc;
}

int octree::closest(punto p, double &d) const{
  // si esta fuera reemplaza por un punto sobre la caja
  for (int i=0;i<dim;i++){
    if (p[i]<pmin[i])  p[i]=pmin[i];
    if (p[i]>=pmax[i]) p[i]=pmax[i]-(pmax[i]-pmin[i])*epsilon;
  }
  caja &t=*(child->cajap(p));
  d=MAXREAL;
  return t.closest(p,d);
}

void octree::add(int ix, int &nc, bool &puesto, double e){
  const punto &p=(*plist)[ix];
  puesto=false; nc=-1;
  if (!check_interior(p)) return;
  caja &t=*(child->cajap(p));
  double d=MAXREAL;
  nc=t.closest(ix,d,e);
  if (d<e) return;
  puesto=true; t.add(ix);
}

void octree::add_no_rep(int ix, int &nc, bool &puesto, double e){
  const punto &p=(*plist)[ix];
  puesto=false; nc=-1;
  if (!check_interior(p)) return;
  caja &t=*(child->cajap(p));
  double d=MAXREAL;
  t.closest_en_caja(p,nc,d,e);
  if (d<e) return;
  if (nc==-1) nc=t.otro(ix);
  puesto=true; t.add(ix);
}

// elimina el nodo ix
void octree::remove(int ix){
  caja *t=child->cajap((*plist)[ix]);
  int i=t->len-1;
  while (i>=0&&t->ip[i]!=ix) i--;
  if (i==-1) return; // no estaba
  for (;i<t->len-1;i++) t->ip[i]=t->ip[i+1];
  t->len--;
  // si len=0 no hay ip
  if (t->len==0) {delete[] t->ip; t->ip=0;}
  // arregla los childs del 1er parent
  if ((t=t->parent)) t->len--; else return;
  if (t->len<maxlen){// chupa los puntos y elimina los child
    t->ip=new int [maxlen];
    t->len=0;
    for (int j=0;j<nchild;j++){
      caja &ch=t->child[j];
      for (i=0;i<ch.len;i++) t->ip[t->len++]=ch.ip[i];}
    delete [] t->child; t->child=0;}
  // arregla len en los parent superiores
  while ((t=t->parent)) t->len--;
}

// intercambia los indices de n1 y n2 (por un swap en la lista de nodos)
bool octree::swap_n(int in1, int in2){
  caja *t1=child->cajap((*plist)[in1]);
  caja *t2=child->cajap((*plist)[in2]);
  if (!t1||!t2) return false;
  if (t1==t2) return true;
  int i;
  for (i=0;i<t1->len;i++){
    if (t1->ip[i]!=in1) continue;
    t1->ip[i]=in2; break;
  } if (i==t1->len) return false; //no estaba
  for (i=0;i<t2->len;i++){
    if (t2->ip[i]!=in2) continue;
    t2->ip[i]=in1; break;
  } if (i==t2->len) return false; //no estaba
  return true;
}

#undef _OADEBUG
