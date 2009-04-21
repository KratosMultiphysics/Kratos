// lee archivos con mallas

// OJO: sscanf da floats y no doubles

#include <fstream>
#include <iomanip> // setw
#include <cstdio> // sscanf
#include <cstring> // memcmp strcat strcpy strtok

#include "case.h" // todo se pasa a minusculas para comparar
//#define temporizar
#include "tiempo.h"
#include "malla.h"

using namespace std;

#define _maxlinea 512

#define _getlinea {if (!getlinea()) return malformed();}
#define _chktoken {if (!token||!token[0]) return malformed();}
#define _gettoken {if (!getlinea()) return malformed();\
  token=strtok(linea,seps); _chktoken;}
#define _nexttoken (token=strtok(NULL,seps))
#define _tokenpoint(start){\
  token=strtok(start,seps);_chktoken;\
  pt[0]=atof(token);\
  token=strtok(NULL,seps);_chktoken;\
  pt[1]=atof(token);\
  if (hayz){\
    token=strtok(NULL,seps);_chktoken;\
    pt[2]=atof(token);\
  }\
}

#define _es(x) !strcmp(linea,x)
#define _noes(x) strcmp(linea,x)

#define _dxfpoint(pt) {if (!dxfpoint(pt)) return false; set=true;}
#define _getdxf {_getlinea; code=atoi(linea); _getlinea; int llen=strlen(linea); if (llen && linea[llen-1]=='\r') linea[llen-1]='\0'; }
#define _getreal(x) {x=atof(linea); _getdxf;}
#define _setreal(c,x) {if (code==c) {x=atof(linea); _getdxf; set=true;}}
#define _setint(c,x) {if (code==c) {x=atoi(linea); _getdxf; set=true;}}
#define _setstring(c,x) {if (code==c) {strcpy(x,linea); _getdxf; set=true;}}

#define _N  const char* arch,array1<nodo> &nl
#define _NE const char* arch,array1<nodo> &nl,array1<elemento> &el

static ifstream *fp=0;
static char
  linea[_maxlinea],
  *token,
  seps[]=" ;,\t\r",
  *stop,
  comment[_maxlinea]="";
static int code;
static size_t comment_len=0;
static malla *m; // para que las funciones static sepan a quien enchufarle el error

static void cierra() {fp->close(); delete fp; fp=0;}

static bool abre
 (const char *filename, const char *trail, bool test=true, bool bin=false) {
  char fn[_max_file_len];strcpy(fn,filename);strcat(fn,".");strcat(fn,trail);
  if (fp)
    cierra();
  if (!bin) fp=new ifstream(fn);
  else fp=new ifstream(fn,ios::in|ios::binary);
  if (!(fp->is_open())) {
    delete fp; fp=0;
    if (test) m->add_error(No_File);
    return false;
  }
  return true;
}

static bool malformed() {
  cierra();
  m->add_error(Bad_Format);
  _savetime(malformed); 
  return false;
}

static bool getlinea() {
  if (fp->eof()) return false;
  fp->getline(linea,255); minusc(linea);
  while (comment_len&&!strncmp(linea,comment,comment_len)) {
    fp->getline(linea,255); minusc(linea);
  }
  return true;
}

static void set_comment(const char* c){
  if (!c) {comment_len=0; comment[0]=0;}
  int len=strlen(c);
  if (len>=_maxlinea) return; else comment_len=len;
  strcpy(comment,c);
}

// busca un target en un archivo, al principio de la linea
// target y limit tienen que venir en minusculas (_getlinea->min)
// en linea vuelve la linea leida que ontiene el string
static bool search(const char* target,const char* limit=0,bool rewind=false)
{
  if (rewind) fp->seekg(0,ios::beg);
  size_t len,tlen=strlen(target),llen=(limit)? strlen(limit) :~0;
  while(1){
    _getlinea;len=strlen(linea);
    if (linea[len-1]=='\r') linea[--len]='\0';
    if (len>=tlen&&!memcmp(linea,target,tlen)) break;
    if (len>=llen&&!memcmp(linea,limit,llen)) return false;
  }
  return true;
}

//====================================
// lee puntos
// archivo .xyz o .xy
// solo una lista de coordenadas hasta el eof (no hay cantidad)
bool malla::lee_xy(_N){
  set_comment("");
  bool hayz=false; // para leer los tokens
  m=this; nl.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;
  nodo pt;

/*
  ifstream &f=*fp;
  // de esta forma no lee separadores
  while(!f.eof()){
    f >> pt[0] >> pt[1] >> ws;
    nl+=pt;
  }
*/
  while (getlinea()&&linea[0]) {
    _tokenpoint(linea);
    if (_nexttoken) pt.f=atoi(token)&fmask;
    nl+=pt;
  }

  cierra();
  _savetime(lee_xy);
  return true;
}

bool malla::lee_xyz(_N){
  set_comment("");
  bool hayz=true; // para leer los tokens
  m=this; nl.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;
  nodo pt;
  while (getlinea()&&linea[0]) {
    _tokenpoint(linea);
    if (_nexttoken) pt.f=atoi(token)&fmask;
    nl+=pt;
  }
  cierra();
  _savetime(lee_xyz);
  return true;
}

// conectividades y nodos por separado (matlab)
bool cascara::agrega_con(const char *stipo, const char *arch){
  set_comment("");
  m=parent;
  char filename[_max_file_len];
  if (arch&&arch[0]) strcpy(filename,arch); else strcpy(filename,parent->nombre);
  char *extaddr=ext_begin(filename); if (extaddr) *(extaddr-1)=0;
  if (!abre(filename,".con",true)) return false;
  _initime;
  elemento e1(elemento::etipo(stipo)); int nv=e1.nv();
  int old_base=elemento::io_base(); elemento::io_base(1); // .con es base 1
  int i;
  // nodos previos (si hay)
  for (i=0;i<n.len;i++) parent->n[n[i]].f.set(flag1); // flag1 -> puesto en n
  // en la primera linea me fijo si hay o no flags  
  _gettoken; int count=1; while(_nexttoken) count++;
  if (count==nv) elemento::flag_mask(0);
  if (count==nv+1) {elemento::flag_mask(fmask); hayfe=true;}
  else return malformed();
  fp->seekg(0,ios::beg); // rewind
  // lee los elementos
  while (!fp->eof()){
    (*fp) >> e1 >> ws;
    e+=e1;
    for (i=0;i<nv;i++) {
      if (e1[i]<0||e1[i]>parent->n.len) return malformed();
      if (parent->n[e1[i]].f.noes(flag1)) n+=e1[i];
    }
  }
  for (i=0;i<n.len;i++) parent->n[n[i]].f.reset(flag1);
  cierra();
  elemento::io_base(old_base); // restaura
  _savetime(agrega_con);
  return true;
}


//====================================
// Lee nodos y elementos
//
// La separacion entre numeros o palabras de una linea puede ser " ;,\t"
//
// Nodos:
// Hay una linea que indica el formato de entrada, por ej:
// 2576 Nodos # x y z h f v
// El primer nro es obligatorio y es la cantidad de nodos, el resto puede
// no estar, en cuyo caso se supone equivalente a x y z
// # indica que hay numeros de nodo, si no esta se supone numerado desde 1
// x e y son obligatorios, z no.
// h indica que hay un double que es h del nodo
// f es un flag binario opcional a gusto del usuario, que se conserva.
// v es un real asiciado al nodo
// Se puede poner base0 base1 base=0 base=1 b0 b1 b=0 o b=1
// Puede faltar la palabra Nodos, o decir cualquier otra cosa
// Las coordenadas vienen separadas por espacios o tabs
//
// Elementos:
// puede no haber elementos
// se dividen en grupos por tipo, hay una linea que dice
// cantidad, tipo y puede tener # y/o f, ej: 123 triangulos # f
// al numero de # no se le da bola
// tipo puede ser cualquier e_tipo definido en elemento.h
// despues vienen los nodos (y caras) de acuerdo al tipo
//

bool malla::lee_dat(_NE){
  set_comment("");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;
  ifstream &f=*fp;

  int qn,qe,i,j,nv,in=-1,ie=-1;

  // formato en que vienen los nodos
  _gettoken;
  sscanf(token,"%d",&qn); // cantidad de nodos
  bool hayn=false,hayz=false,/*hayh=false,*/hayf=false; // hayh y hayv de malla
  hayh=hayv=false;
  int numbase=1; //default
  pline tdat(8); int dat;
  // 0:nro 1:x 2:y 3:z 4:h 5:f 6:v
  while(_nexttoken){
    if (*token=='#') {hayn=true; tdat+=0; continue;}
    if (!strcmp(token,"xy")) {tdat+=1; tdat+=2; continue;}
    if (!strcmp(token,"xyz")) {hayz=true; tdat+=1; tdat+=2; tdat+=3; continue;}
    if (*token=='x') {tdat+=1; continue;}
    if (*token=='y') {tdat+=2; continue;}
    if (*token=='z') {hayz=true; tdat+=3; continue;}
    if (*token=='h') {hayh=true; tdat+=4; continue;}
    if (*token=='f') {hayf=true; tdat+=5; continue;}
    if (*token=='v' && *(token+1)!='e') // no confundir v con vertex
      {hayv=true; tdat+=6; continue;}
    if (*token=='b') {
      if      (!strcmp(token,"base0"))  numbase=0;
      else if (!strcmp(token,"base=0")) numbase=0;
      else if (!strcmp(token,"b0"))     numbase=0;
      else if (!strcmp(token,"b=0"))    numbase=0;
      else if (!strcmp(token,"base1"))  numbase=1;
      else if (!strcmp(token,"base=1")) numbase=1;
      else if (!strcmp(token,"b1"))     numbase=1;
      else if (!strcmp(token,"b=1"))    numbase=1;
      continue;
    }
  }
  if (!tdat.len) {tdat+=1; tdat+=2; tdat+=3;} // default x y z
  // lee los nodos
  pline index(qn+1/*,-1*/);
  nodo pt;
  _initime; // nodos
//#pragma message("<---------------- ***   OJO, CUSTOMIZADO!!!  ***\n")
//double factorp=.0001/RAND_MAX;
  for (i=0;i<qn;i++) {
    pt.ini();
    j=0; _gettoken;
    do {
      dat=tdat[j];
      if (dat==0) {in=atoi(token); if (in<0) {in=-in-1; pt.f.set(n_h);}}
      else if (dat==1) pt[0]=atof(token)/*+factorp*rand()*/;
      else if (dat==2) pt[1]=atof(token)/*+factorp*rand()*/;
      else if (dat==3) pt[2]=atof(token)/*+factorp*rand()*/;
      else if (dat==4) set_min(hmin,(pt.h=atof(token)));
      else if (dat==5) pt.f.set(atoi(token)&fmask); // set por n_h
      else if (dat==6) set_min_max(nvmin,nvmax,(pt.v=atof(token)));
    } while(_nexttoken&&++j<tdat.len);
    if (j<tdat.len-1){
      malformed();
      add_error(Near_Node);
      add_error(Max(i,in));
      return false;
    }
    if (!hayn) in=nl.len+numbase; // si no hay num empieza de 1
    index.resize(in+1); if (index.len<in+1) index.len=in+1;/*for (k=index.len;k<in+1;k++) index+=-1;*/
    index[in]=nl+=pt;
  }
  _savetime(lee_nodos);

  f >> ws; // (witespace) leee el enter por si no hay elementos

  // elementos
  char stipo[20];
  int old_base=elemento::io_base(); elemento::io_base(0); // hay index (va en 0)
  elemento ei; e_tipo tipo;
  _initime;
  while (!f.eof()){
    // formato en que vienen los elementos
    _gettoken; sscanf(token,"%d",&qe); // cantidad de elementos
    _nexttoken; sscanf(token,"%s",stipo); // tipo
    tipo=ei.tipo(stipo);
    if (tipo==e_indefinido){
      char msg[255];
      sprintf(msg,"%s: \"%s\"",Unknown_Element,stipo);
      add_warning(msg);
      for (i=0;i<qe;i++) _getlinea; // lee los desconocidos
    }
    hayn=false; elemento::flag_mask(0);
    while(_nexttoken){
      if (*token=='#') {hayn=true; continue;}
      if (*token=='f') {elemento::flag_mask(fmask); continue;}
    }

    el.resize(el.len+qe);
    for (i=0;i<qe;i++){
      if (f.eof()) {qe=i;break;}
      if (hayn) f >> ie; // de todos modos no le doy bola
      if (f.eof()) {qe=i;break;}
      f >> ei;
      nv=ei.nv(); for (j=0;j<nv;j++) {
        in=ei[j]; if (in<0||in>index.len||index[in]==-1) {
          malformed();
          add_error(Near_Element);
          add_error(Max(ie,i));
          return false;
        }
        ei[j]=index[in];
      }
      el+=ei;
    }
    f >> ws; // (witespace) leee el enter

    // los poliedros de 4 nodos los pone como tetraedros
    if (tipo==e_poliedro)
      for (i=el.len-qe;i<el.len;i++) el[i].p2t();
  }
  _savetime(lee_elementos);
  cierra();
  elemento::io_base(old_base);
  _savetime(lee_dat);
  return true;
}

// formato gid
bool malla::lee_msh(_NE){
  set_comment("");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,false)&&!abre(arch,ext,true)) return false;
  _initime;
  char stipo[20];
  int i,j,dim,nv=-1;
  double coord;
  punto pt;

  // encabezado
  do{
    while (getlinea()&&linea[0]) {
      if (linea[0]=='#') continue;
      token=strtok(linea,seps);
      if (token&&!strcmp("coordinates",token)) break;
      if (token&&strcmp("mesh",token)) continue;
      // hay elementos
      _nexttoken; while (token){
        if (!strcmp("dimension",token)){
          _nexttoken; if (token[0]=='=') _nexttoken; dim=atoi(token);
        }
        else if (!strcmp("elemtype",token)){
          _nexttoken; strcpy(stipo,token);
        }
        else if (!strcmp("nnode",token)){
          _nexttoken; if (token[0]=='=') _nexttoken; nv=atoi(token);
        }
        _nexttoken;
      }
      break;
    }

    // nodos
    _initime;
    //labels de nodos no secuenciales
    pline index(10000); // no se la cantidad, aloco para 10K
    while (!fp->eof()){
      if (token&&!strcmp("coordinates",token)) break;
      _getlinea;
      if (linea[0]=='#') continue;
      token=strtok(linea,seps);
    }
    if (token&&!strcmp("coordinates",token)) while (getlinea()&&linea[0]) {
      if (linea[0]=='#') continue;
      token=strtok(linea,seps); // numero
      if (!strcmp("end",token)) break;
      i=atoi(token); index.resize(i+1); if (index.len<i+1) index.len=i+1;
      j=0; while (j<3&&(token=strtok(NULL,seps))) {
        coord=atof(token);
        pt[j]=coord;
        j++;
      }
      index[i]=nl+=pt;
    }
    _savetime(lee_nodos);

    // elementos
    // segun GID: linear triangle ?????? (averiguar)
    elemento e1((nv==2)? "segmento" : stipo);
    _initime;
    while (!fp->eof()){
      if (token&&!strcmp("elements",token)) break;
      _getlinea;
      if (linea[0]=='#') continue;
      token=strtok(linea,seps);
    }
    if (token&&!strcmp("elements",token))while (getlinea()&&linea[0]) {
      token=strtok(linea,seps);  // numero
      if (!strcmp("end",token)) break;
      j=0; while (j<nv&&(token=strtok(NULL,seps))) {
        i=atoi(token);
        e1[j]=index[i];
        j++;
      }
      if ((token=strtok(NULL,seps))) { // flag (material)
        i=atoi(token); if (i<=_FREE_FLAGS) {
          i=1<<i; // directamente (aunque no haya material 0) para no confundir
          e1.f=i;
          for (j=0;j<nv;j++) nl[e1[j]].f.set(i);
        }
      }
      else e1.f=0;
      el+=e1;
    }
    _savetime(lee_elementos);
  } while (!fp->eof());

  cierra();
/*
  // valores (implementado para escalares)
  if (!abre(arch,ext,false)) return true;

  // lee la linea del encabezado
  while (getlinea()&&linea[0]=='#') {};

  for (i=0;i<nl.len;i++){
    if (fp->eof()) return true;// (mas o menos!!)
    while (getlinea()&&linea[0]=='#') {};
    token=strtok(linea,seps);// numero de nodo
    token=strtok(NULL,seps);// valor
    nl[i].v=atof(token); setear nvmin y nvax
  }

  cierra();
*/
  _savetime(lee_gid);
  return true;
}

// formato samcef
// esto es viejo
bool malla::lee_mai(_NE){
  set_comment("");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;
  char x[20], y[20], z[20];

  int i,ix,qn,qe,base;
  punto pt;
  e_tipo tipo;

  // encabezado
/*
! Number of nodes : 413
! Number of elements : 425
*/
  _getlinea;
  while (memcmp("! number of nodes :",linea,19)) _getlinea;
  sscanf(linea,"! number of nodes : %d",&qn);
  while (memcmp("! number of elements :",linea,22)) _getlinea;
  sscanf(linea,"! number of elements : %d",&qe);

  //labels de nodos no secuenciales
  pline index(qn+1);

  // nodos
/*
.NOEUDS
 I 1 X 10  10  0
 I 2 X 9.44444444  9.44444444  0
*/
  if (!search(".noeuds","return")) return malformed();
  for (i=0;i<qn;i++){
    _getlinea; sscanf(linea," i %d x %s %s %s",&ix, x, y, z);
    pt[0]=strtod(x,&stop);pt[1]=strtod(y,&stop);pt[2]=strtod(z,&stop);
    index.resize(ix+1); if (index.len<ix+1) index.len=ix+1;
    index[ix]=nl+=pt;
  }

  // elementos
/*
.MAILLES ###############################ojo, hay un espacio al final
 I 1 AT 1 0 N 45 3 4 46
 I 2 AT 1 0 N 46 4
*/
  if (!search(".mai","return",true)) {
    _savetime(lee_samcef);
    cierra(); 
    return true;
  }
  pline elnodes(10);
  for (i=0;i<qe;i++){
    _getlinea;
    token=strtok(linea,seps);  // I
    while ((token=strtok(NULL,seps))&&strcmp(token,"n"));
    elnodes.clean(); base=0;
    while ((token=strtok(NULL,seps))) {
      sscanf(token,"%d",&ix);
      if (ix>0) elnodes+=index[ix];
      if (ix==0) base=elnodes.len;
      // negativo => nada (son puntos intermedios)
    }
    if (elnodes.len==3) tipo=e_triangulo;
    else if (elnodes.len==4) {
      if (base==3) tipo=e_tetraedro;
      else tipo=e_cuadrilatero;
    }
    else if (elnodes.len==6&&base==3) tipo=e_wedge;
    else if (elnodes.len==8&&base==4) tipo=e_cubo;
    else tipo=e_poligono; // no hay poliedros.
    el+=elemento(tipo,elnodes.vertex);
  }
  cierra();
  _savetime(lee_samcef);
  return true;
}

// formato OOFELIE
// base 1
// Parcial!!
// El archivo tene nodos con FixationSet
// y elementos con Propelem
// los elementos con mas de una Propelem estan repetidos!!
bool malla::lee_e(_NE){
  set_comment("//");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;

  int i;
  punto pt;

  const char 
    ntarget[]="positionset",
    ptarget[]="propelem",
    etarget[]="elementset";

  // lee los nodos
  if (!search(ntarget)) return malformed();
  getlinea(); // {
  bool hayz=true;
  while (getlinea()&&linea[0]&&linea[0]!='}') {
    token=strtok(linea,seps); _chktoken; // numero
    _tokenpoint(NULL); nl+=pt;
  }
  int qn=nl.len;

  // propiedades
  // busca y cuenta
  if (!search(ptarget,etarget,true)) return malformed();
  int qp=1;
  while (search(ptarget,etarget)) qp++;
  // holder
  class prop{
   public:
    char p[256];
    e_tipo tipo;
    prop() {} 
    const prop &operator=(const char* c) {
      strcpy(p,c);
      if      (strstr(p,"hexa8"))  tipo=e_cubo;
      else if (strstr(p,"quad4"))  tipo=e_cuadrilatero;
      else if (strstr(p,"tetra4")) tipo=e_tetraedro;
      else if (strstr(p,"tri3"))   tipo=e_triangulo;
      else                         tipo=e_indefinido;
      return *this;
    }
    bool operator ==(const char* c) const {
      int lp=strlen(p),lc=strlen(c);
      if (lp!=lc) return false;
      return (!strncmp(p,c,lp));
    }
    bool operator !=(const char* c) const {return (!(*this==c));}
  };
  prop *pr=new prop[qp];
  // lee 
  for (i=0;i<qp;i++) {
    search(ptarget,etarget,!i); // rebobina el 1ro
    token=strtok(linea," ("); _chktoken;
    token=strtok(NULL," ("); _chktoken;
    pr[i]=token; if (pr[i].tipo==e_indefinido) return malformed();
  }

  // elementos
  // si hay elementos frontera deben venir al final
  if (!search(etarget,0,true)) return malformed();
  getlinea(); // {
  pline emap; cpliner comun;
  array1<pline> eldenod(qn); // para verificar elementos frontera
  elemento ei; int nv,dim0=-1,ie,iec;  // -1 para evitar warning de iniciaclizacion
  while (getlinea()&&linea[0]&&linea[0]!='}') {
    token=strtok(linea,seps); _chktoken; // numero
    // propiedad y tipo
    _nexttoken; 
    i=0; while(i<qp&&pr[i]!=token) i++; if (i==qp) return malformed();
    emap+=i;
    ei.tipo(pr[i].tipo); nv=ei.nv();
    // nodos
    for (i=0;i<nv;i++) {_nexttoken; ei.n[i]=atoi(token)-1;}
    ie=el.len;
    // verifica si es base o frontera
    if (ie){
      if (ei.dim()<dim0){
        if (ei.tipo()==e_cuadrilatero){
          comun=eldenod[ei[0]].inters(eldenod[ei[2]]);
        }
        else if (ei.tipo()==e_triangulo){
          comun=(eldenod[ei[0]].inters(eldenod[ei[1]])).inters(eldenod[ei[2]]);
        }
        else if (ei.tipo()==e_segmento){
          comun=eldenod[ei[0]].inters(eldenod[ei[1]]);
        }
        else comun.clean();
        if (comun.len==1) { // repetido
          iec=comun[0]; const elemento &ec=el[iec];
          i=ec.icara(ei); if (i==-1) return malformed();
          emap+=i; emap+=iec; 
          continue; // no agrega el elemento
        }
      }
    }
    else dim0=ei.dim();
    emap+=-1;emap+=ie;
    for (i=0;i<nv;i++) eldenod[ei[i]]+=ie;
    el+=ei;
  }
  int qe=emap.len/3;  
  
  // genera un archivo de traducción
  char tfn[_max_file_len]; strcpy(tfn,arch); strcat(tfn,".tmp");
  ofstream tf(tfn); 
  tf << qp << " Propelem: #,prop" << endl;
  for (i=0;i<qp;i++) tf << i+1 << ',' << pr[i].p << endl;
  tf << qe << " Elements: #,prop,face,from" << endl;
  for (i=0;i<qe;i++) 
    tf << i+1           << ',' 
       << emap[3*i  ]+1   << ',' 
       << emap[3*i+1]+1 << ',' 
       << emap[3*i+2]+1 << endl; 
  tf.close();

  cierra();
  _savetime(lee_oofelie);
  return true;
}

// formato NASTRAN
// los borrados se eliminan en lee()
bool malla::lee_nas(_NE){
  set_comment("");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;

  int i;
  pline index(int(1e6));
  nodo pt; pt.f.set(n_borrado); // no guardar nodos intermedios en cuadraticos
  elemento e3(e_triangulo);

  // encabezado
  do{_gettoken} while (strcmp(token,"grid*"));

  // nodos
  //GRID*,1,0,-0.27189035E-01,-0.11801019E-01,0.36064133E-01,,,
  do{
    _nexttoken; i=atoi(token);
    _nexttoken;
    _nexttoken; pt[0]=atof(token);
    _nexttoken; pt[1]=atof(token);
    _nexttoken; pt[2]=atof(token);
    index.resize(i+1); if (index.len<i+1) index.len=i+1; index[i]=nl+=pt;
    _gettoken;
  }while (!strcmp(token,"grid*"));

  // posible intermedio
  while (memcmp(token,"ctria",5)) _gettoken;

  // elementos
  //CTRIA6*,1,1,1,4,3,111179,111171,111180,
  //CTRIA6*,2,1,3,2,1,111175,111181,111180,
  //CTRIA6*,3,1,11,5,6,111163,111164,111156,
  //CTRIA6*,4,1,12,17,13,111126,111120,111127,
  do{
    _nexttoken; i=atoi(token);
    _nexttoken;
    _nexttoken; e3[0]=i=index[atoi(token)]; nl[i].f.reset(n_borrado);
    _nexttoken; e3[1]=i=index[atoi(token)]; nl[i].f.reset(n_borrado);
    _nexttoken; e3[2]=i=index[atoi(token)]; nl[i].f.reset(n_borrado);
    el+=e3;
    _gettoken;
  }while (!memcmp(token,"ctria",5));

  cierra();
  _savetime(lee_nastran);
  return true;
}

// formato ANSYS
// los borrados se eliminan en lee()
bool malla::lee_ans(_NE){
  set_comment("");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;

  int i;
  pline index(int(1e6));
  nodo pt; pt.f.set(n_borrado); // no guardar nodos intermedios en cuadraticos
  elemento e3(e_triangulo);

  // encabezado
  do{_gettoken} while (strcmp(token,"n"));

  // nodos
  //N,1,-0.271890E-01,-0.118010E-01,0.360641E-01
  //N,2,-0.265244E-01,-0.112842E-01,0.358800E-01
  do{
    _nexttoken; i=atoi(token);
    _nexttoken; pt[0]=atof(token);
    _nexttoken; pt[1]=atof(token);
    _nexttoken; pt[2]=atof(token);
    index.resize(i+1); if (index.len<i+1) index.len=i+1;
    index[i]=nl+=pt;
    _gettoken;
  }while (!strcmp(token,"n"));

  // posible intermedio
  while (strcmp(token,"en")) _gettoken;

  // elementos
  //EN,1,1,4,3,3,111179,111171,3,111180,
  //EN,2,3,2,1,1,111175,111181,1,111180,
  //EN,3,11,5,6,6,111163,111164,6,111156,
  do{
    _nexttoken; i=atoi(token);
    _nexttoken; e3[0]=i=index[atoi(token)]; nl[i].f.reset(n_borrado);
    _nexttoken; e3[1]=i=index[atoi(token)]; nl[i].f.reset(n_borrado);
    _nexttoken; e3[2]=i=index[atoi(token)]; nl[i].f.reset(n_borrado);
    el+=e3;
    _gettoken;
  }while (!strcmp(token,"en"));

  cierra();
  _savetime(lee_ansys);
  return true;
}

// formato ls-dyna
bool malla::lee_key(_NE){
  set_comment("$");
  m=this; nl.clean(); el.clean(); nl.resize(1000); el.resize(1000);
  int maxindex=-1; pline index(1000); // indices de los nodos
  if (!abre(arch,ext,true)) return false;
  _initime;

  // encabezado
  search("*",0,false);
  while (linea[0]=='*'){
    if (_es("*end")) break;
    else if (_es("*element_solid")){
      // tetraedros y wedges y hexaedros
      //*ELEMENT_SOLID
      //$#   eid     pid      n1      n2      n3      n4 ........
      //       1       1      13      12    1503      14
      //       2       1    1337    1318    1319    1338
      int pid,ixe[8];
      e_tipo tipo;
      _gettoken; do{
        //no me interesa el nro
        _nexttoken; pid=atoi(token)-1; // flag
        _nexttoken; ixe[0]=atoi(token);
        _nexttoken; ixe[1]=atoi(token);
        _nexttoken; ixe[2]=atoi(token);
        _nexttoken; ixe[3]=atoi(token);
        _nexttoken; ixe[4]=atoi(token);
        _nexttoken; ixe[5]=atoi(token);
        _nexttoken; ixe[6]=atoi(token);
        _nexttoken; ixe[7]=atoi(token);
             if (ixe[5]==ixe[4]) tipo=e_tetraedro;
        else if (ixe[6]==ixe[5]) tipo=e_wedge;
        else                     tipo=e_cubo;
        el+=elemento(tipo,ixe); if (pid<=_FREE_FLAGS) el.last().f=1<<pid;
        _gettoken;
      }while (token[0]!='*');
    }
    else if (_es("*element_shell")){
      // triangulos y cuadrilateros (facundo lo usa con segmentos)
      //*ELEMENT_SHELL
      //$#   eid     pid      n1      n2      n3      n4
      //       1       1      13      12    1503      14
      //       2       1    1337    1318    1319    1338
      int pid,ixe[4];
      e_tipo tipo;
      _gettoken; do{
        //no me interesa el nro
        _nexttoken; pid=atoi(token)-1; // flag
        _nexttoken; ixe[0]=atoi(token);
        _nexttoken; ixe[1]=atoi(token);
        _nexttoken; ixe[2]=atoi(token);
        _nexttoken; ixe[3]=atoi(token);
             if (ixe[1]==ixe[2]) tipo=e_segmento;
        else if (ixe[2]==ixe[3]) tipo=e_triangulo;
        else                     tipo=e_cuadrilatero;
        el+=elemento(tipo,ixe); if (pid<=_FREE_FLAGS) el.last().f=1<<pid;
        _gettoken;
      }while (token[0]!='*');
    }
    else if (_es("*node")){
      // nodos
      // x, y o z pueden no estar
      // facu pone h despues de z
      // puede venir en formato fijo de 16 o separado por comas
      //*NODE
      //$#   nid               x               y               z      tc      rc
      //0         1         2         3         4         5
      //012345678901234567890123456789012345678901234567890123456789
      //       1     -7.74273157     -6.32851601
      //......
      //      13     -9.94987488     -0.99999905
      //      14     -9.94987488     -0.99999905      0.50000000
      int i,llen;
      nodo pt; if (_es("*node_rigid_surface")) pt.f.set(n_permanente);
      char rd[17];
      _getlinea; do{
        if (strchr(linea,',')){ // formato libre separado por comas
          token=strtok(linea,seps);
          i=atoi(token);
          _nexttoken; if (token&&token[0]) pt[0]=atof(token); else pt[0]=0;
          _nexttoken; if (token&&token[0]) pt[1]=atof(token); else pt[1]=0;
          _nexttoken; if (token&&token[0]) pt[2]=atof(token); else pt[2]=0;
          _nexttoken; if (token&&token[0]) pt.h=atof(token); else pt.h=0;
        }
        else { // formato fijo: puede no estar x y estar y
          llen=strlen(linea); token=&linea[0];
//		  if (llen && linea[llen-1]=='\r') linea[--llen]='\0';
          strncpy(rd,token,8); rd[8]=0; llen-=8; token+=8;
          i=atoi(rd);
          if (llen){
            strncpy(rd,token,16); rd[16]=0; llen-=16; token+=16;
            pt[0]=atof(rd);
          } else pt[0]=0;
          if (llen){
            strncpy(rd,token,16); rd[16]=0; llen-=16; token+=16;
            pt[1]=atof(rd);
          } else pt[1]=0;
          if (llen){
            strncpy(rd,token,16); rd[16]=0; llen-=16; token+=16;
            pt[2]=atof(rd);
          } else pt[2]=0;
          if (llen){
            strncpy(rd,token,16); rd[16]=0; llen-=16; token+=16;
            pt.h=atof(rd);
          } else pt.h=0;
        }
        nl+=pt; index+=i; set_max(maxindex,i);
        _getlinea;
      }while (linea[0]!='*');
	  llen=strlen(linea); if (linea[llen-1]=='\r') linea[llen-1]='\0';
    }
    else search("*",0,false);
  }

  if (el){
    pline anti_index(maxindex+1); anti_index.set(-1); //set asigna len
    int j,k;
    for(j=0;j<index.len;j++) anti_index[index[j]]=j;
    for(j=0;j<el.len;j++) for(k=0;k<el[j].nv();k++) el[j][k]=anti_index[el[j][k]];
  }

  cierra();
  _savetime(lee_key);
  return true;
}


// formato STL
bool malla::lee_stl(_NE){
  set_comment("");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,true,true)) return false;
  _initime;

  nodo pt; punto normal; float vl[12];
  elemento e3(e_triangulo); unsigned short fe;
  ifstream &f=*fp;
  char header[80];

  // encabezado
  f.read(header,80); if (f.eof()) return malformed();
  unsigned long int qe; f.read((char *)&qe,4);
  
  // cada triangulo tiene nodos :-(
  el.resize(qe); nl.resize(3*qe);
  for (unsigned long int i=0;i<qe;i++){
    if (f.eof()) return malformed();
    f.read((char *)vl,48);
//    normal[0]=vl[0]; normal[1]=vl[1]; normal[2]=vl[2];
    pt[0]=vl[3]; pt[1]=vl[4]; pt[2]=vl[5]; e3[0]=nl+=pt;
    pt[0]=vl[6]; pt[1]=vl[7]; pt[2]=vl[8]; e3[1]=nl+=pt;
    pt[0]=vl[9]; pt[1]=vl[10]; pt[2]=vl[11]; e3[2]=nl+=pt;
    f.read((char *)&fe,2); //e3.f=fe;
    el+=e3;
  }

  cierra();
  _savetime(lee_stl);
  return true;
}

//========================================================================
// lee DXF

//entra apuntando a un code y sale apuntando a un code con el ultimo code y linea leidos
static bool dxfpoint(nodo &pt){
  while (1) {
    if (code/10==1) break;
    _getdxf;
  }
  _getreal(pt[0]);
  _getreal(pt[1]);
  if (code/10==3) {_getreal(pt[2]);} else pt[2]=0;
  return true;
}

static int layer_flag(){
  while (1) {
    /////////////// ojo que no se desfase (handle layer o color pueden confundir)
    _getlinea; code=atoi(linea); // code
    _getlinea; // valor
    if (code==8) break;
    if (!code || code/10==1 || code==70) return 0;
  }
  // linea es el nombre de la layer
  if (*linea!='_') return 0;
  char *startf=linea+1;

  if (*startf=='h'||*startf=='H') return n_h; //malla de h

  // F3 f4 F05 6 06 13 ....
  if (*startf=='F'||*startf=='f') {startf++; if (!*startf) return 0;}
  if (*startf>='0'&&*startf<='9'){ // una o dos cifras
    if (startf[1]>='0'&&startf[1]<='9') startf[2]=0;
    else startf[1]=0;
    int k=atoi(startf); // atoi da 0 si no puede convertirlo
    if (k==0&&(*startf!='0'||(startf[1]!=0&&startf[1]!='0'))) return 0;
    if (k>=0) return (1<<k)&fmask;
  }
  return 0;
}

// Lee un DXF
// para los los limites cerrados el nodo inicial es igual al final
bool malla::lee_dxf(_NE){
  set_comment("");
  m=this; nl.clean(); el.clean();
  if (!abre(arch,ext,true)) return false;
  _initime;
  
  code=0; linea[0]=0;
  bool set;
  int i,j,k,lf;
  nodo pelm[4];
  nodo pt;
  char text[255],t0;

  while (_noes("entities")) _getdxf;

  _getdxf; while (1){
         if (_es("eof")) break;
    else if (_es("endsec")) break;
    else if (_es("line")){
      // lee linea (por ahora sin thickness)
      lf=layer_flag();
      elemento elm(e_segmento); elm.f=lf;
      j=0; _getdxf; while (code) {
        set=false;
        if (code/10==1) {_dxfpoint(pt); pelm[j++]=pt;}
        if (!set) _getdxf;
      }
      if (pelm[0].eq_c(pelm[1],1e-8)) continue; // long 0
      elm[0]=nl+=pelm[0]; elm[1]=nl+=pelm[1];
      if (lf) {nl[elm[0]].f.set(lf);nl[elm[1]].f.set(lf);}
      el+=elm;
    }
    else if (_es("3dface")){
      // lee 3DFACE
      lf=layer_flag();
      pt[2]=0;
      j=0; _getdxf; while (code) {
        set=false;
        if (code/10==1) {_dxfpoint(pt); pelm[j++]=pt;}
        if (!set) _getdxf;
      }
      // elimina repetidos
      for (i=j-1;i>0;i--) {
        for (k=i-1;k>=0;k--) if (pelm[i].eq_c(pelm[k],1e-8)) break;
        if (k<0) continue;
        // repetido
        for (k=i;k<j-1;k++) pelm[k]=pelm[k+1];
        j--;
      }
      if (j==1) continue; //todos repetidos
      elemento elm(
        (j==2)? e_segmento :
        (j==3)? e_triangulo :
                e_cuadrilatero);
      for (i=0;i<j;i++) elm[i]=nl+=pelm[i];
      elm.f=lf; el+=elm;
      if (lf) for (i=0;i<j;i++) nl[elm[i]].f.set(lf);
    }
    else if (_es("lwpolyline")){ // le saque el bulge
      // polilinea 2d
      lf=layer_flag();
      int nv=0,flag=-1;
      pt[2]=0;
      while (code!=10){
        set=false;
        _setint(90,nv); // cantidad de vertices
        _setint(70,flag); // flag
        if (!set) _getdxf;
      }
      bool closed=((flag&1)==1);
      int first=nl.len;
      for (i=0;i<nv;i++) {
        _dxfpoint(pt);
        if (!i||!pt.eq_c(nl.last(),1e-8)) nl+=pt;
      }
      elemento elm(e_segmento); elm.f=lf;
      for (i=first;i<nl.len-1;i++) {
        elm[0]=i; elm[1]=i+1; el+=elm;
      }
      if (closed) {// cerrada (duplica el 0)
        elm[0]=nl.len-1; elm[1]=first; el+=elm;
      }
      if (lf) for (i=first;i<nl.len;i++) nl[i].f.set(lf);
    }
    else if (_es("polyline")){
      // lee polilineas de acuerdo al tipo
      lf=layer_flag();
      _getdxf; while (code!=70&&code!=0) _getdxf;
      int flag=0; _setint(70,flag);
      if (!(flag&(16+64))){ // polilinea (2d o 3d)
        pt[2]=0;
        bool closed=((flag&1)==1);
        while (_noes("vertex")) _getdxf;
        int first=nl.len;
        while (_es("vertex")){
          _dxfpoint(pt);
          while (code!=0) _getdxf;
          if (nl.len==first||!pt.eq_c(nl.last(),1e-8)) nl+=pt;
        }
        elemento elm(e_segmento);  elm.f=lf;
        for (i=first;i<nl.len-1;i++) {
          elm[0]=i; elm[1]=i+1; el+=elm;
        }
        if (closed) {// cerrada (duplica el 0)
          elm[0]=nl.len-1; elm[1]=first; el+=elm;
        }
        if (lf) for (i=first;i<nl.len;i++) nl[i].f.set(lf);
      }
      else if (flag&16) { // lee una 3DMesh
        bool mc,nc; int m=-1,n=-1,nv;
        mc=((flag&1)==1);  // mclosed
        nc=((flag&32)==32); // nclosed
        while (_noes("vertex")) {
          set=false;
          _setint(71,m);
          _setint(72,n);
          if (!set) _getdxf;
        }
        nv=m*n;
        pline mapnod(nv);
        // lee los vertices
        for (i=0;i<nv;i++){
          while (_noes("vertex")) _getdxf;
          _dxfpoint(pt); mapnod+=nl+=pt;
          if (lf) nl.last().f.set(lf);
        }
        // genera los cuads
        int i1,j1; elemento elm(e_cuadrilatero); elm.f=lf;
        for (i=0;i<m;i++)
        {
          if (i==m-1) {i1=0;if (!mc) break;} else i1=i+1;
          for (j=0;j<n;j++)
          {
            if (j==n-1) {j1=0;if (!nc) break;} else j1=j+1;
            elm[0]=mapnod[i *n+j ];
            elm[1]=mapnod[i *n+j1];
            elm[2]=mapnod[i1*n+j1];
            elm[3]=mapnod[i1*n+j ];
            el+=elm;
          }
        }
      }
      else if (flag&64){ // lee una polyface
        int nodos=-1,elms=-1;
        while (_noes("vertex")) {
          set=false;
          _setint(71,nodos);
          _setint(72,elms);
          if (!set) _getdxf;
        }
        if (nodos<0) malformed(); // limite de 32767
        if (elms<0) malformed(); // limite de 32767
        pline mapnod(nodos); // porque puede no ser la unica
        // lee los vertices
        for (i=0;i<nodos;i++){
          while (_noes("vertex")) _getdxf;
          _dxfpoint(pt); mapnod+=nl+=pt;
          if (lf) nl.last().f.set(lf);
        }
        if (nodos==4&&elms==4){ // es un tetraedro
          elemento elm(e_tetraedro); elm.f=lf;
          for (i=0;i<4;i++) elm[i]=mapnod[i];
          el+=elm;
          while (_noes("seqend")) _getdxf;
        }
        else if (nodos==8&&elms==6){ // es un cubo
          elemento elm(e_cubo); elm.f=lf;
          for (i=0;i<8;i++) elm[i]=mapnod[i];
          el+=elm;
          while (_noes("seqend")) _getdxf;
        }
        else if (nodos==6&&elms==5){ // es un prisma
          elemento elm(e_wedge); elm.f=lf;
          for (i=0;i<6;i++) elm[i]=mapnod[i];
          el+=elm;
          while (_noes("seqend")) _getdxf;
        }
        else{ // lee los elms triangulos o cuadrilateros (no mas)
          int v1=-1,v2=-1,v3=-1,v4=-1;
          for (i=0;i<elms;i++){
            while (_noes("vertex")) _getdxf;
            v4=-999999;
            _getdxf; while (code) {
              set=false;
              _setint(71,v1);
              _setint(72,v2);
              _setint(73,v3);
              _setint(74,v4);
              if (!set) _getdxf;
            }
            elemento elm((v4==-999999)? e_triangulo : e_cuadrilatero); elm.f=lf;
            elm[0]=mapnod[abs(v1)-1]; elm[1]=mapnod[abs(v2)-1]; elm[2]=mapnod[abs(v3)-1];
            if (v4!=-999999) elm[3]=mapnod[abs(v4)-1];
            el+=elm;
          }
        }
      }
    }
    else if (_es("text")){
      int c72=0,c73=0;
      nodo pi;
      j=0; _getdxf; while (code) {
        set=false;
        if (code/10==1) _dxfpoint(pelm[j++]);
        _setint(72,c72);
        _setint(73,c73);
        _setstring(1,text);
        if (!set) _getdxf;
      }
      if (c73==0&&c72==0) pt=pelm[0]; // default lower left (puede no haber pelm[1])
      else if (c73==0&&(c72==3||c72==5)) pt=(pelm[0]+pelm[1])/2; // ali o fit
      else pi=pelm[1];
      t0=text[0];
      if (t0=='-'||t0=='+'||t0!='.'||(t0>=48&&t0<=57)){
        set_min_max(nvmin,nvmax,(pi.v=atof(text))); // 0 si no puede convertirlo
      }
      else
        continue;
      nl+=pi;
    }
    else if (_es("point")){
      pt.ini();
      pt.f=layer_flag();
      _getdxf; while (code) {
        set=false;
        if (code/10==1) _dxfpoint(pt);
        _setstring(8,text); // layer
        if (!set) _getdxf;
      }
      nl+=pt;
    }
    else _getdxf;
  }
  cierra();
  _savetime(lee_dxf);
  return true;
}


#undef _maxlinea
#undef _getlinea
#undef _chktoken
#undef _gettoken
#undef _nexttoken
#undef _tokenpoint
#undef _es
#undef _noes
#undef _dxfpoint
#undef _getdxf
#undef _getreal
#undef _setreal
#undef _setint
#undef _setstring
#undef _N
#undef _NE
