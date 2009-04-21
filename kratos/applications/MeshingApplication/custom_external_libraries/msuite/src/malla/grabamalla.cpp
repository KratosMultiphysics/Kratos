#include <ostream>
#include <iomanip> // setw
#include <cstring> // strcat strcpy
#include <cstdio> // sprintf remove rename

//#define temporizar
#include "tiempo.h"
#include "malla.h"

using namespace std;

//===============================
//    #defines utiles
//===============================

// endl hace flush, y por ahi para debug es mejor
//#ifdef _DEBUG
//  #define _EOL endl
//#else
  #define _EOL '\n'
//#endif

// width>=precision+7 (+n.ne+xxx)
#define _oi(x) setw(10) << x
#define _spacei "          "
#define _or(x) setprecision(15) << setw(25) << x
#define _spacer "                         "
#define _of(x) setprecision(4) << setw(12) << x
#define _spacef "            "


#define _OF(trail)\
  char fn[_max_file_len]; strcpy(fn,filename);\
  if (!fn[0]) return false;\
  strcat(fn,trail);\
  backup(fn);\
  ofstream f(fn);\
  if (!f||!f.is_open()) {\
    add_error(Not_Open);\
    return false;\
  }

#define _OF2(trail)\
  f.close();\
  strcpy(fn,filename);\
  if (!fn[0]) return false;\
  strcat(fn,trail);\
  backup(fn);\
  f.open(fn);\
  if (!f||!f.is_open()) {\
    add_error(Not_Open);\
    return false;\
  }

#define _close f.close(); return true;

#define _makename(trail) \
  char filename[_max_file_len]; \
  if (arch&&arch[0]) strcpy(filename,arch); \
  else strcpy(filename,nombre);\
  char *extaddr=ext_begin(filename); if (extaddr) *(extaddr-1)=0;\
  strcat(filename,trail);

#define _makename2(trail) \
  if (arch&&arch[0]) strcpy(filename,arch); \
  else {strcpy(filename,nombre);strcat(filename,trail);}


//=================================
// archivos .dat:
// qn # x y [z]
// ....
// qet tipo #                (element sets por tipo)
// nodo1, nodo2,..... nodo#n

bool malla::graba_dat(const char *arch) {
  _makename("");
  _OF(".dat");
  int i,j,k,ix,nc,base=elemento::io_base();
  f << n.len << " Nodes # x y";
  bool hayz=tipo.noes(m_planaxy);

  if (hayz)  f << " z";
  if (hayh)  f << " h";
  if (hayfn)  f << " f";
  if (hayv)  f << " v";
//  if (base==1)  f << " base=1"; no, porque hay numeros
  f << _EOL;
  punto::io_z(hayz); punto::o_precision(15);  // toda
  punto::o_width(25); punto::o_separator(" ");


  for (i=0;i<n.len;i++){
//#pragma message("<---------------- ***   OJO, CUSTOMIZADO!!!  ***\n")
//if (!n[i].e) continue;
    f << _oi(i+base) << n[i]/*.redondea(1e-2)*/;
    if (hayh) f << _or(n[i].h);
    if (hayfn) f << setw(15) << (n[i].f&fmask);
    if (hayv) f << _or(n[i].v);
//for (j=0;j<n[i].e.len;j++) f << ' ' << n[i].e[j];
    f << _EOL;
  }
  punto::io_reset();

  // element sets (debe haber correspondencia con los textos de elemento)
  pline settipo[e_ntipos]; int qet;
  bool print_vecinos=false; // para debug o para quien quiera vecinos
  elemento::flag_mask((hayfe)? fmask : 0); // para que no ponga flags comentar la linea
  if (print_vecinos&&!vecino) mk_vecino();
//mk_esferas(true);
  for (i=0;i<e.len;i++) {
//if (e[i].f.noes(flag3)) continue;
    settipo[e[i].itipo()]+=i;
  }
  for (j=1;j<e_ntipos;j++){ // los indefinidos (0) no salen
    qet=settipo[j].len;
    if (!qet) continue;
    f << qet << " " << e[settipo[j][0]].stipo() << " #";
    if (hayfe) f << " f";
    f << _EOL;
    for (i=0;i<qet;i++){
      ix=settipo[j][i];
//if (!e[ix].have(11949)&&!e[ix].have(12441)) continue;
      f << ix << "\t" << e[ix];
      if (print_vecinos){
        nc=e[ix].nc();
        for (k=0;k<nc;k++) f << "\t" << vecino[ix][k];
      }
//f << ' ' << ce[i] << ' ' << re[i] << ' ' << ve[i];
      f << _EOL;
    }
  }
/*
// nube de vecinos de vecinos
  array1<pline> nube;
  if (vecindad_de_nodo(nube,0,50)){
    _OF2(".nube");
    for (i=0;i<n.len;i++){
      const pline &vi=nube[i];
//      f << vi.len << " ";
      for (j=0;j<vi.len;j++) f << vi[j] << " ";
      f << _EOL;
    }
  }
*/
  _close;
}

// 0=<1e-10; 201=>1e10  (20 decadas)
// intermedio hay 10 por decada
static int hist_index(double d){
  if (d<1e-10)
    return 0;
  if (d>1e10)
    return 201;
  return int(floor(10*log10(d+1e-30)))+100;
}
static double rangestart(int i){
  if (i==0) return 0;
  if (i==201) return 1e10;
  return pow((double)10,double(i-1)/10-10);
}

//=================================
// calidad
//=================================
/*
  #         numero de elemento
  nv        cantidad de vertices
  nf        vertices en la frontera
  cf        caras en la frontera
  v         volumen
  hm        h medio de los nodos
  lhmin     arista/(h medio de sus nodos) minimo
  lhmax     arista/(h medio de sus nodos) maximo
  R3(v)/hm  raiz cubica de v / h medio
  calv      (maximo modulo gradiente funciones de forma * hmin)^-1
  dmin      minimo diedro (grados)
  180-dmax  maximo diedro (grados)
*/
bool malla::graba_calidad(const char *arch) {
  if (!mk_vecino()) return false;

  _initime;

  _makename("_cal");
  _OF(".worse");

  const char encabezado[]=
     "         # nv nf cf           v          hm       lhmin       lhmax    R3(v)/hm        calv        dmin    180-dmax\n";
//   1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
  f << encabezado;

  int elen=e.len;

  mk_esferas(true);
  mk_h_nn_med(true);

  int i,j,k,in0,in1,nv,nf,vf,dim;
  double ffv[64];
  punto Dfv[64];
  punto g_v;
  bool hayff;
  double
    vminm,
    h,hmin_e,hm,hmin_em,
    l,lmin,lminm,lmax,lmaxm,
    amin,aminm,amax,amaxm,
    calv,calvm,
    vh,vhminm,vhmaxm;

  int
    lmin_hist[202],lmax_hist[202],
    amin_hist[202],amax_hist[202],
    calv_hist[202],vh_hist[202];
  memset(lmin_hist,0,202*SZI);memset(lmax_hist,0,202*SZI);
  memset(amin_hist,0,202*SZI);memset(amax_hist,0,202*SZI);
  memset(calv_hist,0,202*SZI);memset(vh_hist,0,202*SZI);

  vminm=amaxm=aminm=calvm=vhminm=hmin_em=lminm=MAXREAL; lmaxm=0; vhmaxm=-MAXREAL;
  for (i=0;i<elen;i++) {
    elemento &ei=e[i]; nv=ei.nv(); dim=ei.dim(); if (dim<2) continue;
    for (vf=j=0;j<e[i].nc();j++) if (vecino[i][j]<0) vf++; // caras de frontera
    g_v=gv(ei);
    angulos(i,amin,amax); set_min(aminm,amin); amax=180-amax; set_min(amaxm,amax);
    hayff=fforma(i,g_v,ffv,Dfv);
    calv=0; hmin_e=lmin=MAXREAL; hm=lmax=0;
    for (nf=j=0;j<nv;j++) {
      in0=ei[j];
      const nodo &nj=n[in0];
      if (nj.f.es(n_frontera)) nf++;
      if (hayff) set_max(calv,Dfv[j].mod2());
      hm+=h=nj.h; set_min(hmin_e,nj.h);
      for (k=j+1;k<nv;k++){
        in1=ei[k];
        l=2.0*nj.distancia(n[in1])/(h+n[in1].h);
        set_min_max(lmin,lmax,l);
      }
    }
    hm/=nv;
    set_min(hmin_em,hmin_e);
    set_min(lminm,lmin); set_max(lmaxm,lmax);
    if (hayff) {
      calv=1.0/(hmin*sqrt(calv)); set_min(calvm,calv); // hmin en vez de lmin(?)
    }
    else calv=-1;
    vh=pow(fabs(ve[i]),1.0/dim)/hm; if (ve[i]<0) vh=-vh;
    set_min_max(vhminm,vhmaxm,vh);
    set_min(vminm,ve[i]);

    lmin_hist[hist_index(lmin)]++; lmax_hist[hist_index(lmax)]++;
    amin_hist[hist_index(amin)]++; amax_hist[hist_index(amax)]++;
    vh_hist[hist_index(vh)]++; if (hayff) calv_hist[hist_index(calv)]++; 

    // salida
    if (
      (hayff&&calv<.1) ||
      amin<5 ||
      amax<5 ||
      vh<.1 || vh>1 ||
      lmin<.5 || lmax>2
      )
      f << _oi(i)
        << setw(3) << nv
        << setw(3) << nf
        << setw(3) << vf
        << _of(ve[i])
        << _of(hm)
        << _of(lmin)
        << _of(lmax)
        << _of(vh)
        << _of(calv)
        << _of(amin)
        << _of(amax)
        << _EOL;
  }
  f << encabezado
    << _spacei // num
    << "         " //nv y nf
    << _of(vminm)
    << _of(hmin_em)
    << _of(lminm)
    << _of(lmaxm)
    << _of(vhminm)
    << _of(calvm)
    << _of(aminm)
    << _of(amaxm)
    << _EOL
    << _spacei // num
    << "         " //nv y nf
    << _spacef // ve
    << _spacef // hm
    << _spacef // lmin
    << _spacef // lmax
    << _of(vhmaxm);


  _OF2(".histogram");
  f << "                  lmin      lmax        vh      calv      amin      amax\n";
//      1234567890123456789012345678901234567890123456789012345678901234567890
  for (i=0;i<202;i++)
    f << _of(rangestart(i))
      << _oi(lmin_hist[i])
      << _oi(lmax_hist[i])
      << _oi(vh_hist[i])
      << _oi(calv_hist[i])
      << _oi(amin_hist[i])
      << _oi(amax_hist[i])
      << _EOL;

  _close;
  hayh=false; // porque es h medio
  _savetime(calidad);
  return true;
}


//=================================
// solo nodos
bool malla::graba_xyz(const char *arch){
  _makename("");

  char fn[_max_file_len]; strcpy(fn,filename); if (!fn[0]) return false;
  if (tipo.es(m_planaxy)){strcat(fn,".xy"); punto::io_z(false);}
  else {strcat(fn,".xyz"); punto::io_z(true);}
  backup(fn); ofstream f(fn);
  if (!f||!f.is_open()) {
    add_error(Not_Open);
    return false;
  }

  punto::o_precision(15);
  punto::o_width(25); punto::o_separator(" ");
  if (hayfn) 
    for (int i=0;i<n.len;i++) f << n[i] << ' ' << (n[i].f&fmask) << _EOL;
  else
    for (int i=0;i<n.len;i++) f << n[i] << _EOL;
  punto::io_reset();
  _close;
}

// archivos para matlab solo elementos
// no mezclado (triangulos segmentos tetras..)
bool malla::graba_con(const char *arch) {
  _makename("");
  _OF(".con");
  int old_base=elemento::io_base(); elemento::io_base(1); // .con es base 1
  elemento::flag_mask((hayfe)? fmask : 0); // para que no ponga flags comentar la linea
  for (int i=0;i<e.len;i++) f << e[i] << _EOL;  
  elemento::io_base(old_base);
  _close;
}


//=================================
// archivos .wrl
// solo superficies y lineas (no mezcla)
bool malla::graba_wrl(const char *arch) {
  if (!e) return false;
  if (tipo.es(m_vol)){
//    return graba_frontera("wrl",true);
    return false;
  }
  _makename("");
  _OF(".wrl");

  int i,j; punto pt; elemento ei;

  f << "#VRML V1.0 ascii" << '\n' << '\n';
  f << "#f File Generated by MESHSUITE: " << '\n' << '\n';
  f << "Separator {" << '\n';
  f << "  DEF BackgroundColor Info { string \"1.0 1.0 1.0\" }" << '\n'; // NO LE DA BOLA
  f << "  Separator {" << _EOL;

  // nodos
  punto::io_z(true); punto::o_precision(15);
  punto::o_width(25); punto::o_separator(" ");
  f << "    Coordinate3 {" << _EOL;
  f << "      point [" << _EOL;
  for (i=0;i<n.len;i++) f << n[i] << "," << _EOL;
  f << "      ] #point" << _EOL;
  f << "    } #Coordinate3" << _EOL;

  // caras
  f << "    Material { diffuseColor .3 .4 .78}" << _EOL;
  if (tipo.es(m_lin)) f << "    IndexedLineSet {" << _EOL;
  else f << "    IndexedFaceSet {" << _EOL;
  f << "      coordIndex [" << _EOL;
  for (i=0;i<e.len;i++){
    ei=e[i];
    for (j=0;j<ei.nv();j++) f  << ei[j] << ", "; // base 0
    f << "-1,\n";
  }
  f << "      ] #coordIndex" << _EOL;
  f << "    } #IndexedFaceSet" << _EOL;
  f << "  } #Separator" << _EOL;
  f << "} #Separator" << _EOL;

  _close;
  return true;
}



//=================================
// archivos .msh o .gid:
// (element sets por tipo)
bool malla::graba_msh(const char *arch) {
  _makename("");
  _OF(".msh");
  bool hayz=tipo.noes(m_planaxy),primertipo=true;
  int i,j;

  // element sets
  pline settipo[e_ntipos]; int qet;
  elemento::flag_mask(0); // para que no ponga flags
  int old_base=elemento::io_base(1); // base1
  for (i=0;i<e.len;i++) settipo[e[i].itipo()]+=i;

  // elementos
  for (j=1;j<e_ntipos;j++){ // los indefinidos (0) no salen
    qet=settipo[j].len;
    if (!qet) continue;
    elemento &e0=e[settipo[j][0]];
    /*http://www.gidhome.com/support_team/gidbeta/gid147.html
         #POSTPROCESS%20DATA%20FILES>Postprocess%20mesh%20
         format:%20ProjectName.post.msh */
    static const char* gid_elm_name[]=
      {"undefined",
       "Point",
       "Linear",
       "Triangle","Quadriateral","Polygon (not defined!!)",
       "Tetrahedra","Hexahedra","Prism","Polyhedron (not defined!!)"};

    f << "MESH dimension " << ((hayz)? "3" : "2");
    f << " ElemType " << gid_elm_name[e0.tipo()];
    f << " Nnode = " << e0.nv();
    f << _EOL;

    // nodos (solo el primer element set
    f << "coordinates" << _EOL;
    if (primertipo){
      primertipo=false;
      punto::io_z(hayz); punto::o_precision(15);
      punto::o_width(25); punto::o_separator(" ");
      for (i=0;i<n.len;i++){
        f << _oi(i+1) << n[i] << _EOL;
      }
      punto::io_reset();
    }
    f << "end coordinates" << _EOL;

    // conectividades
    f << "elements" << _EOL;
    pline &stj=settipo[j];
    for (i=0;i<qet;i++){
      f  << i+1 << '\t' << e[stj[i]] << _EOL; // si se escribe el flag hay que hacer log2
    }
    f << "end elements" << _EOL;
  }
  elemento::io_base(old_base);
  _close;
}


//=================================
// archivos .mai para SAMCEF
bool malla::graba_mai(const char *arch) {
  _makename("");
  _OF(".mai");

  int i,j,nv;
  e_tipo tipo;

  f << "! Number of nodes : " << n.len << _EOL;
  f << "! Number of elements : " << e.len << _EOL;

  // nodos
  punto::io_z(true); punto::o_precision(15);
  punto::o_width(25); punto::o_separator(" ");

  f << "! Nodes Definition\n.NOEUDS\n";
  for (i=0;i<n.len;i++)
    f << "I " << i+1 << " X " << n[i] << _EOL;

  punto::io_reset();

  // elementos

  f << "!\n! Elements Definition\n.MAILLES\n";

  for (i=0;i<e.len;i++){
    elemento &ei=e[i]; tipo=ei.tipo(); nv=ei.nv();
    if (tipo==e_poliedro) continue;
    f << "I " << i+1 << " N ";
    for (j=0;j<nv-1;j++){
      if ((tipo==e_cubo && j==4)||(tipo==e_wedge && j==3)) f << "0 ";
      f << ei[j]+1 << " ";
    }
    if (tipo==e_tetraedro) f << "0 ";
    f << ei[nv-1]+1 << _EOL;
  }

////////////// caras de frontera por grupos (condiciones de borde/flags)

  if (mk_vecino()) {
    f << "!\n! Condiciones de borde por grupos\n.CLT\n";

    // numeracion de caras diferentes entre elementos mios y samcef
    // tetras: cara i opuesta al nodo i
    //        0                                  3
    //        1                                  4
    //        2                                  2
    //        3                                  1
    // wedge:
    //        0: cuadrilatero opuesto al 0       3
    //        1: cuadrilatero opuesto al 1       4
    //        2: cuadrilatero opuesto al 2       2
    //        3: triangulo inferior              1
    //        4: triangulo superior              5
    // cubo:
    //        0: plano xy z=0                    1
    //        1: plano xz y=0                    2
    //        2: plano yz x=0                    5
    //        3: plano xy z=1                    6
    //        4: plano xz y=1                    4
    //        5: plano yz x=1                    3
    int *translate,
      t_tetra[]={3,4,2,1},
      t_wedge[]={3,4,2,1,5},
      t_cubo []={1,2,5,6,4,3};

    static const int ngrupos=8,flag_ini=(1<<ngrupos)-1;
    pline eg[ngrupos],ig[ngrupos];
    for (i=0;i<ngrupos;i++) {eg[i].resize(100);ig[i].resize(100);}

    int nc,k;
    flagtype flag;
    for (i=0;i<e.len;i++){
      elemento &ei=e[i]; nc=ei.nc();
      tipo=ei.tipo();
      if      (tipo==e_tetraedro) translate=t_tetra;
      else if (tipo==e_wedge)     translate=t_wedge;
      else if (tipo==e_cubo)      translate=t_cubo;
      else continue;
      for (j=0;j<nc;j++){
        if (vecino[i][j]>=0) continue;
        // cara de frontera
        elemento c=ei.cara(j); nv=c.nv();
        flag=flag_ini;for (k=1;k<nv;k++) flag.mask(n[c[k]].f);
//?????????        if (flag.es(16)) flag=16; // el 16 (arriba) esta duplicado con otros (ojo, caso particular)
        if (flag==0) continue;
        k=0; while(k<ngrupos-1&&!flag.es(1<<k)) k++;
        eg[k]+=i+1; ig[k]+=translate[j];
      }
    }

    for (i=0;i<ngrupos;i++)
      for (j=0;j<eg[i].len;j++)
        f << "FACE " << eg[i][j] << " " << ig[i][j] << " GRUPO " << i << _EOL;
  }
  else {
    _revienta(1);
  }

//////////////

  f << "RETURN\n!\n!:@#SD:END-OF-FILE\n";
  _close;
}

// archivos ls-dyna (.key)
// enteros en 8 y floating en 15
bool malla::graba_key(const char *arch) {
  _makename("");
  _OF(".key");

  f << "*KEYWORD\n";

  int i,j;

  // elementos
  if (e){
    // cuenta flags de elementos (pid's)
    // hasta free_flags y no valen oreados (usa el mas bajo)
    unsigned int fe;
    ordlist flist(_FREE_FLAGS); pline epid(e.len);
    const int localfmask =(1<<_FREE_FLAGS)-1;
    bool hayshell=false,hayvol=false; //puede haber ambas
    for (i=0;i<e.len;i++){
      if (e[i].pd) {epid+=-1; continue;} // ni poligonos ni poliedros
      if (e[i].dim()==3) hayvol=true; else hayshell=true;
      fe=e[i].f&localfmask;
      if (!fe) {flist+=0; epid+=0; continue;}
      j=0; while (!(fe&(1<<(j++)))); flist+=j; epid+=j;
    }
    for (j=0;j<flist.len;j++)
      f << "*PART\n"
        << "$# title\n"
        << "Elements_with_flag_" << flist[j] << _EOL
        << "$#     pid\n"
        << setw(8) << flist[j] // los pid's puede no ser consecutivos??
        << _EOL;
    int eid=0,nv;
    if (hayshell){
      flist.clean(); // la rearmo solo con pid's de sup para **VOL
      f << "*ELEMENT_SHELL\n"
        << "$#   eid     pid      n1      n2      n3      n4\n";
      for (i=0;i<e.len;i++){
        if (epid[i]<0) continue; //poly
        f << setw(8) << (++eid)
          << setw(8) << epid[i];
        flist+=epid[i];
        nv=e[i].nv();
        for (j=0;j<nv;j++) f << setw(8) << e[i][j]+1;
        for (   ;j< 4;j++) f << setw(8) << e[i][nv-1]+1; // repite el ultimo
        f << _EOL;
      }
      f << "*VOL\n"
        << setw(8) << 1 // unico volumen
        << setw(8) << 1 // part_id??
        << setw(8) << flist.len // cantidad de pid's en este vol
        << setw(8) << 0 // cantidad de pid's de mallas de h (no hay!!)
        << _EOL;
      for (j=0;j<flist.len;j++) f << setw(8) << flist[j] << _EOL;
    }
    if (hayvol){
      f << "*ELEMENT_SOLID\n"
        << "$#   eid     pid      n1      n2      n3      n4"
        << "      n5      n6      n7      n8\n";
      for (i=0;i<e.len;i++){
        if (epid[i]<0) continue; //poly
        f << setw(8) << (++eid)
          << setw(8) << epid[i];
        nv=e[i].nv();
        for (j=0;j<nv;j++) f << setw(8) << e[i][j]+1;
        for (   ;j< 8;j++) f << setw(8) << e[i][nv-1]+1; // repite el ultimo
        f << _EOL;
      }
    }
  }

  // nodos
  // podria sar n_permanente como *NODE_RIGID_SURFACE, pero no lo hago
  // facu pone h despues pero no esta bien con el estandar
  f << "*NODE\n"
    << "$#   nid               x               y               z               h\n";
  char num[17]; num[16]=0;
  if (hayh) for (i=0;i<n.len;i++)
    f << setw(8) << i+1
      << setprecision(8) << setw(16) << n[i][0] // e-001 ????
      << setprecision(8) << setw(16) << n[i][1]
      << setprecision(8) << setw(16) << n[i][2]
      << setprecision(8) << setw(16) << n[i].h
      << _EOL;
  else for (i=0;i<n.len;i++)
    f << setw(8) << i+1
      << setprecision(8) << setw(16) << n[i][0]
      << setprecision(8) << setw(16) << n[i][1]
      << setprecision(8) << setw(16) << n[i][2]
      << _EOL;

  // fin
  f << "*END\n";

  _close;
}


//================================================
// archivos de autocad

static void acadheader(ofstream &f, const char* l1=0, const char* l2=0){
  f << "  0\nSECTION\n  2\nHEADER\n  9\n$ACADVER\n  1\nAC1009\n  0\nENDSEC\n"
    << "  0\nSECTION\n  2\nTABLES\n  0\nTABLE\n  2\nLAYER\n 70\n     3\n";
  if (l1)  {
    f << "  0\nLAYER\n  2\n" << l1
      << "\n 70\n     0\n 62\n     7\n  6\nCONTINUOUS\n";
  if (l2)
    f << "  0\nLAYER\n  2\n" << l2
      << "\n 70\n     0\n 62\n     7\n  6\nCONTINUOUS\n";
  }
  f << "  0\nENDTAB\n  0\nENDSEC\n"
    << "  0\nSECTION\n  2\nENTITIES\n";
}

static void dxfp(ofstream &f, const nodo& p, int n=0){
  f << " " << 10+n << _EOL << setprecision(15) << p[0] << _EOL;
  f << " " << 20+n << _EOL << setprecision(15) << p[1] << _EOL;
  f << " " << 30+n << _EOL << setprecision(15) << p[2] << _EOL;
}

// numero de punto
bool malla::graba_dxf_nnod(const char *arch) {
  double ht=pmin.distancia(pmax)/200;
  _makename("_nnod.dxf");
  _OF("");
  acadheader(f,"nni","nnf");
  double hi=ht;
  for (int i=0;i<n.len;i++){
    f << "  0\nTEXT\n  8\n";
    if (n[i].f.es(n_frontera)) f << "nnf\n"; else f << "nni\n";
    dxfp(f,n[i]);
    if (n[i].h>ERRADM)  hi=n[i].h/5;
    f << " 40\n" << hi << "\n  1\n" << i << _EOL;
   }
  f << "  0\nENDSEC\n  0\nEOF\n";
  _close;
}

// graba una malla como dxf
// falta implementar poliedros
bool malla::graba_dxf(const char *arch, double esize) {

  static const char etetra[]={
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n2\n 72\n3\n 73\n4\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n3\n 72\n1\n 73\n4\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n4\n 72\n1\n 73\n2\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n1\n 72\n3\n 73\n2\n"
    "  0\nSEQEND\n"};
  static const char ewedge[]={
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n2\n 72\n3\n 73\n6\n 74\n5\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n3\n 72\n1\n 73\n4\n 74\n6\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n1\n 72\n2\n 73\n5\n 74\n4\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n1\n 72\n3\n 73\n2\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n4\n 72\n5\n 73\n6\n"
    "  0\nSEQEND\n"};
  static const char ecubo[]={
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n1\n 72\n4\n 73\n3\n 74\n2\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n1\n 72\n2\n 73\n6\n 74\n5\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n1\n 72\n5\n 73\n8\n 74\n4\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n7\n 72\n8\n 73\n5\n 74\n6\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n7\n 72\n3\n 73\n4\n 74\n8\n"
    "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n 71\n7\n 72\n6\n 73\n2\n 74\n3\n"
    "  0\nSEQEND\n"};

  _makename("");
  _OF(".dxf");
  acadheader(f,"nodos","elementos");
  int i,j,nlen=n.len,elen=e.len,nv,nc,dim,k,ix,in;
  punto g,pt[32];

  // nodos sueltos (no tienen elementos)
  for (i=0;i<nlen;i++){
    if (esize==1&&n[i].e) continue;
    f << "  0\nPOINT\n  8\nnodos\n";
    dxfp(f,n[i]);
  }

  if (esize<=ERRADM) esize=1;

  if (tipo.es(m_sup)&&tipo.noes(m_vol|m_lin)&&
       esize==1&&nlen<32767){
    // polyface mesh
    //(32767=maximo permitido por autocad)
    f << "  0\nPOLYLINE\n  8\nelementos\n 66\n1\n 70\n64\n";
    f << " 71\n" << nlen << "\n 72\n" << elen << _EOL;
    // nodos
    for (i=0; i<nlen; i++){
      const punto& ni=n[i];
      f << "  0\nVERTEX\n  8\nelementos\n";
      dxfp(f,ni);
      f << " 70\n192\n";
    }
    // elementos
    for (i=0; i<elen; i++){
      const elemento &ei=e[i]; nv=ei.nv();
      for (ix=1,k=0;k<((nv-1)>>1);k++){
        f << "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n";
        in=ei[0]+1; if (ix!=1) in=-in; // visible si es el primero
        f << " 71\n" << in << _EOL;
        in=ei[ix]+1; ix++; // visible
        f << " 72\n" << in << _EOL;
        in=ei[ix]+1; ix++; // visible
        f << " 73\n" << in << _EOL;
        if (ix<nv){
          in=ei[ix]+1; if (ix!=nv-1) in=-in; // visible si no hay mas
          f << " 74\n" << in << _EOL;
        }
      }
    }
    f << "  0\nSEQEND\n";
  }
  else for (i=0; i<elen; i++) {
    // elementos sueltos
    elemento &ei=e[i]; nv=ei.nv(); nc=ei.nc(); dim=ei.dim();
    if (esize!=1) g=gv(i);
    // nodos
    for (j=0;j<nv;j++) {
      pt[j]=n[ei[j]]; if (esize!=1) pt[j]=pt[j]*esize+g*(1-esize);
    }
    if (dim==3){// una polyface
      f << "  0\nPOLYLINE\n  8\nelementos\n 66\n1\n 70\n64\n 71\n" << nv
        << "\n 72\n" << nc  << _EOL;
      // nodos
      for (j=0;j<nv;j++) {
        f << "  0\nVERTEX\n  8\nelementos\n";
        dxfp(f,pt[j]);
        f << " 70\n192\n";
      }
      // caras
      e_tipo etipo=ei.tipo();
      if (etipo==e_tetraedro) f << etetra;
      else if (etipo==e_wedge) f << ewedge;
      else if (etipo==e_cubo) f << ecubo;
      else if (etipo==e_poliedro){
        int nc=ei.nc(),nvc,ic[3];
        for (j=0;j<nc;j++){
          ei.cara(j,nvc,ic);
          f << "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128"
            << "\n 71\n" << ic[0]+1 << "\n 72\n" << ic[1]+1 << "\n 73\n" << ic[2]+1 << _EOL;
        }
        f << "  0\nSEQEND\n";
      }
    }
    else if (dim==2){
      if (nv==3||nv==4){ // 3dface
        f << "  0\n3DFACE\n  8\nelementos\n";
        for (j=0;j<nv;j++) dxfp(f,pt[j],j);
        if (nv==3) dxfp(f,pt[2],3); //triangulo => repite el tercero
      }
      else { // poligono --> polyface
        f << "  0\nPOLYLINE\n  8\nelementos\n 66\n1\n 70\n64\n";
        f << " 71\n" << nv << "\n 72\n1\n";
        // nodos
        for (j=0;j<nv;j++){
          f << "  0\nVERTEX\n  8\nelementos\n";
          dxfp(f,pt[j]);
          f << " 70\n192\n";
        }
        // elemento
        for (ix=1,k=0;k<=nv/4;k++){
          f << "  0\nVERTEX\n  8\nelementos\n 10\n0\n 20\n0\n 70\n128\n";
          in=1; if (ix!=1) in=-in; // visible si es el primero
          f << " 71\n" << in << _EOL;
          in=ix+1; ix++; // visible
          f << " 72\n" << in << _EOL;
          in=ix+1; ix++; // visible
          f << " 73\n" << in << _EOL;
          if (ix<nv){
            in=ix+1; if (ix!=nv-1) in=-in; // visible si no hay mas
            f << " 74\n" << in << _EOL;
          }
        }
        f << "  0\nSEQEND\n";
      }
    }
    else if (dim==1){
      f << "  0\nLINE\n  8\nelementos\n";
      dxfp(f,pt[0],0);dxfp(f,pt[1],1);
    }
  }
  f << "  0\nENDSEC\n  0\nEOF\n";
  _close;
}

//=================================================================================
// cascara

// solo conectividades
bool cascara::graba_con(const char *arch){
  const char *nombre=parent->nombre;
  _makename("");
//  _OF(".con"); add_error no va en cascaras
  char fn[_max_file_len]; strcpy(fn,filename);
  if (!fn[0]) return false;
  strcat(fn,"_cascara.con"); backup(fn); ofstream f(fn);
  if (!f||!f.is_open()) return false;

  int old_base=elemento::io_base(); elemento::io_base(1); // .con es base 1
  for (int i=0;i<e.len;i++) f  << e[i] << _EOL;
  elemento::io_base(old_base); // restaura
  _close;
}

// nodos y conectividades
bool cascara::graba_dat(const char *arch) {
  const char *nombre=parent->nombre;
  _makename("");
  char fn[_max_file_len]; strcpy(fn,filename);
  if (!fn[0]) return false;
  strcat(fn,".dat"); backup(fn); ofstream f(fn);
  if (!f||!f.is_open()) return false;

  int i,j,ix,base=elemento::io_base();
  f << n.len << " Nodes # x y";
  bool
    hayz=parent->tipo.noes(m_planaxy),
    hayh=parent->hayh,
    hayfn=parent->hayfn,
    hayv=parent->hayv;

  if (hayz)  f << " z";
  if (hayh)  f << " h";
  if (hayfn) f << " f";
  if (hayv)  f << " v";
  if (base==1)  f << " base=1";
  f << _EOL;
  punto::io_z(hayz); punto::o_precision(15);  // toda
  punto::o_width(25); punto::o_separator(" ");

  const array1<nodo> &pn=parent->n;
  for (i=0;i<n.len;i++){
    f << _oi(n[i]+base)
      << pn[n[i]];
    if (hayh) f << _or(pn[n[i]].h);
    if (hayfn) f << setw(15) << (pn[n[i]].f&fmask);
    if (hayv) f << _or(pn[n[i]].v);
    f << _EOL;
  }
  punto::io_reset();

  // element sets (debe haber correspondencia con los textos de elemento)
  pline settipo[e_ntipos]; int qet;
  elemento::flag_mask((parent->hayfe)? fmask : 0); // para que no ponga flags comentar la linea
  for (i=0;i<e.len;i++) settipo[e[i].itipo()]+=i;
  for (j=1;j<e_ntipos;j++){ // los indefinidos (0) no salen
    qet=settipo[j].len;
    if (!qet) continue;
    f << qet << " " << e[settipo[j][0]].stipo() << " #";
    if (parent->hayfe) f << " f";
    f << _EOL;
    for (i=0;i<qet;i++){
      ix=settipo[j][i];
      f << ix << "\t" << e[ix] << _EOL;
    }
  }
  _close;
}

#undef _oi
#undef _spacei
#undef _or
#undef _spacer
#undef _of
#undef _spacef
#undef _OF
#undef _OF2
#undef _close
#undef _makename
#undef _makename2
