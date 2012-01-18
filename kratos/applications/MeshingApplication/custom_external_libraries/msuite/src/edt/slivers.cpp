// Rutinas de voronoi y malla para slivers

//#define temporizar
#include "tiempo.h"

#include "voronoi.h"

//#include <iostream> // salida a cout (odbg) ////////////////sacar!!

using namespace std;

/* OJO: 
  Que pasa con un swap en un sandwich de slivers?
    Muy probablemente se trate de tres slivers en un pentagono circular.
    El resultado (aparentemente) es la eliminacion de uno de ellos sin mas 
      problemas, excepto, quizas, por las vecindades con el resto.
    Parece que pueden quedar no-slivers de vol minimo tb.
  No confundir slivers con elementos anisotropicos, en estos ultimos hay nodos cerca.
  Al testear slivers ordena el sliver 0123 cuadrilatero, siendo 02 el maximo diedro
  Al swappear:
    Se puede generar otro sliver que hay que retestear (sin entrar en loop infinito)
    Tambien puede que sea de frontera, en cuyo caso hay que retestear f_slivers, pero
      cuidado que hay una lista para rehacer, no debe ser peligroso tener uno borrado
*/

static double _volchico=1e-1; // 1/10 del volumen ideal
static pline rehacer_f_slivers; // con flag1
static ordlist rehacer_i_slivers; // con flag2
static const int syf=e_frontera|e_sliver;
static const punto *p00; // tres vectores arista, origen en el 0
static double hc; // h de comparacion
static int imaxd0,imaxd1; // indices ini y fin de la arista con max diedro (imaxd0<imaxd1)


// calcula diedros
#define _testd(ic0,ic1,in0,in1) {cosa=f[ic0]*f[ic1]; \
/*if (odbg) cout << "cosa: " << cosa << endl;*/\
  if (cosa>.98) {\
  ceros++; if (set_max(maxcosa,cosa)) {imaxd0=in0; imaxd1=in1;}\
  }\
  else if (cosa>-.98) return false;\
}

//static bool odbg=false;
static bool _es_sliver() {
  double r,l,p2n,lim=hc/10;// +/- h/10 de la circunferencia
  punto c,av,p2h;
  // para ser sliver n3 debe estar cerca de la circunferencia n0 n1 n2
  // se compara la normal y la proyeccion
  c30(p00[0],p00[1],c,r,av); // circunferencia de 3 puntos
  l=av.mod(); // av=area vector; l=area
  if (l<ERRADM) return false; // no se que corcho es pero no es un sliver
  av/=l; // normal
  // Componente de p2 normal al plano del circulo
  p2n=p00[2]*av;
//if (odbg) cout << endl << "lim: " << lim << endl;
//if (odbg) cout << "p2n: " << p2n << endl;
  if (p2n>lim) return false; // vol+ y p2 separado del plano p0p1 (no va fabs!)
  if (p2n<=0) return true; // lo considero sliver y listo OJO!!!!!
  // Proyeccion de p2 sobre el plano del circulo y centrado
  p2h=(p00[2]-av*p2n)-c; l=p2h.mod(); // l=distancia al centro
//if (odbg) cout << "l-r: " << l-r << endl;
  if (fabs(l-r)>lim) return false; // lejos de la circunferencia
  // El punto esta cerca del circulo, tanto en el plano como normal.
  // Todavia puede ser un tetraedro rectangular chico (de prisma aplastado por ej)
  //   o tb un cap con un nodo cerca de otro.
  // Hay que testear diedros, debe dar 4 de ~0 y dos de ~180 correctamente distribuidos
  // normales a las caras hacia el elemento
  punto f[4]={
    ((p00[2]-p00[0])%(p00[1]-p00[0])).dir(), // cara 0 (opuesta al nodo 0)
     (p00[1]%p00[2]).dir(),                  // cara 1
     (p00[2]%p00[0]).dir(),                  // cara 2
     (p00[0]%p00[1]).dir()                   // cara 3
  };
  // cos >.98 o <-.98 son unos 11 grados de dif con 0 o 180
  double cosa,maxcosa=0; int ceros=0;
  // indice de cara 0,1 indice de vertice de arista 0 1
  _testd(0,1,2,3);
  _testd(0,2,1,3);
  _testd(0,3,1,2);
  _testd(1,2,0,3);
  _testd(1,3,0,2);
  _testd(2,3,0,1);
  if (ceros!=2) return false;
  return true;
}

bool malla::es_sliver(int ie, double h) const {
  elemento &ei=e[ie];
  if (ei.f.es(e_borrado)) return false;
  // primero descarta los de volumen grande
  if (ve&&ve[ie]>h*h*h*_volchico/6) return false; // vol grande
  punto p[3]={n[ei[1]]-n[ei[0]],n[ei[2]]-n[ei[0]],n[ei[3]]-n[ei[0]]}; p00=p;
  if (!ve&&triple(p[0],p[1],p[2])>h*h*h*_volchico) return false; // vol grande
  hc=h;
  return _es_sliver(); // calcula diedros
  // este no marca ni reordena
}

bool voronoi::mark_sliver(int is){
  if (s[is].f.es(e_borrado)) return false;
  esfera &si=s[is]; const array1<nodo> &n=m->n;
//if (is==4506832) odbg=true; else odbg=false;
  hc=Min(n[si[0]].h,Min(n[si[1]].h,Min(n[si[2]].h,n[si[3]].h)));
//if (odbg) cout << "vt: " << si.vt << endl;
//if (odbg) cout << "hc3*_volchico/6: " << hc*hc*hc*_volchico/6 << endl;
  if (si.vt>hc*hc*hc*_volchico/6){si.f.reset(e_sliver); return false;} // vol grande
  punto p[3]={n[si[1]]-n[si[0]],n[si[2]]-n[si[0]],n[si[3]]-n[si[0]]}; p00=p;
  if (!_es_sliver()) {si.f.reset(e_sliver); return false;}
//if (odbg) cout << "ES SLIVER" << endl;
  // es sliver => ordena caras
  si.f.set(e_sliver);
  int nswaps=0;
  if (imaxd0!=0){si.swap(0,imaxd0); nswaps++;}// imaxd0 pasa a ser el nodo 0
  if (imaxd1!=2){si.swap(2,imaxd1); nswaps++;}// imaxd1 pasa a ser el nodo 2
  if (nswaps==1){si.swap(1,3);}// hay que hacer otro para que quede derecho
  return true;
}
bool voronoi::mark_slivers(bool remake){
  _initime;
  if (NV!=4||(!remake&&slivers_marked)) return true;
  if (!mk_h_nn(remake)) return false; // necesita h para no calcular diedros al cuete
  int count=0;
  for (int is=0;is<s.len;is++) if (mark_sliver(is)) count++;
  slivers_marked=true;
  if (m->INFO_CL) cout << count << " slivers" << endl;
  _savetime(markslivers);
  return true;
}

bool malla::mark_sliver(int ie){
  elemento &ei=e[ie];
  if (ei.nv()!=4||ei.tipo()!=e_tetraedro||ei.f.es(e_borrado)) return false;
  hc=Min(n[ei[0]].h,Min(n[ei[1]].h,Min(n[ei[2]].h,n[ei[3]].h)));
  if (ve&&ve[ie]>hc*hc*hc*_volchico/6){ei.f.reset(e_sliver); return false;} // vol grande
  punto p[3]={n[ei[1]]-n[ei[0]],n[ei[2]]-n[ei[0]],n[ei[3]]-n[ei[0]]}; p00=p;
  if (!ve&&triple(p[0],p[1],p[2])>hc*hc*hc*_volchico){
      ei.f.reset(e_sliver); return false; // vol grande
  }
  if (!_es_sliver()) {ei.f.reset(e_sliver); return false;}
  // es sliver => ordena caras
  ei.f.set(e_sliver);
  int swaps=0,*ne=ei.n;
  if (imaxd0!=0){// imaxd0 pasa a ser el nodo 0
    Swap(ne[0],ne[imaxd0]); if (vecino) Swap(vecino[ie][0],vecino[ie][imaxd0]);
    swaps++;
  }
  if (imaxd1!=2){// imaxd1 pasa a ser el nodo 2
    Swap(ne[2],ne[imaxd1]); if (vecino) Swap(vecino[ie][2],vecino[ie][imaxd1]);
    swaps++;
  }
  if (swaps==1){// hay que hacer otro
    Swap(ne[1],ne[3]); if (vecino) Swap(vecino[ie][1],vecino[ie][3]);
  }
  return true;
}

bool malla::mark_slivers(bool remake){
  _initime;
  if (tipo.noes(m_vol)||(!remake&&tipo.es(m_slivers_marked))) return true;
  if (!hayh&&!mk_h_nn()) return false;
  for (int ie=0;ie<e.len;ie++) mark_sliver(ie);
  tipo.set(m_slivers_marked);
  _savetime(markslivers);
  return true;
}

// Elimina slivers interiores
// si el sliver es de frontera, no hace nada.
void voronoi::rm_i_slivers(){
  if (NV!=4) return;
  int is,count=0,laburo=0,oldilen,oldflen;
  static const int bof=e_borrado|e_frontera;
  // asi por slivers perpendiculares
  squeeze();
  if (!mark_slivers()) return;
  _initime;
  ordlist rehacer(s.len/10); 
  for (is=0;is<s.len;is++) {
    flagtype &f=s[is].f;
    if (f.noes(e_sliver)&&f.es_alguno(bof)) continue;
    s[is].f.set(flag2); rehacer+=is;
  }
  rehacer_i_slivers.clean();
  while (rehacer){
    oldilen=rehacer.len; oldflen=rehacer_f_slivers.len;
    while (rehacer){
      is=rehacer.last(); rehacer.len--;
      if (is>=s.len) continue; // si viene por rehacer
      flagtype &f=s[is].f; f.reset(flag2);
      if (f.es(e_sliver)&&f.noes_ninguno(bof)) rm_i_sliver(is);
      if (m->INFO_CL && !((++laburo)%1000)) {
        cout << "\r" << laburo << "done, " << rehacer.len+rehacer_i_slivers.len << " tests remaining                  " << flush;
      }
    }
    rehacer.roba(rehacer_i_slivers); // para reducir perpendiculares
    if (rehacer.len==oldilen&&rehacer_f_slivers.len==oldflen) count++; else count=0;
    if (rehacer_f_slivers) rm_f_slivers(); // para sacar sandwiches
    if (count==4) 
      break;
  }
  if (m->INFO_CL) 
    cout << "\r" << rehacer.len << " slivers remains                                          " << endl;
  _savetime(i_slivers);
}
// si es interior y hay dos tetraedos a un lado, swapea diagonales
bool voronoi::rm_i_sliver(int is){  
  esfera &si=s[is];  
  int 
      v0=si.vecino[0],
      v1=si.vecino[1],
      v2=si.vecino[2],
      v3=si.vecino[3];
  
  //frontera
  if (v0<0||v1<0||v2<0||v3<0) return false;

  // interior (sin vecinos <0)
  // uso la diagonal que tenga dos vecinos, vecinos entre si
  // el mayor diedro es 02 y sigue 13
  if (s[v1].tiene_vecino(v3) && swap3x2(is,v1,v3)) return true; // mayor diedro (0-2)
  if (s[v0].tiene_vecino(v2) && swap3x2(is,v0,v2)) return true; // diagonal opuesta (1-3)
  if (swap4x4(is)) return true;
  si.swap(0,1); si.swap(2,3);// se labura por la arista 0-2
  if (swap4x4(is)) return true;
  rehacer_i_slivers+=is; return false;
}

#define _retest_s(ixs,fl){\
  if(mark_sliver(ixs)){\
    if (fl.es(e_frontera)) {\
      if (fl.noes(flag1)) {\
        fl.set(flag1); rehacer_f_slivers+=ixs;\
      }\
    }\
    else {\
      if (fl.noes(flag2)) {\
        fl.set(flag2); rehacer_i_slivers+=ixs;\
      }\
    }\
  }\
}

// tres tetras compartiendo una arista -> dos compartiendo una cara
bool voronoi::swap3x2(int is0, int is1, int is2){
  array1<nodo>&n=m->n;
  
  // la arista comun de punta tiene dos nodos ni abajo y ns arriba
  // la futura cara tiene tres nodos n0 n1 y n2; ni es opuesto al tetra si

  esfera s0=s[is0],s1=s[is1],s2=s[is2];

  int // cara
    in0=s1[s1.index_vecino(is0)],
    in1=s2[s2.index_vecino(is1)],
    in2=s0[s0.index_vecino(is2)];

  //pone nodos y vecinos ordenados
  // cualquier permutacion de un par de nodos es impar al reordenar (cambia el sentido)
  int ix,perm=0;
  ix=s0.index(in1); if (ix!=0) {perm++; s0.swap(0,ix);}
  ix=s0.index(in2); if (ix!=1) {perm++; s0.swap(1,ix);}
  if (perm&1) s0.swap(2,3);
  
  perm=0;
  ix=s1.index(in2); if (ix!=0) {perm++; s1.swap(0,ix);}
  ix=s1.index(in0); if (ix!=1) {perm++; s1.swap(1,ix);}
  if (perm&1) s1.swap(2,3);
  
  perm=0;
  ix=s2.index(in0); if (ix!=0) {perm++; s2.swap(0,ix);}
  ix=s2.index(in1); if (ix!=1) {perm++; s2.swap(1,ix);}
  if (perm&1) s2.swap(2,3);

  int 
    ini=s0[2],ins=s0[3], // arista
    vs0=s0.vecino[2],vi0=s0.vecino[3],
    vs1=s1.vecino[2],vi1=s1.vecino[3],
    vs2=s2.vecino[2],vi2=s2.vecino[3];

  // dos tretras nuevos, uno con ns y otro con ni
  // s0 vuela, ss reemplaza a s1 y si reemplaza a s2
  esfera ss,si;
  ss[0]=in0; ss[1]=in1; ss[2]=in2; ss[3]=ins; ss.define();
  si[0]=in0; si[1]=in2; si[2]=in1; si[3]=ini; si.define(); // orientacion
  if (ss.vt<ERRADM||si.vt<ERRADM){  // para no hacer macanas
//    s[is0].f.set(flag3); s[is1].f.set(flag3); s[is2].f.set(flag3); 
    return false; // probable sandwich
  }  

  // vecinos
  ss.vecino[0]=vs0; ss.vecino[1]=vs1; ss.vecino[2]=vs2; ss.vecino[3]=is2;
  si.vecino[0]=vi0; si.vecino[1]=vi2; si.vecino[2]=vi1; si.vecino[3]=is1; // orientacion
  if (vs0>=0) s[vs0].replace_vecino(is0,is1); if (vs2>=0) s[vs2].replace_vecino(is2,is1);
  if (vi0>=0) s[vi0].replace_vecino(is0,is2); if (vi1>=0) s[vi1].replace_vecino(is1,is2);

  // elementos de los nodos
  // n0 tenia a s1 y a s2 asi que esta bien  
  n[in1].e.replace1(is0,is1); // n1 pierde a s0 y gana a s1  
  n[in2].e.replace1(is0,is2); // n2 pierde a s0 y gana a s2
  n[ins].e.remove1(is0);n[ins].e.remove1(is2); // ns pierde a s0 y a s2
  n[ini].e.remove1(is0);n[ini].e.remove1(is1); // ni pierde a s0 y a s1

  // por si hay material o algo, pongo los flags oreados (sin sliver ni frontera)
  flagtype &fs=s[is1].f,&fi=s[is2].f;
  ss.f=si.f=((fs|fi)&(~(syf|flag1|flag2)));
  if (vs0<0||vs1<0||vs2<0) ss.f.set(e_frontera);
  if (vi0<0||vi1<0||vi2<0) si.f.set(e_frontera);
  s[is1]=ss; s[is2]=si; // reemplaza  
  _retest_s(is1,fs);  _retest_s(is2,fi); // ojo que mark_sliver puede cambiar la orientacion

  //s[is0].f.set(e_borrado); s.remove(is0);
  // cambia is0 por el ultimo no-sliver (si no es este)
  int isl,j,v;  ordlist &sb=s.borrado;
  if (is0!=s.len-1) { // no es la ultima de s
    isl=s.len-1; while (s[isl].f.es_alguno(e_sliver|e_borrado)&&isl>=0) {isl--;}
    _revienta(isl<0);
    esfera &sl=s[isl]; // la ultima no borrada y no sliver
    const int *nsl=sl.n;
    for (j=0;j<4;j++) {
      n[nsl[j]].e.replace1(isl,is0);
      v=sl.vecino[j]; if (v>=0) s[v].replace_vecino(isl,is0);
    }
    s[is0]=sl;
    sl.f.set(e_borrado); s.remove(isl); // por si isl no es la ultima
  }
  else s.len--; // s0 era la ultima
  while (sb.len&&sb.last()==s.len-1) {sb.len--; s.len--;} // por si es la ultima
  
  return true;
}

// cambia la diagonal en un octaedro
// se supone que ya se intento swap3x2
bool voronoi::swap4x4(int is0){
  static bool primerpaso=true;
  esfera s0=s[is0]; 
  const array1<nodo>&n=m->n;
  cpliner se=n[s0[0]].e.inters(n[s0[2]].e); // elementos de la arista 02
  if (se.len!=4) // array cubico o sliver perpendicualr
    return false;
  
  int is1,is2,is3,in0,in1,in2,in3,in4,in5,iv0,iv1,iv2,iv3,iv4,iv5,iv6,iv7,k;

  // arista comun 02, caras de ss 1 y 3
  // el resto de los elementos, nodos y vecinos
  is2=s0.vecino[3]; is3=s0.vecino[1];
  for (k=0; k<4; k++) {is1=se[k]; if (is1!=is0&&is1!=is2&&is1!=is3) break;}
  esfera s1=s[is1], s2=s[is2], s3=s[is3];
  in0=s0[0]; in1=s0[1]; in2=s0[2]; in3=s0[3];
  iv0=s0.vecino[2]; iv1=s0.vecino[0];
  k=s1.index(in2); iv2=s1.vecino[k]; k=s1.index(in0); iv3=s1.vecino[k];
  k=s2.index_vecino(is0); in4=s2[k]; 
  k=s2.index(in2); iv4=s2.vecino[k]; k=s2.index(in0); iv5=s2.vecino[k];
  k=s3.index_vecino(is0); in5=s3[k];
  k=s3.index(in2); iv6=s3.vecino[k]; k=s3.index(in0); iv7=s3.vecino[k];

  // OJO si hay dos vecinos iguales que no pegue dos tetras por dos caras
  // 0 y 1 no pueden ser iguales porque habria hecho 3x2
  // como 13 es diagonal solo pueden ser iguales 426 o 537
  // pero solo es peligroso 24 o 35 que darian dos contra dos, 
  // 26 o 37 dan slver contra dos que se puede arreglar despues (o no, o loop inifinito)
  if (primerpaso && ((iv2>=0&&iv2==iv4)||(iv3>=0&&iv3==iv5))){
    // se fija si se puede hacer del otro lado
    if ((iv3>=0&&iv3==iv7)||(iv2>=0&&iv2==iv6)) return false;
    s[is0].swap(0,2); s[is0].swap(1,3); // la otra diagonal
    primerpaso=false;
    bool retval=swap4x4(is0);
    primerpaso=true;
    return retval;
  }
  
  // rearma 
  // por si hay material o algo, pongo los flags oreados (sin sliver ni frontera)
  flagtype f=((s0.f|s1.f|s2.f|s3.f)&(~(syf|flag1|flag2)));

  // la arista comun (in1-in5) es siempre del nodo 0 al nodo 3
  s0.n[0]=in1; s0.n[1]=in0; s0.n[2]=in3; s0.n[3]=in5;
  s0.vecino[0]=iv6; s0.vecino[1]=is1; s0.vecino[2]=is2; s0.vecino[3]=iv0;
  s0.define(); s0.f=f; if (iv6<0||iv1<0) s0.f.set(e_frontera);
  
  s1.n[0]=in1; s1.n[1]=in3; s1.n[2]=in2; s1.n[3]=in5;
  s1.vecino[0]=iv7; s1.vecino[1]=is3; s1.vecino[2]=is0; s1.vecino[3]=iv1;
  s1.define(); s1.f=f; if (iv7<0||iv0<0) s1.f.set(e_frontera);
  
  s2.n[0]=in1; s2.n[1]=in4; s2.n[2]=in0; s2.n[3]=in5;
  s2.vecino[0]=iv2; s2.vecino[1]=is0; s2.vecino[2]=is3; s2.vecino[3]=iv4;
  s2.define(); s2.f=f; if (iv2<0||iv4<0) s2.f.set(e_frontera);
  
  s3.n[0]=in1; s3.n[1]=in2; s3.n[2]=in4; s3.n[3]=in5;
  s3.vecino[0]=iv3; s3.vecino[1]=is2; s3.vecino[2]=is1; s3.vecino[3]=iv5;
  s3.define(); s3.f=f; if (iv3<0||iv5<0) s3.f.set(e_frontera);

  if (s0.vt<-ERRADM || s1.vt<-ERRADM || s2.vt<-ERRADM || s3.vt<-ERRADM) {
    if (!primerpaso) return false;
    // se fija si se puede hacer del otro lado
    s0.swap(0,2); s0.swap(1,3); // la otra diagonal
    primerpaso=false;
    bool retval=swap4x4(is0);
    primerpaso=true;
    return retval;
  }

  s[is0]=s0; s[is1]=s1; s[is2]=s2; s[is3]=s3; 

  nodo &n0=n[in0], &n1=n[in1], &n2=n[in2], &n3=n[in3], &n4=n[in4], &n5=n[in5];
  n0.e.remove1(is1); n0.e.remove1(is3); 
  n2.e.remove1(is0); n2.e.remove1(is2);
  n1.e+=is1; n1.e+=is3; 
  n5.e+=is0; n5.e+=is2;
  n3.e.replace1(is3,is1); 
  n4.e.replace1(is1,is3);

  // ojo vecinos iguales (26 o 37) con 3 y 7 hay un problema (is1xis3 e is3xis1)
  if (iv2>=0) s[iv2].replace_vecino(is1,is2); if (iv6>=0) s[iv6].replace_vecino(is3,is0); 
//  if (iv3>=0) s[iv3].replace_vecino(is1,is3); if (iv7>=0) s[iv7].replace_vecino(is3,is1);
  if (iv3>=0) s[iv3].replace_vecino(is1,-s.len); if (iv7>=0) s[iv7].replace_vecino(is3,is1); if (iv3>=0) s[iv3].replace_vecino(-s.len,is3);
  if (iv5>=0) s[iv5].replace_vecino(is2,is3);

  // ojo que mark_sliver puede cambiar la orientacion  
  _retest_s(is0,s0.f);
  _retest_s(is1,s1.f);
  _retest_s(is2,s2.f);
  _retest_s(is3,s3.f);
  
  return true;
}

#undef _retest_s

// slivers de frontera
// para volar caps hay que marcar de otro modo.
void voronoi::rm_f_slivers(){
  if (NV!=4) return;
  _initime;
  array1<nodo> &n=m->n;
  int is,isl,v,k;
  ordlist &sb=s.borrado;
  pline &rehacer=rehacer_f_slivers;
  if (!rehacer){
    squeeze();
    if (!mark_slivers()) return;
    rehacer.ini(); rehacer.natural(s.len); 
    for (is=0;is<s.len;is++) s[is].f.set(flag1);
  }
  while (rehacer){
    is=rehacer.last(); rehacer.len--;
    if (is>=s.len) continue; // si viene de i_slivers por rehacer
    esfera &si=s[is];
    si.f.reset(flag1);
    if (si.f.noes_alguno(syf)) continue;
    if (si.f.es(e_borrado)) continue; // si viene de i_slivers por rehacer
    // verifca que tenga al menos dos caras en la frontera
    // (e_frontera es al menos una en la frontera pero puede haber sandwich)
    // ahora esta ordenado => 02 o 13
    if (
      (si.vecino[0]>=0||si.vecino[2]>=0) // de este lado no son los dos frontera
      &&                                 // y
      (si.vecino[1]>=0||si.vecino[3]>=0) // de este lado tampoco
      ) continue; // no es de frontera
/*    
    // 1:cap o sandwich 2:sliver 3 y 4 rarezas de la naturaleza
    for (k=count=0;k<4&&count<2;k++) {if (si.vecino[k]<0) count++;}
ya no vuela caps!!
    if (!count) continue;
    if (count==1){// puede ser un cap=>vuela
#pragma message("<---------------- ***   verificar si entra aca porque no marca caps!!!  ***\n")
      //cap: area mayor=suma de las otras tres
      punto n0=n[si.n[0]], n1=n[si.n[1]]-n0, n2=n[si.n[2]]-n0, n3=n[si.n[3]]-n0;
      double a0=((n2-n1)%(n3-n1)).mod(), a1=(n2%n3).mod(), a2=(n3%n1).mod(), a3=(n1%n2).mod();
      if (a1>a0) Swap(a0,a1); if (a2>a0) Swap(a0,a2); if (a3>a0) Swap(a0,a3); // la mayor->a0
      if ((a1+a2+a3-a0)>1e-2*a0) continue; // no es cap
    }
    if (count<2) continue; 
*/
    for (k=0;k<4;k++){
      n[si[k]].e.remove1(is);
      v=si.vecino[k]; if (v<0) continue;
      esfera &sv=s[v];
      sv.replace_vecino(is,-1);
      sv.f.set(e_frontera); // => ahora si puede ser sliver de frontera
      if (sv.f.es_todos(syf)&&sv.f.noes(flag1)) 
        {rehacer+=v; sv.f.set(flag1);}
    }
    // si.f.set(e_borrado); s.remove(is);
    // cambia is0 por el ultimo no-sliver (si no es este)
    if (is!=s.len-1) {
      isl=s.len-1; while (s[isl].f.es_alguno(e_sliver|e_borrado)&&isl>=0) {isl--;}
      _revienta(isl<0);
      esfera &sl=s[isl]; // la ultima no borrada y no sliver
      const int *nsl=sl.n;
      for (k=0;k<4;k++) {
        n[nsl[k]].e.replace1(isl,is);
        v=sl.vecino[k]; if (v>=0) s[v].replace_vecino(isl,is);
      }
      s[is]=sl;
      sl.f.set(e_borrado); s.remove(isl); // por si no es la ultima
    }
    else s.len--;
    while (sb.len&&sb.last()==s.len-1) {sb.len--; s.len--;} // por si es la ultima
  }
  _savetime(f_slivers);
}

// vuela elementos chicos de frontera 
// con dos caras en la frontera
void malla::rm_f_slivers(double h3_factor){
  if (tipo.noes(m_vol)) return;
  _initime;
  mk_h_nn(); // si no hay h lo hace de los nn
  if (!ve) volumen();
  if (!vecino) mk_vecino();

  int i,j,nv,nc,vf,vj;
  double hi,vc;
  bool alguno=false;

  ordlist rehacer; rehacer.natural(e.len);
  while (rehacer){
    i=rehacer.last(); rehacer.len--;
    // dos caras en la frontera (ojo sandwiches!!)
    nc=e[i].nc(); const cpline vi=vecino[i];
    for (j=vf=0;j<nc&&vf<2;j++) if (vi[j]<0) vf++;
    if (vf!=2) continue;
    
    // testea volumen vs. h^3
    nv=e[i].nv(); const int *nei=e[i].n;
//  for (hi=n[nei[0]].h,j=1;j<nv;j++) hi+=n[nei[j]].h; hi/=nv; // h medio
    for (hi=n[nei[0]].h,j=1;j<nv;j++) set_min(hi,n[nei[j]].h);// h min
    vc=ve[i]/pown(hi,3);
    if (vc>h3_factor) continue;

    // agrega los vecinos a rehacer pero ojo si uno es el last por el swap de vuelae
    for (j=0;j<nc;j++){
      vj=vi[j];
      if (vj>0) {if (vj==e.len-1) rehacer+=i; else rehacer+=vj;}
    }
    vuelae(i,true); // marca frontera
    alguno=true;
  }
  
  if (alguno) {frontera.ini();nn.ini();}
  _savetime(f_slivers);
}

// engorda slivers
bool malla::engorda_slivers(){
  if ((tipo&(m_vol|m_sup|m_lin))!=m_vol) return false;// no es o no es puro
  if (!mk_vecino()) return false;
  if (!mk_h_nn()) return false;
  if (!mk_nn()) return false;
  if (!ve) volumen(); ce.ini(); re.ini();
  
  _initime;
  mark_slivers();

  int i,j,l,in,iej; int k; double h,h3,hmin_e;
  punto normal,dp,c0;
  pline rehacer(e.len/10);

  for (i=0;i<e.len;i++) if (e[i].f.es(e_sliver)) rehacer+=i;
  pline pasadas(e.len,0);
  while (rehacer){
    i=rehacer.last(); rehacer.len--;
    elemento &ei=e[i];
    if (pasadas[i]>5&&(ve[i]>0||pasadas[i]>10)) continue;    
    pasadas[i]++;ei.f.reset(e_sliver);
    for(k=0;k<4;k++){
      in=ei[k]; nodo &ni=n[in];
      if (ni.f.es(n_frontera)) continue;
      elemento c=ei.cara(k); // hacia afuera del elemento
      c0=n[c[0]];
      normal=((n[c[2]]-c0)%(n[c[1]]-c0)).dir();
      dp=ni-c0; 
      ni+=normal*(ni.h/10-(dp*normal));      
      // rehace h de nn minimo
      hmin_e=MAXREAL; const cpline& nni=nn[in];
      for (j=0;j<nni.len;j++) ni.distancia2_menor(n[nni[j]],hmin_e);
      ni.h=sqrt(hmin_e);
      for (l=0;l<nni.len;l++){
        nodo &nj=n[nni[l]]; hmin_e=MAXREAL; const cpline& nnj=nn[nni[l]];
        for (j=0;j<nnj.len;j++) nj.distancia2_menor(n[nnj[j]],hmin_e);
        nj.h=sqrt(hmin_e);
      }
      // arregla volumen
      h=ni.h; h3=h*h*h;
      const cpline& eni=ni.e;
      for(j=0;j<eni.len;j++) {
        iej=eni[j]; elemento &ej=e[iej];
        ve[iej]=volumen(ej);
        if (ve[iej]>h3*_volchico/6) {// arreglado
          if (ej.f.noes(e_sliver)) continue;
          // lo saca de la lista
          rehacer[rehacer.index(iej)]=rehacer.last(); rehacer.len--;
          ej.f.reset(e_sliver);
        }
        else{
          if (ej.f.es(e_sliver)) continue;
          // lo agraga a la lista
          if (rehacer) {rehacer+=rehacer[0]; rehacer[0]=iej;}
          else rehacer+=iej;
          ej.f.set(e_sliver);
        }
      }
      break;
    }
  }
  _savetime(slivers);
  return true;
}
