#include <cmath> // sqrt
//#define temporizar
#include "tiempo.h"
#include "pline.h"
#include "voronoi.h"

using namespace std;

#ifdef _DEBUG
  static bool dodebugout=false;
  #include "vadebug.h"
#endif

// El principal comentario sobre las rutinas es que los puntos se
//   mueven, alejandose de las esferas, para evitar degeneracion.
// Un nodo se puede mover varias veces, cada vez la mitad de la anterior,
//   por eso empieza siendo grande, para poder dividirlo por 2 unas
//   cuantas veces. Aunque un nodo puede estar en miles de esferas, un
//   nodo nuevo solo puede estar "justo" en muchas esferas si estoy
//   mallando justamente una esfera (o un circulo) con muchos puntos.

///////////////////////////////////////////////////////////////////////////////////////////
//   parametros "tocables"
///////////////////////////////////////////////////////////////////////////////////////////
static bool perturba=false; // (parece que no hace falta) ojo array regular cubico
static bool reponer=true; // reponer posicion de nodos (true)
static bool parcial=false; // pone solo si esta a mas de dnmin*h (false)
static bool test_close_to_boundary=false; // testea y en caso positivo no poner el nodo

// El programa no acepta puntos repetidos ni muy cercanos, dos nodos
//   muy cerca pueden dar tanto esferas minusculas como gigantes, la
//   tolerancia de posicion de nodos cercanos es local y no de la malla.
static const double _pb=1e-8; // epsilon con _pb de la diagonal del bounding box
static const double _ph=1e-4; // epsilon con _ph del h minimo de la malla

static double rmax=.7; //si r>rmax*h pone un punto en el centro de la esfera
static double dmax=1.4; //si dist entre nodos frontera > dmax*h => agrega pto medio
static double dnmin=.5;// d<dnmin*h => no agrega el punto (parcial)
// d_frontera<dfmin*h => no agrega el punto (test_close_to_boundary)
// con mas de .5 los centros de esfera quedan adentro del dominio
// pero eso deja muchos slivers de frontera (??)
static double dfmin=.51;
///////////////////////////////////////////////////////////////////////////////////////////

static int EPN;  // esferas por nodo (para prealocar, 100 o 1000)
static pline rehacer,rhindex; // del generador
static bool gen=false; // avisa que viene del generador
static bool preparado=false; // avisa que existe la lista de posiciones originales (false)
static double epsilon=0; // epsilon fijo, para perturbar posiciones
static punto pmin,pmax; // bbox, static pues lo puede crear uno y usar otro

// copia la posicion original y perturba
// ademas limpia los elementos de nodo y flag de frontera
void voronoi::prepara(){
  if (preparado) return;
  array1<nodo> &n=m->n; nori=n.len;
  if (pori) delete [] pori; pori=new punto[nori];
  punto pr;
  for (int i=0;i<nori;i++) {
    nodo &ni=n[i]; ni.e.clean(); ni.f.reset(n_frontera);
    pori[i]=ni;
    if (!perturba) continue;
    ni+=((NV==4)? pr.randdir(): pr.randdir2d())*epsilon;
  }
  preparado=true;
}

void voronoi::restaura(){
  if (!preparado) return;
  array1<nodo> &n=m->n;
  for (int i=0;i<nori;i++) n[i].setpos(pori[i]);
  delete [] pori; pori=0; nori=0;
  preparado=false;
}

//==============================================================
//                    celda de un nuevo nodo
//==============================================================
// Busca todas las esferas que contengan al nodo
//   y las reemplaza por esferas nuevas con el nodo
// No lo hace si test_close_to_boundary==true y esta muy cerca de la fronera

// Si el nodo in esta justo en la superficie de alguna esfera se mueve.

// Entra con la busqueda pesada resuelta: o se le dice que esfera (is) contiene
// al nuevo nodo (in) o se le pasa un nodo (nc) que "seguro" tenga una esfera
// que lo contiene (por ejemplo: el nodo mas cercano). El dato que no pasa es -1.

static const int flag12=flag1|flag2;
#define _resetflags\
  for (i=0;i<stested.len;i++) s[stested[i]].f.reset(flag12);

bool voronoi::cavidad(
  int in, // indice del nodo
  int is, // si se conoce: una esfera que lo contiene
  int nc, // si no: un nodo con una esfera que lo contiene (ej: el mas cercano)
  bool ini_o_clean  // asigna o devuelve memoria y sale
  ){

  int ixs,i,j,k,l,o;

// temporizadores individuales (comentar)
//static unsigned long int Tacum[5],Tcount=0; clock_t T0;

  //--------------------------------------------------------------------
  // listas prealocadas para no alocar y dealocar a lo loco :-)
  //   el generador de puntos puede dar miles de esferas por nodo
  //   asi que el tamanio de prealocacion es variable segun de donde viene

#ifndef _DEBUG
  static // en version debug no deja verlos => descomentar para que no lo haga static
#endif
          pline
    stested, // esferas testeadas
    sn,      // esferas que contienen al nodo
    tf,vf,   // caras y vecinos de la cavidad
    cn;      // datos sobre caras

  if (ini_o_clean){
//Tacum[0]=Tacum[1]=Tacum[2]=Tacum[3]=Tacum[4]=Tcount=0; // temporizadores (comentar)
    if (stested.size){ // devuelve memoria y sale
      stested.ini(); sn.ini(); tf.ini(); vf.ini(); cn.ini();
    }
    else { // inicializa
      stested.resize(EPN); sn.resize(EPN);
      tf.resize(EPN); vf.resize(EPN); cn.resize(EPN);
    }
    return false;
  }

  // borra el contenido
  stested.clean(); sn.clean(); tf.clean(); vf.clean(); cn.clean();
  //--------------------------------------------------------------------

//if (Tcount==10000){// temporizadores (comentar)
//  cout << '\n' << Tacum[0] << ' ' << Tacum[1] << ' ' << Tacum[2] << ' ' << Tacum[3] << ' ' << Tacum[4] << endl;
//  Tacum[0]=Tacum[1]=Tacum[2]=Tacum[3]=Tacum[4]=Tcount=0;
//} else Tcount++;

  const array1<nodo> &n=m->n;

  nodo ni=n[in];

  // epsilon inicial, despues s_epsilon se va dividiendo por 2
  // empiezo /2 para garantizar movimiento total < epsilon
  // no puede ser muy chico paa poder subdividirlo varias veces
  esfera::s_epsilon=Max(epsilon/2,esfera::s_epsilon);

  // flag1 de esfera => ya fue testeada
  // flag2 de esfera => contiene al nodo
  // flag3 de esfera => esta en rehacer

  // si hay una esfera (is) que contiene al nodo empieza con esa
  // si no, nc debe tener una esfera con el nodo
//T0=clock(); // temporizadores (comentar)
  if (is==nada){ // busca una esfera que contenga al nodo
    _revienta(nc==nada); // debe venir el nodo mas cercano
    // busca entre las esferas del nodo cercano
    const cpline &snc=n[nc].e;
    for (i=0;i<snc.len;i++) {
      is=snc[i];
      esfera &si=s[is];
      if(si.have_moving(ni)) break;
      else {si.f.set(flag1); stested+=is;}
    }
    if (i==snc.len) { // no encontro ninguna (no puede ser)
#ifdef _DEBUG
      debugout(n,s,ni,&snc,0);
      _revienta(1); // debug
#endif
      _resetflags; return false; //release
    }
  }
#ifdef _DEBUG
  else _revienta(!s[is].have_moving(ni)); // testeo si is realmente tiene al nodo
#endif

//Tacum[0]+=clock()-T0; T0=clock(); // temporizadores (comentar)

  // busca el resto de esferas con el nodo
  sn+=is; stested+=is; s[is].f.set(flag1|flag2);
  nodo ns0,ns1,ns2;
  for (i=0;i<sn.len;i++){
    is=sn[i]; const int *v=s[is].vecino;
    for (k=0;k<NV;k++){
      ixs=v[k]; // indice de la esfera vecina
      if (ixs<0&&test_close_to_boundary){
        // verifica si esta cerca de la frontera
        esfera &si=s[is];
        if (NV==3){
          ns0=n[si[(k+1)%3]]; ns1=n[si[(k+2)%3]];
          if (ni.distr(ns0,ns1)<ni.h*dfmin) // distancia con signo por posibles concavidades
            {_resetflags; return false;}
        }
        else{
          ns0=n[si[(k+1)%4]];
          if (k&1) {ns1=n[si[(k+2)%4]]; ns2=n[si[(k+3)%4]];}
          else {ns2=n[si[(k+2)%4]]; ns1=n[si[(k+3)%4]];}
          if (ni.distp(ns0,ns1,ns2)<ni.h*dfmin){// distancia con signo por posibles concavidades
            _resetflags; return false;
          }
        }
        // no esta cerca de la frontera
        ixs=-ixs-2;// indice de esfera exterior
      }
      if (ixs==-1) continue; // frontera
      esfera &sk=s[ixs];
      if (sk.f.es(flag1)) continue; // ya fue testeada
      sk.f.set(flag1); stested+=ixs;
      if (!sk.have_moving(ni)) continue; // no contiene al nodo
      sk.f.set(flag2); sn+=ixs;
    }
  }

#ifdef _DEBUG
  // retestea
  // esto daba mal con hay arrays cuadrados sin perturbar
  // pero ya lo arregle (digo, si falla no debe ser por eso)
  for (i=0;i<stested.len;i++) {
    ixs=stested[i];
    esfera &sk=s[ixs];
    if (sk.have(ni)) {
      if (sk.f.es(flag2)) continue;
    }
    else {
      if (sk.f.noes(flag2)) continue;
    }
    _revienta(1); // verificar epsilon (m->hayh m->hmin)
    if (dodebugout)
      debugout(n,s,ni,&sn,0);
  }
#endif
//Tacum[1]+=clock()-T0; T0=clock(); // temporizadores (comentar)

  // busca las caras exteriores de la cavidad
  // y elimina las esferas de sus nodos
  for (i=0;i<sn.len;i++){
    is=sn[i]; esfera &si=s[is]; ixs=(gen&&si.f.es(e_exterior))? -is-2 :is;
    const int *v=si.vecino; const int *ns=si.n;
    for (j=0;j<NV;j++){
      // n[ns[j]].e.remove1(is); // quita la esfera de sus nodos (pline caro, mejor ordlist?)
      cpline &enj=n[ns[j]].e; if(enj.replace1(is,enj[enj.len-1])) enj.len--; // un poco mas rapido
/*      //mas rapido!!!
      int* enj=n[ns[j]].e.vertex; int &enjlen=n[ns[j]].e.len; k=enjlen;
      while(k){if (enj[--k]==is) {enj[k]=enj[--enjlen]; break;}} //no es mas rapido!!*/
      o=v[j]; if (o<-1&&gen) o=-o-2;  //o=opuesto
      if (o>=0&&s[o].f.es(flag2)) continue; // no es exterior y esta en la cavidad
      vf+=o; // vecino por cara externa de la cavidad
       //indice de la cara en el vecino
      if (o>=0){
        k=s[o].index_vecino(ixs);
        vf+=k;
      } else vf+=-1;
      // caras de frontera ordenadas
      tf+=ns[(j+1)%NV];
      if (NV==4) {
        if (j&1) {tf+=ns[(j+3)%NV]; tf+=ns[(j+2)%NV];} // impar
        else {tf+=ns[(j+2)%NV]; tf+=ns[(j+3)%NV];} // par
      }
      else {tf+=ns[(j+2)%NV]; tf+=0;} // cualquier nro < n.len
    }
  }
//Tacum[2]+=clock()-T0; T0=clock(); // temporizadores (comentar)

#ifdef _DEBUG
  if (dodebugout){
    debugout(n,s,ni,&sn,&tf);
  }
#endif

  //descachea posicion (para define)
  n[in].setpos(ni);

  // hace una esfera con cada cara de la cavidad y la agrega a rehacer
  int c1,c2,iv,ic,ntf=tf.len/3,f1,f2,f3;
  esfera st; st[0]=in; if (gen) st.f.set(flag3);

  for (i=0;i<ntf;i++){
    // indice de la nueva esfera
    if (sn.len) is=sn.last(); // reutiliza (sn.len se updatea al final)
    else is=s.next_index(); // nueva

    ni.e+=is;

    st.f.reset(e_frontera);
    // vecino opuesto
    j=2*i; st.vecino[0]=o=vf[j];
    if (o>=0) {
      esfera &so=s[o];
      so.vecino[vf[j+1]]=is;
      if (gen&&so.f.es(e_exterior)) {
        st.vecino[0]=-o-2;
      }
    } else st.f.set(e_frontera);

    // cara
    j=3*i; f1=tf[j];f2=tf[j+1];f3=tf[j+2];
    st.n[1]=f1; st.n[2]=f2; st.n[3]=f3;
    if (!gen){
      if (f1>=n.len-nvirt||f2>=n.len-nvirt||f3>=n.len-nvirt)
        st.f.set(e_virtual);
      else st.f.reset(e_virtual);
    }
    st.define();

//    _revienta(m->es_sliver(st.n,ni.h));
#ifdef _DEBUG
  if (dodebugout)
    debugout(n,s,ni,0,&tf);
#endif
    _revienta(st.vt<-epsilon); // genero una esfera invertida
    // Esto de arriba jodia cuando se eliminaban las 
    //   esferas exteriores para el generador en f_delaunay (ahora no)
    // pero puede suceder siempre que haya concavidades!! 
    //   un cap lejano puede contener al nodo pero 
    //   quedando fuera de la cara exterior (aun 2D)

    // vecinos adyacentes
    for(k=0;k<NV-1;k++){ // cada cara adyacente (cara k+1)
      c1=st[(k+1)%(NV-1)+1]; // 2 3 1 en 3d o 2 1 en 2d
      n[c1].e+=is; // agrega la esfera al nodo
      if (NV==4) c2=st[(k+2)%3+1]; else c2=-1; // 3 1 2
      for (l=0;l<cn.len;l+=4) {
        if (cn[l]==c1&&cn[l+1]==c2&&cn[l+2]!=is) break;
      }
      if (l==cn.len) { // agrega la cara invertida
        if (NV==3) Swap(c1,c2);
        cn+=c2; cn+=c1; cn+=is; cn+=k+1;
      }
      else {// la cara ya fue puesta
        iv=cn[l+2]; ic=cn[l+3];
        st.vecino[k+1]=iv;
        s[iv].vecino[ic]=is;
        cn.remove(4,l);
      }
    }

    // agrega la esfera
    if (sn.len) { //reutiliza una que va a borrar
      if (gen&&s[is].f.noes(flag3))
        {rhindex[is]=rehacer.len; rehacer+=is;} // no estaba en rehacer
      s[is]=st; sn.len--;
    }
    else { // agrega una a s
      if (gen) {
        if (is==rhindex.len) rhindex+=rehacer.len; // nueva
        else rhindex[is]=rehacer.len; // hueco
        rehacer+=is;
      }
      s+=st;
    }
  }

//Tacum[3]+=clock()-T0; T0=clock(); // temporizadores (comentar)
  _revienta(cn.len);
#ifdef _DEBUG
  if (dodebugout)
    debugout(n,s,ni,0,&tf);
#endif
  // vuela las esferas no huecas que queden y las saca de rehacer
  int ixr,isl;
  for (i=0;i<sn.len;i++){
    is=sn[i];
    if (gen&&s[is].f.es(flag3)) {
      // is estaba en rehacer entonces la saca
      ixr=rhindex[is];
      // si es la ultima la vuela y si no swappea con la ultima
      if (ixr<rehacer.len-1){
        isl=rehacer.last();
        rehacer[ixr]=isl; rhindex[isl]=ixr;
      }
      rehacer.len--;
      s[is].f.reset(flag3);
    }
    s[is].f.set(e_borrado); s.remove(is);
  }

//Tacum[4]+=clock()-T0; // temporizadores (comentar)

  //descachea la lista de esferas
  n[in].e=ni.e;

  _resetflags;
  return true;
}

#undef _resetflags

//==============================================================
//                    esferas de voronoi
//==============================================================
// no agrega un punto si esta "muy cerca" de otro (ver dmin2)
// si parcial=true, no agrega si esta a menos de dnmin*h

bool voronoi::delaunay(){
  array1<nodo> &n=m->n;
  if (!n) return false;
  _initime;

  bool INFO_CL=m->INFO_CL;
  int nlen=n.len;
  int i,j,k,l,in=-1,in1,imin;
  double d2,dmin2,drep=-1.0;
  ordlist nni(1000);
  ordlist nborrado((parcial)? nlen/10 : 10); // nodos no agregados
  nodo nv; nv.h=MAXREAL/4; nv.f.set(n_virtual); // virtual

  EPN=100; cavidad(0,0,0,true); // preasigna listas (100 esferas por nodo)
  _push_1b(test_close_to_boundary,false); // NO HAY FRONTERA!

  // uso un cursor para randomizar la insersion de puntos
  // (cualquier orden encarece mucho, el generador "ordena" por capas y es carisimo)
  // la idea de poner nodos fijos o permanentes al principio es muy mala por la misma razon
  pline index; index.random_index(nlen);

  m->e_ini(); s.clean(); s.resize(Max(7*nlen,100));
  m->bbox(); punto deltap=m->pmax-m->pmin; double rnube=deltap.mod(); // por 2
  pmin=m->pmin-deltap/100; pmax=m->pmax+deltap/100; // para el octree
  punto cnube=(pmax+pmin)/2;
  // para mover nodos cuando estan cerca de la superficie de esferas
  if (m->hayh) epsilon=m->hmin*_ph; // si hay h con h minimo
  else epsilon=rnube*_pb; // si no, con el tamanio de la nube
  double dnmin2=epsilon*epsilon; // (compara distancias al cuadrado)

  // copia la posicion original y perturba (si preturbar=true)
  prepara();

  // octree propio
  octree o(n,pmin,pmax,NV-1,15);

  // primer nodo (no borrado)
  j=0; while (true&&j<n.len){
    in=index[j++];
    if (n[in].f.es(n_borrado)) nborrado+=in;
    else break;
  }
  _revienta(j==n.len);
  o.add(in);

  // Puntos Virtuales:
  // El tema aqui es no ponerlos muy cerca para que la union de
  //   tetraedros que no tengan puntos virtuales se parezcan lo
  //   mas posible al convex hull, y no ponerlos muy lejos para
  //   evitar que una esfera casi tangente a la superficie contenga
  //   otros puntos erroneamente (u obligue a moverlos al pedo)
  // Uso un factor 2 porque no me calienta el "verdadero" convex hull
  // Quizas fuese mejor poner un icosaedro en lugar de un tetraedro.
  {
    // tetraedro cuya esfera inscripta tiene radio 1
    nvirt=NV;
    double tx,ty; punto tetra[4];
    if (NV==3) {
      m->tipo.set(m_planaxy);
      esfera::setclass(2,n);
      tx=SQRT3;
      tetra[0]=punto(  0, 2);
      tetra[1]=punto(-tx,-1);
      tetra[2]=punto( tx,-1);
    }
    else{
      m->tipo.reset(m_planaxy);
      esfera::setclass(3,n);
      tx=SQRT2*SQRT3; ty=SQRT2;
      tetra[0]=punto(  0,   0, 3);
      tetra[1]=punto(  0,2*ty,-1);
      tetra[2]=punto( tx, -ty,-1);
      tetra[3]=punto(-tx, -ty,-1);
    }

    // tetraedro centrado en cnube y escalado con rnube
    // indices virtuales >= nlen
    n.resize(nlen+NV);
    nodo nv; nv.h=MAXREAL/4; nv.f.set(n_virtual|flag1); // no testea esferas
    punto pr;
    for (i=0;i<NV;i++) {
      nv=cnube+(rnube*tetra[i]);
      // muy perturbado (cara grande horizontal del dominio)
      nv+=((NV==3)? pr.randdir2d() : pr.randdir())*(rnube/100);
      n+=nv;
    }
    // esfera virtual
    esfera st; st.f.set(e_virtual); memset(st.vecino,-1,4*SZI);
    for (i=0;i<NV;i++) {st.n[i]=nlen+i; n[nlen+i].e+=0;}
    s+=st.define();
    cavidad(in,0); // primer nodo en la esfera virtual
  }

  /*{
    if (NV==4) {
      m->tipo.reset(m_planaxy);
      esfera::setclass(3,n);
      // icosaedro cuya esfera inscripta tiene radio 1
      nvirt=12;
      n[in].e.ini(); n[in].e.natural(20);
      n.resize(nlen+12);
      nv.e.resize(5); nv.e.len=5; // no testea esferas
      nv[0]=-0.483063277230303;nv[1]= 0.570270527140898;nv[2]=   1.07727266016107;nv.e[0]=14;nv.e[1]=2 ;nv.e[2]=11;nv.e[3]=19;nv.e[4]=18; nv=cnube+(rnube*nv); n+=nv;
      nv[0]= 0.803551549454465;nv[1]= 0.739883687956249;nv[2]=  0.567975491656873;nv.e[0]=10;nv.e[1]=9 ;nv.e[2]=16;nv.e[3]=11;nv.e[4]=14; nv=cnube+(rnube*nv); n+=nv;
      nv[0]= -0.19003186868894;nv[1]=  1.31061874011783;nv[2]=-0.0183265025143908;nv.e[0]=14;nv.e[1]=6 ;nv.e[2]=16;nv.e[3]=5 ;nv.e[4]=18; nv=cnube+(rnube*nv); n+=nv;
      nv[0]= -1.14016239604532;nv[1]= 0.324606026499249;nv[2]=-0.0622524965190791;nv.e[0]=4 ;nv.e[1]=0 ;nv.e[2]=19;nv.e[3]=5 ;nv.e[4]=18; nv=cnube+(rnube*nv); n+=nv;
      nv[0]=-0.707736401013585;nv[1]=-0.667960464361605;nv[2]=  0.647296522501082;nv.e[0]=4 ;nv.e[1]=19;nv.e[2]=2 ;nv.e[3]=1 ;nv.e[4]=15; nv=cnube+(rnube*nv); n+=nv;
      nv[0]=  0.46755242726372;nv[1]=-0.489650040600989;nv[2]=  0.967248243672383;nv.e[0]=9 ;nv.e[1]=2 ;nv.e[2]=11;nv.e[3]=13;nv.e[4]=1 ; nv=cnube+(rnube*nv); n+=nv;
      nv[0]= 0.685533562259055;nv[1]= 0.716434756308743;nv[2]= -0.682582023157354;nv.e[0]=6 ;nv.e[1]=8 ;nv.e[2]=3 ;nv.e[3]=10;nv.e[4]=16; nv=cnube+(rnube*nv); n+=nv;
      nv[0]=   1.1349313208961;nv[1]=-0.389034042720271;nv[2]=0.00966910469191139;nv.e[0]=3 ;nv.e[1]=9 ;nv.e[2]=10;nv.e[3]=13;nv.e[4]=7 ; nv=cnube+(rnube*nv); n+=nv;
      nv[0]= 0.296988455706378;nv[1]= -1.14925911852417;nv[2]= 0.0126500597519787;nv.e[0]=13;nv.e[1]=7 ;nv.e[2]=12;nv.e[3]=1 ;nv.e[4]=15; nv=cnube+(rnube*nv); n+=nv;
      nv[0]= -0.54398749178369;nv[1]= 0.519187388081939;nv[2]=  -1.03370534990735;nv.e[0]=8 ;nv.e[1]=17;nv.e[2]=6 ;nv.e[3]=0 ;nv.e[4]=5 ; nv=cnube+(rnube*nv); n+=nv;
      nv[0]=-0.693693875392516;nv[1]=-0.693177815774587;nv[2]= -0.694554061559543;nv.e[0]=4 ;nv.e[1]=17;nv.e[2]=0 ;nv.e[3]=12;nv.e[4]=15; nv=cnube+(rnube*nv); n+=nv;
      nv[0]=   0.4700821071415;nv[1]=-0.571527261015718;nv[2]=  -1.05302875623334;nv.e[0]=17;nv.e[1]=8 ;nv.e[2]=3 ;nv.e[3]=7 ;nv.e[4]=12; nv=cnube+(rnube*nv); n+=nv;
      rnube*=.8; // estan perturbados .1
      esfera st; int *v=st.vecino; st.f.set(e_virtual); i=nlen-1;
      st[0]=i+ 4;st[1]=  in;st[2]=i+11;st[3]=i+10;v[0]=17;v[1]=-1;v[2]= 5;v[3]= 4;st.define();s+=st;
      st[0]=i+ 9;st[1]=i+ 6;st[2]=  in;st[3]=i+ 5;v[0]= 2;v[1]=15;v[2]=-1;v[3]=13;st.define();s+=st;
      st[0]=i+ 1;st[1]=i+ 6;st[2]=i+ 5;st[3]=  in;v[0]= 1;v[1]=19;v[2]=11;v[3]=-1;st.define();s+=st;
      st[0]=i+12;st[1]=i+ 8;st[2]=i+ 7;st[3]=  in;v[0]=10;v[1]= 8;v[2]= 7;v[3]=-1;st.define();s+=st;
      st[0]=i+ 4;st[1]=i+ 5;st[2]=i+11;st[3]=  in;v[0]=15;v[1]= 0;v[2]=19;v[3]=-1;st.define();s+=st;
      st[0]=i+ 4;st[1]=i+10;st[2]=i+ 3;st[3]=  in;v[0]= 6;v[1]=18;v[2]= 0;v[3]=-1;st.define();s+=st;
      st[0]=i+ 3;st[1]=i+10;st[2]=i+ 7;st[3]=  in;v[0]= 8;v[1]=16;v[2]= 5;v[3]=-1;st.define();s+=st;
      st[0]=i+ 9;st[1]=i+ 8;st[2]=i+12;st[3]=  in;v[0]= 3;v[1]=12;v[2]=13;v[3]=-1;st.define();s+=st;
      st[0]=i+12;st[1]=i+ 7;st[2]=i+10;st[3]=  in;v[0]= 6;v[1]=17;v[2]= 3;v[3]=-1;st.define();s+=st;
      st[0]=i+ 2;st[1]=i+ 8;st[2]=i+ 6;st[3]=  in;v[0]=13;v[1]=11;v[2]=10;v[3]=-1;st.define();s+=st;
      st[0]=i+ 2;st[1]=i+ 7;st[2]=i+ 8;st[3]=  in;v[0]= 3;v[1]= 9;v[2]=16;v[3]=-1;st.define();s+=st;
      st[0]=i+ 2;st[1]=i+ 1;st[2]=  in;st[3]=i+ 6;v[0]= 2;v[1]= 9;v[2]=-1;v[3]=14;st.define();s+=st;
      st[0]=i+ 9;st[1]=i+12;st[2]=i+11;st[3]=  in;v[0]=17;v[1]=15;v[2]= 7;v[3]=-1;st.define();s+=st;
      st[0]=i+ 9;st[1]=  in;st[2]=i+ 6;st[3]=i+ 8;v[0]= 9;v[1]=-1;v[2]= 7;v[3]= 1;st.define();s+=st;
      st[0]=i+ 2;st[1]=i+ 3;st[2]=  in;st[3]=i+ 1;v[0]=18;v[1]=11;v[2]=-1;v[3]=16;st.define();s+=st;
      st[0]=i+ 9;st[1]=  in;st[2]=i+11;st[3]=i+ 5;v[0]= 4;v[1]=-1;v[2]= 1;v[3]=12;st.define();s+=st;
      st[0]=i+ 2;st[1]=  in;st[2]=i+ 3;st[3]=i+ 7;v[0]= 6;v[1]=-1;v[2]=10;v[3]=14;st.define();s+=st;
      st[0]=i+12;st[1]=i+10;st[2]=i+11;st[3]=  in;v[0]= 0;v[1]=12;v[2]= 8;v[3]=-1;st.define();s+=st;
      st[0]=i+ 4;st[1]=i+ 3;st[2]=i+ 1;st[3]=  in;v[0]=14;v[1]=19;v[2]= 5;v[3]=-1;st.define();s+=st;
      st[0]=i+ 4;st[1]=  in;st[2]=i+ 1;st[3]=i+ 5;v[0]= 2;v[1]=-1;v[2]= 4;v[3]=18;st.define();s+=st;
    }
  }*/

    
  // loop principal
  for(;j<nlen;j++){
    if (((j>>10)<<10)==j) {
      if (INFO_CL) cout << "\rnodos: " << j << flush;
    }
    in=index[j];
    if (n[in].f.es(n_borrado)){
      nborrado+=in;
      continue;
    }
    punto pj=n[in];
    // el nodo nuevo debe estar en alguna esfera del nodo mas cercano
    // aun cuando haya mas de un "mas cercano" (o mas de cuatro)
    // debe estar en alguna de c/u (diagrama de Voronoi)

    // No puede haber nodos sin esferas en el octree,
    // si agrega uno que no va => debe sacarlo

    // busca el nodo puesto mas cercano
    i=o.add(in);// busca uno cercano en el octree
    // busca el mas cercano
    if (parcial) dnmin2=pown(n[in].h*dnmin,2);
    while(1){
      const cpline &sni=n[i].e; // ojo si el octree devuelve uno sin esferas!!!!!!!!!!!!!!!!
      // de todos los vecinos naturales del i y el i, busca el mas cercano
      nni.clean();
      for (k=0;k<sni.len;k++) {
        const int *ns=s[sni[k]].n;
        for (l=0;l<NV;l++) nni+=ns[l];
      }
      imin=nni[0]; dmin2=pj.distancia2(n[imin]);
      for (k=1;k<nni.len;k++){
        in1=nni[k];
        d2=pj.distancia2(n[in1]);
        if (d2<dmin2) {dmin2=d2; imin=in1;}
      }
      if (imin==i) break; // este es el mas cercano
      i=imin;
    }

//-------------------------------------------------------------------
// para tocar y experimentar
// el nodo que se va a introducir es muy cercano a otro que ya esta
// �que significa "muy cercano"? �que hay que hacer?
    if (dmin2<dnmin2){ // es muy cercano
      if (n[in].f.es(n_permanente)){ // este nodo no puede volar
        if (!n[i].f.es(n_permanente)){ // si los dos son permanentes lo pone igual
          // el que ya estaba no es permanente => vuela
          // no quiero modificar el nro de orden (in) del nodo
          // => le pongo las esferas y la posicion del no-permanente
          // la posicion es para calcular las esferas nuevas, pero luego se restaura al final
          n[in].e=n[i].e; n[i].e.clean(); const cpline& en=n[in].e; // reemplaza lista de esferas
          for (k=0;k<en.len;k++) s[en[k]].replace(i,in); // reemplazo en las esferas
          n[in].setpos(n[i]); // despues de todo estan muy cerca (si no puede haber vol's negativos)
          n[i].f.set(n_borrado); nborrado+=i; o.remove(i); // borra el i
          if (!parcial) set_max(drep,pj.distancia(n[i]));
          continue;  //no pone el nodo
        }
      }
      else if (n[i].f.es(n_h)&&n[in].f.noes(n_h)){ // vuela el nodo h (si es el i)
        n[in].e=n[i].e; n[i].e.clean(); const cpline& en=n[in].e; // reemplaza lista de esferas
        for (k=0;k<en.len;k++) s[en[k]].replace(i,in); // reemplazo en las esferas
        n[in].setpos(n[i]); // despues de todo estan muy cerca
        n[i].f.set(n_borrado); nborrado+=i; o.remove(i); // borra el i
        if (!parcial) set_max(drep,pj.distancia(n[i]));
        continue;  //no pone el nodo
      }
      else{ // caso estandar
        n[in].f.set(n_borrado); nborrado+=in; o.remove(in); // borra el in
        if (!parcial) set_max(drep,pj.distancia(n[i]));
        continue;  //no pone el nodo
      }
    }
//------------------------------------------------------------------

    // agrega el nodo in a la triangulacion, i es el mas cercano
    // alguna esfera del i contiene al nodo
    if (!cavidad(in,nada,i)){ // reemplaza las esferas por nuevas
      // No hay frontera ni se hace el test_close_to_boundary,
      //   entonces hubo algun error (grave)
      _revienta(1); // debug
      n[in].f.set(n_borrado); nborrado+=in; o.remove(in); // release
    }
  }

  // limpia
  cavidad(0,0,0,true); // inicializa listas
  
  if (INFO_CL){
    cout << "\rnodos: " << nlen << "\tesferas: " << s.len-s.borrado.len << endl;
  }

  if (reponer) restaura(); // devuelve posicion original
/*
pline tempo(n[107119].e); tempo+=n[nlen+2].e;
debugout(n,s,n[107119],&tempo,0);
cout << "out: " << tempo.len << " tetras" << endl;
*/

  rm_virtuales(); // elimina las esferas y nodos virtuales

  // squeeze de los nodos no agregados
  // podria pensarse en devolverlos sin esferas y ya esta (usuario-dependiente)
  // swappeo cada borrado con el ultimo nodo y reduzco la lista de nodos
  if (nborrado&&m->o) {delete m->o; m->o=0;} // no hay borrados en el octree de voronoi
#ifndef _NO_DELETE_NODES
  while (nborrado.len){
    // (virtuales eliminados)
    int ib=nborrado.last(),il=n.len-1,&nh=m->nodosh;
    // los nodos de h complican el proceso
    // si hay y debo eliminar un no-h swappeo antes, el ultimo no-h, por el ultimo (h)
    if (nh){
      if (ib<n.len-nh) { // es un no-h
        m->swap_n(il,il-nh);
        if (ib==il-nh) ib=il; // era justo el ultimo no-h
      }
      else nh--; // es un nodo de h
    }
    if (ib==il) {
      n.len--; nborrado.len--;
      if (nori) nori--;
      if (m->ndir&&m->dir) m->dir.len--; // dir de voronoi es por elementos
      continue;
    }
    // ya no es el ultimo nodo
    // nborado es ordlist=> el ultimo no es borrado
    // swappeo con el ultimo nodo
    nodo &nl=n[il];
    const cpline &snl=nl.e;
    for (i=0;i<snl.len;i++) s[snl[i]].replace(il,ib);
    n[ib]=nl; nborrado.len--; n.len--;
    if (nori) {pori[ib]=pori[il]; nori--;}
    if (m->ndir&&m->dir) {m->dir[ib]=m->dir[il]; m->dir.len--;} // dir de voronoi por elementos
  }
#endif
  if (drep>0) {
    m->add_warning(Repeated_Nodes); m->add_warning(drep);
  }

  if (INFO_CL){
    cout << "\rnodos: " << nlen << "\tesferas: " << s.len-s.borrado.len << endl;
  }
  _savetime(delaunay);
  _pop_1b(test_close_to_boundary); //devuelve el valor original
  return true;
}


//==============================================================
//           delaunay de los nodos de frontera (y de h)
//==============================================================
bool voronoi::f_delaunay(){
  _initime;
  array1<nodo> &n=m->n;

  if (!mk_hydir()) // elementos y datos de la frontera
    {_savetime(error); return false;}
  m->epsilon=m->hmin*_ph; // aca y no en mk_hydir() porque _ph es local

  _push_1b(parcial,false); // no parcial
  _push_1b(reponer,false); // no reponer posicion hasta marcar las exteriorers
  if (!delaunay()){
    _pop_1b(reponer); _pop_1b(parcial);
    restaura();
    _savetime(error); return false;
  }
  _pop_1b(reponer); _pop_1b(parcial);

  // todos los no n_h son frontera
  int i;
  for (i=0;i<n.len-m->nodosh;i++) n[i].f.set(n_frontera);
  frontera();// marca esferas exteriores
  for (;i<n.len;i++) n[i].f.reset(n_h); m->nodosh=0;

  // Si se borran las esferas exteriores para el generador se producen errores de precision
  if (reponer) {
    restaura(); // devuelve posicion
    borra_exteriores(); //esferas exteriores
    rm_f_slivers(); // slivers de frontera
    squeeze();
  }

/*#ifdef _DEBUG
  graba_dat();
#endif*/

  if (m->INFO_CL) {
    cout << "\r\tDelaunay de frontera: "
         << n.len << " nodos y "
         << s.len-s.borrado.len << " elementos" << endl;
  }

  _savetime(f_delaunay);
  return true;
}


//==============================================================
//      a que tetraedro pertenece el punto
//  (debe ser un conjunto convexo de tetraedros)
//==============================================================
// la ordlist tetras fue agregada para ver por que reentraba en
// la lista. Se puede (debe?) mejorar
int voronoi::de_que_tetra(const punto &p, int primer_tetra, double *ff){
  
  if (p[0]<pmin[0]||p[0]>=pmax[0]) return -1;
  if (p[1]<pmin[1]||p[1]>=pmax[1]) return -1;
  if (NV==4&&(p[2]<pmin[2]||p[2]>=pmax[2])) return -1;
//static int maxlen=0;
  int k,kmin,is=primer_tetra;
  double ffmin;
  ordlist tetras(10);
  while (is>=0){
    tetras+=is;
    const esfera &si=s[is];
    si.fforma(p,ff); // centro de si en el tetraedro de st
    for (ffmin=ff[0],kmin=0,k=1;k<NV;k++)
      if (set_min(ffmin,ff[k])) kmin=k;
    if (ffmin>-ERRADM) break; // lo encontro
    is=si.vecino[kmin];
    if (is==-1) break; // esta fuera de todo
    if (tetras.have(is)){
#ifdef _DEBUG
      _revienta(1); // debug      
      if (dodebugout){
        debugout(m->n,s,p,&tetras,0,true);
      }
#endif
      break; // hay algo raro (vol negativo?) pero tomatelas
    }
  }

//if (set_max(maxlen,tetras.len)) cout << maxlen << endl;

  return is;
}



//==============================================================
//                    generador
//==============================================================
bool voronoi::mk_puntos(bool h_radio){
  _initime;

  array1<nodo> &n=m->n;
  bool INFO_CL=m->INFO_CL;

  _push_1b(reponer,false);    // no reponer posicion
  _push_1b(parcial,false);    // no parcial

  int is,ist,in,k;
  double hm,ff[4];
  esfera si;
  punto pt;
  ordlist tetras(16);

  if (h_radio){
    if (!m->hayh){
      static const double hfijo=1;
      for (in=0;in<n.len;in++) {n[in].h=hfijo; n[in].f.set(1);}
      m->hmin=hfijo; m->hayh=true;
    }
    else for (in=0;in<n.len;in++) n[in].f.set(1);
    m->hayfn=true;
  }

  if ((h_radio&&!delaunay())||(!h_radio&&!f_delaunay())) {
    _pop_1b(reponer); _pop_1b(parcial);
    restaura();
    _savetime(error); return false;
  }

/*
#ifdef _DEBUG
  graba_dat();
  return true;
#endif
*/
  // r>rmax*h => pone el nodo
  double oldrmax=rmax; // para mantener el default
  // rmax mide indirectamente el peor angulo
  rmax=1.1*sqrt((float)(NV-1))/2; // cuadrado - cubo
  double h_r=1.5/rmax; // solo se usa si h_radio==true

  EPN=1000; cavidad(0,0,0,true); // preasigna listas

  slivers_marked=false;

  // flag3 => esta en rehacer

  // todas las interiores deben rehacerse
  gen=true;
  _push_1b(test_close_to_boundary,true); // necesario (puede poner nodos justito en la frontera)
  rehacer.resize(10*s.len); rhindex.resize(10*s.len); rhindex.len=s.len;
  static const int e_norehacer=e_exterior|e_borrado;
  for (is=0;is<s.len;is++) {
    if (s[is].f.es_alguno(e_norehacer)) continue;
    rhindex[is]=rehacer.len; rehacer+=is; s[is].f.set(flag3);
  }

  // loop principal
  while (rehacer.len) {
    is=rehacer.last(); rehacer.len--; s[is].f.reset(flag3);
    si=s[is];

    // el h de comparacion "debe ser" el minimo de entre los nodos
    // dice si hay algun nodo a distancia razonable del que tiene el minimo h
    for (hm=n[si[0]].h,k=1;k<NV;k++) set_min(hm,n[si[k]].h); // minimo

    _revienta(hm<ERRADM);

    if (si.r<=rmax*hm) continue;

    // la esfera es grande

    // busca si hay un tetraedro conteniendo el centro (si el centro esta en el dominio)
    pt=si.c;
    ist=de_que_tetra(pt,is,ff);
    if (ist<0)
      continue;
    if (ist!=is) si=s[ist];

    // calcula h para el nodo 
    if (!h_radio){
      // interpolado del tetraedro que lo contiene
      // puede haber un ff negativo (fuera de la malla? elemento invertido?) y el h dar negativo
      for (hm=0,k=0;k<NV;k++) hm+=n[si[k]].h*ff[k];
      _revienta(hm<0);
    }
    else {      
      bool interno=true;
      for (hm=0,k=0;k<NV;k++) {
        hm+=n[si[k]].h*ff[k]; 
        if (interno&&n[si[k]].f.es(1)) interno=false;
      }
      if (!interno) 
        hm=Max(si.r*h_r,hm); // si no, hm
    }
    // agrega el nodo
    in=n+=pt; n[in].h=hm;

    // en cavidad verifica si esta cerca de la superfice
    // y si es asi no lo agrega
    if (!cavidad(in,ist)) // redelauniza
      n.len--;

    if (INFO_CL&&(n.len>>10)<<10==n.len) {cout << "\rnodos: " << n.len; cout.flush();}
  }

/*#ifdef _DEBUG
  graba_dat();
#endif*/

  if (!h_radio) borra_exteriores();
  rm_f_slivers(); // elimina slivers de frontera
  cavidad(0,0,0,true); // limpia listas
  gen=false; rehacer.ini(); rhindex.ini();

  if (INFO_CL){
    cout << "\rnodos: " << n.len << "\tesferas: " << s.len-s.borrado.len << endl;
  }

  squeeze();
  _pop_1b(reponer); _pop_1b(parcial); _pop_1b(test_close_to_boundary);
  restaura(); // devuelve posicion
  rmax=oldrmax;


#ifdef _DEBUG
//  graba_dat();
//  debugout(n,s,n[107119],&(n[107119].e),0);
#endif

  _savetime(puntos);
  return true;
}


// Hace las esferas sacando los puntos muy cercanos y
// agrega puntos en las esferas medianas
// Requiere la malla previa para interpolar y para conocer las aristas largas
bool voronoi::refina_esferas(malla &vieja,double alpha){
  _initime;

  int is,ist,in,k,ie,elen=vieja.e.len;
  int ic,nc,iv,nvc,ix[4],in1,in2,nve;
  double d,hm,ff[4];
  esfera si;

  if (!elen||!vieja.mk_vecino()||(!vieja.hayh&&!vieja.mk_h_nn_min(true)))
  {_savetime(error); return false;}

  array1<nodo> &n=m->n; //nodos de la malla actual (vieja movidos)
  bool INFO_CL=m->INFO_CL;

  // OJO: donde dice frontera no es interfase, es frontera del dominio, de la malla vieja
  // para trabajar con interfases hay que hacer mas trucos

  // elimina nodos cercanos a la frontera
  const array1<elemento> &e=vieja.e; // elementos de la vieja
  const array1<cpline> &vecino=vieja.vecino;
  for (k=-1,ie=0;ie<elen;ie++,k=-1){
    const cpline &vi=vecino[ie]; nc=e[ie].nc(); // cant de caras
    while ((k=vi.index(-1,k+1)) != nc){
      const elemento &ei=e[ie]; nve=ei.nv(); // elemento de frontera
      // verifica si los nodos que no estan en la cara estan cerca de la frontera
      const elemento cf=e[ie].cara(k); // cara de frontera
      for (iv=0;iv<nve;iv++){
        in=ei[iv];
        if (n[in].f.es_alguno(n_permanente|n_frontera)) continue;
        // distancia con signo por posibles concavidades
        if (NV==3) d=n[in].distr(n[cf[0]],n[cf[1]]);
        else d=-(n[in].distp(n[cf[0]],n[cf[1]],n[cf[2]]));
        if (d>n[in].h*dfmin) continue;
        n[in].f.set(n_borrado); // finalmente vuela en delaunay()
      }
    }
  }

  // se agregan nodos nuevos en las aristas largas de frontera
  for (ie=0;ie<elen;ie++){
    const elemento &ei=vieja.e[ie]; nc=ei.nc();
    const cpline &vecino=vieja.vecino[ie];
    for (ic=0;ic<nc;ic++){
      if (vecino[int(ic)]>=0) continue;
      // cara de frontera
      if (NV==4){ //3D
        ei.cara(ic,nvc,ix);
        for (iv=0;iv<nvc;iv++){
          const nodo &n1=n[in1=ei[ix[iv]]],&n2=n[in2=ei[ix[(iv+1)%nvc]]];
          d=n1.distancia2(n2);
          hm=(n1.h+n2.h)/2;
          if (d<pown(dmax*hm,2)) continue;
          // arista larga => agrega el nodo
          in=m->nodefrom(in1,in2); // ojo flags
        }
      }
      else { // 2D
        const nodo &n1=n[in1=ei[ic]], &n2=n[in2=ei.npos(ic)];
        d=n1.distancia2(n2);
        hm=(n1.h+n2.h)/2;
        if (d<pown(dmax*hm,2)) continue;
        // arista larga => agrega el nodo
        in=m->nodefrom(in1,in2); // ojo flags
      }
    }
  }

  // delaunay parcial: no pone nodos cercanos
  parcial=true; reponer=false;
  if (!delaunay()){
    parcial=false; reponer=true; restaura();
    _savetime(error); return false;
  }
  // OJO: en delaunay volaron nodos y hubo swaps
  // asi que los nodos de la nueva ya no son los de la vieja movidos

  EPN=1000; cavidad(0,0,0,true); // preasigna listas

  slivers_marked=false;

  // flag3 => esta en rehacer

  // todas deben rehacerse
  gen=true;
  rehacer.resize(10*s.len); rhindex.resize(10*s.len); rhindex.len=s.len;
  for (is=0;is<s.len;is++) {
    if (s[is].f.es(e_borrado)) continue;
    rhindex[is]=rehacer.len; rehacer+=is; s[is].f.set(flag3);
  }

  // loop principal
  test_close_to_boundary=true;
  while (rehacer.len) {
    is=rehacer.last(); rehacer.len--; si=s[is]; s[is].f.reset(flag3);

    // calcula h
    // si antes se hizo alphashape, se puede utilizar h medio
    // si no, conviene utilizar h minimo.
//    for (k=1,hm=n[si[0]].h;k<NV;k++) hm+=n[si[k]].h; hm/=NV;  // h medio
    for (k=1,hm=n[si[0]].h;k<NV;k++) set_min(hm,n[si[k]].h);  // h min
    _revienta(hm<ERRADM);

    if (si.r<rmax*hm) continue;  //esfera pequenia, sale

    // la esfera es grande
    // si ya paso por alfa, las esferas mayores que alfa no estan,
    // si no hay que preguntar si la esfera es mayor que alfa
    // if(si.r>alfa*hm) continue; //esfera demasiado grande
    // seguir con una grande sirve para rellenar huecos internos

    // agrega el centro si esta dentro del dominio de la vieja
    // pero intenta primero interpolar con los nodos de la nueva

    // verifica si hay algun tetraedro (de la malla actual) que contenga el centro
    nc=si[0];if (nc>=vieja.n.len) nc=-1; ie=-1; vieja.dequien(si.c,ie,nc);
    if ((ist=de_que_tetra(si.c,is,ff))<0){ // no hay
      // verifica si esta dentro de la malla vieja
      if (ie<0) continue; // ni en la nueva ni en la vieja => no lo pone
      // dentro => interpola con la vieja (le agrega el nodo)
      vieja.fforma(ie,si.c,ff);
      in=n+=vieja.n[vieja.nodefrom(ie,ff)];
    }
    else { // dentro de ist en la nueva => interpola con la nueva      
      if (ie<0) continue; // si no estaba dentro de la vieja, no lo agrega
      // habra que agregar el nodo a la vieja en este caso?
      in=m->nodefrom(NV,s[ist].n,ff);
    }

    // en cavidad verifica si esta cerca de la superfice
    // y si es asi lo saca, habria que ponerlo "en" la frontera, pero ya se agregaron ahi
    if (!cavidad(in,is)) // redelauniza, mete el centro de la esfera en la malla
      n.len--; // puede haberlo agregado inutilmente en la vieja pero eso no importa

    if (INFO_CL&&(n.len>>10)<<10==n.len) {cout << "\rnodos: " << n.len; cout.flush();}
  }

  // limpia
  vieja.ini(); //inicializa malla vieja (limpia memoria)
  cavidad(0,0,0,true); // limpia listas
  gen=false; rehacer.ini();

  if (INFO_CL) cout << "\rnodos: " << n.len << endl;

//#ifdef _DEBUG
//  graba_dat();
//#endif

  alpha_shape(alpha); // alpha shape
  rm_f_slivers(); // los de frontera si o si (los otros afuera de aca)
  squeeze();
  parcial=false; reponer=true; restaura(); // devuelve posicion

//#ifdef _DEBUG
//  graba_dat();
//#endif

  _savetime(refina_esferas);
  return true;
}



// Hace las esferas sacando los puntos muy cercanos y
// agrega puntos en las esferas medianas y en los slivers
//nodes_from are the nodes from which we will interpolate
//shape functions are the values to be used in the interpolation
bool voronoi::refina_esferas(double alpha, array1<double>& new_h, pline& nodes_from, array1<double>& shape_functions){
  _initime;

  int k,is,in;
  double hm,ff[4];

  array1<nodo> &n=m->n; //nodos de la malla actual
  bool INFO_CL=m->INFO_CL;

  // delaunay parcial: no pone nodos cercanos
#ifndef _NO_DELETE_NODES
  parcial=true;
#else
  parcial=false;
#endif
  reponer=false;
  delaunay();
  alpha_shape(alpha);
  slivers_marked=false;

  // agerga nodos en el centro de las esferas grandes
  //assign new h
  for(int k=0; k<n.len; k++) n[k].h = new_h[k];

  EPN=1000; cavidad(0,0,0,true); // preasigna listas
  // flag3 => esta en rehacer
  // todas deben rehacerse
  gen=true;
  rehacer.resize(10*s.len); rhindex.resize(10*s.len); rhindex.len=s.len;
  for (is=0;is<s.len;is++) {
    if (s[is].f.es(e_borrado)) continue;
    rhindex[is]=rehacer.len; rehacer+=is; s[is].f.set(flag3);
  }
  // loop principal
  test_close_to_boundary=true;
  while (rehacer.len) {
    is=rehacer.last(); rehacer.len--; esfera &si=s[is]; si.f.reset(flag3);
    // calcula h
    // si antes se hizo alphashape, se puede utilizar h medio
    // si no, conviene utilizar h minimo.
    for (k=1,hm=n[si[0]].h;k<NV;k++) hm+=n[si[k]].h; hm/=NV;  // h medio
//    for (k=1,hm=n[si[0]].h;k<NV;k++) set_min(hm,n[si[k]].h);  // h min
    _revienta(hm<ERRADM);

    if (si.r<rmax*hm) continue;  //esfera pequenia, sale

    // la esfera es grande
    // si ya paso por alfa, las esferas mayores que alfa no estan,
    // si no hay que preguntar si la esfera es mayor que alfa
    // if(si.r>alfa*hm) continue; //esfera demasiado grande
    // seguir con una grande sirve para rellenar huecos internos

    // agrega el centro si esta dentro del dominio

    // verifica si hay algun tetraedro que contenga el centro
    if ((is=de_que_tetra(si.c,is,ff))<0) continue;
    in=m->nodefrom(NV,s[is].n,ff);
    for(k=0;k<NV;k++) {nodes_from+=s[is][k]; shape_functions+=ff[k];}

    // redelauniza, mete el centro de la esfera en la malla
    // en cavidad verifica si esta cerca de la superfice
    // y si es asi, devuelve false y lo saca,
    // habria que ponerlo "en" la frontera...
    if (!cavidad(in,is)) {
      n.len--;
      nodes_from.len-=NV;
      shape_functions.len-=NV;
      continue;
    }


    if (INFO_CL&&(n.len>>10)<<10==n.len) {cout << "\rnodos: " << n.len; cout.flush();}
  }

  // limpia
  cavidad(0,0,0,true); // limpia listas
  gen=false; rehacer.ini();

  if (INFO_CL) cout << "\rnodos: " << n.len << endl;

//#ifdef _DEBUG
//  graba_dat();
//#endif

  rm_f_slivers(); // los de frontera si o si (los otros afuera de aca)



  squeeze();
  parcial=false; reponer=true; restaura(); // devuelve posicion

//#ifdef _DEBUG
//  graba_dat();
//#endif

  _savetime(refina_esferas);
  return true;
}

