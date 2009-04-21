// genera una malla cortando los tetraedros de otra

//#define temporizar
#include "tiempo.h"
#include "malla.h"

static malla *mc, *mm; static pline *_eplane;
static int nodo_inicio, nodo_final;
static int ne[4], qne;
static double nviso;
static double *v=0;
static int step=0;


// funciones estaticas para hacer el corte

// agrega un nodo testeando en el octree
static void agrega_n(double alpha){
  malla &m=*mm;
  malla &corte=*mc;
  array1<nodo> &n=corte.n;

  bool puesto;

  if (alpha<ERRADM) alpha=0; else if(alpha>1-ERRADM) alpha=1;
  int q,nc, in=n+=nodo(m.n[nodo_inicio],m.n[nodo_final],alpha);
  ne[qne]=in; qne++; //avisa al elemento que tiene este nodo
  corte.o->add_no_rep(in,nc,puesto,corte.epsilon);//verifica en el octree si ya estaba
  if(puesto) return;
  //no agrego el nodo debido a que ya existe en la malla (o hay uno muy cerca)
  n.len--;
  //verifico si el nodo cercano era uno de corte de este elm
  //(dos nodos muy cercanos en aristas distintas de este mismo elemento)
  for (q=0;q<qne-1;q++)
    if(ne[q]==nc) break; // el nodo ya existe en el elemento
  if(q==qne-1) ne[qne-1]=nc; // no estaba en el elemento
  else qne--; // ya estaba
}

// ya se sabe que el plano corta a la arista, averigua donde y pone el resultado en n
static void mk_intersection(){
  //cual es punto asociado al valor nodal seteado con la rueda del mouse?
  double alpha,divisor_alpha=(v[nodo_final]-v[nodo_inicio]);
  if (fabs(divisor_alpha)<ERRADM) alpha=0; //todo coincide
  else alpha=(nviso-v[nodo_inicio])/divisor_alpha;
  agrega_n(alpha);
}

// verifica el corte de todos los elementos supuestos tetraedros
#define _nf(ix) nodo_final=ei[ix]
#define _ni(ix) nodo_inicio=ei[ix]

static void cut_tetras(){
  bool n0_arriba,dif1,dif2,dif3;
  int i,j,nt;

  malla &corte=*mc; malla &m=*mm; pline& eplane=*_eplane;
  array1<nodo> &n=corte.n; array1<elemento> &e=corte.e;

  corte.pmin=m.pmin; corte.pmax=m.pmax;
  corte.hayv=m.hayv;corte.nvmin=m.nvmin;corte.nvmax=m.nvmax;
  corte.hayh=m.hayh;corte.hmin=m.hmin;
  corte.epsilon=m.epsilon;
  corte.o=new octree(n,corte.pmin,corte.pmax);

  eplane.clean();
  for(i=0;i<m.e.len;i+=step){
    const elemento &e1=m.e[i]; nt=e1.nt();
    for(j=0;j<nt;j++){ // para cada tetraedro del elemento
      const int *ei=e1.tetra(j);
      qne=0;
      /*   Verifica posicion de nodos respecto al corte
      Si los nodos 1, 2 y 3 tienen el mismo signo que el 0 no hay interseccion.
      Si los nodos 1, 2 y 3 tienen distinto signo que el 0 es el caso A.
      Si solo el nodo 1 es distinto, es el caso B.
      Si solo el nodo 2 es distinto, es el caso C.
      Si solo el nodo 3 es distinto, es el caso D.
      Si solo el nodo 1 es igual, es el caso E.
      Si solo el nodo 2 es igual, es el caso F.
      Si solo el nodo 3 es igual, es el caso G.
      */
      n0_arriba=v[ei[0]]>nviso;
      dif1=(v[ei[1]]>nviso)!=n0_arriba;
      dif2=(v[ei[2]]>nviso)!=n0_arriba;
      dif3=(v[ei[3]]>nviso)!=n0_arriba;

      _ni(0);
      if(dif1){
        _nf(1); mk_intersection(); // 01
        if(dif2){
          _nf(2); mk_intersection(); // 01 y 02
          if(dif3){
            _nf(3); mk_intersection(); // 01, 02 y 03 => A
          }
          else {// 01 y 02 pero no 03 ==> G
            _ni(3); mk_intersection(); // corta 32
            _nf(1); mk_intersection(); // corta 31
          }
        }
        else if(dif3){
          _ni(2); mk_intersection(); // corta 21
          _nf(3); mk_intersection(); // corta 23
          _ni(0); mk_intersection(); // 01 y 03 pero no 02 => F
        }
        else{// no corta 02 ni corta 03 => B
          _ni(2); mk_intersection(); // corta 21
          _ni(3); mk_intersection(); // corta 31
        }
      }
      else if(dif2){
        _nf(2); mk_intersection(); // no corta 01 pero corta 02
        if(dif3){
          _nf(3); mk_intersection(); // corta 03 => E
          _ni(1); mk_intersection(); // corta 13
          _nf(2); mk_intersection(); // corta 12
        }
        else{// no corta 01 ni corta 02 => C
          _ni(3); mk_intersection(); // corta 32
          _ni(1); mk_intersection(); // corta 12
        }
      }
      else if(dif3){
        _nf(3); mk_intersection(); // no corta 01 ni corta 02 pero corta 03 => D
        _ni(1); mk_intersection(); // corta 13
        _ni(2); mk_intersection(); // corta 23
      }

      //agrega y orienta el elemento
      if (qne==3) {
          if(n0_arriba) Swap(ne[1],ne[2]);
          //avisa a cada nodo que tiene el proximo elemento a agregar
          n[ne[0]].e+=e.len;n[ne[1]].e+=e.len;n[ne[2]].e+=e.len;
          e+=elemento(e_triangulo,ne);
          eplane+=i;
      }
      else if (qne==4) {
          if(n0_arriba) Swap(ne[1],ne[3]);
          n[ne[0]].e+=e.len;n[ne[1]].e+=e.len;n[ne[2]].e+=e.len;n[ne[3]].e+=e.len;
          e+=elemento(e_cuadrilatero,ne);
          eplane+=i;
      }
    }
  }//end for

  if(n) corte.tipo.set(m_sup);
}


// construye una malla como corte de otra por un plano (plane_p0,plane_normal)
// plane_p0 es un punto perteneciente al plano de corte
// plane_normal es el vector plane_normal al plano de corte
malla::malla(malla &m, const punto &p0, const punto &pn, pline *mapel, int step_i):_INI_MALLA {
  nombre[0]=ext[0]=0;
  _initime;

  copy_nombre(m); strcat(nombre,"_cut"); // nombre
  if (!m.mk_vecino()) {_savetime(isosup); return;}

  // asigno en v la distancia al plano
  if (mm!=&m) {delete[] v; v=new double[m.n.len];}
  for (int i=0;i<m.n.len;i++) v[i]=m.n[i]*pn;

  mc=this; mm=&m; _eplane=mapel; nviso=p0*pn; step=step_i;
  cut_tetras();

  _savetime(plano);
}

// Isosurface
malla::malla(malla &m, double valor, pline *mapel, int step_i):_INI_MALLA{
  nombre[0]=ext[0]=0;
  _initime;

  sprintf(nombre,"%s_isosup_%.2g",m.nombre, nviso); // nombre
  if (!m.mk_vecino()) {_savetime(isosup); return;}

  // asigno en v los nodevalues
  if (mm!=&m) {delete[] v; v=new double[m.n.len];}
  for (int i=0;i<m.n.len;i++) v[i]=m.n[i].v;

  mc=this; mm=&m; _eplane=mapel; nviso=valor; step=step_i;
  cut_tetras();
  _savetime(isosup);
}
