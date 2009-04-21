// funciones de la malla para hacer delaunay y generar puntos
// (interfaz entre las clases malla y voronoi)

//#define temporizar
#include "tiempo.h"
#include "malla.h"
#include "voronoi.h"

using namespace std;

////////////////////////////////////////////////////
// funciones normales propias
////////////////////////////////////////////////////

// convierte esferas de voronoy en elementos de malla
// por chuncks y con operaciones riesgosas de mantenimiento
// porque en esta rutina es maximo el uso de memoria (esferas "y" elementos)
void voronoi::s2e(){
  _initime;
  m->frontera.ini(); m->tipo.reset(m_cerrada); m->nn.ini(); m->rm_dir();
  delete m->o; m->o=0; // porque los nodos pudieron ser movidos/quitados/agregados
  bool INFO_CL=m->INFO_CL;

  delete [] numf; numf=0;
  delete [] efn; efn=0;
  delete [] gef; gef=0;
  delete [] hef; hef=0;

  // pasa esferas a elementos de la malla
  array1<elemento> &e=m->e;
  array1<cpline>   &v=m->vecino;
  array1<punto>    &ce=m->ce;
  array1<double>   &re=m->re;
  array1<double>   &ve=m->ve;

  squeeze(); // no hay borradas
  s.fit();  // len=size (si no hay que borrar lo que sobra)

  const int chunk=100000,elen=s.len;
  int ie,esize=0;

  esfera **slist=s.list; ////////////ojo
  s.len=0;s.size=0;s.list=0;s.ini();

  // no se debe hacer e_ini
  // porque los elementos de cada nodo son esferas
  e.clean(); v.clean(); ce.clean(); re.clean(); ve.clean();

  e_tipo etipo((NV==3)? e_triangulo : e_tetraedro);

  if (INFO_CL) {cout << endl;}
  for (ie=0;ie<elen;ie++){
    if (!(ie%chunk)){
      // asigna memoria de a chunks
      esize=Min(esize+chunk,elen);
      e.resize(esize);  e.len=esize;
      v.resize(esize);  v.len=esize;
      ce.resize(esize); ce.len=esize;
      re.resize(esize); re.len=esize;
      ve.resize(esize); ve.len=esize;
      if (INFO_CL) {cout << "\resferas->elementos: " << ie; cout.flush();}
    }
    const esfera &si=*(slist[ie]);
    e[ie].set(etipo,si.n); e[ie].f=si.f;
    v[ie].copia(NV,si.vecino);
    ce[ie]=si.c; re[ie]=si.r; ve[ie]=si.vt;
    delete slist[ie];
  }
  free(slist);

  // en 2d los vecinos van como los nodos
  if (NV==3) for (ie=0;ie<elen;ie++) v[ie].first(2);

  // tipo
  if (m->e){
    m->tipo.reset(m_nodos|m_lin|m_sup|m_cerrada);
    m->tipo.set(m_modificada|m_orientada|((NV==3) ? m_sup : m_vol));
    if (slivers_marked) m->tipo.set(m_slivers_marked);
    else m->tipo.reset(m_slivers_marked);
  }

  if (INFO_CL) cout << "\resferas->elementos: " << elen << endl;
  _savetime(s2e);
}

// funciones de malla
bool malla::mk_puntos(bool slivers,bool h_radio){
  if (!e&&!h_radio) {add_error(No_Mesh);return false;}
  _initime;

  voronoi v(this);
  if (!v.mk_puntos(h_radio)) // genera puntos y elimina slivers de frontera
    {_savetime(error);return false;}
  if (slivers) v.rm_i_slivers();
  v.s2e();
//graba_dat();
//graba_calidad();

  _savetime(malla_puntos);
  _infotime("nodos: ");
  _infotime(n.len);
  _infotime("; elementos: ");
  _infotime(e.len);
  _infotime(endl);
  return true;
}

// voronoi de frontera (sin poner puntos interiores)
bool malla::f_delaunay(){
  if (!e) {add_error(No_Mesh);return false;}
  _initime;

  voronoi v(this);
  if (!v.f_delaunay()) // triangulacion
    {_savetime(error);return false;}
  v.s2e();
  _savetime(malla_f_delaunay);
  _infotime("nodos: ");
  _infotime(n.len);
  _infotime("; elementos: ");
  _infotime(e.len);
  _infotime(endl);
/*
// centros
  malla m;
  nodo ni; ni.f=1;
  elemento s("segmento"); s.f=1;
  int i,j,nv=v.NV;
  for (i=0;i<e.len;i++) {
    ni.setpos(ce[i]); m.n+=ni;
    s[0]=i;
    const cpline vi=vecino[i];
    for (j=0;j<nv;j++){
      if (vi[j]>i) {s[1]=vi[j]; m.e+=s;}
    }
  }
  m.copy_nombre(*this); strcat(m.nombre,"_centros"); m.graba_dat();
//
*/
  return true;
}

// triangulacion de una nube de puntos
bool malla::delaunay(bool slivers){
  _initime;
  voronoi v(this);
  if (!v.delaunay()) // triangulacion
    {_savetime(error);return false;}
  if (slivers)
    v.rm_slivers(); // ojo frontera si no paso por alpha

  v.s2e();

  _savetime(malla_delaunay);
  _infotime("nodos: ");
  _infotime(n.len);
  _infotime("; elementos: ");
  _infotime(e.len);
  _infotime(endl);
  return true;
}

bool malla::delaunay_refinado(malla &vieja, double alpha, bool slivers){
  _initime;

  voronoi v(this);

//  graba_dat("../OUTPUT/MSUITE/a_prueba.dat"); 
  if (!v.refina_esferas(vieja,alpha)) // refnada
    {_savetime(error);return false;}
  if (slivers)
    v.rm_i_slivers(); // ojo frontera si no paso por alpha
  v.s2e();

  _savetime(malla_delaunay_refinado);
  _infotime("nodos: ");
  _infotime(n.len);
  _infotime("; elementos: ");
  _infotime(e.len);
  _infotime(endl);
  return true;
}

bool malla::delaunay_refinado(double alpha, array1<double>& new_h, pline& nodes_from, array1<double>& shape_functions){
  _initime;

  voronoi v(this);

//  graba_dat("../OUTPUT/MSUITE/a_prueba.dat");
  if (!v.refina_esferas(alpha, new_h , nodes_from, shape_functions)) // refnada
    {_savetime(error);return false;}
  v.s2e();

  _savetime(malla_delaunay_refinado);
  _infotime("nodos: ");
  _infotime(n.len);
  _infotime("; elementos: ");
  _infotime(e.len);
  _infotime(endl);
  return true;
}

bool malla::interpola(malla &bk){
  // bbox de ambas
  bk.bbox(); bbox();
  bk.pmin.set_min_max(pmin,pmax);
  bk.pmax.set_min_max(pmin,pmax);
  if (!bk.e){
    voronoi v(&bk);
    v.delaunay();
    v.s2e();
//  bk.edt(.1,false); con poligonos no funca bien
  }
  hayv=true;

  int i,ie;
  int nv,j;
  double f[64],sumf;
  for (i=0;i<n.len;i++) {
    nodo &ni=n[i];
    ni.v=0;
    ie=bk.dequien(ni);
    if (ie<0) { // fuera de la malla
      _revienta(1); //en debug paramos a ver
      continue;
    }
    const elemento &ei=bk.e[ie];
    bk.fforma(ei,ni,f);
    nv=ei.nv();
    for (sumf=0,j=0;j<nv;j++){
      const nodo &ne=bk.n[ei[j]];
      if (ne.f.noes(n_virtual)) {ni.v+=f[j]*ne.v; sumf+=f[j];}
    }
    if (fabs(sumf-1)>ERRADM)
      ni.v/=sumf;
    set_min_max(nvmin,nvmax,ni.v);
  }
  bk.ini();
  return true;
}
