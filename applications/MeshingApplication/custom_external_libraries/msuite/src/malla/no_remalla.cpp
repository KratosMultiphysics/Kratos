#include "malla.h"

// En la primera llamada se genera una malla de conectividades fijas
// con los nodos que tienen determinados bits encendidos en el flag.
// flag_fija es el conjunto de bits de los nodos que forman la malla fija
// pueden ser varios si hay varias piezas.
// Un elemento es fijo si tiene los tres nodos con el bit de la pieza en cuestion,
// se averigua enmascarando f1&f2&f3&flag_fija

// En las subsiguientes recupera las conectividades impuestas
// los numeros de nodo de mfija deben coincidir con los actuales
// pero las posiciones no importan

// IMPLEMENTADO SOLO PARA TRIANGULOS PLANOS

// conserva la malla fija y su lista de nodos frontera
static malla* mfija=0; // malla fija
static pline* nf_fijo=0; // nodos de frontera de la malla fija
static int efijos=0;

static int mk_mfija(malla &m, flagtype flag_fija){
  int i, ultimo_fijo=-1;
  // genera una cascara con los elementos que tienen 3 nodos fijos
  cascara cfija; cfija.parent=&m; 

  // busca los nodos fijos y los pone al ppio
  for (i=0;i<m.n.len;i++) {
    if (m.n[i].f.noes_ninguno(flag_fija)) continue;
    if (i>(++ultimo_fijo)) m.swap_n(ultimo_fijo,i);
    cfija.n+=ultimo_fijo;
  }
  // elementos
  pline ix(int(1.2*m.e.len)); cfija.e.resize(int(1.2*m.e.len)); // preasigna espacio
  m.mask_elms(&cfija.e,&ix,flag_fija); 
  // malla fija
  if (mfija) delete mfija; mfija=new malla(cfija); efijos=mfija->e.len;
  // frontera de la fija o intefase
  mfija->mk_vecino();
  mfija->graba_dat();
  if (nf_fijo) delete nf_fijo; nf_fijo=new pline; pline &nff=*nf_fijo;
  for (i=0;i<mfija->n.len;i++) {if (mfija->n[i].f.es(n_frontera)) nff+=i;}
  //pone los elementos fijos al ppio
  for (i=0;i<efijos;i++) m.swap_e(i,ix[i]);
  return efijos;
}

static int recover_elms(malla &m, flagtype flag_fija){
  m.mk_vecino();

  int i,j,k,in,ink,mask,qei,mix,ie,ieo,nflen=nf_fijo->len;

  // recupera la frontera
  // por cada elemento de frontera fija
  // busca elementos que tenga nodos distintos
  // interior de fija y no fija o interior de fija y otra fija
  for (i=0;i<nflen;i++){
    in=(*nf_fijo)[i]; 
    mask=m.n[in].f&flag_fija; // bit del nodo
    const pline &eni=m.n[in].e; qei=eni.len;
    for (j=0;j<qei;j++){
      ie=eni[j]; elemento &ei=m.e[ie];
      for (mix=0,k=0;k<3;k++){
        ink=ei[k];
        if (m.n[ink].f.noes_ninguno(mask)) mix|=1; // no pertenece a esta pieza
        else // pertenece 
          if (mfija->n[ink].f.noes(n_frontera)) mix|=2; //pero no es frontera de fija
      }
      if (mix!=3) continue;
      // hay que corregir
      ieo=m.vecino[ie][ei.index(in)+1];
      if (!m.diagonal_swap(ie,ieo,false))
        return 0;
      i--;break;//testea de nuevo el mismo nodo
    }
  }

  // recupera interior de la malla solida
  pline ix(efijos); m.mask_elms(0,&ix,flag_fija);
  // inicializa listas que se relacionan con elementos
  m.nn.ini(); m.vecino.ini(); m.ve.ini(); m.re.ini(); m.ce.ini(); m.rm_dir();  
  // en la malla actual hay ix.len elms a reemplazar por mfija.e.len
  // reemplaza los elementos por un array nuevo
  array1<elemento>e(mfija->e); e.resize(efijos+m.e.len-ix.len);
  for (i=m.e.len-1;i>=0;i--){ // ix esta en orden creciente
    if (ix.len&&i==ix.last()) ix.len--;
    else e+=m.e[i];
  }
  m.e.roba(e);
  // nueva lista de elementos de cada nodo
  for (i=0;i<m.n.len;i++) m.n[i].e.clean();
  for (i=0;i<m.e.len;i++) for (j=0;j<3;j++) m.n[m.e[i][j]].e+=i;

  return efijos;
}

int malla::mk_mfija(flagtype flag_fija, bool primero){
  delaunay();
  if (primero) return ::mk_mfija(*this,flag_fija);
  else return ::recover_elms(*this,flag_fija);
}

