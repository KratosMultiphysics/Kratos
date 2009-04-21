#include <cstdlib> // alloc
#include <cstring> // memset
//#define temporizar
#include "tiempo.h"
#include "voronoi.h"

using namespace std;

// genera una lista de h a partir del tamanio de elementos
bool voronoi::mk_h_nn(bool remake){
  if (!s) return false;
  if (m->hayh&&!remake) return true;
  array1<nodo> &n=m->n;
  _initime;
  int i,j,k,in,nlen=n.len;
  double d;
  ordlist nni(64);
  double hmin=MAXREAL;
  for (i=0;i<nlen;i++){
    const nodo &p=n[i]; const cpline &eni=p.e;
    nni.clean();
    for (j=0;j<eni.len;j++){
      const esfera &ej=s[eni[j]];
      for (k=0;k<NV;k++){
        in=ej[k]; if (in!=i) nni+=in;
      }
    }
    d=p.distancia2(n[nni[0]]);
    for (j=1;j<nni.len;j++) p.distancia2_menor(n[nni[j]],d);
    n[i].h=sqrt(d);
    set_min(hmin,n[i].h);
  }
  m->hmin=hmin; m->hayh=true;
  _savetime(mk_h_nn);
  return true;
}

// captura la malla como frontera y hace dir de elementos y h de los nodos
bool voronoi::mk_hydir(){
  qef=m->e.len;
  if (!qef) return false;

  //no acepta mezcla de tipos
  flagtype fdim=m->tipo&(m_vol|m_sup|m_lin|m_nodos);
  if (fdim!=m_sup&&fdim!=m_lin) {
    m->add_error(Mixed);
    return false;
  }

  if ((m->tipo).es(m_sup)) NV=4;
  else if ((m->tipo).es(m_lin)) NV=3;

  _initime;

  const array1<nodo> &n=m->n;
  qnf=n.len-m->nodosh;

  delete [] numf; numf=new int[qnf];
  // dir hacia afuera
  if (!m->orienta_cerrada(true,true,numf,true)) {
    delete [] numf; numf=0;
    _savetime(error); return false;
  }

  int i,j,nv;

  // h y h minimo
  double hmin=MAXREAL;
  for (i=0;i<qnf;i++) {
    if (n[i].h<ERRADM) n[i].h=m->h_nn_min(i,false);
    set_min(hmin,n[i].h);
  }

  // elementos de frontera de cada nodo
  if (efn) delete [] efn; efn=new cpline[qnf];
  for (i=0;i<qnf;i++) efn[i].roba(n[i].e);

  // centroide y h de los elementos de frontera
  delete [] gef; gef=new punto [qef];
  delete [] hef; hef=new double [qef];
  delete [] ef; ef=new elemento [qef];
  delete[] dir; dir=new punto[qef];
  for (i=0;i<qef;i++) {
    ef[i]=m->e[i]; dir[i]=m->dir[i]; nv=ef[i].nv();
    gef[i].zero(); hef[i]=0;
    for(j=0;j<nv;j++) {
      const nodo &nj=n[ef[i][j]];
      gef[i]+=nj; hef[i]+=nj.h;
    }
    gef[i]/=nv; hef[i]/=nv;
  }

  m->e_ini(); // aca inicializa los elementos de la malla

  // minh de nodos de h (para epsilon)
  for (;i<n.len;i++) set_min(hmin,n[i].h);
  m->hayh=true;
  m->hmin=hmin;

  _savetime(mk_hydir);
  return true;
}


//=====================================================================

// esferas externas

#define _ext(si,is){\
  s[is].f.set(e_exterior); status[is]=7;\
  for (k=0;k<NV;k++) {\
    vi=si.vecino[k]; \
    if (vi==-1) continue;\
    if (vi<0) vi=-vi-2;\
    s[vi].replace_vecino(is,-is-2);\
  }\
}


// testea con los elementos de frontera del nodo
#define _testdir\
  todosext=todosint=true;\
  for (ien=0;ien<eni.len;ien++){\
    cosn=dirc*dir[eni[ien]];\
    if (cosn<coslim) todosext=false; \
    if (cosn>=-coslim) todosint=false; \
    if (!todosext&&!todosint) break;\
  }\
  if (!todosext&&!todosint) continue;\
  if (todosext) _ext(si,is)\
  else status[is]=6;

#define _swap_last {indefinido[i]=indefinido.last(); i--; indefinido.len--;}

// busca tetraedros interiores/exteriores
// usa flag1 de nodos
void voronoi::frontera(){
  _initime;
  array1<nodo> &n=m->n;

  int i,in,isn,isn1,ino,inn,ien,is,is1,nf0,vi,nc,ie,nv;
  int iv,ixn,vde,vdi,k;
  pline nntested(100);
  esfera si;
  double cosn,coslim,hm,mod;
  bool todosext,todosint;
  punto pt,n0,g,dirc,npos[4];
  pline indefinido(s.len/5); // para no recorrer al pedo
  // 0: nada; 1-4 cant de nodos testeados; 5 indefinido; 6 interior; 7 exterior
  int *status=(int*)calloc(s.len,SZI);

  // elementos con fronteras distintas y malla de h
  for(is=0;is<s.len;is++){
    si=s[is];
    if (si.f.es(e_borrado)) continue;
    in=si[0];
    if (in<qnf) nf0=numf[in]; else {status[is]=6; continue;}
    for (iv=0;iv<NV;iv++) {
      in=si[iv]; if (in<qnf&&numf[in]==nf0) continue;
      // malla de h o fronteras distintas
      status[is]=6; break;
    }
  }
  delete [] numf; numf=0;

  // pone puntos a cada lado de cada cara de frontera a para
  // discriminar los tetraedros que los contengan
  for (ie=0;ie<qef;ie++){
    elemento efi=ef[ie]; nc=efi[0]; nv=efi.nv();
    dirc=dir[ie]*(hef[ie]/10);
    // puntos fuera
    pt=gef[ie]+dirc; is=de_que_tetra(pt,nc,k);
    if (is>=0){
//_revienta(status[is]==6);
      if (status[is]<6) {si=s[is]; _ext(si,is);}
    }
    for (iv=0;iv<nv;iv++){
      is=de_que_tetra((pt+n[efi[iv]])/2,nc,k);
      if (is>=0){
//_revienta(status[is]==6);
        if (status[is]<6) {si=s[is]; _ext(si,is);}
      }
    }
    // puntos dentro
    pt=gef[ie]-dirc; is=de_que_tetra(pt,nc,k);
_revienta(is<-1);
    if (is>=0) status[is]=6;
    for (iv=0;iv<nv;iv++){
      is=de_que_tetra((pt+n[efi[iv]])/2,nc,k);
_revienta(is<-1);
      if (is>=0) status[is]=6;
    }
  }
  delete [] gef; gef=0; delete [] hef; hef=0; delete [] ef; ef=0;

/*
for (is=0;is<s.len;is++) {
if (status[is]==6) _ext(s[is],is)
}
*/

  // aristas exteriores/interiores
  // con coslim generoso porque es muy peligroso pifiarla: una arista casi tangente
  // falsamente clasificada marca todos los elementos que la comparten
  coslim=.25;
  for (in=0;in<qnf;in++){
    nodo &ni=n[in]; n0=(pori) ? pori[in]: (punto&)ni;
    cpline sni=n[in].e; // copia la lista de esferas del nodo i;
    // deja las que no esten hechas
    for (isn=0;isn<sni.len;isn++){
      if (status[sni[isn]]<6) continue; // no esta hecha
      sni[isn]=sni.last();
      sni.len--; isn--;
    }
    // busca aristas interiores o exteriores definibles
    const cpline &eni=efn[in]; // elementos de frontera del nodo
    nntested.clean();
    for (isn=0;isn<sni.len;isn++){
      is=sni[isn]; si=s[is]; ixn=si.index(in); int &stat=status[is];
      for (iv=0;iv<NV-1;iv++){
        ino=si[(ixn+1+iv)%NV]; nodo &no=n[ino];
        if (no.f.es(flag1)) continue; // ya se testeo
        no.f.set(flag1); nntested+=ino;
        // testea con los elementos de frontera del nodo
        dirc=(((pori)? pori[ino]: (punto&)no) - n0).dir(); // saliendo del nodo
        _testdir;
        // definido, avisa a las esferas indefinidas que comparten la arista
        for (isn1=isn+1;isn1<sni.len;isn1++){
          is1=sni[isn1]; esfera &si1=s[is1];
          if (!si1.have(ino)) continue;
          if (todosext) _ext(si1,is1)
          else status[is1]=6;
          sni[isn1]=sni.last(); sni.len--; isn1--; // la saca de sni
        }
        break; // iv
      } // iv: vertices de la esfera
      if (stat>5) continue; // definda => otra esfera del nodo in
      stat++; //un nodo mas para el cual este elemento esta indefinido
      if (stat==NV){stat=5; indefinido+=is;} // indefinido desde todos los nodos
    } // isn: esferas del nodo
    for (inn=0;inn<nntested.len;inn++) n[nntested[inn]].f.reset(flag1);
  }// in: nodo

  // quedan en indefinido:
  //      elementos casi-tangentes
  //      tetraedros vertices o antivertices (todas las aristas en la frontera)
  //      gran variabilidad de normales en todos los nodos (nodos en aristas y vertices)

  // El asunto es no sacar el jamon de un sandwich de slivers.
  coslim=0; // consecuencias???
  for (i=0;i<indefinido.len;i++){
    is=indefinido[i]; si=s[is]; const int *v=si.vecino; int &stat=status[i];

    // direccion del centroide
    g.zero();
    for (iv=0;iv<NV;iv++) g+=npos[iv]=(pori)? pori[si[iv]]: (punto&)n[si[iv]];
    g/=NV;
    for (iv=0;iv<NV;iv++){
      in=si[iv];
      // testea con los elementos de frontera del nodo
      const cpline eni=efn[in]; // elementos de frontera del nodo
      // el centroide puede coincidir con un nodo (segmento con nodo al medio o cap perfecto)
      dirc=g-npos[iv]; mod=dirc.mod2();
      if (mod<ERRADM) _ext(si,is) // cap
      else {
        dirc/=sqrt(mod);
        _testdir;
      }
      // definido
      _swap_last; break; // no mas vertices iv
    }
    if (iv<NV) continue; // i - definido
    // centroide indefinido, no es vertice

    // vecinos
    // Si no es vertice, tres exteriores o interiores no generan dudas.
    // En cuanto a dos y dos, el tema es si matar concavidades o convexidades.
    // Sacar o dejar con dos y dos no es tan grave si no es "gordo" => solo slivers.
    for (vde=vdi=iv=0;iv<NV;iv++) {
      if (v[iv]<0) vde++; // vecino definido exterior
      else if (status[v[iv]]==6) vdi++; // vecino definido interior
    }
    if ((vde==1&&vdi==NV-1)|| // posible cap de frontera
      (NV==4&&vde==2&&vdi==2)){ // posible sliver de frontera
      hm=Min(n[si[0]].h,Min(n[si[1]].h,n[si[2]].h));if (NV==4) hm=Min(hm,n[si[3]].h);
      if (si.vt>pown(hm,NV-1)/100) continue; // 1/100
      _ext(si,is); // sliver o cap de frontera => exterior
      _swap_last; continue; // i - definido
    }
    else if (vde>NV-2){// tres o cuatro exteriores (en 3d)
      _ext(si,is);
      _swap_last; continue; // i - definido
    }
    else if (vdi>NV-2){// tres o cuatro interiores (en 3d)
      stat=6;
      _swap_last; continue; // i - definido
    }

    // centro de esfera
    for (iv=0;iv<NV;iv++){
      in=si[iv];
      // testea con los elementos de frontera del nodo
      const cpline eni=efn[in]; // elementos de frontera del nodo
      dirc=(si.c-npos[iv]).dir(); // saliendo del nodo
      _testdir;
      // definido
      _swap_last; break;  // iv - definido
    }

  }

  delete [] dir; dir=0; delete [] efn; efn=0;


  // quedan en indefinido:
  //      casi-tangentes
  //      todos los nodos con gran variabilidad de normales (aristas y vertices)
  // un primer paso por diferencia de 2 y despues define por mayorita
  int dif=1;
  while (1) {
    int lastlen=0;
    while(lastlen!=indefinido.len){
      lastlen=indefinido.len;
      for (i=0;i<indefinido.len;i++){
        is=indefinido[i]; si=s[is]; const int *v=si.vecino; int &stat=status[is];
        for (vde=vdi=iv=0;iv<NV;iv++) {
          if (v[iv]<0) vde++; // vecino definido exterior
          else if (status[v[iv]]==6) vdi++; // vecino definido interior
        }
        if (vde-vdi>dif){
          _ext(si,is);
          _swap_last; continue; // i - definido
        }
        else if (vdi-vde>dif){
          stat=6;
          _swap_last; continue; // i - definido
        }
      }
    }
    if (!dif) break;
    dif=0;
  }

  // en el ejemplo de la camara algunos casi-slivers conteniendo centros
  // impiden meter el centro si se los considera exteriores.
  // Dejo los slivers y se sacan al final
/*
  // los indefinidos que quedan son quasi-slivers de frontera => exteriores
  for(i=0;i<indefinido.len;i++) {
    is=indefinido[i]; si=s[is]; _ext(si,is);
  }
*/
// array1<esfera> st(indefinido.len);             ????????????
//  for(i=0;i<indefinido.len;i++) st+=s[indefinido[i]];
//  s.roba(st);

  free(status);
  _savetime(frontera);
}

// Esta implementado asi para minimizar el laburo de la ordlist de s.borrado
void voronoi::borra_exteriores(){
  _initime;
  array1<nodo> &n=m->n;
  int is,v;
  int k;

  // primero un squeeze para borrar los pocos que hay en s.borrado
  squeeze();
  // ahora agrega en serie
  for (is=0;is<s.len;is++){
    esfera &si=s[is];
    if (si.f.noes(e_exterior)) continue;
    const int *nsi=si.n;
    for (k=0;k<NV;k++) {
      v=si.vecino[k];
      n[nsi[k]].e.remove1(is);
      if (v==-1) continue;
      if (v<0) v=-v-2; else s[v].f.set(e_frontera);
      s[v].replace_vecino(-is-2,-1);
    }
    si.f.set(e_borrado); s.remove(is);
  }
  _savetime(borra_exteriores);
}

//===========================================================================
//          a que tetraedro pertenece un punto
//
//    esto esta aca porque no se donde meterlo y es para averiguar
//    en frontera() si un tetra es exterior o interior
//===========================================================================
// usa flag2 de esferas
// puede rebotar en slivers
int voronoi::de_que_tetra(const punto &p, int &nc, int &vlejano) {
  const array1<nodo> &n=m->n;
  int ie=n[nc].e[0]; if (ie<0) ie=-ie-2;
  pline en(100); // setea flag2;

  while (true){
    esfera &ei=s[ie]; nc=ei[0];
    if (ei.f.es(flag2)) {
//_revienta(1);
      ie=vlejano=-1; break;
    } // ya paso por este
    // testeo si es interior y lo agrego a la lista de pasados
    if (ei.have_t(p,&vlejano)) break;
    en+=ie; ei.f.set(flag2);
    ie=ei.vecino[vlejano];
    if (ie>=0) continue;
    // esfera exterior, si pasa a una convcavidad puede reentrar, pero si
    // pasa a un tetra virtual no (porque no esta la info y no porque no pueda pasar)
    if (ie==-1) {vlejano=-1; break;}
    ie=-ie-2; // esfera exterior no-virtual
  }

  for (int i=0;i<en.len;i++) s[en[i]].f.reset(flag2);

  return ie;
}
