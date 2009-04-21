// rutinas de suavizado
//#define temporizar
#include "tiempo.h"
#include "malla.h"
#include "cone.h"
using namespace std;

// empareja el h real en superficies
// es una implementacion muy naive, para un patch muy suave

// Primero vuela nodos muy cercanos a otros 
//   (genera deformaciones porque mueve vertices de elementos)
// Despues agrega nodos en los huecos grandes
// Despues agrregla (es un decir) las deformaciones por swap de diagonales
// OJO: NO USAR CUANDO HAY MUCHA DIFERENCIA DE H
// SI SACA VARIOS PUNTOS CONECTADOS SIN SUAVIZADO, NECESARIAMENTE VA A HABER LIOS

bool malla::smooth_surface(){
  if (!e||e[0].dim()!=2) return false;
  _initime;
  
  int i,j,k,l,nc,ivuela,iqueda,ilast,ie,ie1,ie2,ixq,ixv,count=0;
  int in,in1,in2,ivk,iek,nci,ncj,jmax=-1,iv1,iv2,nv;
  double hm,vm,h21,h22,c,cmin,cmax;
  bool nohacerlo, inversion, unohecho=true, retval=false;
  elemento e3(e_triangulo);
  cpline vec3(3); vec3.len=3;
  
  mk_vecino();
  if (!hayh) mk_h_nn(true); else mk_nn();
  orienta(false,true); //para delaunay y edges (con dir)
  n_arista(); // marca para no mover los nodos de aristas
  mk_esferas(false,true); //para delaunay
  if (o) delete o; o=0;
  
  // captura los nodos de h aparte
  array1<nodo>nh(nodosh);
  for (i=0;i<nodosh;i++) nh[i]=n[n.len-nodosh+i];
  n.len-=nodosh; if (nn) nn.len-=nodosh;
  
  static const int n_intocable=n_permanente+n_frontera+n_edge;
  
  while (unohecho&&(++count)<5){
    unohecho=false;
    // vuela nodos cercanos
    for (i=0;i<n.len-1;i++) { // compara con uno mayor
      nodo &ni=n[i]; 
      h21=pown(n[i].h,2); if (h21<ERRADM) continue; // si h~0 refina ad-infinitum
      // nodos vecinos
      cpline &nni=nn[i];
      // el mas cercano de indice mayor
      for (cmin=MAXREAL,nc=-1,j=0;j<nni.len;j++) {
        if (nni[j]<i) continue;
        //c(comparador)=dist^2/((h1^2+h2^2)/2) el minimo es 1/4 (esta al cuadrado => 1/2h)
        h22=pown(n[nni[j]].h,2); if (h22<ERRADM) continue; // si h~0 refina ad-infinitum
        c=(8*ni.distancia2(n[nni[j]]))/(h21+h22);
        if (c<cmin) {cmin=c; nc=nni[j];}
      } if (cmin>1) continue; 
      // dist<h/2
      ivuela=nc; iqueda=i;
      if (n[ivuela].f.es_alguno(n_intocable)) {
        Swap(ivuela,iqueda);
        if (n[ivuela].f.es_alguno(n_intocable)) continue;
      }
      // busca elementos con ambos nodos, uno o dos, uno si son diagonal de cuadrado
      ie1=-1;ie2=-1;
      cpline &enq=n[iqueda].e; // elem del nodo que queda
      nohacerlo=false;
      for (j=0;j<enq.len;j++) {
        if (!e[enq[j]].have(ivuela)) continue;
        const elemento &ei=e[enq[j]]; // el elemento ei contiene los dos nodos
        // verifico si el tercer nodo no tiene solo tres elementos, 
        // en cuyo caso no puedo volar este elemento
        //    salvo que sean dos nodos sucesivos en un cuadrilatero
        nv=ei.nv();
        // si no es cuad con los dos conscutivos busco 
        // que el resto de nodos no tenga solo 3 elementos
        // AGREGADO: no permito que queden tres nodos, puede hacer un monio 
        //    (monio: dos triangulos vecinos con tres elementos en los dos nodos interfase)
        //    (tres elementos=> angulos>90 => un tri no puede tener dos nodos con tres elms)
        ///    REVISAR PARA QUADS!!
        ixq=ei.index(iqueda);
        if(
          nv!=4 // no es cuadrilatero
          ||
          !(ei[(ixq+1)%4]==ivuela||ei[(ixq+3)%4]==ivuela)//no son consecutivos (??????????????)
          )
        {// busca si el otro tiene (3) no, 4!
          for (k=1;k<nv;k++){
            l=(ixq+k)%nv;
            if (ei[l]==ivuela) continue;
            //          if (n[ei[l]].e.len==3) {nohacerlo=true;break;}
            if (n[ei[l]].e.len<5) {nohacerlo=true;break;}
          }
          if(nohacerlo) break;
        }
        if (ie1==-1) ie1=enq[j]; else {ie2=enq[j]; break;}
      }
      if (nohacerlo) continue;
      // elimina los elementos con ambos nodos (ie1 e ie2)
      ie=ie1; while(1){
        // reduce o elimina el elemento
        elemento &evuela=e[ie]; cpline &vie=vecino[ie];
        e_tipo t=evuela.tipo();
        if (t==e_triangulo) {
          //arregla los vecinos
          for(k=0;k<3;k++){//se posiciona antes de la arista que va a volar
            if (
              (evuela[k]==ivuela&&evuela[(k+1)%3]==iqueda) ||
              (evuela[k]==iqueda&&evuela[(k+1)%3]==ivuela)
              ) break;
          }
          //los otros dos vecinos son vecinos entre si
          iv1=vie[(k+1)%3]; iv2=vie[(k+2)%3];
          if (iv1>=0) vecino[iv1].replace1(ie,iv2);
          if (iv2>=0) vecino[iv2].replace1(ie,iv1);
          // vuela el elemento (swappeando con el ultimo)
          vuelae(ie,false); 
        }
        // es mas complicado que un simplice
        // los nn pueden ser diagonales de un cuadrado => un solo elemento con ambos
        else if (t==e_cuadrilatero){
          ixv=evuela.index(ivuela),ixq=evuela.index(iqueda);
          if (abs(ixq-ixv)==2) {//  diagonal => vuela
            //arregla los vecinos
            iv1=vie[ixv]; iv2=vie[(ixv+1)%4];
            if (iv1>=0) vecino[iv1].replace1(ie,iv2);
            if (iv2>=0) vecino[iv2].replace1(ie,iv1);
            iv1=vie[ixq]; iv2=vie[(ixq+1)%4];
            if (iv1>=0) vecino[iv1].replace1(ie,iv2);
            if (iv2>=0) vecino[iv2].replace1(ie,iv1);
            // vuela el elemento
            vuelae(ie,false);
          }
          else {
            // transforma en triangulo 
            // la arista ixq-ixv o viceversa vuela y su vecino
            if (((ixq+1)%4)==ixv) l=ixv; else l=ixq; // punto de partida
            //arma un triangulo y arregla los vecinos
            e3[0]=iqueda; vec3[0]=vie[l];
            e3[1]=evuela[(l+1)%4]; vec3[1]=vie[(l+1)%4];
            e3[2]=evuela[(l+2)%4]; vec3[2]=vie[(l+2)%4];
            evuela=e3; vecino[ie]=vec3; 
            n[ivuela].e.remove1(ie);
          }
        }
        else{
          _revienta(true); ////////////////////implementar poliedros
        }
        if(ie==ie1&&ie2!=-1) { // pasa ie a ie2 para volarlo
          // puede haber pasado que ie2 era el ultimo, entonces al volar ie1 swappeo
          // y ahora el ie2 es el que era ie1
          if (ie2!=e.len) ie=ie2; // por si justo era el ultimo (swappeado)
          else ie1=-1; // conserva ie pero hace esto para que no vuelva a entrar aca
        }
        else break;
      }

      // esto mueve nodos y puede hacer lios en las superficies
//      combine(ivuela,iqueda);// no usar combine
      if (n[iqueda].f.noes_ninguno(n_intocable)) {
        nodo &nvuela=n[ivuela];
        nodo &nqueda=n[iqueda];
        nqueda.f|=nvuela.f;
        nqueda.setpos((nqueda+nvuela)/2);
        nqueda.h=Min(nqueda.h,nvuela.h);
        nqueda.v=(nqueda.v+nvuela.v)/2;
      }
    
      retval=unohecho=true;
      // reemplaza en los elementos
      const cpline &env=n[ivuela].e;
      for (j=0;j<env.len;j++) e[env[j]].replace(ivuela,iqueda);
      enq+=env;
      // recalcula esferas
      inversion=false;
      for (j=0;j<enq.len;j++) {
        ie=enq[j];
        punto olddir=dir[ie];
        esfera_e(ie,ce[ie],re[ie],ve[ie],&dir[ie]);
        if (olddir*dir[ie]<0) inversion=true; // se invirtio un elemento
      }
      // recalcula nn
      cpline &nnq=nn[iqueda]; nn1(iqueda,&nnq,true);
      for (j=0;j<nnq.len;j++) nn1(nnq[j],&nn[nnq[j]],true);
/*      if (inversion&&n[iqueda].f.noes_ninguno(n_intocable)) {// laplaciano 
//        laplace_smooth(iqueda,1); // no llamar porque nn es trucho
        punto g(n[nnq[0]]); 
        for (int j=1;j<nnq.len;j++) g+=n[nnq[j]]; g/=nnq.len;
        n[iqueda].setpos(g);
      }*/
      // elimina el nodo swappeando con el ultimo
      ilast=n.len-1;
      if (ivuela!=ilast){
        const cpline &enl=n[ilast].e;
        for (j=0;j<enl.len;j++) e[enl[j]].replace(ilast,ivuela);
        const cpline nnl=nn[ilast];
        for (j=0;j<nnl.len;j++) nn[nnl[j]].replace1(ilast,ivuela);
        n.swap(ivuela,ilast); nn.swap(ivuela,ilast);
      }    
      n.len--; nn.len--; i--;
    }
//count=1000; continue;
    // agrega nodos intermedios
    cpline nodc,vecc,newe,newv(3); newv.len=3;
    for (i=0;i<e.len-1;i++){
      elemento &ei=e[i]; nci=ei.nc();
      cpline &vi=vecino[i]; 
      in1=in2=-1;
      for (j=0;j<nci;j++) {
        if (vi[j]<i) continue; // compara solo con vecino mayor
        in1=ei[j]; in2=ei[(j+1)%nci];
        h21=pown(n[in1].h,2); if (h21<ERRADM) continue;
        h22=pown(n[in2].h,2); if (h22<ERRADM) continue;
        if (n[in1].f.es_alguno(n_intocable)&&n[in2].f.es_alguno(n_intocable)) continue;
        //d(comparador)=dist^2/((h1^2+h2^2)/2) el maximo es 4 (esta al cuadrado => 2h)
        c=n[in1].distancia2(n[in2])/(2*(h21+h22));
        if (c>1) break;
      }
      if (j==nci) continue;
    
      // la arista n1 n2 es larga (>2h)
      nodo &n1=n[in1], &n2=n[in2];    
      nodo nm=(n1+n2)/2;// punto medio
      elemento &ej=e[vi[j]]; cpline &vj=vecino[vi[j]];
      // inicializa
      nodc.clean(); vecc.clean(); newe.clean();
      hm=vm=0;
      //lista nodos y vecinos y verifica que no sea muy cercano a alguno
      // elemento i
      in=in2; k=(j+1)%nci; do{
        nodo &nin=n[in];
        nodc+=in; hm+=nin.h; vm+=nin.v;
        if (nm.distancia2(nin)<pown(nin.h,2)/4) break;
        vecc+=vi[k];
        in=ei.npos(k); k=(k+1)%nci;
      }while(in!=in1);
      if (in!=in1) continue;
      // elemento opuesto
      j=vi[j]; ncj=ej.nc();
      k=ej.index(in1); in=in1; do{
        nodo &nin=n[in];
        nodc+=in; hm+=nin.h; vm+=nin.v;
        if (nm.distancia2(nin)<pown(nin.h,2)/4) break;
        vecc+=vj[k];
        in=ej.npos(k); k=(k+1)%ncj;
      }while(in!=in2);
      if (in!=in2) continue;
    
      // agrega el nodo y hace los nuevos elementos
      nc=nodc.len; // total
      nm.h=hm/=nc; nm.v=vm/=nc;
      for (k=0;k<nc;k++) {
        cpline &enk=n[nodc[k]].e;
        enk.remove1(i); enk.remove1(j);
      }
      e3[0]=n.len;
      // los primeros dos reemplazan a i y j
      e3[1]=nodc[0]; e3[2]=nodc[1]; newe+=i; ei=e3;
      nm.e+=i; n[nodc[0]].e+=i; n[nodc[1]].e+=i; vi.len=3;
      e3[1]=nodc[1]; e3[2]=nodc[2]; newe+=j; ej=e3;
      nm.e+=j; n[nodc[1]].e+=j; n[nodc[2]].e+=j; vj.len=3;
      for (k=2;k<nc;k++) {
        e3[1]=nodc[k]; e3[2]=nodc[k+1];
        newe+=iek=e+=e3;
        nm.e+=iek; n[nodc[k]].e+=iek; n[nodc[k+1]].e+=iek;
        vecino+=newv; re+=0; ce+=pzero; ve+=0; dir+=pzero;
      }
      n+=nm;
      // vecinos,esferas,dir
      for (k=0;k<nc;k++) {
        iek=newe[k];
        esfera_e(iek,ce[iek],re[iek],ve[iek],&dir[iek]);
        dir[iek]=((n[e[iek][1]]-n[e[iek][0]])%(n[e[iek][2]]-n[e[iek][0]])).dir();
        cpline &vk=vecino[iek]; vk[0]=newe[k-1]; vk[1]=ivk=vecc[k];
        if (ivk>=0) vecino[ivk].replaceall((k<nci-1)? i: j,iek);
        vk[2]=newe[k+1];
      }
      // nn
      nn+=nodc;
      for (k=0;k<nc;k++) nn1(nodc[k],&nn[nodc[k]],true);
      i--; retval=unohecho=true;
    }

//count=1000;continue;
    
    // "delaunay" (swap de diagonales)
    int lasti=-1,d_count=0; // para evitar loops infinitos
    for (i=0;i<e.len;i++){
      if (i!=lasti) {d_count=0; lasti=i;} else d_count++;
      elemento &ei=e[i]; 
      if (ei.tipo()!=e_triangulo) continue;
      cpline &vi=vecino[i];
      for (hm=0,j=0;j<3;j++) hm+=n[ei[j]].h; hm/=3;
      if (ve[i]<hm*hm/16){ // elemento aplastado (1/8 del bonito)
        //busca la arista mas larga;
        for (cmax=0,jmax=-1,j=0;j<3;j++) {
          c=n[ei[j]].distancia(n[ei[(j+1)%3]]);
          if (set_max(cmax,c)) jmax=j;
        }
        if (vi[jmax]<0) continue; // la mas larga es de frontera
        if (e[vi[jmax]].tipo()!=e_triangulo&&!q2t(vi[jmax])) continue;
        if (n[ei[jmax]].f.es_alguno(n_intocable)&& // ambos intocables
          n[ei[(jmax+1)%3]].f.es_alguno(n_intocable)) continue;
        if (diagonal_swap(i,vi[jmax],1)) {
          retval=unohecho=true; if (d_count<4) i--; continue;
        }
      }
      // testea delaunay
      for (j=0;j<3;j++) {
        if (vi[j]<i) continue;
        if (e[vi[j]].tipo()!=e_triangulo) continue;
        if (diagonal_swap(i,vi[j],1)) {
          retval=unohecho=true; if (d_count<4) i--; break;
        }
      }
    }
//count=1000;
  }

  // busca cuadrados con angulo muy obtuso para dividirlos en triangulos
//  q2t_all(150);

  // restaura los nodos de h
  if (nh){
    n+=nh;
    for (i=0;i<nh.len;i++) nn+=pline();
  }

  _savetime(smooth_surface)
  return retval;
}


/*     swap de diagonales en 2d  (muy parecido a q2t en meshelm.cpp)

                   3                            3
                  / \                          / \
            v3  /     \   v2                 / 2|1 \
              /    2    \                  /    |    \
            /             \              /      |      \
         0 ----------------- 2   ->   0 < 0  1  |  2   0> 2
            \             /              \      |      /
              \    1    /                  \    |    /
            v0  \     /   v1                 \ 1|2 /
                  \ /                          \ /
                   1                            1
   
  metodo 1=delaunay  2=maximo angulo  3=maxima minima altura
  metodo==0 es un swap incondicional
*/
bool malla::diagonal_swap(int ie1, int ie2, int metodo){
  elemento &e1=e[ie1],&e2=e[ie2];
  if (e1.tipo()!=e_triangulo||e2.tipo()!=e_triangulo) return false;
  mk_vecino();
  if (!metodo||tipo.es(m_planaxy)) orienta(false,false); // solo orientada
  else orienta(false,true); // necesito dir
  if (metodo==1&&!(ve&&re&&ce)) mk_esferas(); // necesito volumen centros y radios

  int i21,i12,in0,in1,in2,in3,v0,v1,v2,v3;

  cpline &vc1=vecino[ie1],&vc2=vecino[ie2];

  i21=vc1.index(ie2); if (i21==3) return false;
  i12=vc2.index(ie1); if (i12==3) return false;

  in0=e1[(i21+1)%3]; v0=vc1[(i21+1)%3];
  in1=e1[(i21+2)%3]; v1=vc1[(i21+2)%3];
  in2=e1[i21];       v2=vc2[(i12+1)%3];
  in3=e2[(i12+2)%3]; v3=vc2[(i12+2)%3];

  // en casos muy retorcidos puede hacer macana en nodos de valencia 3
  if ((v0>=0&&v0==v3)||(v1>=0&&v1==v2)) return false;

  if (metodo){
    // no se swapea si la interfase es una arista geometrica 
    // considero suficientemente coplanar entre 140º y 220º (180+/-40)
    // pero si los triangulos estan opuestos a menos de 5º se considera inversion (180+/-175)
    static const double coslim1=cos(g2r(40)),coslim2=cos(g2r(175));
    double cosa=1; // si es plana
    if (tipo.noes(m_planaxy)){
      cosa=dir[ie1]*dir[ie2];
      if (cosa<coslim1&&cosa>coslim2) return false; // arista y no-inversion
    }
    // es llano o casi llano
    if (cosa<coslim2) return false; // lo considero inversion
    // se puede swapear
    if (metodo==1){ // delaunay
      // verifica inclusion en el circulo
      if (fabs(ve[ie1])>ERRADM){
        if ((n[in3]-ce[ie1]).mod()>re[ie1]) return false;
      }
      else if (fabs(ve[ie2])>ERRADM){
        if ((n[in1]-ce[ie2]).mod()>re[ie2]) return false;
      }
      else return false; // ehhh, los dos son muy finitos, no vale!!
    }
    else { // no-delaunay
      // versores arista
      punto ar[4],ntmp=ar[3]=n[in0];
      ar[3]-=(ar[2]=n[in3]); // 3->0
      ar[2]-=(ar[1]=n[in2]); // 2->3
      ar[1]-=(ar[0]=n[in1]); // 1->2
      ar[0]-=ntmp;           // 0->1
      
      if (metodo==2){ // máximo angulo
        // uso angulo y no coseno por si hay invertidos
        int j,jmax=-1; double a,amax=0;
        for (j=0;j<4;j++) ar[j].dir();
        for (j=0;j<4;j++){
          // la segunda arista debe ir opuesta ==> tangente no cambia
          a=atan2((ar[j]%ar[(j+3)%4]).mod(),ar[j]*ar[(j+3)%4]); // -180 a 180
          if (a<0) a+=DOSPI; // <0 -> >180
          if (set_max(amax,a)) jmax=j;
        }
        if (jmax&1) return false; // si es 1 o 3 no swappea
      }
      else if (metodo==3){ // maximiza minima altura
        punto normal=(n[in2]-n[in0])%(n[in3]-n[in1]).dir();//para ver el signo del area
			  double
				  d0=ar[0].mod(),
				  d1=ar[1].mod(),
				  d2=ar[2].mod(),
				  d3=ar[3].mod(),
				  d02=(n[in2]-n[in0]).mod(),
				  d13=(n[in3]-n[in1]).mod(),
				  // areas
				  aa1=(ar[0]%ar[1])*normal,
				  ad2=(ar[1]%ar[2])*normal,
				  aa2=(ar[2]%ar[3])*normal,
				  ad1=(ar[3]%ar[0])*normal,
				  				  
				  mha1 = aa1/Max3(d0,d1,d02),
				  mha2 = aa2/Max3(d2,d3,d02),
				  mhd1 = ad1/Max3(d3,d0,d13),
				  mhd2 = ad2/Max3(d1,d2,d13);
						  
        if (Min(mha1,mha2)>Min(mhd1,mhd2)) return false;
      }
    }
  }

  n[in0].e.remove1(ie2); n[in2].e.remove1(ie1);
  n[in1].e+=ie2; n[in3].e+=ie1;

  if (v1>=0) vecino[v1].replace1(ie1,ie2);
  if (v3>=0) vecino[v3].replace1(ie2,ie1);

  e1[0]=in0; e1[1]=in1; e1[2]=in3; vc1[0]=v0; vc1[1]=ie2; vc1[2]=v3;
  e2[0]=in2; e2[1]=in3; e2[2]=in1; vc2[0]=v2; vc2[1]=ie1; vc2[2]=v1;

  // esferas y direccion
  if (ve||re||ce||(edir&&dir)){
    punto dir1,dir2,cen1,cen2; double rad1,rad2,vol1,vol2;
    esfera_e(ie1,cen1,rad1,vol1,&dir1);
    esfera_e(ie2,cen2,rad2,vol2,&dir2);
    if (ce) {ce[ie1]=cen1; ce[ie2]=cen2;}
    if (re) {re[ie1]=rad1; re[ie2]=rad2;}
    if (ve) {ve[ie1]=vol1; ve[ie2]=vol2;}
    if (edir&&dir) {dir[ie1]=dir1; dir[ie2]=dir2;}
  }

  if (nn) {
    nn[in0].remove1(in2); nn[in2].remove1(in0);
    nn[in1]+=in3; nn[in3]+=in1;
  }
  
  return true;
}

bool malla::diagonal_swap(int metodo) {
	bool retval=false;
  mk_vecino();
  if (!metodo||tipo.es(m_planaxy)) orienta(false,false); // solo orientada
  else {
    orienta(false,true); // necesito dir
    if (metodo==1&&!(ve&&re&&ce)) mk_esferas(); // necesito volumen centros y radios
  }
	int ie1,ie2,j;
  pline rehacer(e.len);
  for (j=0;j<e.len-1;j++) {rehacer+=j;e[j].f.set(flag1);} // el ultimo no
  array1<pline> pareja(e.len); // evita circulacion (orden e.len*max nn[i].len) (gracias a Pooyan)
  while (rehacer){
    ie1=rehacer[--rehacer.len];
    e[ie1].f.reset(flag1);
		if (e[ie1].tipo()!=e_triangulo) continue;
    pline &pasado=pareja[ie1];
		for (j=0;j<3;j++) {
			ie2=vecino[ie1][j];
      if (pasado.have(ie2)) continue; // ya se testeo esta pareja
      pasado+=ie2;
			if (ie1>ie2 || e[ie2].tipo()!=e_triangulo) continue; //(> una sola vez y evita vecino -1)
			if (!diagonal_swap(ie1,ie2,metodo)) continue;
      retval=true;
      if (e[ie2].f.noes(flag1)){rehacer+=ie2; e[ie2].f.set(flag1);}
      if (e[ie1].f.noes(flag1)){rehacer+=ie1; e[ie1].f.set(flag1);}
    }
  }
  return retval;
}


// suavizado laplaciano (mallas euclideas)
// factor=1/pasos
bool malla::laplace_smooth(int i, double factor){
  if ((!nn)||(!mk_nn())) return false;
  static const int n_intocable=n_frontera|n_edge|n_permanente|n_h;
  if (n[i].f.es_alguno(n_intocable)) return false;
  cpline &nni=nn[i]; punto g(0,0,0);
  for (int j=0;j<nni.len;j++) g+=n[nni[j]]; g/=nni.len;
  nodo &ni=n[i];
  if (g.eq_c(ni)) return false;
  if (factor==1) ni.setpos(g);
  else ni.setpos(ni*(1-factor)+g*factor);
  return true;
}

bool malla::laplace_smooth(double factor){
  if (!mk_nn()) return false;
  if (!mk_vecino()) return false;
  if (tipo.noes(m_vol)&&tipo.noes(m_planaxy)&&!n_arista()) return false;
  _initime;
  ce.ini(); re.ini(); ve.ini(); rm_dir();
  bool cambio_alguno=false;
  do{
    cambio_alguno=false;
    for (int i=0;i<n.len;i++){
      if (laplace_smooth(i,factor)) cambio_alguno=true;
    }
  }while (cambio_alguno);
    if (cambio_alguno) tipo.set(m_modificada);
  _savetime(laplace);
  return true;
}

//marca nodos de aristas geometricas
//setea n_edge
bool malla::n_arista(double grados){
  
   //solo superficies
  if (tipo.es(m_planaxy) ||
      ((tipo&(m_vol|m_sup|m_lin))!=m_sup)) return false;

  if (!mk_vecino()) return false;

  if (dir&&ndir) {dir.clean(); ndir=false;}
  if (!dir) {dir.resize(e.len); edir=true;}
  if (!orienta(false,true)) return false;

  const double coslim=cos(g2r(grados));
  int i,v,j,nv;

  //resetea n_edge
  for (i=0;i<n.len;i++) n[i].f.reset(n_edge);

  for (i=0;i<e.len-1;i++){
    punto diri=dir[i];
    const elemento &ei=e[i];
    const cpline &vi=vecino[i]; nv=vi.len;
    for (j=0;j<nv;j++){
      v=vi[j];
      if (v>=0&&(v<i||diri*dir[v]>coslim)) continue;
      // frontera o angulo chico
      n[ei[j]].f.set(n_edge); n[ei[(j+1)%nv]].f.set(n_edge);
    }
  }

  return true;
}
//graba un archivo de aristas geometricas
//setea n_edge y usa flag1 en nodos
bool malla::graba_aristas(double grados,const char *exti){
  
   //solo superficies
  if (tipo.es(m_planaxy) ||
      ((tipo&(m_vol|m_sup|m_lin))!=m_sup)) return false;

  if (!mk_vecino()) return false;

  if (dir&&ndir) {dir.clean(); ndir=false;}
//????  if (!dir) {dir.resize(e.len); edir=true;}
  if (!orienta(false,true)) return false;

  const double coslim=cos(g2r(grados));
  int i,v,j,nv;
  elemento s(e_segmento);
  
  //resetea n_edge
  for (i=0;i<n.len;i++) n[i].f.reset(n_edge);

  // genera los segmentos
  cascara a; a.parent=this;
  for (i=0;i<e.len-1;i++){
    punto diri=dir[i];
    const elemento &ei=e[i];
    const cpline &vi=vecino[i]; nv=vi.len;
    for (j=0;j<nv;j++){
      v=vi[j];
      if (v>=0&&(v<i||diri*dir[v]>coslim)) continue;
      // frontera o angulo chico
      s[0]=ei[j]; s[1]=ei[(j+1)%nv]; a.e+=s;
      if (n[s[0]].f.noes(flag1)){n[s[0]].f.set(flag1|n_edge);a.n+=s[0];}
      if (n[s[1]].f.noes(flag1)){n[s[1]].f.set(flag1|n_edge);a.n+=s[1];}
    }
  }
  for (i=0;i<a.n.len;i++) n[a.n[i]].f.reset(flag1);

  // graba
  malla m(a);
  // nombre y extension
  sprintf(m.nombre,"%s_edges",nombre);
  if (exti&&exti[0]) sprintf(m.nombre+strlen(m.nombre),".%s",exti);
  if(!(m.graba())) return false;
  return true;
}


// Elimina aristas interiores muy chichas y todos sus tetraedros 
// Esto fue armado para recomponer el desastre metrico que hace la capa
//   limite al meterse aplastando tetraedros
// Desaparecen algunos nodos
// Solo se aplastan simplices
bool malla::colapsa_aristas_chicas(){
  e_tipo etipo;
  if (tipo.es(m_vol)) etipo=e_tetraedro;
  else if (tipo.es(m_sup)) etipo=e_triangulo;
  else return false;
  
  // necesidades
  if (!mk_vecino()) return false;
  if (!mk_nn()) return false;
  if (!mk_h_nn_med()) return false; // con remake forzado

  rm_esferas();

  int in0,in1,i,j,k,iej,v,v0,v1,ilast,nv=(etipo==e_triangulo)? 3 : 4,nvl;

  for(in0=0;in0<n.len-1;in0++){
    nodo &n0=n[in0]; if (n0.f.es(n_frontera)) continue;
    pline &nn0=nn[in0];
    for(i=0;i<nn0.len;i++){
      in1=nn0[i]; if (in1<in0) continue;
      nodo &n1=n[in1]; if (n1.f.es(n_frontera)) continue;
      if (n0.distancia2(n1)>pown((n0.h+n1.h)/10,2)) continue; //h/5
//      if (n0.distancia2(n1)>pown((n0.h+n1.h)/4,2)) continue; //h/2
      ordlist elms=n0.e.inters(n1.e); // elementos comunes
      // verifica que sean tetras
      for (j=0;j<elms.len&&e[elms[j]].tipo()==etipo;j++); // no hace nada
      if (j<elms.len) 
        continue; // hay algun no simplice
      // verifica que no haya solo tres en ninguna arista exterior
      for (j=0;j<elms.len;j++) {
        iej=elms[j]; elemento &ej=e[iej]; cpline &vj=vecino[iej];
        v0=vj[ej.index(in0)]; v1=vj[ej.index(in1)]; 
        if (v0>=0&&v1>=0&&vecino[v0].have(v1)) break;
      }
      if (j<elms.len) 
        continue; // hay alguna arista con tres
      // vuela todos los elementos comunes
      // swapea con el ultimo ==> de mayor a menor
      for (j=elms.len-1;j>=0;j--) {
        iej=elms[j]; elemento &ej=e[iej]; cpline &vj=vecino[iej];
        v0=vj[ej.index(in0)]; v1=vj[ej.index(in1)]; 
        if (v0>=0) vecino[v0].replace1(iej,v1);
        if (v1>=0) vecino[v1].replace1(iej,v0);
        for (k=0;k<nv;k++) n[ej[k]].e.remove1(iej);
        // swap
        ilast=e.len-1;
        if (iej!=ilast){// si es el ultimo, simplemente lo elimina
          // no es el ultimo  
          // swap con el ultimo
          const elemento &el=e[ilast]; nvl=el.nv();
          for (k=0;k<nvl;k++) n[el[k]].e.replace1(ilast,iej);
          cpline &vl=vecino[ilast];
          for (k=0;k<vl.len;k++) {v=vl[k]; if(v>=0) vecino[v].replace1(ilast,iej);}
          vecino.swap(iej,ilast);
          if (ve) ve.swap(iej,ilast);
          if (edir&&dir){
            if (dir.len==e.len) dir.swap(iej,ilast);
            else rm_dir();
          }
          e.swap(iej,ilast);
        }
        if (edir&&dir) dir.len--;
        if (ve) ve.len--;
        vecino.len--;
        e.len--;
      }
      // vuela n1
      // combina datos
      n0.f|=n1.f;
      n0.setpos((n0+n1)/2);
      n0.h=(n0.h+n1.h)/2; // h medio
      n0.v=(n0.v+n1.v)/2;
      // elementos
      pline &en1=n1.e;
      for (j=0;j<en1.len;j++) e[en1[j]].replace(in1,in0);
      n0.e+=n1.e; // no hay comunes
      // recalcula nn
      nn1(in0,&nn0,true); for (j=0;j<nn0.len;j++) nn1(nn0[j],&nn[nn0[j]],true);
      // swappea con el ultimo
      ilast=n.len-1;
      if (in1!=ilast){
        const cpline &enl=n[ilast].e;
        for (j=0;j<enl.len;j++) e[enl[j]].replace(ilast,in1);
        const cpline nnl=nn[ilast];
        for (j=0;j<nnl.len;j++) nn[nnl[j]].replace1(ilast,in1);
        n.swap(in1,ilast); nn.swap(in1,ilast);
      }    
      if (dir&&ndir) {if (dir.len==n.len) dir.len--; else rm_dir();}
      n.len--; nn.len--;
      i--;
    }
  }
  return true;
}

void malla::j_min_med(int in,double &jmin, double &jmed) const{
  const pline &en=n[in].e;
  double ji=ve[en[0]]; jmin=jmed=ji;
  for (int i=1;i<en.len;i++) {
    ji=ve[en[i]]; set_min(jmin,ji); jmed+=ji;
  } jmed/=en.len;
}

#define _agregar(in)\
if (n[in].f.noes_ninguno(noincluir)&&\
    (counter[in]<150||(ji[2*in]<0&&counter[in]<1000))){\
  n[in].f.set(flag1); rehacer+=in; counter[in]++;\
}

#define _recalc_v {\
  const pline &en=ni.e;\
  for (i=0;i<en.len;i++) {\
    ie=en[i]; ve[ie]=volumen(ie); const int* ne=e[ie].n; nv=e[ie].nv(); hm=0;\
    for (k=0;k<nv;k++) hm+=n[ne[k]].h;\
    ve[ie]*=pown(nv/hm,3);\
  }\
}

// empareja el jacobiano
// ps=p0 y n(unitario) del plano de simetria
// si hay plano de simetria los inmoviles deben venir con n_permanente
/*
static int outnode=-1;
bool malla::suavej(double clim_i,const punto ps[2]){
  if (!mk_vecino()) return false;  //n_frontera
  if (!mk_h_nn_min(true)) return false; // movimientos en relacion a h min

  _initime;
  if (INFO_CL) cout << "\nSuavizado" << endl;

  ofstream* mov=0;
  if (outnode>-1) {
    char fn[255]; sprintf(fn,"%s_mov_nod_%d.xyz",nombre,outnode);
    mov=new ofstream(fn);
    (*mov) << n[outnode] << endl;
  }

  double jc,newjc,hi,hm,
    *ji=(double*)malloc(2*n.len*SZD);
  int ie,nv,in,i,k,ink,maxcounter=0,maxcnode,
    *counter=(int*)calloc(n.len,SZI);

  pline rehacer(n.len/2);
  punto jdir,lastpos;

  int noincluir=n_permanente|flag1; if (!ps) noincluir|=n_frontera;

  if (!ve) volumen(); // para los tetraedros
  for (ie=0;ie<e.len;ie++){// volumen de cada elemento
    const int* ne=e[ie].n; nv=e[ie].nv(); hm=0;
    for (k=0;k<nv;k++) hm+=n[ne[k]].h;
    ve[ie]*=pown(nv/hm,3);
  }

  // si hay nodos offset se guardan los nodos para no alejarlos
  array1<punto> nori; int oldnlen=n.len; // son los ultimos
  if (n.last().f.es(n_offset)){
    for (in=n.len-1;in>=0;in--) 
      if (n[in].f.noes(n_offset)||n[in].f.es(n_permanente)) break;
    oldnlen=in+1;
    nori.resize(n.len-oldnlen);
    for (in=oldnlen;in<n.len;in++) nori+=n[in];
  }

  // se guardan jmin y jmed de cada nodo  
  for (in=0;in<n.len;in++) {
    j_min_med(in,ji[2*in],ji[2*in+1]);
  }

  double clim=Min(clim_i,1e-3);
  while (1){
    // rejunta los que hay que tocar
    for (in=0;in<n.len;in++) {
      if (ji[2*in]/ji[2*in+1]<clim) _agregar(in);
    }

    pline hacer(rehacer.len); 
    while (rehacer){
      if (INFO_CL) cout << "\rrehacer: " << rehacer.len << "\t\t." << flush;
      rehacer.swap(hacer);
      while (hacer){
        in=hacer.last(); hacer.len--;
        nodo &ni=n[in]; ni.f.reset(flag1); hi=ni.h/500; lastpos=ni;
        const pline &nni=nn[in];
        if (set_max(maxcounter,counter[in])) maxcnode=in;

        jc=ji[2*in]/ji[2*in+1];
                   ni[0]+=hi; _recalc_v // x
        j_min_med(in,ji[2*in],ji[2*in+1]); jdir[0]=ji[2*in]/ji[2*in+1]-jc;
        ni[0]-=hi; ni[1]+=hi; _recalc_v // y
        j_min_med(in,ji[2*in],ji[2*in+1]); jdir[1]=ji[2*in]/ji[2*in+1]-jc;
        ni[1]-=hi; ni[2]+=hi; _recalc_v // z
        j_min_med(in,ji[2*in],ji[2*in+1]); jdir[2]=ji[2*in]/ji[2*in+1]-jc;
        ni[2]-=hi;

        if (jdir.mod()<ERRADM) {// no cambia
          _recalc_v
          j_min_med(in,ji[2*in],ji[2*in+1]);// recalcula
          if (ji[2*in]/ji[2*in+1]<clim) _agregar(in); // puede haberlo arreglado otro
          continue; 
        }

        jdir*=(ni.h/jdir.mod()/100);
        ni+=jdir; 
        if (ps&&ni.f.es(n_simetria)) ni-=ps[1]*((ni-ps[0])*ps[1]);
        if (ni.f.es(n_offset)){ // si es permanente no entro aca
          // no mas alla de h/3
          punto &pori=nori[in-oldnlen]; hi=(ni-pori).mod();
          if (hi>ni.h/3) ni.setpos(pori+(ni-pori)*(ni.h/3/hi));
        }
        _recalc_v
        j_min_med(imin,ji[2*imin],ji[2*imin+1]); newjc=ji[2*imin]/ji[2*imin+1];
        if (newjc<jc) {// empeoro
          ni.setpos(lastpos); newjc=jc; _recalc_v
          // recalcula
          j_min_med(in,ji[2*in],ji[2*in+1]);
          continue;
        }
        // queda modificado
        if (newjc<clim) _agregar(in);
        for (k=0;k<nni.len;k++){
          ink=nni[k];j_min_med(ink,ji[2*ink],ji[2*ink+1]);
          if (ji[2*ink]/ji[2*ink+1]<clim) _agregar(ink);
        }
        if (in==outnode) (*mov) << ni << endl;
      }
    }
    if (INFO_CL) {
      cout << "\rSuavizado Completo :" << clim << endl;
      cout << "\rMaximo contador: " << maxcounter << " nodo: " << maxcnode << endl;
    }
    if (clim==clim_i) break;
    clim=Min(clim*=10,clim_i);
    memset(counter,0,n.len*SZI);
  }
  free(counter); free(ji);
  
  if (outnode>-1) mov->close();

  ve.ini();

  _savetime(suavej);
  return true;
}
*/
