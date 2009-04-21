// Poliedrizacion

#include<cstdio>
//#define temporizar
#include "tiempo.h"
#include "malla.h"

using namespace std;

// todos los poliedros a tetraedros
void malla::p2t_all(){
  for (int i=0;i<e.len;i++) {
    e_tipo t=e[i].tipo();
    if (t==e_poliedro||t==e_poligono) {p2t(i); i--;}
  }
}


// conviene? hacer primero la union y despues la captura
// en caso contrario propaga por captura y genera mas concavidades
// aparte para la union se compara solo una vez con indice mayor
// la captura la decide el sliver, para dejar la menor oncavidad

// siempre volar el de indice mas grande juntae(chico,grande)
// porque swapea con el ultimo

//    en 3d ei.nc()==2*ei.nv()-4 si el poliedro tiene topologia esferica
bool malla::edt(double delta) {
  _initime;

  if (!e||!mk_vecino()||!mk_esferas()) {_savetime(error);return false;}

  int i,nexti,iej,io,count;
  int j,k,nc=-1,ncj,dim=e[0].dim();
  double dc,ri,rj,dmax;
  punto ci;
  bool added;

  delta/=2; // dividido para promediar

  // primera pasada para esferas cercanas
  count=0; nexti=0; while (nexti<e.len){ //  cada elemento
    i=nexti;nexti=i+1;
    if (INFO_CL&&(i>>10)<<10==i)
      {cout << "\rEDT (union): " << i << "\t" << count; cout.flush();}
    do{
      added=false;
      ri=re[i]; ci=ce[i]; cpline &vi=vecino[i]; nc=vi.len;
      for (j=0;j<nc;j++) { // cada cara
        iej=vi[j];
        if (!added&&iej<i) continue; // anterior o frontera
        else if (iej<0) continue; // frontera
        dc=ci.distancia(ce[iej]); rj=re[iej];
        // dos esferas se unen si dc<delta*(ri+rj)/2
        if (dc>delta*(ri+rj)) continue;
        // esferas cercanas
        if (iej<i) {nexti=i; Swap(iej,i);} // debe volar el de mayor indice
        juntae(i,iej); // lo agrega
        added=true; count++; break;
      }
    }while(added); //si junta vuelve a revisar
  }
  if (INFO_CL) cout << "\rEDT (union): " << e.len << "\t" << count << endl;

  // capturas y desarme de concavidades
  if (dim==3) {
    tipo.reset(m_slivers_marked);
    // captura
    nexti=0; count=0; while (nexti<e.len){ //  cada elemento
      i=nexti;nexti=i+1;
      e[i].f.reset(e_sliver);
      if (INFO_CL&&(i>>10)<<10==i)
        {cout << "\rEDT (captura): " << i << "\t" << count; cout.flush();}
      do{
        cpline &vi=vecino[i]; nc=vi.len;
        if (nc==4) break; // no es poliedro
        // poliedro
        added=false;

        // captura totalmente contenidos
        for (j=0;j<nc;j++) { // cada cara
          if (
               (iej=vi[j])<0 // frontera
            || !e[i].have(e[iej],true) // no esta contenido
            ) continue;
          // contenido
          cpline &vj=vecino[iej]; ncj=vj.len;
          // verifica si este es el mas cercano que lo contiene
          for (k=0;k<ncj-1;k++){
            if (
                 (io=vj[k])!=i  // no es este
              && io>=0          // no es frontera
              && vj.index(io,k+1)<ncj // es vecino por dos caras
              && ce[iej].distancia(ce[io])<=ce[iej].distancia(ce[i]) // es mas cercano
            ) break;
          }
          if (k==ncj) continue;
          ri=re[i]; ci=ce[i];
          if (iej<i) {nexti=i; Swap(iej,i);} // debe volar el de mayor indice
          juntae(i,iej); // lo agrega 
          re[i]=ri; ce[i]=ci; // sin modificar centro y radio
          added=true; count++; break; // retestea solo porque cambio la estructura de caras
        }
      }while(added); //si junta vuelve a revisar
    }
    if (INFO_CL) cout << "\rEDT (captura): " << e.len << "\t" << count << endl;
    // desarme de poliedros muy cocavos
    for (count=0,i=0;i<e.len;i++){ //  cada elemento
      if (INFO_CL&&(i>>10)<<10==i)
        {cout << "\rEDT (concavos): " << i << "\t" << count; cout.flush();}
      if (e[i].nc()>4&&maxd(i,dmax,j,k)>185)  // 5º de exceso
        {p2t(i);count++;}
    }
    if (INFO_CL) cout << "\rEDT (concavos): " << e.len << "\t" << count << endl;
  }

  _savetime(edt);
  _infotime("nodos: ");
  _infotime("; elementos: ");
  _infotime(e.len);
  _infotime(endl);
  return true;
}

