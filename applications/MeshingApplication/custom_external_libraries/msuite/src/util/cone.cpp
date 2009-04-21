#include <cmath> // fabs
#include "utiles.h"
#include "punto.h"
#include "cone.h"

using namespace std;

// ================================================================================================
// cone returns versor vc minimizing maximum angle [0,90) with the given versors
// cc is the maximum cosine (aperture)
// The resulting cone is invalid if it does not contains all the points (spans more than a hemisphere)
static punto mdir; // mean normal, must never be null
static punto vc; static double cc;

static void cone2(const punto& p0, const punto& p1){
	vc=(p0+p1); double m=vc.mod();
  if (m<ERRADM) {cc=0; return;}
	vc/=m; cc=p0*vc;
	if (cc<0){vc=-vc; cc=-cc;}
}

static void cone3(const punto &p0,const punto &p1,const punto &p2){  
	vc=(p1-p0)%(p2-p0); double m=vc.mod();
	if (m>ERRADM) {
		vc/=m; cc=p0*vc;
		if (cc<0){vc=-vc; cc=-cc;}
	}
	else {
		// coincident or "colinear" points => maximum distance (minimum cos)
		double c01=p0*p1,c12=p1*p2,c20=p2*p0;
		if      (c01<=c12&&c01<=c20) cone2(p0,p1);
		else if (c12<=c20&&c12<=c01) cone2(p1,p2);
		else                         cone2(p2,p0); 
	}
}

// Enveloping cone of l with f as defining boundary set
// entran dos conjuntos: l y f, calcula el minimo casquete con f en la frontera
// puntos iguales joden cuando ninguno entra en el cono del otro
// por ser unitarios, ERRADM es absoluto (asi que no hay problemas)
static bool cone(const punto* const l, int nl, int f[3], int nf) {
	if (nf && cc<=0) return false;
	if (nf==3){cone3(l[f[0]],l[f[1]],l[f[2]]); return (cc>0);}
	if (nf==2){cone2(l[f[0]],l[f[1]]);}
	int nll=0;
	if (nf==1) {nll=1; cone2(l[f[0]],l[0]);}
	if (nf==0) {nll=2; cone2(l[0],l[1]);}
	while (nl-nll) {
		if (vc*l[nll]<cc-ERRADM) { // si esta justito fuera, pasa de largo
			f[nf]=nll; 
			cone(l,nll,f,nf+1);
		}
		nll++;
	}
	return (cc>0);
}

bool cone(const punto* const l, int nl, punto &_vc, double &_cc){
  int f[3];
  bool retval=cone(l,nl,f,0);
  for(int i=0;i<nl;i++) if (vc*l[i]<cc-ERRADM) {_cc=0; return false;}
  _vc=vc; _cc=cc;
  return retval;
};
