// PUNTOS O VECTORES

#include <cstring> // strcpy (i/o)
#include <iomanip> // manipuladores (i/o)
#include <cstdlib> // random

#include "punto.h"

using namespace std;

/*
static inline double _Max3
  (const double &x,const double &y,const double &z)
    {return ((x>y) ? ((x>z) ? x : z) : ((y>z) ? y : z));}
*/

static inline double _p2(const double &x){return (x*x);}
// #define _p2(x) ((x)*(x)) // de este modo opera dos veces

//compara y asigna componente a componente
void punto::set_min_max(punto& pmin, punto &pmax) const{
  if(x[0]<pmin[0]) pmin[0]=x[0];
  if(x[1]<pmin[1]) pmin[1]=x[1];
  if(x[2]<pmin[2]) pmin[2]=x[2];

  if(x[0]>pmax[0]) pmax[0]=x[0];
  if(x[1]>pmax[1]) pmax[1]=x[1];
  if(x[2]>pmax[2]) pmax[2]=x[2];
}
// comparacion
bool punto::eq_b(const double* b,const double &r) const
{
  double dx[3];
  return (
    ((dx[0]=fabs(x[0]-b[0]))<=r) &&
    ((dx[1]=fabs(x[1]-b[1]))<=r) &&
    ((dx[2]=fabs(x[2]-b[2]))<=r) &&
    (_p2(dx[0])+_p2(dx[1])+_p2(dx[2])<=_p2(r))
          );
}
bool punto::eq_b2(const double* b, const double &r2) const
{
  double dx,dy,dz;
  dx=x[0]-b[0];  dy=x[1]-b[1];  dz=x[2]-b[2];
  return (dx*dx+dy*dy+dz*dz<=r2);
}
bool punto::eq_b2(const punto &b, const double &r2) const
{
  double dx,dy,dz;
  dx=x[0]-b.x[0];  dy=x[1]-b.x[1];  dz=x[2]-b.x[2];
  return (dx*dx+dy*dy+dz*dz<=r2);
}
bool punto::eq_c(const double* b,const double &r) const
  {return (fabs(x[0]-b[0])<=r && fabs(x[1]-b[1])<=r && fabs(x[2]-b[2])<=r);}

// +=
punto& punto::operator+=(const punto& b)
  {x[0]+=b.x[0];x[1]+=b.x[1];x[2]+=b.x[2];return *this;}

// -=
punto& punto::operator-=(const punto& b)
  {x[0]-=b.x[0];x[1]-=b.x[1];x[2]-=b.x[2];return *this;}

// *=
punto& punto::operator*=(double f)
  {x[0]*=f;x[1]*=f;x[2]*=f;return *this;}

// /=
punto& punto::operator/=(double f)
  {x[0]/=f;x[1]/=f;x[2]/=f;return *this;}

// opuesto
punto punto::operator-() const
  {return punto(-x[0],-x[1],-x[2]);}

// suma
punto punto::operator+(const punto& b) const
  {return punto(x[0]+b.x[0],x[1]+b.x[1],x[2]+b.x[2]);}

// resta
punto punto::operator-(const punto& b) const
  {return punto(x[0]-b.x[0],x[1]-b.x[1],x[2]-b.x[2]);}

// producto escalar
double punto::operator*(const punto& b) const
  {return (x[0]*b.x[0]+x[1]*b.x[1]+x[2]*b.x[2]);}

// escalar * punto
punto punto::operator*(double f) const
  {return punto(x[0]*f,x[1]*f,x[2]*f);}

// cociente de un punto por un escalar (sin test)
punto punto::operator/(double f) const
  {return punto(x[0]/f,x[1]/f,x[2]/f);}

// producto vectorial, se usa el % para tener la misma precedencia que el *
punto punto::operator%(const punto& b) const
  {return punto(x[1]*b.x[2]-x[2]*b.x[1],
                 x[2]*b.x[0]-x[0]*b.x[2],
                 x[0]*b.x[1]-x[1]*b.x[0]);}
// triple producto
double punto::triple(const punto& p1,const punto& p2) const
   {return x[0]*(p1.x[1]*p2.x[2]-p1.x[2]*p2.x[1])+
           x[1]*(p1.x[2]*p2.x[0]-p1.x[0]*p2.x[2])+
           x[2]*(p1.x[0]*p2.x[1]-p1.x[1]*p2.x[0]);}

//modulo
double punto::mod() const
  {return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);}
//modulo al cuadrado
double punto::mod2() const
  {return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];}
//maxima componente en valor abs
double punto::modc() const
//  {return _Max3(fabs(x[0]),fabs(x[1]),fabs(x[2]));}
//  {return Max(Max(fabs(x[0]),fabs(x[1])),fabs(x[2]));}
{
  double M=fabs(x[0]),m=fabs(x[1]);
  if (M<m) M=m;
  m=fabs(x[2]);
  if (M<m) return m;
  return M;
}

// distancia
double punto::distancia(const punto& b) const
  {return (sqrt(_p2(x[0]-b.x[0])+_p2(x[1]-b.x[1])+_p2(x[2]-b.x[2])));}
// distancia al cuadrado
double punto::distancia2(const punto& b) const
  {return (_p2(x[0]-b.x[0])+_p2(x[1]-b.x[1])+_p2(x[2]-b.x[2]));}
// metrica del box
double punto::distanciac(const punto& b) const
//  {return _Max3(fabs(x[0]-b.x[0]),fabs(x[1]-b.x[1]),fabs(x[2]-b.x[2]));} compiler error en release??
//  {return Max(Max(fabs(x[0]-b.x[0]),fabs(x[1]-b.x[1])),fabs(x[2]-b.x[2]));}
{
  double M=fabs(x[0]-b.x[0]),m=fabs(x[1]-b.x[1]);
  if (M<m) M=m;
  m=fabs(x[2]-b.x[2]);
  if (M<m) return m;
  return M;
}

// calcula distancia solo si es <
bool punto::distancia_menor
(const punto& b, double &dmin) const{
  double d,
  di=fabs(x[0]-b.x[0]); if (di>dmin) return false; d =di*di;
  di=fabs(x[1]-b.x[1]); if (di>dmin) return false; d+=di*di;
  di=fabs(x[2]-b.x[2]); if (di>dmin) return false; d+=di*di;
  dmin=sqrt(d); return true;
}
bool punto::distancia2_menor
(const punto& b, double &dmin2) const{
  double
  d2= _p2(x[0]-b.x[0]); if (d2>dmin2) return false;
  d2+=_p2(x[1]-b.x[1]); if (d2>dmin2) return false;
  d2+=_p2(x[2]-b.x[2]); if (d2>dmin2) return false;
  dmin2=d2; return true;
}
bool punto::distanciac_menor
(const punto& b, double &dminc) const{
  double dc,
  delta=fabs(x[0]-b.x[0]); if (delta>dminc) return false;
  dc=delta;
  delta=fabs(x[1]-b.x[1]); if (delta>dminc) return false;
  if (delta>dc) dc=delta;
  delta=fabs(x[2]-b.x[2]); if (delta>dminc) return false;
  if (delta>dc) dc=delta;
  dminc=dc; return true;
}

// angulos, giros y 2d
punto& punto::giro(const punto &n,double ang){
  punto
    q=n*(operator*(n)/n.mod2()), // vector en la dir de n
    d=operator-(q); // diferencia
  double md=d.mod();
  if (md<ERRADM) return (*this); // no gira
  punto norm=(operator%(n)).dir(); // versor normal a n y este
  return operator=(q+d*cos(ang)+norm*(d.mod()*sin(ang)));
}

// interseccion en el plano
bool inters2d(
    const punto &p00,const punto &p01,
    const punto &p10,const punto &p11,
    punto &p){
  punto d0=p01-p00, d1=p11-p10, d00=p10-p00;
  double det=d1.pv2d(d0),x0/*,x1*/;
  if (fabs(det)<ERRADM) {p.zero(); return false;}
  x0=(d00[1]*d1[0]-d00[0]*d1[1])/det;
//  x1=(d00[1]*d0[0]-d00[0]*d0[1])/det;
  p=p00+x0*d0;
  return true;
}

// punto de minima distancia a dos rectas
punto& punto::prrmedio
  (const punto &p00,const punto &p01,
   const punto &p10,const punto &p11){
  punto
    d0=p01-p00,d1=p11-p10,d01=p10-p00,
    m[3]={d0,-d1,d0%d1},f=d01,a;
  transpose(m);
  a.solve(m,f);
  *this=(p00+d0*a[0]+p10+d1*a[1])/2;
  return *this;
}
punto prrmedio
  (const punto &p00,const punto &p01,
   const punto &p10,const punto &p11,
   punto &p02, punto &p12){
  punto
    d0=p01-p00,d1=p11-p10,d01=p10-p00,
    m[3]={d0,-d1,d0%d1},f=d01,a;
  transpose(m);
  a.solve(m,f);
  p02=p00+d0*a[0];
  p12=p10+d1*a[1];
  return (p02+p12)/2;
}

// distribucion esferica
// direccion aleatoria
// arquimedes: el area de la esfera es = a la del cilindro tangente (d(teta)*dz)
// cada slice en z tiene la misma cant de puntos (aunque parezca mentira!)
punto& punto::randdir()
{
  static const double factort=DOSPI/((unsigned int)RAND_MAX+1); // [0,2pi)
  static const double factorz=2.0/RAND_MAX; // [0,2]
  x[2]=rand()*factorz-1; // [-1,1]
  double 
    proy=sqrt(1-(x[2]*x[2])), // z proyectado en xy
    t=rand()*factort; // angulo equiprobable en xy
  x[0]=proy*cos(t);
  x[1]=proy*sin(t);
  return *this;
}

punto& punto::randdir2d()
{
  static const double factort=DOSPI/((unsigned int)RAND_MAX+1); // [0,2pi)
  double t=rand()*factort; // angulo equiprobable en xy
  x[0]=cos(t); x[1]=sin(t); x[2]=0;
  return *this;
}

// distribucion cubica
// punto aleatorio  [0,1)
punto& punto::randpt()
{
  static const double factor=double(1.0/((unsigned int)RAND_MAX+1));
  x[0]=rand()*factor;
  x[1]=rand()*factor;
  x[2]=rand()*factor;
  return *this;
}
punto& punto::randpt2d()
{
  static const double factor=double(1.0/((unsigned int)RAND_MAX+1));
  x[0]=rand()*factor;
  x[1]=rand()*factor;
  x[2]=0;
  return *this;
}

// solucion de m[ij]*x[j]=f[i]
// (se pueden eficientizar los calls)
punto& punto::solve(
 const punto m[3],
 const punto &f){
  punto adj[3]; double d;
  return solve(m,f,adj,d);
}
punto& punto::solve(
 const punto m[3],
 const punto &f,
 double &d){
  punto adj[3];
  // podria llamar al otro solve, pero asi es mas barato
  x[0]=x[1]=x[2]=0;
  d=m[0].operator*(adj[0]=m[1].operator%(m[2]));
//  if (fabs(d)<ERRADM) {d=0; return *this;}
  if (d==0) {
    return *this;
  }
  adj[1]=m[2].operator%(m[0]);
  adj[2]=m[0].operator%(m[1]);
  transpose(adj);
  x[0]=(adj[0].operator*(f))/d;
  x[1]=(adj[1].operator*(f))/d;
  x[2]=(adj[2].operator*(f))/d;
  return *this;
}
punto& punto::solve(
 const punto m[3],
 const punto &f,
 punto adj[3],
 double &d){
  x[0]=x[1]=x[2]=0;
  d=m[0].operator*(adj[0]=m[1].operator%(m[2]));
//  if (fabs(d)<ERRADM) {d=0; return *this;}
  if (d==0) {
    return *this;
  }
  adj[1]=m[2].operator%(m[0]);
  adj[2]=m[0].operator%(m[1]);
  transpose(adj);
  x[0]=(adj[0].operator*(f))/d;
  x[1]=(adj[1].operator*(f))/d;
  x[2]=(adj[2].operator*(f))/d;
  return *this;
}

void transpose(punto m[3]){
  double
  t=m[0].x[1]; m[0].x[1]=m[1].x[0]; m[1].x[0]=t;
  t=m[0].x[2]; m[0].x[2]=m[2].x[0]; m[2].x[0]=t;
  t=m[1].x[2]; m[1].x[2]=m[2].x[1]; m[2].x[1]=t;
}

// circunferencia y esfera
// (Pi-C)^2 = r^2 = C^2 = Pi^2-2Pi*C+C^2 => Pi*C=Pi^2/2

// circunferencia de tres puntos
void c3(
    const punto &p0,
    const punto &p1,
    const punto &p2,
    punto &c,
    double &r,
    punto &a){
  punto p[3]={p1-p0,p2-p0,p[0]%p[1]};
  double xv=p[2].mod();
  if (xv<ERRADM) {r=0; return;}
  p[2]/=xv; // el cuarto punto perpendicular modulo 1
  punto l(p[0].mod2()/2,p[1].mod2()/2,0.5);
  double a2;
  c.solve(p,l,a2);
  if (a2==0) {r=0; return;}
  a=p[2]*(a2/2); c-=p[2]*(c*p[2]);  r=c.mod(); c+=p0;
}
void c30(// el primero es el origen
    const punto &p1,
    const punto &p2,
    punto &c,
    double &r,
    punto &a){
  punto p[3]={p1,p2,p1%p2};
  double xv=p[2].mod();
  if (xv<ERRADM) {r=0; return;}
  p[2]/=xv; // el cuarto punto perpendicular modulo 1
  punto l(p[0].mod2()/2,p[1].mod2()/2,0.5);
  double a2;
  c.solve(p,l,a2);
  if (a2==0) {r=0; return;}
  a=p[2]*(a2/2); c-=p[2]*(c*p[2]);  r=c.mod();
}

// arco por tres puntos
void a3(
    const punto &p0,
    const punto &p1,
    const punto &p2,
    punto &c,
    double &r,
    punto &b){
  punto p[3]={p1-p0,p2-p0,(p[0]%p[1]).dir()};
  punto l(p[0].mod2()/2,p[1].mod2()/2,0.5);
  double a2;
  c.solve(p,l,a2);
  if (a2==0) {r=0; return;}
  c-=p[2]*(c*p[2]); r=c.mod(); c+=p0;
  double k=r/sqrt(l[1]/2);
  b=p[2]*(k-sqrt(k*k-1));
}

// esfera de cuatro puntos
void s4(
    const punto &p0,
    const punto &p1,
    const punto &p2,
    const punto &p3,
    punto &c,
    double &r,
    double &v){
  punto p[3]={p1-p0,p2-p0,p3-p0};
  punto l(p[0].mod2()/2,p[1].mod2()/2,p[2].mod2()/2);
  c.solve(p,l,v);
  if (v==0) {
    r=0; return;
  }
  v/=6; r=c.mod(); c+=p0;
}
void s40(// el primero es el origen
    const punto &p1,
    const punto &p2,
    const punto &p3,
    punto &c,
    double &r,
    double &v){
  punto p[3]={p1,p2,p3};
  punto l(p[0].mod2()/2,p[1].mod2()/2,p[2].mod2()/2);
  c.solve(p,l,v);
  if (v==0) {
    r=0; return;
  }
  v/=6; r=c.mod();
}

// centro de la circunsfera de n puntos
punto center(int n,const punto *p){
  if (n==2) {return (p[0]+p[1])/2;}
  if (n==3) {
    double r; punto a,c;
    c3(p[0],p[1],p[2],c,r,a);
    return c;
  }
  if (n==4) {
    double r,v; punto c;
    s4(p[0],p[1],p[2],p[3],c,r,v);
    return c;
  }
  return punto(0,0,0);
}

// diedro entre pderecho arista y pizquierdo
double diedro(
   const punto &pderecho,
   const punto &arista,
   const punto &pizquierdo,
   int modo){
  double m;
  punto xi,xd,a;
  m=arista.mod(); if (m<ERRADM) return 360; a=arista/m;
  xd=pderecho%arista;
  m=(xd).mod(); if (m<ERRADM) return 360; xd/=m;
  xi=arista%pizquierdo;
  m=(xi).mod(); if (m<ERRADM) return 360; xi/=m;
  if (modo<=1){// angulo
    double
    s=triple(xi,xd,a), // sin(180-d)
    c=-xi*xd,          // cos  (cos(180-d)=-cos(d))
    ang=atan2(s,c);      // -180 a 180
    if (ang<0) ang+=DOSPI;  // 0 - 360
    if (modo==0) return ang; // en radianes
    return r2g(ang); // en grados
  }
  else if (modo==2) return triple(xi,xd,a);  // sin
  else if (modo==3) return -(xi*xd);  // cos
  return 360;
}

/////////////////
// io

static bool hayz=true;
static int width=0, prec=16;
static char sep[10]="\t";

bool punto::io_z(int i)
  {if (i!=-1) hayz=(!i) ? false : true; return hayz;}
int punto::o_width(int i)
  {if (i!=-1) width=i; return width;}
int punto::o_precision(int i)
  {if (i!=-1) prec=i; return prec;}
const char* punto::o_separator(const char* s)
  {if (s) strcpy(sep,s); return sep;}

ostream& punto::graba(ostream &a) const{
  if (width) a << setw(width);
  a << setprecision(prec);
//  a << setiosflags(ios::scientific);
  a << x[0] << ((sep[0])? sep : "\t");
  if (width) a << setw(width);
  a << setprecision(prec);
//  a << setiosflags(ios::scientific);
  if (!hayz) return a << x[1];
  // hay z
  a << x[1] << ((sep[0])? sep : "\t");
  if (width) a << setw(width);
  a << setprecision(prec);
//  a << setiosflags(ios::scientific);
  a << x[2];
  return a;
}

istream& punto::lee(istream &a)
  {a >> x[0] >> x[1]; if (hayz) a >> x[2]; return a;}
