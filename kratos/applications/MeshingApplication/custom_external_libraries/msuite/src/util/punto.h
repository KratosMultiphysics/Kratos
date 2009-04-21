//              PUNTOS O VECTORES (prealocado)

#ifndef PUNTO_H
#define PUNTO_H

#include <iostream> // overload de <<
#include <cmath> // fabs, atan...
#include "utiles.h" //erradm DOSPI

class punto{
protected:
  double x[3];

  std::istream& lee(std::istream &);
  std::ostream& graba(std::ostream &) const;

public:
  // constructores
  punto() {x[0]=x[1]=x[2]=0;} // prealeocado e inicializado
  punto(double xx, double yy, double zz=0) {x[0]=xx; x[1]=yy; x[2]=zz;}
  punto(const punto& p) {x[0]=p.x[0]; x[1]=p.x[1]; x[2]=p.x[2];}
  punto(const double *xx) {x[0]=xx[0]; x[1]=xx[1]; x[2]=xx[2];}
  // x, y, z de otros
  punto(const punto &px, const punto &py) {x[0]=px.x[0]; x[1]=py.x[1]; x[2]=0;}
  punto(const punto &px, const punto &py, const punto &pz)
    {x[0]=px.x[0]; x[1]=py.x[1]; x[2]=pz.x[2];}

  // destructor
  ~punto() {}

  //componenete
  double& operator [](int n) {return x[n];}
  const double& operator [](int n) const {return x[n];}

  //referencia
  operator const double*() const {return x;}

  // copia
  punto& operator=(const double* p) 
    {x[0]=p[0]; x[1]=p[1]; x[2]=p[2]; return *this;}
  punto& operator=(const punto& p) 
    {x[0]=p.x[0]; x[1]=p.x[1]; x[2]=p.x[2]; return *this;}
  // copia parcial
  punto& set_x(const punto &p){x[0]=p.x[0]; return *this;}
  punto& set_y(const punto &p){x[1]=p.x[1]; return *this;}
  punto& set_z(const punto &p){x[2]=p.x[2]; return *this;}
  punto& set_xy(const punto &p){x[0]=p.x[0];x[1]=p.x[1]; return *this;}
  punto& set_xz(const punto &p){x[0]=p.x[0];x[2]=p.x[2]; return *this;}
  punto& set_yz(const punto &p){x[1]=p.x[1];x[2]=p.x[2]; return *this;}

  punto& zero(){x[0]=x[1]=x[2]=0; return *this;}

  // comparacion
  bool eq_b(const double*,const double &r=ERRADM) const; // bola
  bool eq_b(const punto& p,const double &r=ERRADM) const {return eq_b(p.x,r);}
  bool eq_b2(const double*,const double &r2=ERRADM) const; // bola dado r^2
  bool eq_b2(const punto& p,const double &r2=ERRADM) const;
  bool eq_c(const double*,const double &r=ERRADM) const; // cubo
  bool eq_c(const punto& p,const double &r=ERRADM) const {return eq_c(p.x,r);}
  bool operator ==(const double* p) const {return eq_c(p);}
  bool operator ==(const punto& p) const {return eq_c(p.x);}
  bool operator !=(const double* p) const {return !eq_c(p);}
  bool operator !=(const punto& p) const {return !eq_c(p.x);}

  //compara y asigna componente a componente
  void set_min_max(punto& pmin_guardado,punto& pmax_guardado) const;
/*  friend void set_min_max(punto& pmin, punto& pmax, const punto& p)
    {p.set_min_max(pmin,pmax);}*/

  // operadores y operaciones
  punto& operator+=(const punto&);  // +=
  punto& operator-=(const punto&);  // -=
  punto& operator*=(double);        // *=
  punto& operator/=(double);        // /=

  punto operator-() const;              // opuesto
  punto operator+(const punto&)  const; // suma
  punto operator-(const punto&)  const; // resta
  punto operator%(const punto&)  const; // producto vectorial (% =precedencia *)
  double operator*(const punto&) const; // producto escalar
  punto operator*(double)        const; // escalar * punto y viceversa
  friend punto operator*(double f,const punto &p) {return (p*f);}
  punto operator/(double)        const; // cociente por un escalar (sin test)
  double triple(const punto&,const punto&) const; // triple producto
  friend double triple (const punto &p0,const punto &p1 ,const punto &p2)
    {return p0.triple(p1,p2);}

  double pv2d(const punto &p) const // producto vectorial 2d
    {return (x[0]*p.x[1]-x[1]*p.x[0]);}
  friend double pv2d(const punto &p1,const punto &p2)
    {return p1.pv2d(p2);}

  // modulo
  double mod() const;    // modulo
  friend double mod(const punto &p) {return p.mod();}
  double mod2() const;   // modulo al cuadrado
  friend double mod2(const punto &p) {return p.mod2();}   
  double modc() const;   // maxima componente (abs)
  friend double modc(const punto &p) {return p.modc();}

  punto &dir() {return operator/=(mod());} // transforma este en versor (sin test de 0)

  // distancia
  double distancia(const punto&) const;  // distancia
  friend double distancia(const punto &a,const punto &b)
    {return a.distancia(b);}
  double distancia2(const punto&) const; // distancia al cuadrado
  friend double distancia2(const punto &a,const punto &b)
    {return a.distancia2(b);}
  double distanciac(const punto&) const; // maximo abs delta x y o z
  friend double distanciac(const punto &a, const punto &b)
    {return a.distanciac(b);}
  // calcula rapido si es < dmin (si en una componente es mayor vuelve)
  bool distancia_menor(const punto&, double &dmin) const;
  bool distancia2_menor(const punto&, double &dmin2) const;
  bool distanciac_menor(const punto&, double &dminc) const;

  // angulos, giros y 2d
  double angle(const punto& ref) const // angulo respecto a ref [0,pi]
    {return fabs(atan2(operator %(ref).mod() , operator *(ref)));}
  double angle() const  // angulo respecto al eje x [0,pi]
    {return fabs(atan2(sqrt(x[1]*x[1]+x[2]*x[2]) , x[0]));}
  double anglex() const // angulo  respecto al eje x [0,pi]
    {return fabs(atan2(sqrt(x[1]*x[1]+x[2]*x[2]) , x[0]));}
  double angley() const // angulo  respecto al eje y [0,pi]
    {return fabs(atan2(sqrt(x[2]*x[2]+x[0]*x[0]) , x[1]));}
  double anglez() const // angulo  respecto al eje z [0,pi]
    {return fabs(atan2(sqrt(x[0]*x[0]+x[1]*x[1]) , x[2]));}

  double angle2d(const punto& ref) const // angulo 2d respecto a ref (-pi,pi]
    {return atan2(ref.x[0]*x[1]-ref.x[1]*x[0] , ref.x[0]*x[0]+ref.x[1]*x[1]);}
  double anglex2d() const {return atan2(x[1],x[0]);} // respecto al eje x (-pi,pi]

  punto &giro(const punto &n,double ang); // ang en radianes alrededor de n
  punto &giro90() // gira la proyeccion 90 grados a la izq
    {double t=x[0];x[0]=-x[1];x[1]=t; return *this;}

  punto &polar2d(double a,double l) // polar a partir del punto
    {x[0]+=l*cos(a);x[1]+=l*sin(a); return *this;}

  // diedro [0 - 360) entre pderecho arista y pizquierdo 360=error
  // modo: 0=rad 1=grados 2=sin 3=cos
  friend double diedro(
   const punto &pderecho,
   const punto &arista,
   const punto &pizquierdo,
   int modo=0);

  // miscelaneo
  friend bool inters2d(
    const punto &p00,const punto &p01,
    const punto &p10,const punto &p11,
    punto &p);
  friend double distrr // distancia entre rectas
    (const punto &p00,const punto &p01,
     const punto &p10,const punto &p11)
       {return fabs(((p01-p00)%(p11-p10)).dir()*(p10-p00));}
  punto& prrmedio // punto de minima distancia a dos rectas
    (const punto &p00,const punto &p01,
     const punto &p10,const punto &p11); // this -> pmedio
  friend punto prrmedio // punto de minima distancia a dos rectas
    (const punto &p00,const punto &p01,
     const punto &p10,const punto &p11,
     punto &p02, punto &p12); // calcula los puntos en las rectas y devuelve pmedio
  double distr // distancia a recta (con signo de acuerdo a la orientacion del segmento)
    (const punto &p0, const punto &p1) const
    {return (*this-p0)*((p1-p0).dir().giro90());}
  double distp // distancia a plano (con signo de acuerdo a la orientacion del triangulo)
    (const punto &p0, const punto &p1, const punto &p2) const
    {return (*this-p0)*(((p1-p0)%(p2-p0)).dir());}

  // triangulos y tetraedros
  friend double area2d(const punto &p1,const punto &p2) //p0 es el origen
    {return p1.pv2d(p2)/2;}
  friend double area2d(const punto &p0,const punto &p1,const punto &p2)
    {return (p1-p0).pv2d(p2-p0)/2;}
  friend double area(const punto &p1,const punto &p2) //p0 es el origen
    {return (p1%p2).mod()/2;}
  friend double area(const punto &p0,const punto &p1,const punto &p2)
    {return ((p1-p0)%(p2-p0)).mod()/2;}
  friend double vol(const punto &p1,const punto &p2,const punto &p3) //p0 es el origen
    {return ((p1%p2)*p3)/6;}
  friend double vol(const punto &p0,const punto &p1,const punto &p2,const punto &p3)
    {return (((p1-p0)%(p2-p0))*(p3-p0))/6;}


  //solucion de M[ij]*x[j]=F[i] (3x3)
  //cada punto de M son los coeficientes de una ecuacion
  punto& solve
    (const punto M[3], const punto &f);
  punto& solve // devuelve el denominador
    (const punto M[3], const punto &f, double &d);
  punto& solve // devuelve el denominador y el adjunto
    (const punto M[3], const punto &f, punto Adj[3], double &d);
  friend void transpose(punto M[3]);
  double trace() const {return x[0]+x[1]+x[2];}

  // testear error con r==0
  // v: triple poroducto (6*volumen)  a: prodcto vectorial (modulo=2*area)
  friend void s4 // esfera de 4 puntos
    (const punto &,const punto &, const punto &,const punto &,
       punto &c, double &r, double &v);
  friend void s4 (const punto *p,punto &c, double &r, double &v)    
    {s4(p[0],p[1],p[2],p[3],c,r,v);}
  friend void s40 // esfera de 4 puntos, el primero es el origen
    (const punto &, const punto &,const punto &, punto &c, double &r, double &v);
  friend void c3 // circunferencia de 3 puntos
    (const punto &,const punto &, const punto &, punto &c, double &r, punto &a);
  friend void c3 (const punto *p,punto &c, double &r, punto &a)
    {c3(p[0],p[1],p[2],c,r,a);}
  friend void c30 // circunferencia de 3 puntos, el primero es el origen
    (const punto &, const punto &, punto &c, double &r, punto &a);
  friend punto center(int n,const punto *); // centro de la circunsfera de n puntos

  friend void a3 // arco de 3 puntos (b es el bulge vector!)
    (const punto &,const punto &, const punto &, punto &c, double &r, punto &b);
  friend void a3 (const punto *p,punto &c, double &r, punto &b)
    {a3(p[0],p[1],p[2],c,r,b);}

  // random 
  punto &randdir();   // direccion aleatoria
  punto &randdir2d(); // direccion aleatoria 2D
  punto &randpt();    // componentes aleatorias en cubo [0,1)
  punto &randpt2d();  // componentes aleatorias en cuadrado [0,1)

  //enteros
  punto &entero(){x[0]=int(x[0]);x[1]=int(x[1]);x[2]=int(x[2]);return *this;}
  punto &floor(){x[0]=::floor(x[0]);x[1]=::floor(x[1]);x[2]=::floor(x[2]);return *this;}
  punto &ceil(){x[0]=::ceil(x[0]);x[1]=::ceil(x[1]);x[2]=::ceil(x[2]);return *this;}
  punto &redondea(){x[0]=::floor(x[0]+.5);x[1]=::floor(x[1]+.5);x[2]=::floor(x[2]+.5);return *this;}
  punto redondea(double r) // r es p.ej: .001
    {x[0]=::floor(x[0]/r+.5)*r; x[1]=::floor(x[1]/r+.5)*r; x[2]=::floor(x[2]/r+.5)*r; return *this;}
  friend punto entero(const punto &p){return punto(int(p[0]),int(p[1]),int(p[2]));}
  friend punto floor(const punto &p){return punto(::floor(p[0]),::floor(p[1]),::floor(p[2]));}
  friend punto ceil(const punto &p){return punto(::ceil(p[0]),::ceil(p[1]),::ceil(p[2]));}
  friend punto redondea(const punto &p){return punto(::floor(p[0]+.5),::floor(p[1]+.5),::floor(p[2]+.5));}
  friend punto redondea(const punto &p,double r) // r es p.ej: .001
    {return punto(::floor(p.x[0]/r+.5)*r,::floor(p.x[1]/r+.5)*r,::floor(p.x[2]/r+.5)*r);}

  // overload de >> y <<
  // el default es sin ancho ni precision fijos y separado por tabs
  // el -1 es para que la misma funcion sirva para setear y preguntar
  static bool io_z(int i=-1); // lee/graba z
  static int o_width(int i=-1);  // ancho
  static int o_precision(int i=-1);  // precision
  static const char* o_separator(const char* s=0); // separador
  static void io_reset()
    {io_z(true); o_width(0); o_precision(0); o_separator("\t");}
  friend std::istream& operator >>(std::istream &a, punto &p) {return p.lee(a);}
  friend std::ostream& operator <<(std::ostream &a, const punto &p) {return p.graba(a);}
};

// constantes globales
static const punto pzero(0,0,0);
static const punto _ex(1,0,0),_ey(0,1,0),_ez(0,0,1);
static const punto _e[3]={_ex,_ey,_ez};
static const size_t SZPt=sizeof(punto);

#endif
