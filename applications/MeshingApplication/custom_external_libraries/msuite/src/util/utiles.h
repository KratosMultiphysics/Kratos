#ifndef _UTIL_H
#define _UTIL_H

#include <cstddef> // size_t

#ifdef _DEBUG
  #ifdef _MSC_VER
    // windows
    #define _revienta(cond) if (cond) __asm int 3;
  #else
    // no - windows
    #include <cassert>
    #define _revienta(cond) assert(!(cond));
  #endif
#else
  // release
  #define _revienta(cond) ((void)0)
#endif

#ifdef _WIN32
#define _erase "deltree /y "
#define _bar '\\'
#else
#define _erase "rm -f -r "
#define _bar '/'
#endif

typedef unsigned int uint;

#ifndef REAL
// float        4 3.4E +/- 38   (7 digits)
// double       8 1.7E +/- 308  (15 digits)
// long double 10 1.2E +/- 4932 (19 digits)
  #define REAL double
  #define real double
#endif
static const double ERRADM=1e-12;
static const double MINREAL=1e-300;
static const double MAXREAL=1e+300;

static const double PI       =  3.1415926535897932384626433832795;
static const double DOSPI    =  6.28318530717958647692528676655901;
static const double PIDOS    =  1.57079632679489661923132169163975;
static const double PICUATRO =  0.785398163397448309615660845819876;
static const double G2R      =  0.0174532925199432957692369076848861;
static const double R2G       = 57.2957795130823208767981548141052;
static const double SQRT2    =  1.4142135623730950488016887242097;
static const double SQRT2_2  =  0.707106781186547524400844362104849;
static const double SQRT3    =  1.73205080756887729352744634150587;

#define nada (~0) // flag entero, depende de si es signed o unsigned

static const size_t SZB=sizeof(bool);
static const size_t SZS=sizeof(short);
static const size_t SZI=sizeof(int);
static const size_t SZD=sizeof(double);
static const size_t SZP=sizeof(void *);

static inline int ciclo(int indice, int longitud)
  {int j=(indice%longitud); if (j<0) return j+longitud; else return j;}

//================================
// bits el primero es el 0 (b es numero (1,2,4,8...), no bit)
#define _mk1(n,b) n|=b
#define _mk0(n,b) n&=~b
#define _is1(n,b) n&b // tambien sirve como _is1_any
#define _is0(n,b) ~n&b// tambien sirve como _is0_any
#define _is1_all(n,b) n&b==b
#define _is0_all(n,b) ~n&b==b

//================================
static inline REAL g2r(REAL x) {return (G2R*x);}
static inline REAL r2g(REAL x) {return (R2G*x);}

//================================
// indice vectorial de matriz triangular superior con diagonal (c>=f)
static inline int TSD(int f,int c, int n) {return ((f*((n<<1)-f+1))>>1)+c;}
// indice vectorial de matriz triangular superior sin diagonal (c>f)
static inline int TS(int f,int c, int n) {return ((f*(((n-1)<<1)-f+1))>>1)+c-1;}

// cambia el contenido y no la direccion (lento pero seguro)
template <class T> static inline void Swap(T&t1, T&t2) //swap en conflicto con stl
  {T t0=t1; t1=t2; t2=t0;}

//#define min(a,b)  (((a) < (b)) ? (a) : (b)) // en stdlib, pero...
// minimo
//#define Min(t1, t2)  {(t1<t2) ? t1 : t2;}
template <class T> static inline const T &Min(const T &t1,const T &t2)
{if (t1<t2) return t1; else return t2;}
template <class T> static inline const T &Min3(const T &t1,const T &t2,const T &t3)
{if (t1<t2) {if (t1<t3) return t1; else return t3;} else {if (t2<t3) return t2; else return t3;}}
// maximo
template <class T> static inline const T &Max(const T &t1,const T &t2)
{if (t1>t2) return t1; else return t2;}
template <class T> static inline const T &Max3(const T &t1,const T &t2,const T &t3)
{if (t1>t2) {if (t1>t3) return t1; else return t3;} else {if (t2>t3) return t2; else return t3;}}


// potencia >0
template <class T> static inline T pown(T inp, int p)
{T ret=1;  while (p--) ret*=inp; return ret;} // ojo el 1

// combinatorio (solo un test de validez)
static inline long unsigned int Comb(int qty, int total)
{
  if (!qty||(qty==total)) return 1;
  return Comb(qty-1,total-1)+Comb(qty,total-1);
}

template <class T> inline bool set_min(T& tminguardado,const T& tcompare)
{if (tcompare<tminguardado) {tminguardado=tcompare; return true;} else return false;}
template <class T> inline bool set_max(T& tmaxguardado,const T& tcompare)
{if (tcompare>tmaxguardado) {tmaxguardado=tcompare; return true;} else return false;}
template <class T> inline void set_min_max(T& tminguardado,T& tmaxguardado,const T& tcompare){
  if (tcompare>tmaxguardado) tmaxguardado=tcompare;
  if (tcompare<tminguardado) tminguardado=tcompare;
}

static inline double redondea(double value,double r)
{return (int(value/r+.5))*r;}


// guarda temporariamente el valor de un bool
#define _push_1b(nombre,valor) bool nombre##temp=nombre; nombre=valor;
#define _pop_1b(nombre) nombre=nombre##temp;

#define _push_1i(nombre,valor) int nombre##temp=nombre; nombre=valor;
#define _pop_1i(nombre) nombre=nombre##temp;

#define _push_1r(nombre,valor) double nombre##temp=nombre; nombre=valor;
#define _pop_1r(nombre) nombre=nombre##temp;
#endif
