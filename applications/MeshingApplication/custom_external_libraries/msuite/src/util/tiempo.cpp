#include "tiempo.h"
#include <cstdio> // sprintf
#include <cstdlib> // alloc
#if defined _MSC_VER && _MSC_VER>=1400
  #define sprintf sprintf_s
#endif
#if defined(_WIN32) || defined(_MAC)
  #include <conio.h> // _kbhit
#endif
#include <iostream> // cout
#include <ostream> // <<

using namespace std;

// temporizador con formato
const char *take_time_f(bool reset){
  static char T[128];
//  static const char *format="CPU  Time: %dhs %dmin %ds\nREAL Time: %dhs %dmin %ds\n";
  static const char *format="CPU: %dhs %dmin %ds\n";
  static clock_t c0;
//  static time_t t0;
  if (reset)
    {c0=clock(); /*t0=time(NULL);*/ return T;}

  clock_t c1=clock();
//  time_t t1=time(NULL);
  int tc,sc,mc/*,tr,sr,mr*/;

  tc=int((c1-c0)*1.0/CLOCKS_PER_SEC); // segundos CPU
  sc=tc%60; tc/=60; mc=tc%60; tc/=60;

//  tr=int(difftime(t1,t0)); // segundos
//  sr=tr%60; tr/=60; mr=tr%60; tr/=60;

  sprintf(T,format,tc,mc,sc/*,tr,mr,sr*/);

  return T;
}

const char *dia_y_hora() {
  static char T[128];
  time_t t; time(&t);        /* Get time in seconds */
#if defined _MSC_VER && _MSC_VER>=1400
  struct tm today; localtime_s(&today,&t);
  strftime(T,128,"%d/%m/%Y %H:%M:%S",&today);
#else
  struct tm *today=localtime(&t);        /* Convert time to struct tm form */
  strftime(T,128,"%d/%m/%Y %H:%M:%S",today);
#endif
  return T;
}

const char *dia() {
  static char T[128];
  time_t t; time(&t);        /* Get time in seconds */
#if defined _MSC_VER && _MSC_VER>=1400
  struct tm today; localtime_s(&today,&t);
  strftime(T,128,"%d/%m/%Y",&today);
#else
  struct tm *today=localtime(&t);        /* Convert time to struct tm form */
  strftime(T,128,"%d/%m/%Y",today);
#endif
  return T;
}

// pausa en segundos
void pausa(int tp, const char *message)
{
  if (message&&message[0])
    cout << message; cout.flush();
  long int t0=take_time_t();
#if !defined(_WIN32) && !defined(_MAC)
  while (take_time_t()-t0<tp) ;
#else
  int t=0,k=0;
  while (!(k=_kbhit()) && t<tp) {
    t=take_time_t()-t0;
  }
  if (k) _getch(); // captura el input
# endif
}

//===================
// salida a archivo

ofstream* _F_T_::file_tiempo=0;
long int *_F_T_::initime_stack=0;
int _F_T_::stack_length=0;
int _F_T_::stack_size=0;

_F_T_::_F_T_(){
  if (!file_tiempo){
    file_tiempo=new ofstream("tiempos.txt",ios::out|ios::app);
    (*file_tiempo) << dia_y_hora() << endl;
  }
  if (!initime_stack) {
    stack_size=32;
    initime_stack=(long int *)malloc(stack_size*sizeof(long int));    
  }
}

_F_T_::~_F_T_(){
  if (file_tiempo){
    if (file_tiempo->is_open()){
      (*file_tiempo) << endl;
      file_tiempo->close();
    }
    delete file_tiempo; file_tiempo=0;
  }
  delete [] initime_stack; initime_stack=0;
  stack_length=0; stack_size=0;
}

void _F_T_::push_time(){
  if (stack_length==stack_size-1){
    stack_size+=32;
    initime_stack=(long int *)realloc(initime_stack,stack_size*sizeof(long int));
  }
  initime_stack[stack_length++]=clock();
}

void _F_T_::pop_time(const char *rutina){
  (*file_tiempo) << rutina << ": " 
     << float(((clock()-initime_stack[--stack_length])*100)/CLOCKS_PER_SEC)/100 
     << endl;
}
