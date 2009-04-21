// temporizador

// Para no andar cambiando rutina a rutina las llamadas a las funciones
// y para que sean "gratis" si no se usan, las acemos con el preprocesador 
// las funciones son _initime()  _savetime(rutina) e _infotime(texto)

// Incluirlo no basta, hay que hacer #define temporizar "ANTES" de incluirlo

// time_t y clock_t son long int (time.h)

#ifndef _TIEMPO_H_
#define _TIEMPO_H_

#include <ctime> // time
#include <fstream> // file

// con formato
const char *take_time_f(bool reset=false);
static inline void reset_time_f() {take_time_f(true);}

// cpu en decimas de segundo
inline long int take_time_c(){
  return clock();
}
// tiempo en segundos
inline long int take_time_t(){
  return (long int) (time(NULL));
}

// texto con fecha y hora
const char *dia_y_hora();
const char *dia();

// pausa en segundos
void pausa(int tp=5, const char* message="\nhit any key\n");


///////////////// salida a archivo

class _F_T_{
public:
  static std::ofstream* file_tiempo;
  static long int *initime_stack;
  static int stack_length;
  static int stack_size;

  _F_T_();
  ~_F_T_();

  
  void push_time(); // agrega un tiempo al stack
  void pop_time(const char *);  // graba el tiempo transcurrido y el texto
};

#ifdef temporizar
#pragma message("<---------------- ***   TEMPORIZADO -> tiempos.txt  ***\n")

  static _F_T_ _FT_;

  #define _initime _FT_.push_time();
  #define _savetime(rutina) _FT_.pop_time(#rutina);
  #define _infotime(texto) (*_FT_.file_tiempo) << texto;

#else

  #define _initime
  #define _savetime(rutina)
  #define _infotime(texto)

#endif


#endif
