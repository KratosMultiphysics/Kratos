#pragma once

#ifdef _WIN32
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <unistd.h>
#include <sys/times.h>
#include <sys/timeb.h>
#include <time.h> // CLK_TCK
#  ifndef CLK_TCK
#   define CLK_TCK      CLOCKS_PER_SEC
#  endif
#endif

inline clock_t currentTime() {
// #ifndef _WIN32
//   struct tms time;
//   times(&time);
//   // return time.tms_utime; // CPU time
//   // return time.tms_stime; // sys time
//   return time.tms_utime + time.tms_stime; // CPU + sys time
// #else
   return clock();
// #endif
}

// inline time_t currentWallTime() {
//   return time();
// }
inline int currentWallTime( struct timeb *tp) {
#ifdef _WIN32
  // MS Wind de los co...... que no respeta los estandares.
  ftime( tp);
  return 0;
#elif __linux__
  // ftime is being deprecated in linux
  struct timespec tmp;
  int res = clock_gettime( CLOCK_REALTIME, &tmp);
  tp->time = tmp.tv_sec;
  tp->millitm = ( unsigned short)( tmp.tv_nsec / 1000);
  return res;
#else
  // macOS
  return ftime( tp);
#endif
}

inline bool equalWallTime( const struct timeb &a, const struct timeb &b ) {
  return ( a.millitm == b.millitm ) && ( a.time == b.time );
}

typedef enum {
  CRONO_SYS_TIME = 0, CRONO_WALL_TIME = 1
} _t_crono_measure_type;

class Crono {
public:
  void ini( void) {
    if ( _time_type == CRONO_SYS_TIME)
      _clock_start = currentTime();
    else
      currentWallTime( &_time_start);
  }
  float get( void) {
    float res = 0.0;
    if ( _time_type == CRONO_SYS_TIME)
      res = ( float)( currentTime() - _clock_start) / ( float)CLK_TCK;
    else {
      struct timeb tmp;
      currentWallTime( &tmp);
      time_t dif_sec = tmp.time - _time_start.time;
      int dif_msec = tmp.millitm - _time_start.millitm;
      res = ( float)( ( double)( dif_sec * 1000 + dif_msec) / 1000.0);
    }
    return res;
  }
  float end( void) {
    float res = get();
    ini();
    return res;
  }
  Crono( _t_crono_measure_type tt = CRONO_SYS_TIME): _time_type( tt) {
    ini();
  };

private:
  _t_crono_measure_type _time_type;
  clock_t _clock_start;
  struct timeb _time_start;
};
