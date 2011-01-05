#ifdef _OPENMP
#include <omp.h>
#else
  #ifndef NO_TIMER
    #include "ctime"
  #endif
#endif
double SuperLU_timer_() {
	    #ifndef _OPENMP
	      #ifndef NO_TIMER
		  return std::clock()/static_cast<double>(CLOCKS_PER_SEC);
	      #else
		  return 0.0;
	      #endif
	    #else
		  return  omp_get_wtime();
	    #endif
	  }

