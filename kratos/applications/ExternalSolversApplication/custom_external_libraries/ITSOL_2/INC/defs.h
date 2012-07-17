
#define qsplit qsplit_
#define dgemv dgemv_
#define readmtc readmtc_
#define csrcsc csrcsc_
#define coicsr coicsr_
#define roscal roscal_
#define coscal coscal_

#ifdef _SGI
#define dnrm2 dnrm2_
#define ddot ddot_
#define daxpy daxpy_
#else
#ifdef _LINUX
#define dnrm2 dnrm2_
#define ddot ddot_
#define daxpy daxpy_
#else
#ifdef _IBM
#define dnrm2 dnrm2
#define ddot ddot
#define daxpy daxpy
#else
#define dnrm2 dnrm2_
#define ddot ddot_
#define daxpy daxpy_
#endif
#endif
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

#define MAX_LINE        256
#define MAX_HBNAME      64
#define MAX_MAT	   100

