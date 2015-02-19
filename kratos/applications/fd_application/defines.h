#ifndef DEFINES_H
#define DEFINES_H

#define FETCHTIME ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

#ifdef USE_NOVEC
  typedef double  PrecisionType;
  const uint BW              = 2;        //  Boundary width   ( needed if we want to use SIMD )
#endif

#ifdef USE_DOUBLE
  typedef double  PrecisionType;
  const uint BW              = 8;        //  Boundary width   ( needed if we want to use SIMD )
#endif

#ifdef USE_FLOAT
  typedef float  PrecisionType;
  const uint BW              = 16;       //  Boundary width   ( needed if we want to use SIMD )
#endif

const uint BWP             = BW / 2;   //  Boundary padding ( needed if we want to use SIMD )
const PrecisionType ONESIX = 1.0/6.0;  //  1/6

#endif