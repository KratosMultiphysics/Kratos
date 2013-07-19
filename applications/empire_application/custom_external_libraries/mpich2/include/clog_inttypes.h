/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
/*
   MPE Logging's generated header file for fixed size integer types
   and their printf format specifiers used in CLOG2.
*/
#if !defined( _CLOG_INTTYPES )
#define  _CLOG_INTTYPES

#include <stdint.h>
#include <inttypes.h>

typedef int8_t         CLOG_int8_t;
typedef int16_t        CLOG_int16_t;
typedef int32_t        CLOG_int32_t;
typedef int64_t        CLOG_int64_t;

/* 
   Define address-sized integer for the use of MPI_Comm_xxx_attr
   in clog2_commset.c.
*/
typedef CLOG_int64_t          CLOG_Pint;

#define i8fmt        "%"PRId8
#define i16fmt       "%"PRId16
#define i32fmt       "%"PRId32
#define i64fmt       "%"PRId64

#endif
