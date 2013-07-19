/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
/* src/include/mpi.h.  Generated from mpi.h.in by configure. */
#ifndef MPI_INCLUDED
#define MPI_INCLUDED

/* user include file for MPI programs */

/* Keep C++ compilers from getting confused */
#if defined(__cplusplus)
extern "C" {
#endif

#if defined(__has_attribute)
#  if __has_attribute(pointer_with_type_tag) && \
      __has_attribute(type_tag_for_datatype) && \
      !defined(MPICH_NO_ATTR_TYPE_TAGS)
#    define MPICH_ATTR_POINTER_WITH_TYPE_TAG(buffer_idx, type_idx)  __attribute__((pointer_with_type_tag(MPI,buffer_idx,type_idx)))
#    define MPICH_ATTR_TYPE_TAG(type)                               __attribute__((type_tag_for_datatype(MPI,type)))
#    define MPICH_ATTR_TYPE_TAG_LAYOUT_COMPATIBLE(type)             __attribute__((type_tag_for_datatype(MPI,type,layout_compatible)))
#    define MPICH_ATTR_TYPE_TAG_MUST_BE_NULL()                      __attribute__((type_tag_for_datatype(MPI,void,must_be_null)))
#    include <stddef.h>
#  endif
#endif

#if !defined(MPICH_ATTR_POINTER_WITH_TYPE_TAG)
#  define MPICH_ATTR_POINTER_WITH_TYPE_TAG(buffer_idx, type_idx)
#  define MPICH_ATTR_TYPE_TAG(type)
#  define MPICH_ATTR_TYPE_TAG_LAYOUT_COMPATIBLE(type)
#  define MPICH_ATTR_TYPE_TAG_MUST_BE_NULL()
#endif

#if !defined(INT8_C)
/* stdint.h was not included, see if we can get it */
#  if defined(__cplusplus)
#    if __cplusplus >= 201103
#      include <cstdint>
#    endif
#  endif
#endif

#if !defined(INT8_C)
/* stdint.h was not included, see if we can get it */
#  if defined(__STDC_VERSION__)
#    if __STDC_VERSION__ >= 199901
#      include <stdint.h>
#    endif
#  endif
#endif

#if defined(INT8_C)
/* stdint.h was included, so we can annotate these types */
#  define MPICH_ATTR_TYPE_TAG_STDINT(type) MPICH_ATTR_TYPE_TAG(type)
#else
#  define MPICH_ATTR_TYPE_TAG_STDINT(type)
#endif

#ifdef __STDC_VERSION__ 
#if __STDC_VERSION__ >= 199901
#  define MPICH_ATTR_TYPE_TAG_C99(type) MPICH_ATTR_TYPE_TAG(type)
#else
#  define MPICH_ATTR_TYPE_TAG_C99(type)
#endif
#else 
#  define MPICH_ATTR_TYPE_TAG_C99(type)
#endif

#if defined(__cplusplus)
#  define MPICH_ATTR_TYPE_TAG_CXX(type) MPICH_ATTR_TYPE_TAG(type)
#else
#  define MPICH_ATTR_TYPE_TAG_CXX(type)
#endif

/* Results of the compare operations. */
#define MPI_IDENT     0
#define MPI_CONGRUENT 1
#define MPI_SIMILAR   2
#define MPI_UNEQUAL   3

typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)0x4c000101)
#define MPI_SIGNED_CHAR    ((MPI_Datatype)0x4c000118)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)0x4c000102)
#define MPI_BYTE           ((MPI_Datatype)0x4c00010d)
#define MPI_WCHAR          ((MPI_Datatype)0x4c00040e)
#define MPI_SHORT          ((MPI_Datatype)0x4c000203)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)0x4c000204)
#define MPI_INT            ((MPI_Datatype)0x4c000405)
#define MPI_UNSIGNED       ((MPI_Datatype)0x4c000406)
#define MPI_LONG           ((MPI_Datatype)0x4c000807)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)0x4c000808)
#define MPI_FLOAT          ((MPI_Datatype)0x4c00040a)
#define MPI_DOUBLE         ((MPI_Datatype)0x4c00080b)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)0x4c00100c)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)0x4c000809)
#define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_LONG_LONG      MPI_LONG_LONG_INT

static const MPI_Datatype mpich_mpi_char               MPICH_ATTR_TYPE_TAG(char)               = MPI_CHAR;
static const MPI_Datatype mpich_mpi_signed_char        MPICH_ATTR_TYPE_TAG(signed char)        = MPI_SIGNED_CHAR;
static const MPI_Datatype mpich_mpi_unsigned_char      MPICH_ATTR_TYPE_TAG(unsigned char)      = MPI_UNSIGNED_CHAR;
/*static const MPI_Datatype mpich_mpi_byte               MPICH_ATTR_TYPE_TAG(char)               = MPI_BYTE;*/
static const MPI_Datatype mpich_mpi_wchar              MPICH_ATTR_TYPE_TAG(wchar_t)            = MPI_WCHAR;
static const MPI_Datatype mpich_mpi_short              MPICH_ATTR_TYPE_TAG(short)              = MPI_SHORT;
static const MPI_Datatype mpich_mpi_unsigned_short     MPICH_ATTR_TYPE_TAG(unsigned short)     = MPI_UNSIGNED_SHORT;
static const MPI_Datatype mpich_mpi_int                MPICH_ATTR_TYPE_TAG(int)                = MPI_INT;
static const MPI_Datatype mpich_mpi_unsigned           MPICH_ATTR_TYPE_TAG(unsigned)           = MPI_UNSIGNED;
static const MPI_Datatype mpich_mpi_long               MPICH_ATTR_TYPE_TAG(long)               = MPI_LONG;
static const MPI_Datatype mpich_mpi_unsigned_long      MPICH_ATTR_TYPE_TAG(unsigned long)      = MPI_UNSIGNED_LONG;
static const MPI_Datatype mpich_mpi_float              MPICH_ATTR_TYPE_TAG(float)              = MPI_FLOAT;
static const MPI_Datatype mpich_mpi_double             MPICH_ATTR_TYPE_TAG(double)             = MPI_DOUBLE;
static const MPI_Datatype mpich_mpi_long_double        MPICH_ATTR_TYPE_TAG(long double)        = MPI_LONG_DOUBLE;
static const MPI_Datatype mpich_mpi_long_long_int      MPICH_ATTR_TYPE_TAG(long long int)      = MPI_LONG_LONG_INT;
static const MPI_Datatype mpich_mpi_unsigned_long_long MPICH_ATTR_TYPE_TAG(unsigned long long) = MPI_UNSIGNED_LONG_LONG;

#define MPI_PACKED         ((MPI_Datatype)0x4c00010f)
#define MPI_LB             ((MPI_Datatype)0x4c000010)
#define MPI_UB             ((MPI_Datatype)0x4c000011)

/* 
   The layouts for the types MPI_DOUBLE_INT etc are simply
   struct { 
       double var;
       int    loc;
   }
   This is documented in the man pages on the various datatypes.   
 */
#define MPI_FLOAT_INT         ((MPI_Datatype)0x8c000000)
#define MPI_DOUBLE_INT        ((MPI_Datatype)0x8c000001)
#define MPI_LONG_INT          ((MPI_Datatype)0x8c000002)
#define MPI_SHORT_INT         ((MPI_Datatype)0x8c000003)
#define MPI_2INT              ((MPI_Datatype)0x4c000816)
#define MPI_LONG_DOUBLE_INT   ((MPI_Datatype)0x8c000004)

struct mpich_struct_mpi_float_int       { float f; int i; };
struct mpich_struct_mpi_double_int      { double d; int i; };
struct mpich_struct_mpi_long_int        { long l; int i; };
struct mpich_struct_mpi_short_int       { short s; int i; };
struct mpich_struct_mpi_2int            { int i1; int i2; };
struct mpich_struct_mpi_long_double_int { long double ld; int i; };

static const MPI_Datatype mpich_mpi_float_int       MPICH_ATTR_TYPE_TAG(struct mpich_struct_mpi_float_int)       = MPI_FLOAT_INT;
static const MPI_Datatype mpich_mpi_double_int      MPICH_ATTR_TYPE_TAG(struct mpich_struct_mpi_double_int)      = MPI_DOUBLE_INT;
static const MPI_Datatype mpich_mpi_long_int        MPICH_ATTR_TYPE_TAG(struct mpich_struct_mpi_long_int)        = MPI_LONG_INT;
static const MPI_Datatype mpich_mpi_short_int       MPICH_ATTR_TYPE_TAG(struct mpich_struct_mpi_short_int)       = MPI_SHORT_INT;
static const MPI_Datatype mpich_mpi_2int            MPICH_ATTR_TYPE_TAG(struct mpich_struct_mpi_2int)            = MPI_2INT;
static const MPI_Datatype mpich_mpi_long_double_int MPICH_ATTR_TYPE_TAG(struct mpich_struct_mpi_long_double_int) = MPI_LONG_DOUBLE_INT;

/* Fortran types */
#define MPI_COMPLEX           ((MPI_Datatype)1275070494)
#define MPI_DOUBLE_COMPLEX    ((MPI_Datatype)1275072546)
#define MPI_LOGICAL           ((MPI_Datatype)1275069469)
#define MPI_REAL              ((MPI_Datatype)1275069468)
#define MPI_DOUBLE_PRECISION  ((MPI_Datatype)1275070495)
#define MPI_INTEGER           ((MPI_Datatype)1275069467)
#define MPI_2INTEGER          ((MPI_Datatype)1275070496)
/* 
 * MPI_2COMPLEX and MPI_2DOUBLE_COMPLEX were defined by accident in 
 * MPI 1.0 and removed in MPI 1.1.  
 *
 * This definition provides backward compatibility.  These definitions
 * will be removed in a subsequent MPICH release
 */
#ifdef MPICH_DEFINE_2COMPLEX
#define MPI_2COMPLEX          ((MPI_Datatype)1275072548)
#define MPI_2DOUBLE_COMPLEX   ((MPI_Datatype)1275076645)
#endif 
#define MPI_2REAL             ((MPI_Datatype)1275070497)
#define MPI_2DOUBLE_PRECISION ((MPI_Datatype)1275072547)
#define MPI_CHARACTER         ((MPI_Datatype)1275068698)

/* Size-specific types (see MPI-2, 10.2.5) */
#define MPI_REAL4             ((MPI_Datatype)0x4c000427)
#define MPI_REAL8             ((MPI_Datatype)0x4c000829)
#define MPI_REAL16            ((MPI_Datatype)0x4c00102b)
#define MPI_COMPLEX8          ((MPI_Datatype)0x4c000828)
#define MPI_COMPLEX16         ((MPI_Datatype)0x4c00102a)
#define MPI_COMPLEX32         ((MPI_Datatype)0x4c00202c)
#define MPI_INTEGER1          ((MPI_Datatype)0x4c00012d)
#define MPI_INTEGER2          ((MPI_Datatype)0x4c00022f)
#define MPI_INTEGER4          ((MPI_Datatype)0x4c000430)
#define MPI_INTEGER8          ((MPI_Datatype)0x4c000831)
#define MPI_INTEGER16         ((MPI_Datatype)MPI_DATATYPE_NULL)

/* C99 fixed-width datatypes */
#define MPI_INT8_T            ((MPI_Datatype)0x4c000137)
#define MPI_INT16_T           ((MPI_Datatype)0x4c000238)
#define MPI_INT32_T           ((MPI_Datatype)0x4c000439)
#define MPI_INT64_T           ((MPI_Datatype)0x4c00083a)
#define MPI_UINT8_T           ((MPI_Datatype)0x4c00013b)
#define MPI_UINT16_T          ((MPI_Datatype)0x4c00023c)
#define MPI_UINT32_T          ((MPI_Datatype)0x4c00043d)
#define MPI_UINT64_T          ((MPI_Datatype)0x4c00083e)

static const MPI_Datatype mpich_mpi_int8_t   MPICH_ATTR_TYPE_TAG_STDINT(int8_t)   = MPI_INT8_T;
static const MPI_Datatype mpich_mpi_int16_t  MPICH_ATTR_TYPE_TAG_STDINT(int16_t)  = MPI_INT16_T;
static const MPI_Datatype mpich_mpi_int32_t  MPICH_ATTR_TYPE_TAG_STDINT(int32_t)  = MPI_INT32_T;
static const MPI_Datatype mpich_mpi_int64_t  MPICH_ATTR_TYPE_TAG_STDINT(int64_t)  = MPI_INT64_T;
static const MPI_Datatype mpich_mpi_uint8_t  MPICH_ATTR_TYPE_TAG_STDINT(uint8_t)  = MPI_UINT8_T;
static const MPI_Datatype mpich_mpi_uint16_t MPICH_ATTR_TYPE_TAG_STDINT(uint16_t) = MPI_UINT16_T;
static const MPI_Datatype mpich_mpi_uint32_t MPICH_ATTR_TYPE_TAG_STDINT(uint32_t) = MPI_UINT32_T;
static const MPI_Datatype mpich_mpi_uint64_t MPICH_ATTR_TYPE_TAG_STDINT(uint64_t) = MPI_UINT64_T;

/* other C99 types */
#define MPI_C_BOOL                 ((MPI_Datatype)0x4c00013f)
#define MPI_C_FLOAT_COMPLEX        ((MPI_Datatype)0x4c000840)
#define MPI_C_COMPLEX              MPI_C_FLOAT_COMPLEX
#define MPI_C_DOUBLE_COMPLEX       ((MPI_Datatype)0x4c001041)
#define MPI_C_LONG_DOUBLE_COMPLEX  ((MPI_Datatype)0x4c002042)

static const MPI_Datatype mpich_mpi_c_bool                MPICH_ATTR_TYPE_TAG_C99(_Bool)           = MPI_C_BOOL;
static const MPI_Datatype mpich_mpi_c_float_complex       MPICH_ATTR_TYPE_TAG_C99(float _Complex)  = MPI_C_FLOAT_COMPLEX;
static const MPI_Datatype mpich_mpi_c_double_complex      MPICH_ATTR_TYPE_TAG_C99(double _Complex) = MPI_C_DOUBLE_COMPLEX;
static const MPI_Datatype mpich_mpi_c_long_double_complex MPICH_ATTR_TYPE_TAG_C99(double _Complex) = MPI_C_LONG_DOUBLE_COMPLEX;

/* address/offset types */
#define MPI_AINT          ((MPI_Datatype)0x4c000843)
#define MPI_OFFSET        ((MPI_Datatype)0x4c000844)

/* typeclasses */
#define MPI_TYPECLASS_REAL 1
#define MPI_TYPECLASS_INTEGER 2
#define MPI_TYPECLASS_COMPLEX 3

/* Communicators */
typedef int MPI_Comm;
#define MPI_COMM_WORLD ((MPI_Comm)0x44000000)
#define MPI_COMM_SELF  ((MPI_Comm)0x44000001)

/* Groups */
typedef int MPI_Group;
#define MPI_GROUP_EMPTY ((MPI_Group)0x48000000)

/* RMA and Windows */
typedef int MPI_Win;
#define MPI_WIN_NULL ((MPI_Win)0x20000000)

/* File and IO */
/* This define lets ROMIO know that MPI_File has been defined */
#define MPI_FILE_DEFINED
/* ROMIO uses a pointer for MPI_File objects.  This must be the same definition
   as in src/mpi/romio/include/mpio.h.in  */
typedef struct ADIOI_FileD *MPI_File;
#define MPI_FILE_NULL ((MPI_File)0)

/* Collective operations */
typedef int MPI_Op;

#define MPI_MAX     (MPI_Op)(0x58000001)
#define MPI_MIN     (MPI_Op)(0x58000002)
#define MPI_SUM     (MPI_Op)(0x58000003)
#define MPI_PROD    (MPI_Op)(0x58000004)
#define MPI_LAND    (MPI_Op)(0x58000005)
#define MPI_BAND    (MPI_Op)(0x58000006)
#define MPI_LOR     (MPI_Op)(0x58000007)
#define MPI_BOR     (MPI_Op)(0x58000008)
#define MPI_LXOR    (MPI_Op)(0x58000009)
#define MPI_BXOR    (MPI_Op)(0x5800000a)
#define MPI_MINLOC  (MPI_Op)(0x5800000b)
#define MPI_MAXLOC  (MPI_Op)(0x5800000c)
#define MPI_REPLACE (MPI_Op)(0x5800000d)
#define MPIX_NO_OP  (MPI_Op)(0x5800000e)

/* Permanent key values */
/* C Versions (return pointer to value),
   Fortran Versions (return integer value).
   Handled directly by the attribute value routine
   
   DO NOT CHANGE THESE.  The values encode:
   builtin kind (0x1 in bit 30-31)
   Keyval object (0x9 in bits 26-29)
   for communicator (0x1 in bits 22-25)
   
   Fortran versions of the attributes are formed by adding one to
   the C version.
 */
#define MPI_TAG_UB           0x64400001
#define MPI_HOST             0x64400003
#define MPI_IO               0x64400005
#define MPI_WTIME_IS_GLOBAL  0x64400007
#define MPI_UNIVERSE_SIZE    0x64400009
#define MPI_LASTUSEDCODE     0x6440000b
#define MPI_APPNUM           0x6440000d

/* In addition, there are 5 predefined window attributes that are
   defined for every window */
#define MPI_WIN_BASE          0x66000001
#define MPI_WIN_SIZE          0x66000003
#define MPI_WIN_DISP_UNIT     0x66000005
#define MPIX_WIN_CREATE_FLAVOR 0x66000007
#define MPIX_WIN_MODEL         0x66000009

/* Define some null objects */
#define MPI_COMM_NULL      ((MPI_Comm)0x04000000)
#define MPI_OP_NULL        ((MPI_Op)0x18000000)
#define MPI_GROUP_NULL     ((MPI_Group)0x08000000)
#define MPI_DATATYPE_NULL  ((MPI_Datatype)0x0c000000)
#define MPI_REQUEST_NULL   ((MPI_Request)0x2c000000)
#define MPI_ERRHANDLER_NULL ((MPI_Errhandler)0x14000000)
#define MPIX_MESSAGE_NULL   ((MPIX_Message)MPI_REQUEST_NULL)
#define MPIX_MESSAGE_NO_PROC ((MPIX_Message)0x6c000000)

static const MPI_Datatype mpich_mpi_datatype_null MPICH_ATTR_TYPE_TAG_MUST_BE_NULL() = MPI_DATATYPE_NULL;

/* These are only guesses; make sure you change them in mpif.h as well */
#define MPI_MAX_PROCESSOR_NAME 128
#define MPI_MAX_ERROR_STRING   1024
#define MPI_MAX_PORT_NAME      256
#define MPI_MAX_OBJECT_NAME    128

/* Pre-defined constants */
#define MPI_UNDEFINED      (-32766)
#define MPI_KEYVAL_INVALID 0x24000000

/* MPI-3 window flavors */
typedef enum MPIR_Win_flavor_e {
    MPIX_WIN_FLAVOR_CREATE      = 1,
    MPIX_WIN_FLAVOR_ALLOCATE    = 2,
    MPIX_WIN_FLAVOR_DYNAMIC     = 3,
    MPIX_WIN_FLAVOR_SHARED      = 4
} MPIR_Win_flavor_t;

/* MPI-3 window consistency models */
typedef enum MPIR_Win_model_e {
    MPIX_WIN_SEPARATE   = 1,
    MPIX_WIN_UNIFIED    = 2
} MPIR_Win_model_t;

/* Upper bound on the overhead in bsend for each message buffer */
#define MPI_BSEND_OVERHEAD 88

/* Topology types */
typedef enum MPIR_Topo_type { MPI_GRAPH=1, MPI_CART=2, MPI_DIST_GRAPH=3 } MPIR_Topo_type;

#define MPI_BOTTOM      (void *)0
extern int * const MPI_UNWEIGHTED;

#define MPI_PROC_NULL   (-1)
#define MPI_ANY_SOURCE 	(-2)
#define MPI_ROOT        (-3)
#define MPI_ANY_TAG     (-1)

#define MPI_LOCK_EXCLUSIVE  234
#define MPI_LOCK_SHARED     235

#ifndef MPICH2_CONST
/* #define MPICH2_CONST    const */
#define MPICH2_CONST
#endif

/* C functions */
typedef void (MPI_Handler_function) ( MPI_Comm *, int *, ... );
typedef int (MPI_Comm_copy_attr_function)(MPI_Comm, int, void *, void *, 
					  void *, int *);
typedef int (MPI_Comm_delete_attr_function)(MPI_Comm, int, void *, void *);
typedef int (MPI_Type_copy_attr_function)(MPI_Datatype, int, void *, void *, 
					  void *, int *);
typedef int (MPI_Type_delete_attr_function)(MPI_Datatype, int, void *, void *);
typedef int (MPI_Win_copy_attr_function)(MPI_Win, int, void *, void *, void *,
					 int *);
typedef int (MPI_Win_delete_attr_function)(MPI_Win, int, void *, void *);
/* added in MPI-2.2 */
typedef void (MPI_Comm_errhandler_function)(MPI_Comm *, int *, ...);
typedef void (MPI_File_errhandler_function)(MPI_File *, int *, ...);
typedef void (MPI_Win_errhandler_function)(MPI_Win *, int *, ...);
/* names that were added in MPI-2.0 and deprecated in MPI-2.2 */
typedef MPI_Comm_errhandler_function MPI_Comm_errhandler_fn;
typedef MPI_File_errhandler_function MPI_File_errhandler_fn;
typedef MPI_Win_errhandler_function MPI_Win_errhandler_fn;

/* Built in (0x1 in 30-31), errhandler (0x5 in bits 26-29, allkind (0
   in 22-25), index in the low bits */
#define MPI_ERRORS_ARE_FATAL ((MPI_Errhandler)0x54000000)
#define MPI_ERRORS_RETURN    ((MPI_Errhandler)0x54000001)
/* MPIR_ERRORS_THROW_EXCEPTIONS is not part of the MPI standard, it is here to
   facilitate the c++ binding which has MPI::ERRORS_THROW_EXCEPTIONS. 
   Using the MPIR prefix preserved the MPI_ names for objects defined by
   the standard. */
#define MPIR_ERRORS_THROW_EXCEPTIONS ((MPI_Errhandler)0x54000002)
typedef int MPI_Errhandler;

/* Make the C names for the dup function mixed case.
   This is required for systems that use all uppercase names for Fortran 
   externals.  */
/* MPI 1 names */
#define MPI_NULL_COPY_FN   ((MPI_Copy_function *)0)
#define MPI_NULL_DELETE_FN ((MPI_Delete_function *)0)
#define MPI_DUP_FN         MPIR_Dup_fn
/* MPI 2 names */
#define MPI_COMM_NULL_COPY_FN ((MPI_Comm_copy_attr_function*)0)
#define MPI_COMM_NULL_DELETE_FN ((MPI_Comm_delete_attr_function*)0)
#define MPI_COMM_DUP_FN  ((MPI_Comm_copy_attr_function *)MPI_DUP_FN)
#define MPI_WIN_NULL_COPY_FN ((MPI_Win_copy_attr_function*)0)
#define MPI_WIN_NULL_DELETE_FN ((MPI_Win_delete_attr_function*)0)
#define MPI_WIN_DUP_FN   ((MPI_Win_copy_attr_function*)MPI_DUP_FN)
#define MPI_TYPE_NULL_COPY_FN ((MPI_Type_copy_attr_function*)0)
#define MPI_TYPE_NULL_DELETE_FN ((MPI_Type_delete_attr_function*)0)
#define MPI_TYPE_DUP_FN ((MPI_Type_copy_attr_function*)MPI_DUP_FN)

/* MPI request opjects */
typedef int MPI_Request;

/* MPI message objects for Mprobe and related functions */
typedef int MPIX_Message;

/* User combination function */
typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * ); 

/* MPI Attribute copy and delete functions */
typedef int (MPI_Copy_function) ( MPI_Comm, int, void *, void *, void *, int * );
typedef int (MPI_Delete_function) ( MPI_Comm, int, void *, void * );

#define MPI_VERSION    2
#define MPI_SUBVERSION 2
#define MPICH_NAME     2
#define MPICH2         1
#define MPICH_HAS_C2F  1


/* MPICH2_VERSION is the version string. MPICH2_NUMVERSION is the
 * numeric version that can be used in numeric comparisons.
 *
 * MPICH2_VERSION uses the following format:
 * Version: [MAJ].[MIN].[REV][EXT][EXT_NUMBER]
 * Example: 1.0.7rc1 has
 *          MAJ = 1
 *          MIN = 0
 *          REV = 7
 *          EXT = rc
 *          EXT_NUMBER = 1
 *
 * MPICH2_NUMVERSION will convert EXT to a format number:
 *          ALPHA (a) = 0
 *          BETA (b)  = 1
 *          RC (rc)   = 2
 *          PATCH (p) = 3
 * Regular releases are treated as patch 0
 *
 * Numeric version will have 1 digit for MAJ, 2 digits for MIN, 2
 * digits for REV, 1 digit for EXT and 2 digits for EXT_NUMBER. So,
 * 1.0.7rc1 will have the numeric version 10007201.
 */
#define MPICH2_VERSION "1.5"
#define MPICH2_NUMVERSION 10500300

#define MPICH2_RELEASE_TYPE_ALPHA  0
#define MPICH2_RELEASE_TYPE_BETA   1
#define MPICH2_RELEASE_TYPE_RC     2
#define MPICH2_RELEASE_TYPE_PATCH  3

#define MPICH2_CALC_VERSION(MAJOR, MINOR, REVISION, TYPE, PATCH) \
    (((MAJOR) * 10000000) + ((MINOR) * 100000) + ((REVISION) * 1000) + ((TYPE) * 100) + (PATCH))

/* for the datatype decoders */
enum MPIR_Combiner_enum {
    MPI_COMBINER_NAMED            = 1,
    MPI_COMBINER_DUP              = 2,
    MPI_COMBINER_CONTIGUOUS       = 3, 
    MPI_COMBINER_VECTOR           = 4,
    MPI_COMBINER_HVECTOR_INTEGER  = 5,
    MPI_COMBINER_HVECTOR          = 6,
    MPI_COMBINER_INDEXED          = 7,
    MPI_COMBINER_HINDEXED_INTEGER = 8, 
    MPI_COMBINER_HINDEXED         = 9, 
    MPI_COMBINER_INDEXED_BLOCK    = 10, 
    MPIX_COMBINER_HINDEXED_BLOCK  = 11,
    MPI_COMBINER_STRUCT_INTEGER   = 12,
    MPI_COMBINER_STRUCT           = 13,
    MPI_COMBINER_SUBARRAY         = 14,
    MPI_COMBINER_DARRAY           = 15,
    MPI_COMBINER_F90_REAL         = 16,
    MPI_COMBINER_F90_COMPLEX      = 17,
    MPI_COMBINER_F90_INTEGER      = 18,
    MPI_COMBINER_RESIZED          = 19
};

/* for info */
typedef int MPI_Info;
#define MPI_INFO_NULL         ((MPI_Info)0x1c000000)
#define MPI_MAX_INFO_KEY       255
#define MPI_MAX_INFO_VAL      1024

/* for subarray and darray constructors */
#define MPI_ORDER_C              56
#define MPI_ORDER_FORTRAN        57
#define MPI_DISTRIBUTE_BLOCK    121
#define MPI_DISTRIBUTE_CYCLIC   122
#define MPI_DISTRIBUTE_NONE     123
#define MPI_DISTRIBUTE_DFLT_DARG -49767

#define MPI_IN_PLACE  (void *) -1

/* asserts for one-sided communication */
#define MPI_MODE_NOCHECK      1024
#define MPI_MODE_NOSTORE      2048
#define MPI_MODE_NOPUT        4096
#define MPI_MODE_NOPRECEDE    8192
#define MPI_MODE_NOSUCCEED   16384 

/* predefined types for MPIX_Comm_split_type */
#define MPIX_COMM_TYPE_SHARED    1

/* Definitions that are determined by configure. */
typedef long MPI_Aint;
typedef int MPI_Fint;

static const MPI_Datatype mpich_mpi_aint   MPICH_ATTR_TYPE_TAG(MPI_Aint)   = MPI_AINT;

/* FIXME: The following two definition are not defined by MPI and must not be
   included in the mpi.h file, as the MPI namespace is reserved to the MPI 
   standard */
#define MPI_AINT_FMT_DEC_SPEC "%ld"
#define MPI_AINT_FMT_HEX_SPEC "%lx"

/* Let ROMIO know that MPI_Offset is already defined */
#define HAVE_MPI_OFFSET
/* MPI_OFFSET_TYPEDEF is set in configure and is 
      typedef $MPI_OFFSET MPI_Offset;
   where $MPI_OFFSET is the correct C type */
typedef long long MPI_Offset;

static const MPI_Datatype mpich_mpi_offset MPICH_ATTR_TYPE_TAG(MPI_Offset) = MPI_OFFSET;

/* The order of these elements must match that in mpif.h */
typedef struct MPI_Status {
    int count;
    int cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
    
} MPI_Status;

/* types for the MPIX_T_ interface */
struct MPIR_T_enum;
typedef struct MPIR_T_enum * MPIX_T_enum;
struct MPIR_T_cvar_handle;
typedef struct MPIR_T_cvar_handle * MPIX_T_cvar_handle;
struct MPIR_T_pvar_handle;
typedef struct MPIR_T_pvar_handle * MPIX_T_pvar_handle;
struct MPIR_T_pvar_session;
typedef struct MPIR_T_pvar_session * MPIX_T_pvar_session;

/* extra const at front would be safer, but is incompatible with MPIX_T_ prototypes */
extern struct MPIR_T_pvar_handle * const MPIX_T_PVAR_ALL_HANDLES;

#define MPIX_T_ENUM_NULL         ((MPIX_T_enum)NULL)
#define MPIX_T_CVAR_HANDLE_NULL  ((MPIX_T_cvar_handle)NULL)
#define MPIX_T_PVAR_HANDLE_NULL  ((MPIX_T_pvar_handle)NULL)
#define MPIX_T_PVAR_SESSION_NULL ((MPIX_T_pvar_session)NULL)

/* the MPI_T_ interface requires that these VERBOSITY constants occur in this
 * relative order with increasing values */
enum MPIR_T_verbosity_t {
    /* don't name-shift this if/when MPI_T_ is accepted, this is an MPICH2-only
     * extension */
    MPIX_T_VERBOSITY_INVALID = 0,

    /* arbitrarily shift values to aid debugging and reduce accidental errors */
    MPIX_T_VERBOSITY_USER_BASIC = 221,
    MPIX_T_VERBOSITY_USER_DETAIL,
    MPIX_T_VERBOSITY_USER_ALL,

    MPIX_T_VERBOSITY_TUNER_BASIC,
    MPIX_T_VERBOSITY_TUNER_DETAIL,
    MPIX_T_VERBOSITY_TUNER_ALL,

    MPIX_T_VERBOSITY_MPIDEV_BASIC,
    MPIX_T_VERBOSITY_MPIDEV_DETAIL,
    MPIX_T_VERBOSITY_MPIDEV_ALL
};

enum MPIR_T_bind_t {
    /* don't name-shift this if/when MPI_T_ is accepted, this is an MPICH2-only
     * extension */
    MPIX_T_BIND_INVALID = 0,

    /* arbitrarily shift values to aid debugging and reduce accidental errors */
    MPIX_T_BIND_NO_OBJECT = 9700,
    MPIX_T_BIND_MPI_COMM,
    MPIX_T_BIND_MPI_DATATYPE,
    MPIX_T_BIND_MPI_ERRHANDLER,
    MPIX_T_BIND_MPI_FILE,
    MPIX_T_BIND_MPI_GROUP,
    MPIX_T_BIND_MPI_OP,
    MPIX_T_BIND_MPI_REQUEST,
    MPIX_T_BIND_MPI_WIN,
    MPIX_T_BIND_MPI_MESSAGE,
    MPIX_T_BIND_MPI_INFO
};

enum MPIR_T_scope_t {
    /* don't name-shift this if/when MPI_T_ is accepted, this is an MPICH2-only
     * extension */
    MPIX_T_SCOPE_INVALID = 0,

    /* arbitrarily shift values to aid debugging and reduce accidental errors */
    MPIX_T_SCOPE_READONLY = 60439,
    MPIX_T_SCOPE_LOCAL,
    MPIX_T_SCOPE_GROUP,
    MPIX_T_SCOPE_GROUP_EQ,
    MPIX_T_SCOPE_ALL,
    MPIX_T_SCOPE_ALL_EQ
};

enum MPIR_T_pvar_class_t {
    /* don't name-shift this if/when MPI_T_ is accepted, this is an MPICH2-only
     * extension */
    MPIX_T_PVAR_CLASS_INVALID = 0,

    /* arbitrarily shift values to aid debugging and reduce accidental errors */
    MPIX_T_PVAR_CLASS_STATE = 240,
    MPIX_T_PVAR_CLASS_LEVEL,
    MPIX_T_PVAR_CLASS_SIZE,
    MPIX_T_PVAR_CLASS_PERCENTAGE,
    MPIX_T_PVAR_CLASS_HIGHWATERMARK,
    MPIX_T_PVAR_CLASS_LOWWATERMARK,
    MPIX_T_PVAR_CLASS_COUNTER,
    MPIX_T_PVAR_CLASS_AGGREGATE,
    MPIX_T_PVAR_CLASS_TIMER,
    MPIX_T_PVAR_CLASS_GENERIC
};

/* Handle conversion types/functions */

/* Programs that need to convert types used in MPICH should use these */
#define MPI_Comm_c2f(comm) (MPI_Fint)(comm)
#define MPI_Comm_f2c(comm) (MPI_Comm)(comm)
#define MPI_Type_c2f(datatype) (MPI_Fint)(datatype)
#define MPI_Type_f2c(datatype) (MPI_Datatype)(datatype)
#define MPI_Group_c2f(group) (MPI_Fint)(group)
#define MPI_Group_f2c(group) (MPI_Group)(group)
#define MPI_Info_c2f(info) (MPI_Fint)(info)
#define MPI_Info_f2c(info) (MPI_Info)(info)
#define MPI_Request_f2c(request) (MPI_Request)(request)
#define MPI_Request_c2f(request) (MPI_Fint)(request)
#define MPI_Op_c2f(op) (MPI_Fint)(op)
#define MPI_Op_f2c(op) (MPI_Op)(op)
#define MPI_Errhandler_c2f(errhandler) (MPI_Fint)(errhandler)
#define MPI_Errhandler_f2c(errhandler) (MPI_Errhandler)(errhandler)
#define MPI_Win_c2f(win)   (MPI_Fint)(win)
#define MPI_Win_f2c(win)   (MPI_Win)(win)
#define MPIX_Message_c2f(msg) ((MPI_Fint)(msg))
#define MPIX_Message_f2c(msg) ((MPIX_Message)(msg))

/* PMPI versions of the handle transfer functions.  See section 4.17 */
#define PMPI_Comm_c2f(comm) (MPI_Fint)(comm)
#define PMPI_Comm_f2c(comm) (MPI_Comm)(comm)
#define PMPI_Type_c2f(datatype) (MPI_Fint)(datatype)
#define PMPI_Type_f2c(datatype) (MPI_Datatype)(datatype)
#define PMPI_Group_c2f(group) (MPI_Fint)(group)
#define PMPI_Group_f2c(group) (MPI_Group)(group)
#define PMPI_Info_c2f(info) (MPI_Fint)(info)
#define PMPI_Info_f2c(info) (MPI_Info)(info)
#define PMPI_Request_f2c(request) (MPI_Request)(request)
#define PMPI_Request_c2f(request) (MPI_Fint)(request)
#define PMPI_Op_c2f(op) (MPI_Fint)(op)
#define PMPI_Op_f2c(op) (MPI_Op)(op)
#define PMPI_Errhandler_c2f(errhandler) (MPI_Fint)(errhandler)
#define PMPI_Errhandler_f2c(errhandler) (MPI_Errhandler)(errhandler)
#define PMPI_Win_c2f(win)   (MPI_Fint)(win)
#define PMPI_Win_f2c(win)   (MPI_Win)(win)
#define PMPIX_Message_c2f(msg) ((MPI_Fint)(msg))
#define PMPIX_Message_f2c(msg) ((MPIX_Message)(msg))

#define MPI_STATUS_IGNORE (MPI_Status *)1
#define MPI_STATUSES_IGNORE (MPI_Status *)1
#define MPI_ERRCODES_IGNORE (int *)0

/* See 4.12.5 for MPI_F_STATUS(ES)_IGNORE */
#define MPIU_DLL_SPEC
extern MPIU_DLL_SPEC MPI_Fint * MPI_F_STATUS_IGNORE;
extern MPIU_DLL_SPEC MPI_Fint * MPI_F_STATUSES_IGNORE;
/* The annotation MPIU_DLL_SPEC to the extern statements is used 
   as a hook for systems that require C extensions to correctly construct
   DLLs, and is defined as an empty string otherwise
 */

/* The MPI standard requires that the ARGV_NULL values be the same as
   NULL (see 5.3.2) */
#define MPI_ARGV_NULL (char **)0
#define MPI_ARGVS_NULL (char ***)0

/* For supported thread levels */
#define MPI_THREAD_SINGLE 0
#define MPI_THREAD_FUNNELED 1
#define MPI_THREAD_SERIALIZED 2
#define MPI_THREAD_MULTIPLE 3

/* Typedefs for generalized requests */
typedef int (MPI_Grequest_cancel_function)(void *, int); 
typedef int (MPI_Grequest_free_function)(void *); 
typedef int (MPI_Grequest_query_function)(void *, MPI_Status *); 

/* MPI's error classes */
#define MPI_SUCCESS          0      /* Successful return code */
/* Communication argument parameters */
#define MPI_ERR_BUFFER       1      /* Invalid buffer pointer */
#define MPI_ERR_COUNT        2      /* Invalid count argument */
#define MPI_ERR_TYPE         3      /* Invalid datatype argument */
#define MPI_ERR_TAG          4      /* Invalid tag argument */
#define MPI_ERR_COMM         5      /* Invalid communicator */
#define MPI_ERR_RANK         6      /* Invalid rank */
#define MPI_ERR_ROOT         7      /* Invalid root */
#define MPI_ERR_TRUNCATE    14      /* Message truncated on receive */

/* MPI Objects (other than COMM) */
#define MPI_ERR_GROUP        8      /* Invalid group */
#define MPI_ERR_OP           9      /* Invalid operation */
#define MPI_ERR_REQUEST     19      /* Invalid mpi_request handle */

/* Special topology argument parameters */
#define MPI_ERR_TOPOLOGY    10      /* Invalid topology */
#define MPI_ERR_DIMS        11      /* Invalid dimension argument */

/* All other arguments.  This is a class with many kinds */
#define MPI_ERR_ARG         12      /* Invalid argument */

/* Other errors that are not simply an invalid argument */
#define MPI_ERR_OTHER       15      /* Other error; use Error_string */

#define MPI_ERR_UNKNOWN     13      /* Unknown error */
#define MPI_ERR_INTERN      16      /* Internal error code    */

/* Multiple completion has two special error classes */
#define MPI_ERR_IN_STATUS   17      /* Look in status for error value */
#define MPI_ERR_PENDING     18      /* Pending request */

/* New MPI-2 Error classes */
#define MPI_ERR_FILE        27      /* */
#define MPI_ERR_ACCESS      20      /* */
#define MPI_ERR_AMODE       21      /* */
#define MPI_ERR_BAD_FILE    22      /* */
#define MPI_ERR_FILE_EXISTS 25      /* */
#define MPI_ERR_FILE_IN_USE 26      /* */
#define MPI_ERR_NO_SPACE    36      /* */
#define MPI_ERR_NO_SUCH_FILE 37     /* */
#define MPI_ERR_IO          32      /* */
#define MPI_ERR_READ_ONLY   40      /* */
#define MPI_ERR_CONVERSION  23      /* */
#define MPI_ERR_DUP_DATAREP 24      /* */
#define MPI_ERR_UNSUPPORTED_DATAREP   43  /* */

/* MPI_ERR_INFO is NOT defined in the MPI-2 standard.  I believe that
   this is an oversight */
#define MPI_ERR_INFO        28      /* */
#define MPI_ERR_INFO_KEY    29      /* */
#define MPI_ERR_INFO_VALUE  30      /* */
#define MPI_ERR_INFO_NOKEY  31      /* */

#define MPI_ERR_NAME        33      /* */
#define MPI_ERR_NO_MEM      34      /* Alloc_mem could not allocate memory */
#define MPI_ERR_NOT_SAME    35      /* */
#define MPI_ERR_PORT        38      /* */
#define MPI_ERR_QUOTA       39      /* */
#define MPI_ERR_SERVICE     41      /* */
#define MPI_ERR_SPAWN       42      /* */
#define MPI_ERR_UNSUPPORTED_OPERATION 44 /* */
#define MPI_ERR_WIN         45      /* */

#define MPI_ERR_BASE        46      /* */
#define MPI_ERR_LOCKTYPE    47      /* */
#define MPI_ERR_KEYVAL      48      /* Erroneous attribute key */
#define MPI_ERR_RMA_CONFLICT 49     /* */
#define MPI_ERR_RMA_SYNC    50      /* */ 
#define MPI_ERR_SIZE        51      /* */
#define MPI_ERR_DISP        52      /* */
#define MPI_ERR_ASSERT      53      /* */

#define MPIX_ERR_PROC_FAIL_STOP 54   /* Process failure */

#define MPIX_ERR_RMA_RANGE  55      /* */
#define MPIX_ERR_RMA_ATTACH 56      /* */
#define MPIX_ERR_RMA_SHARED 57      /* */
#define MPIX_ERR_RMA_WRONG_FLAVOR 58    /* */

#define MPI_ERR_LASTCODE    0x3fffffff  /* Last valid error code for a 
					   predefined error class */
/* WARNING: this is also defined in mpishared.h.  Update both locations */
#define MPICH_ERR_LAST_CLASS 58     /* It is also helpful to know the
				       last valid class */
/* End of MPI's error classes */

/* Function type defs */
typedef int (MPI_Datarep_conversion_function)(void *, MPI_Datatype, int, 
             void *, MPI_Offset, void *);
typedef int (MPI_Datarep_extent_function)(MPI_Datatype datatype, MPI_Aint *,
					  void *);
#define MPI_CONVERSION_FN_NULL ((MPI_Datarep_conversion_function *)0)

/* 
   For systems that may need to add additional definitions to support
   different declaration styles and options (e.g., different calling 
   conventions or DLL import/export controls).  
*/
/* --Insert Additional Definitions Here-- */

/* Include any device specific definitions */
/* ... no device specific definitions ... */

/*
 * Normally, we provide prototypes for all MPI routines.  In a few weird
 * cases, we need to suppress the prototypes.
 */
#ifndef MPICH_SUPPRESS_PROTOTYPES
/* We require that the C compiler support prototypes */
/* Begin Prototypes */
int MPI_Send(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Get_count(MPICH2_CONST MPI_Status *, MPI_Datatype, int *);
int MPI_Bsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Ssend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Rsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Buffer_attach( void*, int);
int MPI_Buffer_detach( void*, int *);
int MPI_Isend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Ibsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Issend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Irsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Wait(MPI_Request *, MPI_Status *);
int MPI_Test(MPI_Request *, int *, MPI_Status *);
int MPI_Request_free(MPI_Request *);
int MPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
int MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Waitall(int, MPI_Request *, MPI_Status *);
int MPI_Testall(int, MPI_Request *, int *, MPI_Status *);
int MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Iprobe(int, int, MPI_Comm, int *, MPI_Status *);
int MPI_Probe(int, int, MPI_Comm, MPI_Status *);
int MPI_Cancel(MPI_Request *);
int MPI_Test_cancelled(MPICH2_CONST MPI_Status *, int *);
int MPI_Send_init(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Bsend_init(MPICH2_CONST void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Ssend_init(MPICH2_CONST void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Rsend_init(MPICH2_CONST void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Recv_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Start(MPI_Request *);
int MPI_Startall(int, MPI_Request *);
int MPI_Sendrecv(MPICH2_CONST void *, int, MPI_Datatype,int, int, void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8);
int MPI_Sendrecv_replace(void*, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int MPI_Type_indexed(int, MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hindexed(int, MPICH2_CONST int *, MPICH2_CONST MPI_Aint *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_struct(int, MPICH2_CONST int *, MPICH2_CONST MPI_Aint *, MPICH2_CONST MPI_Datatype *, MPI_Datatype *);
int MPI_Address(MPICH2_CONST void*, MPI_Aint *);
/* We could add __attribute__((deprecated)) to routines like MPI_Type_extent */
int MPI_Type_extent(MPI_Datatype, MPI_Aint *);
/* See the 1.1 version of the Standard.  The standard made an 
   unfortunate choice here, however, it is the standard.  The size returned 
   by MPI_Type_size is specified as an int, not an MPI_Aint */
int MPI_Type_size(MPI_Datatype, int *);
/* MPI_Type_count was withdrawn in MPI 1.1 */
int MPI_Type_lb(MPI_Datatype, MPI_Aint *);
int MPI_Type_ub(MPI_Datatype, MPI_Aint *);
int MPI_Type_commit(MPI_Datatype *);
int MPI_Type_free(MPI_Datatype *);
int MPI_Get_elements(MPICH2_CONST MPI_Status *, MPI_Datatype, int *);
int MPI_Pack(MPICH2_CONST void*, int, MPI_Datatype, void *, int, int *,  MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Unpack(MPICH2_CONST void*, int, int *, void *, int, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPI_Pack_size(int, MPI_Datatype, MPI_Comm, int *);
int MPI_Barrier(MPI_Comm );
int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Gather(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPI_Gatherv(MPICH2_CONST void* , int, MPI_Datatype, void*, MPICH2_CONST int *,
                MPICH2_CONST int *, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int MPI_Scatter(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPI_Scatterv(MPICH2_CONST void* , MPICH2_CONST int *, MPICH2_CONST int *,  MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7);
int MPI_Allgather(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPI_Allgatherv(MPICH2_CONST void* , int, MPI_Datatype, void*, MPICH2_CONST int *,
                   MPICH2_CONST int *, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int MPI_Alltoall(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPI_Alltoallv(MPICH2_CONST void* , MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype,
                  void*, MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int MPI_Reduce(MPICH2_CONST void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPI_Op_create(MPI_User_function *, int, MPI_Op *);
int MPI_Op_free( MPI_Op *);
int MPI_Allreduce(MPICH2_CONST void*, void*, int, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPI_Reduce_scatter(MPICH2_CONST void* , void*, MPICH2_CONST int *, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPI_Scan(MPICH2_CONST void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPI_Group_size(MPI_Group, int *);
int MPI_Group_rank(MPI_Group, int *);
int MPI_Group_translate_ranks (MPI_Group, int, MPICH2_CONST int *, MPI_Group, int *);
int MPI_Group_compare(MPI_Group, MPI_Group, int *);
int MPI_Comm_group(MPI_Comm, MPI_Group *);
int MPI_Group_union(MPI_Group, MPI_Group, MPI_Group *);
int MPI_Group_intersection(MPI_Group, MPI_Group, MPI_Group *);
int MPI_Group_difference(MPI_Group, MPI_Group, MPI_Group *);
int MPI_Group_incl(MPI_Group, int, MPICH2_CONST int *, MPI_Group *);
int MPI_Group_excl(MPI_Group, int, MPICH2_CONST int *, MPI_Group *);
int MPI_Group_range_incl(MPI_Group, int, int [][3], MPI_Group *);
int MPI_Group_range_excl(MPI_Group, int, int [][3], MPI_Group *);
int MPI_Group_free(MPI_Group *);
int MPI_Comm_size(MPI_Comm, int *);
int MPI_Comm_rank(MPI_Comm, int *);
int MPI_Comm_compare(MPI_Comm, MPI_Comm, int *);
int MPI_Comm_dup(MPI_Comm, MPI_Comm *);
int MPI_Comm_create(MPI_Comm, MPI_Group, MPI_Comm *);
int MPI_Comm_split(MPI_Comm, int, int, MPI_Comm *);
int MPI_Comm_free(MPI_Comm *);
int MPI_Comm_test_inter(MPI_Comm, int *);
int MPI_Comm_remote_size(MPI_Comm, int *);
int MPI_Comm_remote_group(MPI_Comm, MPI_Group *);
int MPI_Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm * );
int MPI_Intercomm_merge(MPI_Comm, int, MPI_Comm *);
int MPI_Keyval_create(MPI_Copy_function *, MPI_Delete_function *, int *, void*);
int MPI_Keyval_free(int *);
int MPI_Attr_put(MPI_Comm, int, void*);
int MPI_Attr_get(MPI_Comm, int, void *, int *);
int MPI_Attr_delete(MPI_Comm, int);
int MPI_Topo_test(MPI_Comm, int *);
int MPI_Cart_create(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int, MPI_Comm *);
int MPI_Dims_create(int, int, int *);
int MPI_Graph_create(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int, MPI_Comm *);
int MPI_Graphdims_get(MPI_Comm, int *, int *);
int MPI_Graph_get(MPI_Comm, int, int, int *, int *);
int MPI_Cartdim_get(MPI_Comm, int *);
int MPI_Cart_get(MPI_Comm, int, int *, int *, int *);
int MPI_Cart_rank(MPI_Comm, MPICH2_CONST int *, int *);
int MPI_Cart_coords(MPI_Comm, int, int, int *);
int MPI_Graph_neighbors_count(MPI_Comm, int, int *);
int MPI_Graph_neighbors(MPI_Comm, int, int, int *);
int MPI_Cart_shift(MPI_Comm, int, int, int *, int *);
int MPI_Cart_sub(MPI_Comm, MPICH2_CONST int *, MPI_Comm *);
int MPI_Cart_map(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int *);
int MPI_Graph_map(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int *);
int MPI_Get_processor_name(char *, int *);
int MPI_Get_version(int *, int *);
int MPI_Errhandler_create(MPI_Handler_function *, MPI_Errhandler *);
int MPI_Errhandler_set(MPI_Comm, MPI_Errhandler);
int MPI_Errhandler_get(MPI_Comm, MPI_Errhandler *);
int MPI_Errhandler_free(MPI_Errhandler *);
int MPI_Error_string(int, char *, int *);
int MPI_Error_class(int, int *);
double MPI_Wtime(void);
double MPI_Wtick(void);
#ifndef MPI_Wtime
double PMPI_Wtime(void);
double PMPI_Wtick(void);
#endif
int MPI_Init(int *, char ***);
int MPI_Finalize(void);
int MPI_Initialized(int *);
int MPI_Abort(MPI_Comm, int);


/* Note that we may need to define a @PCONTROL_LIST@ depending on whether 
   stdargs are supported */
int MPI_Pcontrol(const int, ...);

int MPI_DUP_FN ( MPI_Comm, int, void *, void *, void *, int * );


/* MPI-2 functions */

/* Process Creation and Management */
int MPI_Close_port(MPICH2_CONST char *);
int MPI_Comm_accept(MPICH2_CONST char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int MPI_Comm_connect(MPICH2_CONST char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int MPI_Comm_disconnect(MPI_Comm *);
int MPI_Comm_get_parent(MPI_Comm *);
int MPI_Comm_join(int, MPI_Comm *);
int MPI_Comm_spawn(MPICH2_CONST char *, char *[], int, MPI_Info, int, MPI_Comm, MPI_Comm *,
                   int []);
int MPI_Comm_spawn_multiple(int, char *[], char **[], MPICH2_CONST int [], MPICH2_CONST MPI_Info [], int,
			    MPI_Comm, MPI_Comm *, int []); 
int MPI_Lookup_name(MPICH2_CONST char *, MPI_Info, char *);
int MPI_Open_port(MPI_Info, char *);
int MPI_Publish_name(MPICH2_CONST char *, MPI_Info, MPICH2_CONST char *);
int MPI_Unpublish_name(MPICH2_CONST char *, MPI_Info, MPICH2_CONST char *);

/* One-Sided Communications */
int MPI_Accumulate(MPICH2_CONST void *, int, MPI_Datatype, int, MPI_Aint, int,
                   MPI_Datatype,  MPI_Op, MPI_Win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
            MPI_Win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Put(MPICH2_CONST void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
            MPI_Win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPI_Win_complete(MPI_Win);
int MPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *);
int MPI_Win_fence(int, MPI_Win);
int MPI_Win_free(MPI_Win *);
int MPI_Win_get_group(MPI_Win, MPI_Group *);
int MPI_Win_lock(int, int, int, MPI_Win);
int MPI_Win_post(MPI_Group, int, MPI_Win);
int MPI_Win_start(MPI_Group, int, MPI_Win);
int MPI_Win_test(MPI_Win, int *);
int MPI_Win_unlock(int, MPI_Win);
int MPI_Win_wait(MPI_Win);

/* MPI-3 One-Sided Communication Routines */
int MPIX_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                      MPI_Comm comm, void *baseptr, MPI_Win *win);
int MPIX_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                             void *baseptr, MPI_Win *win);
int MPIX_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
int MPIX_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win);
int MPIX_Win_attach(MPI_Win win, void *base, MPI_Aint size);
int MPIX_Win_detach(MPI_Win win, const void *base);

int MPIX_Get_accumulate(const void *origin_addr, int origin_count,
                        MPI_Datatype origin_datatype, void *result_addr, int result_count,
                        MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                        int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
int MPIX_Fetch_and_op(const void *origin_addr, void *result_addr,
                      MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                      MPI_Op op, MPI_Win win);
int MPIX_Compare_and_swap(const void *origin_addr, const void *compare_addr,
                          void *result_addr, MPI_Datatype datatype, int target_rank,
                          MPI_Aint target_disp, MPI_Win win);
int MPIX_Rput(const void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
              int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request);
int MPIX_Rget(void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
              int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request);
int MPIX_Raccumulate(const void *origin_addr, int origin_count,
                     MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                     int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                     MPI_Request *request);
int MPIX_Rget_accumulate(const void *origin_addr, int origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, int result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                         MPI_Request *request);

int MPIX_Win_lock_all(int assert, MPI_Win win);
int MPIX_Win_unlock_all(MPI_Win win);
int MPIX_Win_flush(int rank, MPI_Win win);
int MPIX_Win_flush_all(MPI_Win win);
int MPIX_Win_flush_local(int rank, MPI_Win win);
int MPIX_Win_flush_local_all(MPI_Win win);
int MPIX_Win_sync(MPI_Win win);


/* Extended Collective Operations */
int MPI_Alltoallw(MPICH2_CONST void *, MPICH2_CONST int [], MPICH2_CONST int [],
                  MPICH2_CONST MPI_Datatype [], void *, MPICH2_CONST int [],
		  MPICH2_CONST int [], MPICH2_CONST MPI_Datatype [], MPI_Comm);
int MPI_Exscan(MPICH2_CONST void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
 
/* External Interfaces */
int MPI_Add_error_class(int *);
int MPI_Add_error_code(int, int *);
int MPI_Add_error_string(int, MPICH2_CONST char *);
int MPI_Comm_call_errhandler(MPI_Comm, int);
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *, 
                           MPI_Comm_delete_attr_function *, int *, void *);
int MPI_Comm_delete_attr(MPI_Comm, int);
int MPI_Comm_free_keyval(int *);
int MPI_Comm_get_attr(MPI_Comm, int, void *, int *);
int MPI_Comm_get_name(MPI_Comm, char *, int *);
int MPI_Comm_set_attr(MPI_Comm, int, void *);
int MPI_Comm_set_name(MPI_Comm, MPICH2_CONST char *);
int MPI_File_call_errhandler(MPI_File, int);
int MPI_Grequest_complete(MPI_Request);
int MPI_Grequest_start(MPI_Grequest_query_function *, 
                       MPI_Grequest_free_function *, 
                       MPI_Grequest_cancel_function *, void *, MPI_Request *);
int MPI_Init_thread(int *, char ***, int, int *);
int MPI_Is_thread_main(int *);
int MPI_Query_thread(int *);
int MPI_Status_set_cancelled(MPI_Status *, int);
int MPI_Status_set_elements(MPI_Status *, MPI_Datatype, int);
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *, 
                           MPI_Type_delete_attr_function *, int *, void *);
int MPI_Type_delete_attr(MPI_Datatype, int);
int MPI_Type_dup(MPI_Datatype, MPI_Datatype *);
int MPI_Type_free_keyval(int *);
int MPI_Type_get_attr(MPI_Datatype, int, void *, int *);
int MPI_Type_get_contents(MPI_Datatype, int, int, int, int [], MPI_Aint [], 
                          MPI_Datatype []);
int MPI_Type_get_envelope(MPI_Datatype, int *, int *, int *, int *);
int MPI_Type_get_name(MPI_Datatype, char *, int *);
int MPI_Type_set_attr(MPI_Datatype, int, void *);
int MPI_Type_set_name(MPI_Datatype, MPICH2_CONST char *);
int MPI_Type_match_size( int, int, MPI_Datatype *);
int MPI_Win_call_errhandler(MPI_Win, int);
int MPI_Win_create_keyval(MPI_Win_copy_attr_function *, 
                         MPI_Win_delete_attr_function *, int *, void *);
int MPI_Win_delete_attr(MPI_Win, int);
int MPI_Win_free_keyval(int *);
int MPI_Win_get_attr(MPI_Win, int, void *, int *);
int MPI_Win_get_name(MPI_Win, char *, int *);
int MPI_Win_set_attr(MPI_Win, int, void *);
int MPI_Win_set_name(MPI_Win, MPICH2_CONST char *);

/* Miscellany */
#ifdef FOO
MPI_Comm MPI_Comm_f2c(MPI_Fint);
MPI_Datatype MPI_Type_f2c(MPI_Fint);
MPI_File MPI_File_f2c(MPI_Fint);
MPI_Fint MPI_Comm_c2f(MPI_Comm);
MPI_Fint MPI_File_c2f(MPI_File);
MPI_Fint MPI_Group_c2f(MPI_Group);
MPI_Fint MPI_Info_c2f(MPI_Info);
MPI_Fint MPI_Op_c2f(MPI_Op);
MPI_Fint MPI_Request_c2f(MPI_Request);
MPI_Fint MPI_Type_c2f(MPI_Datatype);
MPI_Fint MPI_Win_c2f(MPI_Win);
MPI_Group MPI_Group_f2c(MPI_Fint);
MPI_Info MPI_Info_f2c(MPI_Fint);
MPI_Op MPI_Op_f2c(MPI_Fint);
MPI_Request MPI_Request_f2c(MPI_Fint);
MPI_Win MPI_Win_f2c(MPI_Fint);
#endif

int MPI_Alloc_mem(MPI_Aint, MPI_Info info, void *baseptr);
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function *, MPI_Errhandler *);
int MPI_Comm_get_errhandler(MPI_Comm, MPI_Errhandler *);
int MPI_Comm_set_errhandler(MPI_Comm, MPI_Errhandler);
int MPI_File_create_errhandler(MPI_File_errhandler_function *, MPI_Errhandler *);
int MPI_File_get_errhandler(MPI_File, MPI_Errhandler *);
int MPI_File_set_errhandler(MPI_File, MPI_Errhandler);
int MPI_Finalized(int *);
int MPI_Free_mem(void *);
int MPI_Get_address(MPICH2_CONST void *, MPI_Aint *);
int MPI_Info_create(MPI_Info *);
int MPI_Info_delete(MPI_Info, MPICH2_CONST char *);
int MPI_Info_dup(MPI_Info, MPI_Info *);
int MPI_Info_free(MPI_Info *info);
int MPI_Info_get(MPI_Info, MPICH2_CONST char *, int, char *, int *);
int MPI_Info_get_nkeys(MPI_Info, int *);
int MPI_Info_get_nthkey(MPI_Info, int, char *);
int MPI_Info_get_valuelen(MPI_Info, MPICH2_CONST char *, int *, int *);
int MPI_Info_set(MPI_Info, MPICH2_CONST char *, MPICH2_CONST char *);
int MPI_Pack_external(MPICH2_CONST char *, MPICH2_CONST void *, int, MPI_Datatype, void *,
                      MPI_Aint, MPI_Aint *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPI_Pack_external_size(MPICH2_CONST char *, int, MPI_Datatype, MPI_Aint *);
int MPI_Request_get_status(MPI_Request, int *, MPI_Status *);
int MPI_Status_c2f(MPICH2_CONST MPI_Status *, MPI_Fint *);
int MPI_Status_f2c(MPICH2_CONST MPI_Fint *, MPI_Status *);
int MPI_Type_create_darray(int, int, int, MPICH2_CONST int [], MPICH2_CONST int [],
                           MPICH2_CONST int [], MPICH2_CONST int [], int,
                           MPI_Datatype, MPI_Datatype *);
int MPI_Type_create_hindexed(int, MPICH2_CONST int [], MPICH2_CONST MPI_Aint [], MPI_Datatype,
                             MPI_Datatype *);
int MPI_Type_create_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int MPI_Type_create_indexed_block(int, int, MPICH2_CONST int [], MPI_Datatype,
                                  MPI_Datatype *);
int MPIX_Type_create_hindexed_block(int, int, const MPI_Aint [], MPI_Datatype, MPI_Datatype *);
int MPI_Type_create_resized(MPI_Datatype, MPI_Aint, MPI_Aint, MPI_Datatype *);
int MPI_Type_create_struct(int, MPICH2_CONST int [], MPICH2_CONST MPI_Aint [],
                           MPICH2_CONST MPI_Datatype [], MPI_Datatype *);
int MPI_Type_create_subarray(int, MPICH2_CONST int [], MPICH2_CONST int [], MPICH2_CONST int [],
                             int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_get_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int MPI_Type_get_true_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int MPI_Unpack_external(MPICH2_CONST char *, MPICH2_CONST void *, MPI_Aint, MPI_Aint *, void *,
                        int, MPI_Datatype) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7);
int MPI_Win_create_errhandler(MPI_Win_errhandler_function *, MPI_Errhandler *);
int MPI_Win_get_errhandler(MPI_Win, MPI_Errhandler *);
int MPI_Win_set_errhandler(MPI_Win, MPI_Errhandler);

/* Fortran 90-related functions.  These routines are available only if
   Fortran 90 support is enabled 
*/
int MPI_Type_create_f90_integer( int, MPI_Datatype * );
int MPI_Type_create_f90_real( int, int, MPI_Datatype * );
int MPI_Type_create_f90_complex( int, int, MPI_Datatype * );

/* MPI-2.2 functions */
int MPI_Reduce_local(MPICH2_CONST void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype, MPI_Op op) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPI_Op_commutative(MPI_Op op, int *commute);
int MPI_Reduce_scatter_block(MPICH2_CONST void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, MPICH2_CONST int [], MPICH2_CONST int [], int outdegree, MPICH2_CONST int [], MPICH2_CONST int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, MPICH2_CONST int [], MPICH2_CONST int [], MPICH2_CONST int [], MPICH2_CONST int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted);
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int [], int [], int maxoutdegree, int [], int []);

/* MPI-3 matched probe functionality, currently MPIX_ extensions */
int MPIX_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPIX_Message *message, MPI_Status *status);
int MPIX_Imrecv(void *buf, int count, MPI_Datatype datatype, MPIX_Message *message, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPIX_Mprobe(int source, int tag, MPI_Comm comm, MPIX_Message *message, MPI_Status *status);
int MPIX_Mrecv(void *buf, int count, MPI_Datatype datatype, MPIX_Message *message, MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);

/* MPI-3 nonblocking comm dup */
int MPIX_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request);

/* MPI-3 nonblocking collectives, currently MPIX_ extensions */
int MPIX_Ibarrier(MPI_Comm comm, MPI_Request *request);
int MPIX_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int MPIX_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int MPIX_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Iscatterv(const void *sendbuf, const int *sendcounts, const int *displs, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7);
int MPIX_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *displs, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int MPIX_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Ialltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int MPIX_Ialltoallw(const void *sendbuf, const int *sendcounts, const int *sdispls, const MPI_Datatype *sendtypes, void *recvbuf, const int *recvcounts, const int *rdispls, const MPI_Datatype *recvtypes, MPI_Comm comm, MPI_Request *request);
int MPIX_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPIX_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPIX_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int *recvcounts, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPIX_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPIX_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int MPIX_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);

/* MPI-3 neighborhood collectives */
int MPIX_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int MPIX_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int MPIX_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
int MPIX_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int MPIX_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int MPIX_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int MPIX_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm);

/* MPI-3 shared memory */
int MPIX_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);

/* MPI-3 noncollective communicator creation */
int MPIX_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm * newcomm);

/* MPI-3 FT */
int MPIX_Comm_group_failed(MPI_Comm, MPI_Group *);
int MPIX_Comm_remote_group_failed(MPI_Comm, MPI_Group *);
int MPIX_Comm_reenable_anysource(MPI_Comm, MPI_Group *);

/* RMA Mutexes extension declarations: */
struct mpixi_mutex_s;
typedef struct mpixi_mutex_s * MPIX_Mutex;

int MPIX_Mutex_create(int count, MPI_Comm comm, MPIX_Mutex *hdl);
int MPIX_Mutex_free(MPIX_Mutex *hdl);
int MPIX_Mutex_lock(MPIX_Mutex hdl, int mutex, int proc);
int MPIX_Mutex_unlock(MPIX_Mutex hdl, int mutex, int proc);

/* MPI-3 MPI_T interface, currently available as MPIX_T extensions */

/* The MPI_T routines are available only in C bindings - tell tools that they
   can skip these prototypes */
/* Begin Skip Prototypes */
int MPIX_T_init_thread(int required, int *provided);
int MPIX_T_finalize(void);
int MPIX_T_enum_get_info(MPIX_T_enum enumtype, int num, char *name, int *name_len);
int MPIX_T_enum_get_item(MPIX_T_enum enumtype, int num, int *value, char *name, int *name_len);
int MPIX_T_cvar_get_num(int *num_cvar);
int MPIX_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity, MPI_Datatype *datatype, MPIX_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *scope);
int MPIX_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPIX_T_cvar_handle *handle, int *count);
int MPIX_T_cvar_handle_free(MPIX_T_cvar_handle *handle);
int MPIX_T_cvar_read(MPIX_T_cvar_handle handle, void *buf);
int MPIX_T_cvar_write(MPIX_T_cvar_handle handle, void *buf);
int MPIX_T_pvar_get_num(int *num_pvar);
int MPIX_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class, MPI_Datatype *datatype, MPIX_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *readonly, int *continuous, int *atomic);
int MPIX_T_pvar_session_create(MPIX_T_pvar_session *session);
int MPIX_T_pvar_session_free(MPIX_T_pvar_session *session);
int MPIX_T_pvar_handle_alloc(MPIX_T_pvar_session session, int pvar_index, void *obj_handle, MPIX_T_pvar_handle *handle, int *count);
int MPIX_T_pvar_handle_free(MPIX_T_pvar_session session, MPIX_T_pvar_handle *handle);
int MPIX_T_pvar_start(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle);
int MPIX_T_pvar_stop(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle);
int MPIX_T_pvar_read(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle, void *buf);
int MPIX_T_pvar_write(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle, void *buf);
int MPIX_T_pvar_reset(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle);
int MPIX_T_pvar_readreset(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle, void *buf);
int MPIX_T_category_get_num(int *num_cat);
int MPIX_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len, int *num_controlvars, int *num_pvars, int *num_categories);
int MPIX_T_category_get_cvars(int cat_index, int len, int indices[]);
int MPIX_T_category_get_pvars(int cat_index[], int len, int indices[]);
int MPIX_T_category_get_categories(int cat_index, int len, int indices[]);
int MPIX_T_category_changed(int *stamp);
/* End Skip Prototypes */
/* End Prototypes */
#endif /* MPICH_SUPPRESS_PROTOTYPES */



/* Here are the bindings of the profiling routines */
#if !defined(MPI_BUILD_PROFILING)
int PMPI_Send(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Get_count(MPICH2_CONST MPI_Status *, MPI_Datatype, int *);
int PMPI_Bsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Ssend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Rsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Buffer_attach( void* buffer, int);
int PMPI_Buffer_detach( void* buffer, int *);
int PMPI_Isend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Ibsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Issend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Irsend(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Wait(MPI_Request *, MPI_Status *);
int PMPI_Test(MPI_Request *, int *, MPI_Status *);
int PMPI_Request_free(MPI_Request *);
int PMPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
int PMPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
int PMPI_Waitall(int, MPI_Request *, MPI_Status *);
int PMPI_Testall(int, MPI_Request *, int *, MPI_Status *);
int PMPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);
int PMPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);
int PMPI_Iprobe(int, int, MPI_Comm, int *, MPI_Status *);
int PMPI_Probe(int, int, MPI_Comm, MPI_Status *);
int PMPI_Cancel(MPI_Request *);
int PMPI_Test_cancelled(MPICH2_CONST MPI_Status *, int *);
int PMPI_Send_init(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Bsend_init(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Ssend_init(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Rsend_init(MPICH2_CONST void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Recv_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Start(MPI_Request *);
int PMPI_Startall(int, MPI_Request *);
int PMPI_Sendrecv(MPICH2_CONST void *, int, MPI_Datatype, int, int, void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(6,8);
int PMPI_Sendrecv_replace(void*, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_indexed(int, MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_hindexed(int, MPICH2_CONST int *, MPICH2_CONST MPI_Aint *, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_struct(int, MPICH2_CONST int *, MPICH2_CONST MPI_Aint *, MPICH2_CONST MPI_Datatype *, MPI_Datatype *);
int PMPI_Address(MPICH2_CONST void*, MPI_Aint *);
int PMPI_Type_extent(MPI_Datatype, MPI_Aint *);
int PMPI_Type_size(MPI_Datatype, int *);
int PMPI_Type_lb(MPI_Datatype, MPI_Aint *);
int PMPI_Type_ub(MPI_Datatype, MPI_Aint *);
int PMPI_Type_commit(MPI_Datatype *);
int PMPI_Type_free(MPI_Datatype *);
int PMPI_Get_elements(MPICH2_CONST MPI_Status *, MPI_Datatype, int *);
int PMPI_Pack(MPICH2_CONST void*, int, MPI_Datatype, void *, int, int *,  MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Unpack(MPICH2_CONST void*, int, int *, void *, int, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPI_Pack_size(int, MPI_Datatype, MPI_Comm, int *);
int PMPI_Barrier(MPI_Comm );
int PMPI_Bcast(void* buffer, int, MPI_Datatype, int, MPI_Comm );
int PMPI_Gather(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPI_Gatherv(MPICH2_CONST void* , int, MPI_Datatype, void*, MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int PMPI_Scatter(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPI_Scatterv(MPICH2_CONST void* , MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7);
int PMPI_Allgather(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPI_Allgatherv(MPICH2_CONST void* , int, MPI_Datatype, void*, MPICH2_CONST int *,
                    MPICH2_CONST int *, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int PMPI_Alltoall(MPICH2_CONST void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPI_Alltoallv(MPICH2_CONST void* , MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype,
                   void*, MPICH2_CONST int *, MPICH2_CONST int *, MPI_Datatype, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int PMPI_Reduce(MPICH2_CONST void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPI_Op_create(MPI_User_function *, int, MPI_Op *);
int PMPI_Op_free( MPI_Op *);
int PMPI_Allreduce(MPICH2_CONST void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPI_Reduce_scatter(MPICH2_CONST void* , void*, MPICH2_CONST int *, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPI_Scan(MPICH2_CONST void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPI_Group_size(MPI_Group, int *);
int PMPI_Group_rank(MPI_Group, int *);
int PMPI_Group_translate_ranks(MPI_Group, int, MPICH2_CONST int *, MPI_Group, int *);
int PMPI_Group_compare(MPI_Group, MPI_Group, int *);
int PMPIX_Comm_group_failed(MPI_Comm, MPI_Group *);
int PMPIX_Comm_reenable_anysource(MPI_Comm, MPI_Group *);
int PMPIX_Comm_remote_group_failed(MPI_Comm, MPI_Group *);
int PMPIX_Mutex_create(int count, MPI_Comm comm, MPIX_Mutex *hdl);
int PMPIX_Mutex_free(MPIX_Mutex *hdl);
int PMPIX_Mutex_lock(MPIX_Mutex hdl, int mutex, int proc);
int PMPIX_Mutex_unlock(MPIX_Mutex hdl, int mutex, int proc);
int PMPI_Comm_group(MPI_Comm, MPI_Group *);
int PMPI_Group_union(MPI_Group, MPI_Group, MPI_Group *);
int PMPI_Group_intersection(MPI_Group, MPI_Group, MPI_Group *);
int PMPI_Group_difference(MPI_Group, MPI_Group, MPI_Group *);
int PMPI_Group_incl(MPI_Group, int, MPICH2_CONST int *, MPI_Group *);
int PMPI_Group_excl(MPI_Group, int, MPICH2_CONST int *, MPI_Group *);
int PMPI_Group_range_incl(MPI_Group, int, int [][3], MPI_Group *);
int PMPI_Group_range_excl(MPI_Group, int, int [][3], MPI_Group *);
int PMPI_Group_free(MPI_Group *);
int PMPI_Comm_size(MPI_Comm, int *);
int PMPI_Comm_rank(MPI_Comm, int *);
int PMPI_Comm_compare(MPI_Comm, MPI_Comm, int *);
int PMPI_Comm_dup(MPI_Comm, MPI_Comm *);
int PMPI_Comm_create(MPI_Comm, MPI_Group, MPI_Comm *);
int PMPI_Comm_split(MPI_Comm, int, int, MPI_Comm *);
int PMPI_Comm_free(MPI_Comm *);
int PMPI_Comm_test_inter(MPI_Comm, int *);
int PMPI_Comm_remote_size(MPI_Comm, int *);
int PMPI_Comm_remote_group(MPI_Comm, MPI_Group *);
int PMPI_Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm *);
int PMPI_Intercomm_merge(MPI_Comm, int, MPI_Comm *);
int PMPI_Keyval_create(MPI_Copy_function *, MPI_Delete_function *, int *, void*);
int PMPI_Keyval_free(int *);
int PMPI_Attr_put(MPI_Comm, int, void*);
int PMPI_Attr_get(MPI_Comm, int, void *, int *);
int PMPI_Attr_delete(MPI_Comm, int);
int PMPI_Topo_test(MPI_Comm, int *);
int PMPI_Cart_create(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int, MPI_Comm *);
int PMPI_Dims_create(int, int, int *);
int PMPI_Graph_create(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int, MPI_Comm *);
int PMPI_Graphdims_get(MPI_Comm, int *, int *);
int PMPI_Graph_get(MPI_Comm, int, int, int *, int *);
int PMPI_Cartdim_get(MPI_Comm, int *);
int PMPI_Cart_get(MPI_Comm, int, int *, int *, int *);
int PMPI_Cart_rank(MPI_Comm, MPICH2_CONST int *, int *);
int PMPI_Cart_coords(MPI_Comm, int, int, int *);
int PMPI_Graph_neighbors_count(MPI_Comm, int, int *);
int PMPI_Graph_neighbors(MPI_Comm, int, int, int *);
int PMPI_Cart_shift(MPI_Comm, int, int, int *, int *);
int PMPI_Cart_sub(MPI_Comm, MPICH2_CONST int *, MPI_Comm *);
int PMPI_Cart_map(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int *);
int PMPI_Graph_map(MPI_Comm, int, MPICH2_CONST int *, MPICH2_CONST int *, int *);
int PMPI_Get_processor_name(char *, int *);
int PMPI_Get_version(int *, int *);
int PMPI_Errhandler_create(MPI_Handler_function *, MPI_Errhandler *);
int PMPI_Errhandler_set(MPI_Comm, MPI_Errhandler);
int PMPI_Errhandler_get(MPI_Comm, MPI_Errhandler *);
int PMPI_Errhandler_free(MPI_Errhandler *);
int PMPI_Error_string(int, char *, int *);
int PMPI_Error_class(int, int *);

/* Wtime done above */
int PMPI_Init(int *, char ***);
int PMPI_Finalize(void);
int PMPI_Initialized(int *);
int PMPI_Abort(MPI_Comm, int);

int PMPI_Pcontrol(const int, ...);

/* MPI-2 functions */

/* Process Creation and Management */
int PMPI_Close_port(MPICH2_CONST char *);
int PMPI_Comm_accept(MPICH2_CONST char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int PMPI_Comm_connect(MPICH2_CONST char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int PMPI_Comm_disconnect(MPI_Comm *);
int PMPI_Comm_get_parent(MPI_Comm *);
int PMPI_Comm_join(int, MPI_Comm *);
int PMPI_Comm_spawn(MPICH2_CONST char *, char *[], int, MPI_Info, int, MPI_Comm, MPI_Comm *,
                   int []);
int PMPI_Comm_spawn_multiple(int, char *[], char **[], MPICH2_CONST int [], MPICH2_CONST MPI_Info [], int,
			    MPI_Comm, MPI_Comm *, int []); 
int PMPI_Lookup_name(MPICH2_CONST char *, MPI_Info, char *);
int PMPI_Open_port(MPI_Info, char *);
int PMPI_Publish_name(MPICH2_CONST char *, MPI_Info, MPICH2_CONST char *);
int PMPI_Unpublish_name(MPICH2_CONST char *, MPI_Info, MPICH2_CONST char *);

/* One-Sided Communications */
int PMPI_Accumulate(MPICH2_CONST void *, int, MPI_Datatype, int, MPI_Aint, int,
		   MPI_Datatype,  MPI_Op, MPI_Win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype, 
	    MPI_Win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Put(MPICH2_CONST void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
	    MPI_Win) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPI_Win_complete(MPI_Win);
int PMPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *);
int PMPI_Win_fence(int, MPI_Win);
int PMPI_Win_free(MPI_Win *);
int PMPI_Win_get_group(MPI_Win, MPI_Group *);
int PMPI_Win_lock(int, int, int, MPI_Win);
int PMPI_Win_post(MPI_Group, int, MPI_Win);
int PMPI_Win_start(MPI_Group, int, MPI_Win);
int PMPI_Win_test(MPI_Win, int *);
int PMPI_Win_unlock(int, MPI_Win);
int PMPI_Win_wait(MPI_Win);
 
/* MPI-3 One-Sided Communication Routines */
int PMPIX_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                       MPI_Comm comm, void *baseptr, MPI_Win *win);
int PMPIX_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                              void *baseptr, MPI_Win *win);
int PMPIX_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
int PMPIX_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win);
int PMPIX_Win_attach(MPI_Win win, void *base, MPI_Aint size);
int PMPIX_Win_detach(MPI_Win win, const void *base);

int PMPIX_Get_accumulate(const void *origin_addr, int origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, int result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
int PMPIX_Fetch_and_op(const void *origin_addr, void *result_addr,
                       MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                       MPI_Op op, MPI_Win win);
int PMPIX_Compare_and_swap(const void *origin_addr, const void *compare_addr,
                           void *result_addr, MPI_Datatype datatype, int target_rank,
                           MPI_Aint target_disp, MPI_Win win);
int PMPIX_Rput(const void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               int target_count, MPI_Datatype target_datatype, MPI_Win win,
               MPI_Request *request);
int PMPIX_Rget(void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               int target_count, MPI_Datatype target_datatype, MPI_Win win,
               MPI_Request *request);
int PMPIX_Raccumulate(const void *origin_addr, int origin_count,
                      MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                      int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                      MPI_Request *request);
int PMPIX_Rget_accumulate(const void *origin_addr, int origin_count,
                          MPI_Datatype origin_datatype, void *result_addr, int result_count,
                          MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                          int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                          MPI_Request *request);

int PMPIX_Win_lock_all(int assert, MPI_Win win);
int PMPIX_Win_unlock_all(MPI_Win win);
int PMPIX_Win_flush(int rank, MPI_Win win);
int PMPIX_Win_flush_all(MPI_Win win);
int PMPIX_Win_flush_local(int rank, MPI_Win win);
int PMPIX_Win_flush_local_all(MPI_Win win);
int PMPIX_Win_sync(MPI_Win win);

/* Extended Collective Operations */
int PMPI_Alltoallw(MPICH2_CONST void *, MPICH2_CONST int [], MPICH2_CONST int [],
                   MPICH2_CONST MPI_Datatype [], void *, MPICH2_CONST int [],
                   MPICH2_CONST int [], MPICH2_CONST MPI_Datatype [], MPI_Comm);
int PMPI_Exscan(MPICH2_CONST void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
 
/* External Interfaces */
int PMPI_Add_error_class(int *);
int PMPI_Add_error_code(int, int *);
int PMPI_Add_error_string(int, MPICH2_CONST char *);
int PMPI_Comm_call_errhandler(MPI_Comm, int);
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function *, 
                           MPI_Comm_delete_attr_function *, int *, void *);
int PMPI_Comm_delete_attr(MPI_Comm, int);
int PMPI_Comm_free_keyval(int *);
int PMPI_Comm_get_attr(MPI_Comm, int, void *, int *);
int PMPI_Comm_get_name(MPI_Comm, char *, int *);
int PMPI_Comm_set_attr(MPI_Comm, int, void *);
int PMPI_Comm_set_name(MPI_Comm, MPICH2_CONST char *);
int PMPI_File_call_errhandler(MPI_File, int);
int PMPI_Grequest_complete(MPI_Request);
int PMPI_Grequest_start(MPI_Grequest_query_function *, 
                       MPI_Grequest_free_function *, 
                       MPI_Grequest_cancel_function *, void *, MPI_Request *);
int PMPI_Init_thread(int *, char ***, int, int *);
int PMPI_Is_thread_main(int *);
int PMPI_Query_thread(int *);
int PMPI_Status_set_cancelled(MPI_Status *, int);
int PMPI_Status_set_elements(MPI_Status *, MPI_Datatype, int);
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function *, 
                           MPI_Type_delete_attr_function *, int *, void *);
int PMPI_Type_delete_attr(MPI_Datatype, int);
int PMPI_Type_dup(MPI_Datatype, MPI_Datatype *);
int PMPI_Type_free_keyval(int *);
int PMPI_Type_get_attr(MPI_Datatype, int, void *, int *);
int PMPI_Type_get_contents(MPI_Datatype, int, int, int, int [], MPI_Aint [], 
                          MPI_Datatype []);
int PMPI_Type_get_envelope(MPI_Datatype, int *, int *, int *, int *);
int PMPI_Type_get_name(MPI_Datatype, char *, int *);
int PMPI_Type_set_attr(MPI_Datatype, int, void *);
int PMPI_Type_set_name(MPI_Datatype, MPICH2_CONST char *);
int PMPI_Type_match_size( int, int, MPI_Datatype *);
int PMPI_Win_call_errhandler(MPI_Win, int);
int PMPI_Win_create_keyval(MPI_Win_copy_attr_function *, 
                         MPI_Win_delete_attr_function *, int *, void *);
int PMPI_Win_delete_attr(MPI_Win, int);
int PMPI_Win_free_keyval(int *);
int PMPI_Win_get_attr(MPI_Win, int, void *, int *);
int PMPI_Win_get_name(MPI_Win, char *, int *);
int PMPI_Win_set_attr(MPI_Win, int, void *);
int PMPI_Win_set_name(MPI_Win, MPICH2_CONST char *);

/* Fortran 90-related functions.  These routines are available only if
   Fortran 90 support is enabled 
*/
int PMPI_Type_create_f90_integer( int, MPI_Datatype * );
int PMPI_Type_create_f90_real( int, int, MPI_Datatype * );
int PMPI_Type_create_f90_complex( int, int, MPI_Datatype * );

/* Miscellany */
int PMPI_Alloc_mem(MPI_Aint, MPI_Info info, void *baseptr);
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_function *, MPI_Errhandler *);
int PMPI_Comm_get_errhandler(MPI_Comm, MPI_Errhandler *);
int PMPI_Comm_set_errhandler(MPI_Comm, MPI_Errhandler);
int PMPI_File_create_errhandler(MPI_File_errhandler_function *, MPI_Errhandler *);
int PMPI_File_get_errhandler(MPI_File, MPI_Errhandler *);
int PMPI_File_set_errhandler(MPI_File, MPI_Errhandler);
int PMPI_Finalized(int *);
int PMPI_Free_mem(void *);
int PMPI_Get_address(MPICH2_CONST void *, MPI_Aint *);
int PMPI_Info_create(MPI_Info *);
int PMPI_Info_delete(MPI_Info, MPICH2_CONST char *);
int PMPI_Info_dup(MPI_Info, MPI_Info *);
int PMPI_Info_free(MPI_Info *info);
int PMPI_Info_get(MPI_Info, MPICH2_CONST char *, int, char *, int *);
int PMPI_Info_get_nkeys(MPI_Info, int *);
int PMPI_Info_get_nthkey(MPI_Info, int, char *);
int PMPI_Info_get_valuelen(MPI_Info, MPICH2_CONST char *, int *, int *);
int PMPI_Info_set(MPI_Info, MPICH2_CONST char *, MPICH2_CONST char *);
int PMPI_Pack_external(MPICH2_CONST char *, MPICH2_CONST void *, int, MPI_Datatype, void *,
                       MPI_Aint, MPI_Aint *) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPI_Pack_external_size(MPICH2_CONST char *, int, MPI_Datatype, MPI_Aint *);
int PMPI_Request_get_status(MPI_Request, int *, MPI_Status *);
int PMPI_Status_c2f(MPICH2_CONST MPI_Status *, MPI_Fint *);
int PMPI_Status_f2c(MPICH2_CONST MPI_Fint *, MPI_Status *);
int PMPI_Type_create_darray(int, int, int, MPICH2_CONST int [], MPICH2_CONST int [],
                            MPICH2_CONST int [], MPICH2_CONST int [], int,
                           MPI_Datatype, MPI_Datatype *);
int PMPI_Type_create_hindexed(int, MPICH2_CONST int [], MPICH2_CONST MPI_Aint [], MPI_Datatype,
                             MPI_Datatype *);
int PMPI_Type_create_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_create_indexed_block(int, int, MPICH2_CONST int [], MPI_Datatype,
                                  MPI_Datatype *);
int PMPIX_Type_create_hindexed_block(int, int, const MPI_Aint [], MPI_Datatype, MPI_Datatype *);
int PMPI_Type_create_resized(MPI_Datatype, MPI_Aint, MPI_Aint, MPI_Datatype *);
int PMPI_Type_create_struct(int, MPICH2_CONST int [], MPICH2_CONST MPI_Aint [],
                            MPICH2_CONST MPI_Datatype [], MPI_Datatype *);
int PMPI_Type_create_subarray(int, MPICH2_CONST int [], MPICH2_CONST int [], MPICH2_CONST int [],
                              int, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_get_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int PMPI_Type_get_true_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int PMPI_Unpack_external(MPICH2_CONST char *, MPICH2_CONST void *, MPI_Aint, MPI_Aint *, void *,
                         int, MPI_Datatype) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7);
int PMPI_Win_create_errhandler(MPI_Win_errhandler_function *, MPI_Errhandler *);
int PMPI_Win_get_errhandler(MPI_Win, MPI_Errhandler *);
int PMPI_Win_set_errhandler(MPI_Win, MPI_Errhandler);
int PMPI_Reduce_local(MPICH2_CONST void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype, MPI_Op op) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPI_Op_commutative(MPI_Op op, int *commute);
int PMPI_Reduce_scatter_block(MPICH2_CONST void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, MPICH2_CONST int [], MPICH2_CONST int [], int outdegree, MPICH2_CONST int [], MPICH2_CONST int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int PMPI_Dist_graph_create(MPI_Comm comm_old, int n, MPICH2_CONST int [], MPICH2_CONST int [], MPICH2_CONST int [], MPICH2_CONST int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int PMPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted);
int PMPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int [], int [], int maxoutdegree, int [], int []);
int PMPIX_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request);
int PMPIX_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPIX_Message *message, MPI_Status *status);
int PMPIX_Imrecv(void *buf, int count, MPI_Datatype datatype, MPIX_Message *message, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPIX_Mprobe(int source, int tag, MPI_Comm comm, MPIX_Message *message, MPI_Status *status);
int PMPIX_Mrecv(void *buf, int count, MPI_Datatype datatype, MPIX_Message *message, MPI_Status *status) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPIX_Ibarrier(MPI_Comm comm, MPI_Request *request);
int PMPIX_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3);
int PMPIX_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int PMPIX_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Iscatterv(const void *sendbuf, const int *sendcounts, const int *displs, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,7);
int PMPIX_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *displs, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int PMPIX_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Ialltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int PMPIX_Ialltoallw(const void *sendbuf, const int *sendcounts, const int *sdispls, const MPI_Datatype *sendtypes, void *recvbuf, const int *recvcounts, const int *rdispls, const MPI_Datatype *recvtypes, MPI_Comm comm, MPI_Request *request);
int PMPIX_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPIX_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPIX_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int *recvcounts, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPIX_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPIX_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPIX_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(2,4);
int PMPIX_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int PMPIX_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int PMPIX_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
int PMPIX_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,7);
int PMPIX_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,3) MPICH_ATTR_POINTER_WITH_TYPE_TAG(4,6);
int PMPIX_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm) MPICH_ATTR_POINTER_WITH_TYPE_TAG(1,4) MPICH_ATTR_POINTER_WITH_TYPE_TAG(5,8);
int PMPIX_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm);
int PMPIX_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);
int PMPIX_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm * newcomm);
int PMPIX_T_init_thread(int required, int *provided);
int PMPIX_T_finalize(void);
int PMPIX_T_enum_get_info(MPIX_T_enum enumtype, int num, char *name, int *name_len);
int PMPIX_T_enum_get_item(MPIX_T_enum enumtype, int num, int *value, char *name, int *name_len);
int PMPIX_T_cvar_get_num(int *num_cvar);
int PMPIX_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity, MPI_Datatype *datatype, MPIX_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *scope);
int PMPIX_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPIX_T_cvar_handle *handle, int *count);
int PMPIX_T_cvar_handle_free(MPIX_T_cvar_handle *handle);
int PMPIX_T_cvar_read(MPIX_T_cvar_handle handle, void *buf);
int PMPIX_T_cvar_write(MPIX_T_cvar_handle handle, void *buf);
int PMPIX_T_pvar_get_num(int *num_pvar);
int PMPIX_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class, MPI_Datatype *datatype, MPIX_T_enum *enumtype, char *desc, int *desc_len, int *bind, int *readonly, int *continuous, int *atomic);
int PMPIX_T_pvar_session_create(MPIX_T_pvar_session *session);
int PMPIX_T_pvar_session_free(MPIX_T_pvar_session *session);
int PMPIX_T_pvar_handle_alloc(MPIX_T_pvar_session session, int pvar_index, void *obj_handle, MPIX_T_pvar_handle *handle, int *count);
int PMPIX_T_pvar_handle_free(MPIX_T_pvar_session session, MPIX_T_pvar_handle *handle);
int PMPIX_T_pvar_start(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle);
int PMPIX_T_pvar_stop(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle);
int PMPIX_T_pvar_read(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle, void *buf);
int PMPIX_T_pvar_write(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle, void *buf);
int PMPIX_T_pvar_reset(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle);
int PMPIX_T_pvar_readreset(MPIX_T_pvar_session session, MPIX_T_pvar_handle handle, void *buf);
int PMPIX_T_category_get_num(int *num_cat);
int PMPIX_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len, int *num_controlvars, int *num_pvars, int *num_categories);
int PMPIX_T_category_get_cvars(int cat_index, int len, int indices[]);
int PMPIX_T_category_get_pvars(int cat_index[], int len, int indices[]);
int PMPIX_T_category_get_categories(int cat_index, int len, int indices[]);
int PMPIX_T_category_changed(int *stamp);
#endif  /* MPI_BUILD_PROFILING */
/* End of MPI bindings */

/* feature advertisement */
#define MPIIMPL_ADVERTISES_FEATURES 1
#define MPIIMPL_HAVE_MPI_INFO 1                                                 
#define MPIIMPL_HAVE_MPI_COMBINER_DARRAY 1                                      
#define MPIIMPL_HAVE_MPI_TYPE_CREATE_DARRAY 1
#define MPIIMPL_HAVE_MPI_COMBINER_SUBARRAY 1                                    
#define MPIIMPL_HAVE_MPI_TYPE_CREATE_DARRAY 1
#define MPIIMPL_HAVE_MPI_COMBINER_DUP 1                                         
#define MPIIMPL_HAVE_MPI_GREQUEST 1      
#define MPIIMPL_HAVE_STATUS_SET_BYTES 1
#define MPIIMPL_HAVE_STATUS_SET_INFO 1

#include "mpio.h"

#if defined(__cplusplus)
}
/* Add the C++ bindings */
/* 
   If MPICH_SKIP_MPICXX is defined, the mpicxx.h file will *not* be included.
   This is necessary, for example, when building the C++ interfaces.  It
   can also be used when you want to use a C++ compiler to compile C code,
   and do not want to load the C++ bindings.  These definitions can
   be made by the C++ compilation script
 */
#if !defined(MPICH_SKIP_MPICXX)
/* mpicxx.h contains the MPI C++ binding.  In the mpi.h.in file, this 
   include is in an autoconf variable in case the compiler is a C++ 
   compiler but MPI was built without the C++ bindings */
#include "mpicxx.h"
#endif 
#endif

/* BEGIN: non-standard but public extensions to MPI */
/* Generalized requests extensions as proposed in "Extending the MPI-2
 * Generalized Request Interface" */
typedef int (MPIX_Grequest_poll_function)(void *, MPI_Status *);
typedef int (MPIX_Grequest_wait_function)(int, void **, double, MPI_Status *);

typedef int MPIX_Grequest_class;
int MPIX_Grequest_class_create(MPI_Grequest_query_function *, 
                       MPI_Grequest_free_function *, 
                       MPI_Grequest_cancel_function *, 
		       MPIX_Grequest_poll_function *,
		       MPIX_Grequest_wait_function *, 
		       MPIX_Grequest_class *);

int MPIX_Grequest_class_allocate(MPIX_Grequest_class,
		       void *,
		       MPI_Request *);

int MPIX_Grequest_start(MPI_Grequest_query_function *, 
                       MPI_Grequest_free_function *, 
                       MPI_Grequest_cancel_function *, 
		       MPIX_Grequest_poll_function *,
		       MPIX_Grequest_wait_function *, void *, MPI_Request *);
#if !defined(MPI_BUILD_PROFILING)
int PMPIX_Grequest_class_create(MPI_Grequest_query_function *, 
                       MPI_Grequest_free_function *, 
                       MPI_Grequest_cancel_function *, 
		       MPIX_Grequest_poll_function *,
		       MPIX_Grequest_wait_function *, 
		       MPIX_Grequest_class *);

int PMPIX_Grequest_class_allocate(MPIX_Grequest_class,
		       void *,
		       MPI_Request *);
int PMPIX_Grequest_start(MPI_Grequest_query_function *, 
                       MPI_Grequest_free_function *, 
                       MPI_Grequest_cancel_function *, 
		       MPIX_Grequest_poll_function *,
		       MPIX_Grequest_wait_function *, void *, MPI_Request *);
#endif

/* END: non-standard but public extensions to MPI */

#endif
