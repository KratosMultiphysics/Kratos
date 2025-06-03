#include "gidpost.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>    // strdup
#ifndef WIN32
#include <strings.h>   // strcasecmp
#else // WIN32
#define strcasecmp  _stricmp
#endif // WIN32
#include <ctype.h>     // tolower

int OpenExampleFiles( GiD_FILE *dst_fdm, GiD_FILE *dst_fdr, char *test_filename, size_t buf_size, 
                    const char *prefix, const char *format) {
  int fail = 0;
  GiD_FILE fdm = 0, fdr = 0;
  if ( !strcasecmp( format, "ascii" ) ) {
    // GiD_PostSetFormatReal( "%.8g" );
    snprintf( test_filename, buf_size, "%s.post.msh", prefix );
    fdm = GiD_fOpenPostMeshFile( test_filename, GiD_PostAscii );
    fail = ( fdm == 0 ) ? 1 : 0;
    snprintf( test_filename, buf_size, "%s.post.res", prefix );
    // fdr = GiD_fOpenPostResultFile( buf, GiD_PostAsciiZipped );
    fdr = GiD_fOpenPostResultFile( test_filename, GiD_PostAscii );
    fail += ( fdr == 0 ) ? 1 : 0;
    printf( "Creating ASCII gid post files '%s.post.msh' and '%s.post.res'.\n", prefix, prefix );
  } else if ( !strcasecmp( format, "bin" ) ) {
    snprintf( test_filename, buf_size, "%s.post.bin", prefix );
    fdm = fdr = GiD_fOpenPostResultFile( test_filename, GiD_PostBinary );
    fail = ( fdm == 0 ) ? 1 : 0;
    printf( "Creating BINary gid post file '%s'.\n", test_filename );
  } else if ( !strcasecmp( format, "hdf5" ) ) {
    snprintf( test_filename, buf_size, "%s.post.h5", prefix );
    printf( "Creating HDF5 gid post file '%s'.\n", test_filename );
    fdm = fdr = GiD_fOpenPostResultFile( test_filename, GiD_PostHDF5 );
    fail = ( fdm == 0 ) ? 1 : 0;
  } else {
    printf( "Unkown format '%s'.\n", format );
    fail = 1;
  }
  if ( !fail) {
    *dst_fdm = fdm;
    *dst_fdr = fdr;
  }
  return fail;
}

void print_help_and_exit( const char *full_prog_name) {
  // get only the name of the executable
  const char *prog_name = &full_prog_name[ strlen( full_prog_name ) - 1 ];
  for ( ; prog_name > full_prog_name; prog_name-- ) {
    if ( ( *prog_name == '/' ) || ( *prog_name == '\\' ) ) {
      break;
    }
  }
  prog_name++;
  printf( "Usage: %s [ -h] [ -f ascii|bin|hdf5] filename_prefix\n", prog_name);
  printf( "    will create filename_prefix.post.{ msh,res | bin | h5} depending on the choosen format.\n");
  printf( "    Use: %s -- filename_prefix          if filename_prefix begins with a '-'\n", prog_name);
  printf( "            default filename_prefix is 'ex1_sphere_fd'\n" );
  exit( 0);
}

int main( int argc, char *argv[]) {
#define BUF_SIZE 1024
  char buf[ BUF_SIZE ];
  char test_filename[ BUF_SIZE ];
  char *prefix = NULL;
  char *format = NULL;
  int skip_options = 0;
  buf[ 0] = '\0';
  test_filename[ 0] = '\0';

  for ( int ia = 1; ia < argc; ia++) {
    if ( !skip_options && ( argv[ ia][ 0] == '-')) {
      char opt = ( char)tolower( argv[ ia][ 1]);
      if ( ( opt == 'h') || !strcasecmp( argv[ ia], "--help") ) {
        print_help_and_exit( argv[ 0 ] );
      } else if ( opt == 'f' ) {
        ia++;
        if ( ia < argc ) {
          format = strdup( argv[ ia ] );
        } else {
          printf( "Missing arguments.\n" );
          print_help_and_exit( argv[ 0 ] );
        }
      } else if ( opt == '-' ) {
        skip_options = 1;
      } else {
        printf( "Unknown option '%s'.\n", argv[ ia]);
        print_help_and_exit( argv[ 0]);
      }
    } else {
      prefix = strdup( argv[ ia]);
      break;
    }
  }

  // default values
  if ( !prefix )
    prefix = strdup( "ex1_sphere_fd" );
  if ( !format )
    format = strdup( "ascii" );

  printf( "version = %s\n", GiD_PostGetVersion());

  GiD_PostInit();

  GiD_FILE fdm = 0, fdr = 0;
  int fail = OpenExampleFiles( &fdm, &fdr, test_filename, BUF_SIZE, prefix, format);

  if ( fail ) {
    printf( "Error opening file '%s' with format '%s'.\n", test_filename, format );
    print_help_and_exit( argv[ 0]);
    exit( 1 );
  }


  GiD_fBeginMeshColor( fdm, "esferitas", GiD_3D, GiD_Sphere, 1, 0.7, 0.7, 0.7 );
  {
    GiD_fBeginCoordinates( fdm );
    {
      GiD_fWriteCoordinates( fdm, 1, 0, 0, 0 );
      GiD_fWriteCoordinates( fdm, 2, 100, 0, 0 );
    }
    GiD_fEndCoordinates( fdm );

    GiD_fBeginElements( fdm );
    {
      GiD_fWriteSphere( fdm, 1, 1, 1.0 );
      GiD_fWriteSphere( fdm, 2, 2, 10.0 );
    }
    GiD_fEndElements( fdm );
  }
  GiD_fEndMesh( fdm );

  GiD_fBeginMeshColor( fdm, "circulitos", GiD_3D, GiD_Circle, 1, 0.7, 0.7, 0.0 );
  {
    GiD_fBeginCoordinates( fdm );
    {
      GiD_fWriteCoordinates( fdm, 3, 30, 0, 0 );
      GiD_fWriteCoordinates( fdm, 4, 70, 0, 0 );
    }
    GiD_fEndCoordinates( fdm );

    GiD_fBeginElements( fdm );
    {
      GiD_fWriteCircle( fdm, 11, 3, 3.0, 1.0, 0.0, 0.0 );
      GiD_fWriteCircle( fdm, 12, 4, 2.0, 0.77, 0.77, 0.0 );
    }
    GiD_fEndElements( fdm );
  }
  GiD_fEndMesh( fdm );

  if ( !strcasecmp( format, "ascii" ) ) {
    GiD_fClosePostMeshFile( fdm );
    // we already did open the results file for writing above
    fdm = 0;
  }

  /* some dummy results */
  {
    double step = 1.0;
    for ( step = 1.0; step < 11.0; step += 1.0 ) {
      /* a scalar */
      int err = GiD_fBeginResultHeader( fdr, "a scalar", "test", step, GiD_Scalar, GiD_OnNodes, NULL );
      {
        int in = 0;
        GiD_fResultUnit( fdr, "m" );
        for ( in = 0; in < 4; in++ ) {
          GiD_fWriteScalar( fdr, in + 1, ( double )in * step );
        }
      }
      GiD_fEndResult( fdr);
      /* a vector */
      GiD_fBeginResultHeader( fdr, "a vector", "test", step, GiD_Vector, GiD_OnNodes, NULL );
      {
        /*
          const char *comp_names[] = { "sin()", "cos()", "sin() * cos()"};
        */
        int in = 0;
        const char *comp_names[] = { "in*step", "2in*step", "x+y" };
        GiD_fResultUnit( fdr, "m" );
        GiD_fResultComponents( fdr, 3, comp_names );
        for ( in = 0; in < 4; in++ ) {
          /*
            double ang = ( 1.0 + in * step) / 3.14159265355897932384626;
            double sa = sin( ang);
            double ca = cos( ang);
          */
          double sa = in * step;
          double ca = in * 2 * step;
          GiD_fWriteVector( fdr, in + 1, sa, ca, sa + ca );
        }
      }
      GiD_fEndResult( fdr);

      /* complex scalar */
      GiD_fBeginResultHeader( fdr, "complex scalar", "test", step, GiD_ComplexScalar, GiD_OnNodes, NULL );
      {
        int in = 0;
        const char *comp_names[] = { "scalar-real", "scalar-imag" };
        GiD_fResultUnit( fdr, "m" );
        GiD_fResultComponents( fdr, 2, comp_names );
        for ( in = 0; in < 4; in++ ) {
          GiD_fWriteComplexScalar( fdr, in + 1, ( double )in * step, ( double )( in * step + 0.5 ) );
        }
      }
      GiD_fEndResult( fdr);
      /* complex vector */
      GiD_fBeginResultHeader( fdr, "complex vector", "test", step, GiD_ComplexVector, GiD_OnNodes, NULL );
      {
        /*
           const char *comp_names[] = { "sin().r", "sin().i", "cos().r", "cos().i", "sin*cos.r",
           "sin*cos.i"};
        */
        int in = 0;
        const char *comp_names[] = { "in*step",      "in*step+0.5", "2in*step",
                                     "2in*step+0.5", "x_r+y_r",     "x_i+y_i" };
        GiD_fResultUnit( fdr, "m" );
        GiD_fResultComponents( fdr, 6, comp_names );
        for ( in = 0; in < 4; in++ ) {
          /*
            double ang_r = ( 1.0 + in * step) / 3.14159265355897932384626;
            double sa_r = sin( ang_r);
            double ca_r = cos( ang_r);
            double ang_i = ( 1.0 - in * step) / 3.14159265355897932384626;
            double sa_i = sin( ang_i);
            double ca_i = cos( ang_i);
          */
          double sa_r = in * step;
          double sa_i = in * step + 0.5;
          double ca_r = in * 2 * step;
          double ca_i = in * 2 * step + 0.5;
          GiD_fWriteComplexVector( fdr, in + 1, sa_r, sa_i, ca_r, ca_i, sa_r + ca_r, sa_i + ca_i );
        }
      }
      GiD_fEndResult( fdr);
    }
  }

  GiD_fClosePostResultFile( fdr);

  GiD_PostDone();

  printf( "%s written with format '%s'\n", test_filename, format );

  return 0;
}
