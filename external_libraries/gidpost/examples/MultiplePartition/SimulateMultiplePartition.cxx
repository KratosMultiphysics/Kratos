#include <stdio.h>
#include <string.h>
#include "gidpost.h"
#include <iostream>
#include <sstream>
#if defined(HAVE_BOOST_THREAD)
#include <boost/thread/thread.hpp>
#endif

#ifdef WIN32
#define strcasecmp  _stricmp
#define strdup  _strdup
#else
#include <strings.h>
#endif

#include "my_crono.h"

bool G_thread_verbose = false;
GiD_PostMode G_format = GiD_PostAscii;
char *G_prefix = NULL;

double Random()
{
  return rand()/double(RAND_MAX);
}

class  PartitionInfo;

class DomainInfo 
{
public:
 
  DomainInfo()
  {
    this->m_SizeX = this->m_SizeY = 1.0;
    this->m_NumberOfPartitionX = 1;
    this->m_NumberOfPartitionY = 1;
  }

  void SetSizeX( double x )
  {
    this->m_SizeX = x;
  }

  double GetSizeX() const
  {
    return this->m_SizeX;
  }

  void SetSizeY( double y )
  {
    this->m_SizeY = y;
  }

  double GetSizeY() const
  {
    return this->m_SizeY;
  }

  void SetNumberOfPartitionX( int v )
  {
    this->m_NumberOfPartitionX = v;
  }

  int GetNumberOfPartitionX() const
  {
    return this->m_NumberOfPartitionX;
  }

  void SetNumberOfPartitionY( int v )
  {
    this->m_NumberOfPartitionY = v;
  }

  int GetNumberOfPartitionY() const
  {
    return this->m_NumberOfPartitionY;
  }

  void Print( ) const
  {
    std::cout << "SizeX: " << this->m_SizeX << std::endl;
    std::cout << "SizeY: " << this->m_SizeY << std::endl;
    std::cout << "NumberOfPartitionX: " << this->m_NumberOfPartitionX << std::endl;
    std::cout << "NumberOfPartitionY: " << this->m_NumberOfPartitionY << std::endl;
  }
  
protected:

  double m_SizeX;
  double m_SizeY;
  int m_NumberOfPartitionX;
  int m_NumberOfPartitionY;
};

class PartitionInfo : public DomainInfo
{
public:

  PartitionInfo() : DomainInfo()
  {
    this->m_IndexPartitionX = 0;
    this->m_IndexPartitionY = 0;
    this->m_LocalNumberOfElementsX = 1;
    this->m_LocalNumberOfElementsY = 1;
    this->m_file = 0;
  }

  GiD_FILE getFileID() const { return m_file;}
  void setFileID( GiD_FILE fd) { m_file = fd;}

  int GetPartitionID( ) const
  {
    return this->GetNumberOfPartitionX() * this->GetIndexPartitionY() + this->GetIndexPartitionX();
  }

  int GetIndexPartitionX( ) const
  {
    return this->m_IndexPartitionX;
  }

  void SetIndexPartitionX( int v )
  {
    this->m_IndexPartitionX = v;
  }
  
  int GetIndexPartitionY( ) const
  {
    return this->m_IndexPartitionY;
  }

  void SetIndexPartitionY( int v )
  {
    this->m_IndexPartitionY = v;
  }
  
  int GetLocalNumberOfElementsX() const
  {
    return this->m_LocalNumberOfElementsX;
  }
  
  void SetLocalNumberOfElementsX( int v )
  {
    this->m_LocalNumberOfElementsX = v;
  }
  
  int GetLocalNumberOfElementsY() const
  {
    return this->m_LocalNumberOfElementsY;
  }
  
  void SetLocalNumberOfElementsY( int v )
  {
    this->m_LocalNumberOfElementsY = v;
  }
  
  int GetGlobalIndexX( int k ) const
  {
    return this->GetIndexPartitionX() * this->GetLocalNumberOfElementsX() + k;
  }

  int GetGlobalIndexY( int l ) const
  {
    return this->GetIndexPartitionY() * this->GetLocalNumberOfElementsY() + l;
  }

  int GetGlobalElementId( int ix, int iy ) const
  {
    int i = this->GetGlobalIndexX( ix );
    int j = this->GetGlobalIndexY( iy );
    return ( j * this->GetNumberOfPartitionX() * this->GetLocalNumberOfElementsX() + i ) + 1;
  }

  int GetGlobalNodeId( int ix, int iy ) const
  {
    int i = this->GetGlobalIndexX( ix );
    int j = this->GetGlobalIndexY( iy );
    return ( j * this->GetNumberOfPartitionX() * ( this->GetLocalNumberOfElementsX() + 1 ) + i ) + 1;
  }

  double GetGlobalCoordinateX( int i ) const
  {
    double delta = this->GetSizeX() / ( this->GetNumberOfPartitionX() * this->GetLocalNumberOfElementsX() );
    return ( this->GetGlobalIndexX( i ) - 1 ) * delta;
  }

  double GetGlobalCoordinateY( int i ) const
  {
    double delta = this->GetSizeY() / ( this->GetNumberOfPartitionY() * this->GetLocalNumberOfElementsY() );
    return ( this->GetGlobalIndexY( i ) - 1 ) * delta;
  }

  std::string GetOutputFileName( const std::string& prefix ) const
  {
    std::stringstream ss;

    ss << prefix << "_" << this->GetPartitionID();
    return ss.str();
  }

  void Print( ) const
  {
    DomainInfo::Print();
    std::cout << "IndexPartitionX: " << this->m_IndexPartitionX << std::endl;
    std::cout << "IndexPartitionY: " << this->m_IndexPartitionY << std::endl;
    std::cout << "LocalNumberOfElementsX: " << this->m_LocalNumberOfElementsX << std::endl;
    std::cout << "LocalNumberOfElementsY: " << this->m_LocalNumberOfElementsY << std::endl;
    std::cout << "File ID: " << this->m_file << std::endl;
  }

protected:
  int m_IndexPartitionX;
  int m_IndexPartitionY;
  int m_LocalNumberOfElementsX;
  int m_LocalNumberOfElementsY;
  GiD_FILE m_file;
};

GiD_FILE GeneratePartitionMesh( const PartitionInfo &partition )
{
  //partition.Print();
  std::string filename = partition.GetOutputFileName( G_prefix );
  GiD_FILE fd = 0;
  if ( G_format == GiD_PostAscii) {
    filename += ".post.msh";
    fd = GiD_fOpenPostMeshFile( filename.c_str(), G_format );
  } else if ( G_format == GiD_PostBinary) {
    filename += ".post.bin";
    fd = GiD_fOpenPostResultFile( filename.c_str(), G_format );
  } else if ( G_format == GiD_PostHDF5) {
    filename += ".post.h5";
    fd = GiD_fOpenPostResultFile( filename.c_str(), G_format );
  }
  if ( fd == 0 ) {
    std::cout << "Can not open mesh file for writing: " << filename << std::endl;
    return 0;
  }
  GiD_fBeginMesh( fd, "TestMsh", GiD_2D, GiD_Quadrilateral, 4 );
  GiD_fBeginCoordinates( fd );
  for( int j = 0; j < partition.GetLocalNumberOfElementsY() + 1; j++ )
    {
    for( int i = 0; i < partition.GetLocalNumberOfElementsX() + 1; i++ )
      {
      int id = partition.GetGlobalNodeId( i, j );
      double x = partition.GetGlobalCoordinateX( i );
      double y = partition.GetGlobalCoordinateY( j );
      GiD_fWriteCoordinates2D( fd, id, x, y );
      }
    }
  GiD_fEndCoordinates( fd );
  GiD_fBeginElements( fd );
  int element[4];
  for( int j = 0; j < partition.GetLocalNumberOfElementsY(); j++ )
    {
    for( int i = 0; i < partition.GetLocalNumberOfElementsX(); i++ )
      {
      int id = partition.GetGlobalElementId( i, j );
      element[0] = partition.GetGlobalNodeId( i, j );
      element[1] = partition.GetGlobalNodeId( i + 1, j );
      element[2] = partition.GetGlobalNodeId( i + 1, j + 1);
      element[3] = partition.GetGlobalNodeId( i, j + 1);
      GiD_fWriteElement( fd, id, element );
      }
    }
  GiD_fEndElements( fd );
  GiD_fEndMesh( fd );
  if ( G_format == GiD_PostAscii) {
    GiD_fClosePostMeshFile( fd );
    fd = 0;
  }
  return fd;
}

void SimulateResultsForPartition( const PartitionInfo &partition )
{
  std::string filename = partition.GetOutputFileName( G_prefix );
  GiD_FILE fd = partition.getFileID();
  if ( G_format == GiD_PostAscii) {
    filename += ".post.res";
    fd = GiD_fOpenPostResultFile( filename.c_str(), GiD_PostAscii );
  }
  if ( fd == 0 ) {
    std::cout << "Can not open results file for writing: " << filename << std::endl;
    return;
  }
  GiD_fBeginResult( fd,
                    "EscalarNodos", "Analysis", 1.0, GiD_Scalar, GiD_OnNodes,
                    NULL, NULL, 0, NULL);
  for( int j = 0; j < partition.GetLocalNumberOfElementsY() + 1; j++ )
    {
    for( int i = 0; i < partition.GetLocalNumberOfElementsX() + 1; i++ )
      {
      int id = partition.GetGlobalNodeId( i, j );
      GiD_fWriteScalar( fd, id, Random() );
      }
    }
  GiD_fEndResult(fd );
  GiD_fClosePostResultFile( fd );
}

void VisitPartitions( const PartitionInfo &partition, int i0, int j0, int n, int m ) {
  PartitionInfo partitionCopy( partition );
  for ( int jp = j0, j = 0; j < m; jp++, ++j ) {
    for ( int ip = i0, i = 0; i < n; ip++, ++i ) {
      partitionCopy.SetIndexPartitionX( ip );
      partitionCopy.SetIndexPartitionY( jp );
      int p_id = partitionCopy.GetPartitionID();
      std::cout << "-- partition " << p_id << "\n";
      if ( G_thread_verbose ) {
        std::cout << "-- partition " << p_id << " initializing\n";
      }
      partitionCopy.setFileID( 0 );
      if ( G_thread_verbose ) {
        std::cout << "-- partition " << p_id << " mesh generate and write\n";
      }
      GiD_FILE fd = GeneratePartitionMesh( partitionCopy );
      partitionCopy.setFileID( fd );
      if ( G_thread_verbose ) {
        std::cout << "-- partition " << p_id << " results generate and write\n";
      }
      SimulateResultsForPartition( partitionCopy );
      partitionCopy.setFileID( 0 );
      if ( G_thread_verbose ) {
        std::cout << "-- partition " << p_id << " done\n";
      }
    }
  }
}


void print_help_and_exit( const char *full_prog_name ) {
  // get only the name of the executable
  const char *prog_name = &full_prog_name[ strlen( full_prog_name ) - 1 ];
  for ( ; prog_name > full_prog_name; prog_name-- ) {
    if ( ( *prog_name == '/' ) || ( *prog_name == '\\' ) ) {
      break;
    }
  }
  prog_name++;
  printf( "Usage: %s [ -h] [ -v] [ -f ascii | bin | hdf5] filename_prefix\n", prog_name );
  printf( "    will create filename_prefix_PartitionID.post.{ msh,res | bin | h5}\n" );
  printf( "    Use: %s -- filename_prefix          if filename_prefix begins with a '-'\n", prog_name );
  exit( 0 );
}

// #define PART_N 10
// #define PART_M 10
#define PART_N 4
#define PART_M 4

int main( int argc, const char* argv[] )
{
#define BUF_SIZE 1024
  char buf[ BUF_SIZE ];
  char *prefix = NULL;
  char *format = NULL;
  int skip_options = 0;

  for ( int ia = 1; ia < argc; ia++ ) {
    if ( !skip_options && ( argv[ ia ][ 0 ] == '-' ) ) {
      char opt = ( char )tolower( argv[ ia ][ 1 ] );
      if ( ( opt == 'h' ) || !strcasecmp( argv[ ia ], "--help" ) ) {
        print_help_and_exit( argv[ 0 ] );
      } else if ( opt == 'v' ) {
        G_thread_verbose = true;
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
        printf( "Unknown option '%s'.\n", argv[ ia ] );
        print_help_and_exit( argv[ 0 ] );
      }
    } else {
      prefix = strdup( argv[ ia ] );
      break;
    }
  }
  if ( !prefix ) {
    prefix = strdup( "ExamplePartitioned");
  }

  if ( !format) {
    format = strdup( "ascii");
  }
  if ( !strcasecmp( format, "ascii")) {
      G_format = GiD_PostAscii;
      strcpy( buf, "_*.post.msh | _*.post.res");
  } else if ( !strcasecmp( format, "bin")) {
      G_format = GiD_PostBinary;
      strcpy( buf, "_*.post.bin");
  } else if ( !strcasecmp( format, "hdf5")) {
      G_format = GiD_PostHDF5;
      strcpy( buf, "_*.post.h5");
  }
  free( format);
  G_prefix = prefix;

  printf( "version = %s\n", GiD_PostGetVersion() );

  Crono clk;
  PartitionInfo partition;

  GiD_PostInit();

  partition.SetNumberOfPartitionX( PART_N );
  partition.SetNumberOfPartitionY( PART_M );
  partition.SetLocalNumberOfElementsX( 250 );
  partition.SetLocalNumberOfElementsY( 250 );
  const int n = PART_N / 2;
  const int m = PART_M / 2;

  std::cout << "Number of partitions to write = " << partition.GetNumberOfPartitionX() * partition.GetNumberOfPartitionY() << std::endl;
#if defined(HAVE_BOOST_THREAD)
  std::cout << "Using boost::thread\n";
  boost::thread **ptrThreads;
  
  ptrThreads = new boost::thread*[4];

  ptrThreads[0] = new boost::thread( VisitPartitions, partition, 0, 0, n, m );
  ptrThreads[1] = new boost::thread( VisitPartitions, partition, n, 0, n, m );
  ptrThreads[2] = new boost::thread( VisitPartitions, partition, n, m, n, m );
  ptrThreads[3] = new boost::thread( VisitPartitions, partition, 0, m, n, m );
  for ( int i = 0; i < 4; i++ ) {
    std::cout << "- joining boost::thread " << i + 1 << " / 4" << std::endl;
    if ( ptrThreads[ i ]->joinable() ) {
      ptrThreads[ i ]->join();
    }
  }
  std::cout << "- end\n";

#else
  std::cout << "Not using boost::thread\n";
  VisitPartitions( partition, 0, 0, n, m );
  VisitPartitions( partition, n, 0, n, m );
  VisitPartitions( partition, n, m, n, m );
  VisitPartitions( partition, 0, m, n, m );
#endif
  
  GiD_PostDone();
  std::cout << "Files created of the form : " << G_prefix << buf << std::endl;
  std::cout << "  duration (s.) = " << clk.end() << std::endl;

#if defined(HAVE_BOOST_THREAD)
  for( int i = 0; i < 4; i++ )
    {
    delete ptrThreads[i];
    }
  delete[] ptrThreads;
#endif

  free( G_prefix );
  G_prefix = NULL;

  return 0;
}
