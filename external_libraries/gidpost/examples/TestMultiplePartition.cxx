#include <stdio.h>
#include <string.h>
#include "gidpost.h"
#include <iostream>
#include <sstream>

#ifdef WIN32
#define strcasecmp _stricmp
#define strdup _strdup
#else
#include <strings.h>
#endif

#include "my_crono.h"

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
  }

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
  }

protected:
  int m_IndexPartitionX;
  int m_IndexPartitionY;
  int m_LocalNumberOfElementsX;
  int m_LocalNumberOfElementsY;
};

void GeneratePartitionMesh( const PartitionInfo &partition )
{
  //partition.Print();
  std::string filename = partition.GetOutputFileName( G_prefix );
  filename += ".post.msh";
  GiD_FILE fd = GiD_fOpenPostMeshFile( filename.c_str(), GiD_PostAscii );
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
  GiD_fClosePostMeshFile( fd );
}

void SimulateResultsForPartition( const PartitionInfo &partition )
{
  std::string filename = partition.GetOutputFileName( "ExamplePartitioned" );
  filename += ".post.res";
  GiD_FILE fd = GiD_fOpenPostResultFile( filename.c_str(), GiD_PostAscii );
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

void VisitPartitions( const PartitionInfo &partition, int i0, int j0, int n, int m )
{
  PartitionInfo partitionCopy( partition );
  for ( int jp = j0, j = 0; j < m; jp++, ++j )
    {
    for ( int ip = i0, i = 0; i < n; ip++, ++i )
      {
      partitionCopy.SetIndexPartitionX( ip );
      partitionCopy.SetIndexPartitionY( jp );
      int p_id = partitionCopy.GetPartitionID();
      std::cout << "-- partition " << p_id << "\n";
      GeneratePartitionMesh( partitionCopy );
      SimulateResultsForPartition( partitionCopy );
      }
    }
}

// #define PART_N 10
// #define PART_M 10
#define PART_N 4
#define PART_M 4

int main( int argc, const char *argv[] ) {
#define BUF_SIZE 1024
    char buf[ BUF_SIZE ];
    char *format = NULL;

    G_prefix = strdup( "ExamplePartitioned" );
    if ( !format ) {
      format = strdup( "ascii" );
    }
    if ( !strcasecmp( format, "ascii" ) ) {
      G_format = GiD_PostAscii;
      strcpy( buf, "_*.post.msh | _*.post.res" );
    } else if ( !strcasecmp( format, "bin" ) ) {
      G_format = GiD_PostBinary;
      strcpy( buf, "_*.post.bin" );
    } else if ( !strcasecmp( format, "hdf5" ) ) {
      G_format = GiD_PostHDF5;
      strcpy( buf, "_*.post.h5" );
    }

    printf( "version = %s\n", GiD_PostGetVersion() );

    Crono clk;
    GiD_PostInit();

    PartitionInfo partition;
    partition.SetNumberOfPartitionX( PART_N );
    partition.SetNumberOfPartitionY( PART_M );
    partition.SetLocalNumberOfElementsX( 10 );
    partition.SetLocalNumberOfElementsY( 10 );
    const int n = PART_N / 2;
    const int m = PART_M / 2;
    std::cout << "Number of partitions to write = "
              << partition.GetNumberOfPartitionX() * partition.GetNumberOfPartitionY() << std::endl;
    VisitPartitions( partition, 0, 0, n, m );
    VisitPartitions( partition, n, 0, n, m );
    VisitPartitions( partition, n, m, n, m );
    VisitPartitions( partition, 0, m, n, m );

    GiD_PostDone();
    std::cout << "Files created of the form : " << G_prefix << buf << std::endl;
    std::cout << "  duration (s.) = " << clk.end() << std::endl;

    free( format);
    format = NULL;

    return 0;
}
