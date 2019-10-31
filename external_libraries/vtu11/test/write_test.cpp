#include "external/catch2/catch.hpp"
#include "vtu11.hpp"

#include <sstream>
#include <fstream>
#include <iostream> // remove

namespace vtu11
{

TEST_CASE("writeAscii_test")
{

  std::vector<double> points
  {
      0.0, 0.0, 0.5,    0.0, 0.3, 0.5,    0.0, 0.7, 0.5,    0.0, 1.0, 0.5, // 0,  1,  2,  3
      0.5, 0.0, 0.5,    0.5, 0.3, 0.5,    0.5, 0.7, 0.5,    0.5, 1.0, 0.5, // 4,  5,  6,  7
      1.0, 0.0, 0.5,    1.0, 0.3, 0.5,    1.0, 0.7, 0.5,    1.0, 1.0, 0.5  // 8,  9, 10, 11
  };

  std::vector<size_t> connectivity
  {
     0,  4,  5,  1, // 0
     1,  5,  6,  2, // 1
     2,  6,  7,  3, // 2
     4,  8,  9,  5, // 3
     5,  9, 10,  6, // 4
     6, 10, 11,  7  // 5
  };

  std::vector<size_t> offsets { 4, 8, 12, 16, 20, 24 };
  std::vector<VtkCellType> types { 9, 9, 9, 9, 9, 9 };

  Vtu11UnstructuredMesh mesh{ points, connectivity, offsets, types };

  std::vector<double> pointData1 { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 };
  std::vector<double> pointData2 { 0.1, -0.2, 0.3, -0.4, 0.5, 0.6, -0.7, 0.8, 0.9, 1.0, 1.1, -1.2 };
  std::vector<double> cellData1 { 3.2, 4.3, 5.4, 6.5, 7.6, 8.7 };
  std::vector<double> cellData2 { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0  };
  std::vector<double> cellData3 = cellData1;

  std::vector<DataSet> pointData
  {
    DataSet { std::string( "pointData1" ), 1, pointData1 },
    DataSet { std::string( "pointData2" ), 1, pointData2 },
  };

  std::vector<DataSet> cellData
  {
    DataSet { std::string( "cellData1" ), 1, cellData1 },
    DataSet { std::string( "cellData2" ), 1, cellData2 },
    DataSet { std::string( "cellData3" ), 1, cellData3 }
  };

  auto readFile = []( const std::string& filename )
  {
    std::ifstream file( filename );

    std::string contents, str;

    while( std::getline( file, str ) )
    {
      contents += str + "\n";
    }

    file.close( );

    return contents;
  };

  std::string filename = "2x3_test.vtu";

  SECTION( "ascii" )
  {
    REQUIRE_NOTHROW( write( filename, mesh, pointData, cellData ) );

    auto written = readFile( filename );
    auto expected = readFile( "testfiles/2x3_ascii.vtu" );

    CHECK( written == expected );
  }

  // The files assume that, we need to add a big endian version
  REQUIRE( endianness( ) == "LittleEndian" );

  SECTION( "base64" )
  {
    Base64BinaryWriter writer;

    REQUIRE_NOTHROW( write( filename, mesh, pointData, cellData, writer ) );

    auto written = readFile( filename );
    auto expected = readFile( "testfiles/2x3_base64.vtu" );

    CHECK( written == expected );
  }

  SECTION( "base64appended" )
  {
    Base64BinaryAppendedWriter writer;

    REQUIRE_NOTHROW( write( filename, mesh, pointData, cellData, writer ) );

    auto written = readFile( filename );
    auto expected = readFile( "testfiles/2x3_base64appended.vtu" );

    CHECK( written == expected );
  }

  SECTION( "raw" )
  {
    RawBinaryAppendedWriter writer;

    REQUIRE_NOTHROW( write( filename, mesh, pointData, cellData, writer ) );

    auto written = readFile( filename );
    auto expected = readFile( "testfiles/2x3_raw.vtu" );

    CHECK( written == expected );
  }

#ifdef VTU11_ENABLE_ZLIB
  SECTION( "raw_compressed" )
  {
    CompressedRawBinaryAppendedWriter writer;

    REQUIRE_NOTHROW( write( filename, mesh, pointData, cellData, writer ) );

    auto written = readFile( filename );
    auto expected = readFile( "testfiles/2x3_compressed.vtu" );

    CHECK( written == expected );
  }
#endif

}

} // namespace vtu11

