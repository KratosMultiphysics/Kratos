//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <map>
#include <numeric>

// External includes

// Project includes
#include "vtu_output.h"
#include "containers/model.h"
#include "includes/kratos_filesystem.h"
#include "processes/fast_transfer_between_model_parts_process.h"
#include "utilities/parallel_utilities.h"

#define VTU11_ENABLE_ZLIB

//          __        ____ ____
// ___  ___/  |_ __ _/_   /_   |
// \  \/ /\   __\  |  \   ||   |
//  \   /  |  | |  |  /   ||   |
//   \_/   |__| |____/|___||___|
//
//  License: BSD License ; see LICENSE
//

// AUTOMATICALLY GENERATED SINGLE HEADER FILE.


#ifndef VTU11_ALIAS_HPP
#define VTU11_ALIAS_HPP

#include <string>
#include <map>
#include <utility>
#include <vector>

namespace vtu11
{

using StringStringMap = std::map<std::string, std::string>;

enum class DataSetType : int
{
    PointData = 0, CellData = 1
};

using DataSetInfo = std::tuple<std::string, DataSetType, size_t>;
using DataSetData = std::vector<double>;

using VtkCellType = std::int8_t;
using VtkIndexType = std::int64_t;

using HeaderType = size_t;
using Byte = unsigned char;

} // namespace vtu11

#ifndef VTU11_ASCII_FLOATING_POINT_FORMAT
    #define VTU11_ASCII_FLOATING_POINT_FORMAT "%.6g"
#endif

#if defined(__cplusplus) && __cplusplus >= 201703L
    #if __has_include(<filesystem>) // has_include is C++17
        #include <filesystem>
        namespace vtu11fs = std::filesystem;
    #elif __has_include(<experimental/filesystem>)
        #include <experimental/filesystem>
        namespace vtu11fs = std::experimental::filesystem;
    #else
        #include "ghc/filesystem.hpp"
        namespace vtu11fs = ghc::filesystem;
    #endif
#else
    #include "ghc/filesystem.hpp"
    namespace vtu11fs = ghc::filesystem;
#endif

#endif // VTU11_ALIAS_HPP

#ifndef VTU11_WRITER_HPP
#define VTU11_WRITER_HPP


namespace vtu11
{

struct AsciiWriter
{
  template<typename T>
  void writeData( std::ostream& output,
                  const std::vector<T>& data );

  void writeAppended( std::ostream& output );

  void addHeaderAttributes( StringStringMap& attributes );
  void addDataAttributes( StringStringMap& attributes );

  StringStringMap appendedAttributes( );
};

struct Base64BinaryWriter
{
  template<typename T>
  void writeData( std::ostream& output,
                  const std::vector<T>& data );

  void writeAppended( std::ostream& output );

  void addHeaderAttributes( StringStringMap& attributes );
  void addDataAttributes( StringStringMap& attributes );

  StringStringMap appendedAttributes( );
};

struct Base64BinaryAppendedWriter
{
  template<typename T>
  void writeData( std::ostream& output,
                  const std::vector<T>& data );

  void writeAppended( std::ostream& output );

  void addHeaderAttributes( StringStringMap& attributes );
  void addDataAttributes( StringStringMap& attributes );

  StringStringMap appendedAttributes( );

  size_t offset = 0;

  std::vector<std::pair<const char*, HeaderType>> appendedData;
};

struct RawBinaryAppendedWriter
{
  template<typename T>
  void writeData( std::ostream& output,
                  const std::vector<T>& data );

  void writeAppended( std::ostream& output );

  void addHeaderAttributes( StringStringMap& attributes );
  void addDataAttributes( StringStringMap& attributes );

  StringStringMap appendedAttributes( );

  size_t offset = 0;

  std::vector<std::pair<const char*, HeaderType>> appendedData;
};

} // namespace vtu11


#endif // VTU11_WRITER_HPP

#ifndef VTU11_UTILITIES_HPP
#define VTU11_UTILITIES_HPP


#include <functional>
#include <limits>
#include <type_traits>

namespace vtu11
{

#define VTU11_THROW( message ) throw std::runtime_error( message )
#define VTU11_CHECK( expr, message ) if( !( expr ) ) VTU11_THROW ( message )

std::string endianness( );

template<typename Iterator>
std::string base64Encode( Iterator begin, Iterator end );

size_t encodedNumberOfBytes( size_t rawNumberOfBytes );

class ScopedXmlTag final
{
public:
    ScopedXmlTag( std::ostream& output,
                  const std::string& name,
                  const StringStringMap& attributes );

    ~ScopedXmlTag( );

private:
    std::function<void( )> closeTag;
};

void writeEmptyTag( std::ostream& output,
                    const std::string& name,
                    const StringStringMap& attributes );

template<typename T> inline
std::string appendSizeInBits( const char* str )
{
    return str + std::to_string( sizeof( T ) * 8 );
}

template<typename T> inline
typename std::enable_if<std::numeric_limits<T>::is_integer &&
                        std::numeric_limits<T>::is_signed, std::string>::type
    dataTypeString( )
{
    return appendSizeInBits<T>( "Int" );
}

template<typename T> inline
typename std::enable_if<std::numeric_limits<T>::is_integer &&
                       !std::numeric_limits<T>::is_signed, std::string>::type
    dataTypeString( )
{
    return appendSizeInBits<T>( "UInt" );
}

template<typename T> inline
typename std::enable_if<std::is_same<T, double>::value ||
                        std::is_same<T, float>::value, std::string>::type
    dataTypeString( )
{
    return appendSizeInBits<T>( "Float" );
}

} // namespace vtu11


#endif // VTU11_UTILITIES_HPP

#ifndef VTU11_ZLIBWRITER_HPP
#define VTU11_ZLIBWRITER_HPP


#ifdef VTU11_ENABLE_ZLIB

namespace vtu11
{

struct CompressedRawBinaryAppendedWriter
{
  template<typename T>
  void writeData( std::ostream& output,
                  const std::vector<T>& data );

  void writeAppended( std::ostream& output );

  void addHeaderAttributes( StringStringMap& attributes );
  void addDataAttributes( StringStringMap& attributes );

  StringStringMap appendedAttributes( );

  size_t offset = 0;

  std::vector<std::vector<std::vector<std::uint8_t>>> appendedData;
  std::vector<std::vector<HeaderType>> headers;
};

} // namespace vtu11


#endif // VTU11_ENABLE_ZLIB

#endif // VTU11_ZLIBWRITER_HPP

#ifndef VTU11_VTU11_HPP
#define VTU11_VTU11_HPP


namespace vtu11
{

struct Vtu11UnstructuredMesh
{
  std::vector<double>& points_;
  std::vector<VtkIndexType>& connectivity_;
  std::vector<VtkIndexType>& offsets_;
  std::vector<VtkCellType>& types_;

  std::vector<double>& points( ){ return points_; }
  std::vector<VtkIndexType>& connectivity( ){ return connectivity_; }
  std::vector<VtkIndexType>& offsets( ){ return offsets_; }
  std::vector<VtkCellType>& types( ){ return types_; }

  size_t numberOfPoints( ){ return points_.size( ) / 3; }
  size_t numberOfCells( ){ return types_.size( ); }
};

/*! Write modes (not case sensitive):
 *
 *  - Ascii
 *  - Base64Inline
 *  - Base64Appended
 *  - RawBinary
 *  - RawBinaryCompressed
 *
 *  Comments:
 *  - RawCompressedBinary needs zlib to be present. If VTU11_ENABLE_ZLIB
 *    is not defined, the uncompressed version is used instead.
 *  - Compressing data takes more time than writing more data uncompressed
 *  - Ascii produces surprisingly small files, is nice to debug, but
 *    is rather slow to read in Paraview. Archiving ascii .vtu files using
 *    a standard zip tool (for example) produces decently small file sizes.
 *  - Writing raw binary data breakes the xml standard. To still produce
 *    valid xml files you can use base64 encoding, at the cost of having
 *    about 30% times larger files.
 *  - Both raw binary modes use appended format
 */

template<typename MeshGenerator>
void writeVtu( const std::string& filename,
               MeshGenerator& mesh,
               const std::vector<DataSetInfo>& dataSetInfo,
               const std::vector<DataSetData>& dataSetData,
               const std::string& writeMode = "RawBinaryCompressed" );

void writePVtu( const std::string& path,
                const std::string& baseName,
                const std::vector<DataSetInfo>& dataSetInfo,
                size_t numberOfFiles );

template<typename MeshGenerator>
void writePartition( const std::string& path,
                     const std::string& baseName,
                     MeshGenerator& mesh,
                     const std::vector<DataSetInfo>& dataSetInfo,
                     const std::vector<DataSetData>& dataSetData,
                     size_t fileId,
                     const std::string& writeMode = "RawBinaryCompressed" );

} // namespace vtu11


#endif // VTU11_VTU11_HPP

#ifndef VTU11_UTILITIES_IMPL_HPP
#define VTU11_UTILITIES_IMPL_HPP

#include <array>

namespace vtu11
{
namespace detail
{

inline void writeTag( std::ostream& output,
                      const std::string& name,
                      const StringStringMap& attributes,
                      const std::string& tagEnd )
{
    output << "<" << name;

    for( const auto& attribute : attributes )
    {
        output << " " << attribute.first << "=\"" << attribute.second << "\"";
    }

    output << tagEnd << "\n";
}

} // namespace detail

inline ScopedXmlTag::ScopedXmlTag( std::ostream& output,
                                   const std::string& name,
                                   const StringStringMap& attributes ) :
   closeTag( [ &output, name ]( ){ output << "</" << name << ">\n"; } )
{
    detail::writeTag( output, name, attributes, ">" );
}

inline ScopedXmlTag::~ScopedXmlTag( )
{
  closeTag( );
}

inline void writeEmptyTag( std::ostream& output,
                           const std::string& name,
                           const StringStringMap& attributes )
{
  detail::writeTag( output, name, attributes, "/>" );
}

inline std::string endianness( )
{
   int i = 0x0001;

   if( *reinterpret_cast<char*>( &i ) != 0 )
   {
     return "LittleEndian";
   }
   else
   {
     return "BigEndian";
   }
}

constexpr char base64Map[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

template<typename Iterator>
inline std::string base64Encode( Iterator begin, Iterator end )
{
    constexpr size_t size = sizeof( decltype( *begin ) );

    size_t length = static_cast<size_t>( std::distance( begin, end ) );
    size_t rawBytes = length * size;
    size_t encodedBytes = ( rawBytes / 3 + 1 ) * 4;

    std::string result;

    result.reserve( encodedBytes );

    auto it = begin;
    size_t byteIndex = 0;

    auto next = [&]( )
    {
        char byte = *( reinterpret_cast<const char*>( &( *it ) ) + byteIndex++ );

        if( byteIndex == size )
        {
            it++;
            byteIndex = 0;
        }

        return byte;
    };

    auto encodeTriplet = [&]( std::array<char, 3> bytes, size_t padding )
    {
        char tmp[5] = { base64Map[(   bytes[0] & 0xfc ) >> 2],
                        base64Map[( ( bytes[0] & 0x03 ) << 4 ) + ( ( bytes[1] & 0xf0 ) >> 4 )],
                        base64Map[( ( bytes[1] & 0x0f ) << 2 ) + ( ( bytes[2] & 0xc0 ) >> 6 )],
                        base64Map[bytes[2] & 0x3f],
                        '\0' };

        std::fill( tmp + 4 - padding, tmp + 4, '=' );

        result += tmp;
    };

    for( size_t i = 0; i < rawBytes / 3; ++i )
    {
        encodeTriplet( { next( ), next( ), next( ) }, 0 );
    }

    if( it != end )
    {
        std::array<char, 3> bytes { '\0', '\0', '\0' };

        size_t remainder = static_cast<size_t>( std::distance( it, end ) ) * size - static_cast<size_t>( byteIndex );

        for( size_t i = 0; i < remainder; ++i )
        {
            bytes[i] = next( );
        }

        encodeTriplet( bytes, 3 - remainder );
    }

    return result;
}

inline size_t encodedNumberOfBytes( size_t rawNumberOfBytes )
{
    if( rawNumberOfBytes != 0 )
    {
        return ( ( rawNumberOfBytes - 1 ) / 3 + 1 ) * 4;
    }
    else
    {
        return 0;
    }
}

} // namespace vtu11

#endif // VTU11_UTILITIES_IMPL_HPP

#ifndef VTU11_WRITER_IMPL_HPP
#define VTU11_WRITER_IMPL_HPP


#include <fstream>

namespace vtu11
{
namespace detail
{

template<typename T> inline
void writeNumber( char (&buffer)[64], T value )
{
    VTU11_THROW( "Invalid data type." );
}

#define __VTU11_WRITE_NUMBER_SPECIALIZATION( string, type )    \
template<> inline                                             \
void writeNumber<type>( char (&buffer)[64], type value )         \
{                                                             \
    std::snprintf( buffer, sizeof( buffer ), string, value ); \
}

__VTU11_WRITE_NUMBER_SPECIALIZATION( VTU11_ASCII_FLOATING_POINT_FORMAT, double )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%lld", long long int )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%ld" , long int )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%d"  , int )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%hd" , short )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%hhd", char )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%llu", unsigned long long int )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%ld" , unsigned long int )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%d"  , unsigned int )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%hd" , unsigned short )
__VTU11_WRITE_NUMBER_SPECIALIZATION( "%hhd", unsigned char )

} // namespace detail

template<typename T>
inline void AsciiWriter::writeData( std::ostream& output,
                                    const std::vector<T>& data )
{
    char buffer[64];

    for( auto value : data )
    {
        detail::writeNumber( buffer, value );

        output << buffer << " ";
    }

    output << "\n";
}

template<>
inline void AsciiWriter::writeData( std::ostream& output,
                                    const std::vector<std::int8_t>& data )
{
  for( auto value : data )
  {
      output << static_cast<int>( value ) << " ";
  }

  output << "\n";
}

inline void AsciiWriter::writeAppended( std::ostream& )
{

}

inline void AsciiWriter::addHeaderAttributes( StringStringMap& )
{
}

inline void AsciiWriter::addDataAttributes( StringStringMap& attributes )
{
  attributes["format"] = "ascii";
}

inline StringStringMap AsciiWriter::appendedAttributes( )
{
  return { };
}


template<typename T>
inline void Base64BinaryWriter::writeData( std::ostream& output,
                                           const std::vector<T>& data )
{
  HeaderType numberOfBytes = data.size( ) * sizeof( T );

  output << base64Encode( &numberOfBytes, &numberOfBytes + 1 );
  output << base64Encode( data.begin( ), data.end( ) );

  output << "\n";
}

inline void Base64BinaryWriter::writeAppended( std::ostream& )
{

}

inline void Base64BinaryWriter::addHeaderAttributes( StringStringMap& attributes )
{
  attributes["header_type"] = dataTypeString<HeaderType>( );
}

inline void Base64BinaryWriter::addDataAttributes( StringStringMap& attributes )
{
  attributes["format"] = "binary";
}

inline StringStringMap Base64BinaryWriter::appendedAttributes( )
{
  return { };
}


template<typename T>
inline void Base64BinaryAppendedWriter::writeData( std::ostream&,
                                                   const std::vector<T>& data )
{
  HeaderType rawBytes = data.size( ) * sizeof( T );

  appendedData.emplace_back( reinterpret_cast<const char*>( &data[0] ), rawBytes );

  offset += encodedNumberOfBytes( rawBytes + sizeof( HeaderType ) );
}

inline void Base64BinaryAppendedWriter::writeAppended( std::ostream& output )
{
  for( auto dataSet : appendedData )
  {
    std::vector<char> data( dataSet.second + sizeof( HeaderType ) );

    *reinterpret_cast<HeaderType*>( &data[0] ) = dataSet.second;

    std::copy( dataSet.first, dataSet.first + dataSet.second, &data[ sizeof( HeaderType ) ] );

    output << base64Encode( data.begin( ), data.end( ) );
  }

  output << "\n";
}

inline void Base64BinaryAppendedWriter::addHeaderAttributes( StringStringMap& attributes )
{
  attributes["header_type"] = dataTypeString<HeaderType>( );
}

inline void Base64BinaryAppendedWriter::addDataAttributes( StringStringMap& attributes )
{
  attributes["format"] = "appended";
  attributes["offset"] = std::to_string( offset );
}

inline StringStringMap Base64BinaryAppendedWriter::appendedAttributes( )
{
  return { { "encoding", "base64" } };
}


template<typename T>
inline void RawBinaryAppendedWriter::writeData( std::ostream&,
                                                const std::vector<T>& data )
{
  HeaderType rawBytes = data.size( ) * sizeof( T );

  appendedData.emplace_back( reinterpret_cast<const char*>( &data[0] ), rawBytes );

  offset += sizeof( HeaderType ) + rawBytes;
}

inline void RawBinaryAppendedWriter::writeAppended( std::ostream& output )
{
  for( auto dataSet : appendedData )
  {
    const char* headerBegin = reinterpret_cast<const char*>( &dataSet.second );

    for( const char* ptr = headerBegin; ptr < headerBegin + sizeof( HeaderType ); ++ptr )
    {
      output << *ptr;
    }

    for( const char* ptr = dataSet.first; ptr < dataSet.first + dataSet.second; ++ptr )
    {
      output << *ptr;
    }
  }

  output << "\n";
}

inline void RawBinaryAppendedWriter::addHeaderAttributes( StringStringMap& attributes )
{
  attributes["header_type"] = dataTypeString<HeaderType>( );
}

inline void RawBinaryAppendedWriter::addDataAttributes( StringStringMap& attributes )
{
  attributes["format"] = "appended";
  attributes["offset"] = std::to_string( offset );
}

inline StringStringMap RawBinaryAppendedWriter::appendedAttributes( )
{
  return { { "encoding", "raw" } };
}

} // namespace vtu11

#endif // VTU11_WRITER_IMPL_HPP

#ifndef VTU11_ZLIBWRITER_IMPL_HPP

#include "zlib.h"

namespace vtu11
{
namespace detail
{

template<typename T>
std::vector<HeaderType> zlibCompressData( const std::vector<T>& data,
                                          std::vector<std::vector<Byte>>& targetBlocks,
                                          size_t blocksize = 32768 ) // 2^15
{
  std::vector<HeaderType> header( 3, 0 );

  if( data.empty( ) )
  {
    return header;
  }

  auto compressedBuffersize = compressBound( blocksize );

  Byte* buffer = new Byte[compressedBuffersize];
  Byte* currentByte = const_cast<Byte*>( reinterpret_cast<const Byte*>( &data[0] ) );

  size_t numberOfBytes = data.size( ) * sizeof( T );
  size_t numberOfBlocks = ( numberOfBytes - 1 ) / blocksize + 1;

  using ZlibSizeType = decltype( compressedBuffersize );

  auto compressBlock = [&]( ZlibSizeType numberOfBytesInBlock )
  {
    ZlibSizeType compressedLength = compressedBuffersize;

    int errorCode = compress( buffer, &compressedLength, currentByte, numberOfBytesInBlock );

    if( errorCode != Z_OK )
    {
      delete[] buffer;

      throw std::runtime_error( "Error in zlib compression (code " + std::to_string( errorCode ) + ")." );
    }

    targetBlocks.emplace_back( buffer, buffer + compressedLength );
    header.push_back( compressedLength );

    currentByte += numberOfBytesInBlock;
  };

  for( size_t iBlock = 0; iBlock < numberOfBlocks - 1; ++iBlock )
  {
    compressBlock( blocksize );
  }

  size_t remainder = numberOfBytes - ( numberOfBlocks - 1 ) * blocksize;

  compressBlock( remainder );

  delete[] buffer;

  header[0] = header.size( ) - 3;
  header[1] = blocksize;
  header[2] = remainder;

  return header;
}

} // detail

template<typename T>
inline void CompressedRawBinaryAppendedWriter::writeData( std::ostream&,
                                                          const std::vector<T>& data )
{
  std::vector<std::vector<Byte>> compressedBlocks;

  auto header = detail::zlibCompressData( data, compressedBlocks );

  offset += sizeof( HeaderType ) * header.size( );

  for( const auto& compressedBlock : compressedBlocks )
  {
    offset += compressedBlock.size( );
  }

  appendedData.push_back( std::move( compressedBlocks ) );
  headers.push_back( std::move( header ) );
}

inline void CompressedRawBinaryAppendedWriter::writeAppended( std::ostream& output )
{
  for( size_t iDataSet = 0; iDataSet < appendedData.size( ); ++iDataSet )
  {
    const char* headerBegin = reinterpret_cast<const char*>( &headers[iDataSet][0] );
    size_t numberOfHeaderBytes = headers[iDataSet].size( ) * sizeof( HeaderType );

    for( const char* ptr = headerBegin; ptr < headerBegin + numberOfHeaderBytes; ++ptr )
    {
      output << *ptr;
    }

    for( const auto& compressedBlock : appendedData[iDataSet] )
    {
      for( auto ptr = compressedBlock.begin( ); ptr < compressedBlock.end( ); ++ptr )
      {
        output << *ptr;
      }
    } // for compressedBLock
  } // for iDataSet

  output << "\n";
}

inline void CompressedRawBinaryAppendedWriter::addHeaderAttributes( StringStringMap& attributes )
{
  attributes["header_type"] = dataTypeString<HeaderType>( );
  attributes["compressor"] = "vtkZLibDataCompressor";
}

inline void CompressedRawBinaryAppendedWriter::addDataAttributes( StringStringMap& attributes )
{
  attributes["format"] = "appended";
  attributes["offset"] = std::to_string( offset );
}

inline StringStringMap CompressedRawBinaryAppendedWriter::appendedAttributes( )
{
  return { { "encoding", "raw" } };
}

} // namespace vtu11

#endif // VTU11_ZLIBWRITER_IMPL_HPP

#ifndef VTU11_VTU11_IMPL_HPP
#define VTU11_VTU11_IMPL_HPP


#include <limits>

namespace vtu11
{
namespace detail
{

template<typename DataType, typename Writer> inline
StringStringMap writeDataSetHeader( Writer&& writer,
                                    const std::string& name,
                                    size_t ncomponents )
{
    StringStringMap attributes = { { "type", dataTypeString<DataType>( ) } };

    if( name != "" )
    {
        attributes["Name"] = name;
    }

    if( ncomponents > 1 )
    {
        attributes["NumberOfComponents"] = std::to_string( ncomponents );
    }

    writer.addDataAttributes( attributes );

    return attributes;
}

template<typename Writer, typename DataType> inline
void writeDataSet( Writer& writer,
                   std::ostream& output,
                   const std::string& name,
                   size_t ncomponents,
                   const std::vector<DataType>& data )
{
    auto attributes = writeDataSetHeader<DataType>( writer, name, ncomponents );

    if( attributes["format"] != "appended" )
    {
        ScopedXmlTag dataArrayTag( output, "DataArray", attributes );

        writer.writeData( output, data );
    }
    else
    {
        writeEmptyTag( output, "DataArray", attributes );

        writer.writeData( output, data );
    }
}

template<typename Writer> inline
void writeDataSets( const std::vector<DataSetInfo>& dataSetInfo,
                    const std::vector<DataSetData>& dataSetData,
                    std::ostream& output, Writer& writer, DataSetType type )
{
    for( size_t iDataset = 0; iDataset < dataSetInfo.size( ); ++iDataset )
    {
        const auto& metadata = dataSetInfo[iDataset];

        if( std::get<1>( metadata ) == type )
        {
            detail::writeDataSet( writer, output, std::get<0>( metadata ),
                std::get<2>( metadata ), dataSetData[iDataset] );
        }
    }
}

template<typename Writer> inline
void writeDataSetPVtuHeaders( const std::vector<DataSetInfo>& dataSetInfo,
                              std::ostream& output, Writer& writer, DataSetType type )
{
    for( size_t iDataset = 0; iDataset < dataSetInfo.size( ); ++iDataset )
    {
        const auto& metadata = dataSetInfo[iDataset];

        if( std::get<1>( metadata ) == type )
        {
            auto attributes = detail::writeDataSetHeader<double>( writer,
               std::get<0>( metadata ), std::get<2>( metadata ) );

            writeEmptyTag( output, "PDataArray", attributes );
        }
    }
}

template<typename Writer, typename Content> inline
void writeVTUFile( const std::string& filename,
                   const char* type,
                   Writer&& writer,
                   Content&& writeContent )
{
    std::ofstream output( filename, std::ios::binary );

    VTU11_CHECK( output.is_open( ), "Failed to open file \"" + filename + "\"" );

    std::vector<char> buffer( 32 * 1024 );

    output.rdbuf( )->pubsetbuf( buffer.data( ), static_cast<std::streamsize>( buffer.size( ) ) );

    output << "<?xml version=\"1.0\"?>\n";

    StringStringMap headerAttributes { { "byte_order",  endianness( ) },
                                       { "type"      ,  type          },
                                       { "version"   ,  "0.1"         } };

    writer.addHeaderAttributes( headerAttributes );

    {
        ScopedXmlTag vtkFileTag( output, "VTKFile", headerAttributes );

        writeContent( output );

    } // VTKFile

    output.close( );
}

template<typename MeshGenerator, typename Writer> inline
void writeVtu( const std::string& filename,
               MeshGenerator& mesh,
               const std::vector<DataSetInfo>& dataSetInfo,
               const std::vector<DataSetData>& dataSetData,
               Writer&& writer )
{
    detail::writeVTUFile( filename, "UnstructuredGrid", writer, [&]( std::ostream& output )
    {
        {
            ScopedXmlTag unstructuredGridFileTag( output, "UnstructuredGrid", { } );
            {
                ScopedXmlTag pieceTag( output, "Piece",
                {
                    { "NumberOfPoints", std::to_string( mesh.numberOfPoints( ) ) },
                    { "NumberOfCells" , std::to_string( mesh.numberOfCells( )  ) }

                } );

                {
                    ScopedXmlTag pointDataTag( output, "PointData", { } );

                    detail::writeDataSets( dataSetInfo, dataSetData,
                        output, writer, DataSetType::PointData );

                } // PointData

                {
                    ScopedXmlTag cellDataTag( output, "CellData", { } );

                    detail::writeDataSets( dataSetInfo, dataSetData,
                        output, writer, DataSetType::CellData );

                } // CellData

                {
                    ScopedXmlTag pointsTag( output, "Points", { } );

                    detail::writeDataSet( writer, output, "", 3, mesh.points( ) );

                } // Points

                {
                    ScopedXmlTag pointsTag( output, "Cells", { } );

                    detail::writeDataSet( writer, output, "connectivity", 1, mesh.connectivity( ) );
                    detail::writeDataSet( writer, output, "offsets", 1, mesh.offsets( ) );
                    detail::writeDataSet( writer, output, "types", 1, mesh.types( ) );

                } // Cells

            } // Piece
        } // UnstructuredGrid

        auto appendedAttributes = writer.appendedAttributes( );

        if( !appendedAttributes.empty( ) )
        {
            ScopedXmlTag appendedDataTag( output, "AppendedData", appendedAttributes );

            output << "_";

            writer.writeAppended( output );

        } // AppendedData

    } ); // writeVTUFile

} // writeVtu

} // namespace detail

template<typename MeshGenerator> inline
void writeVtu( const std::string& filename,
               MeshGenerator& mesh,
               const std::vector<DataSetInfo>& dataSetInfo,
               const std::vector<DataSetData>& dataSetData,
               const std::string& writeMode )
{
    auto mode = writeMode;

    std::transform( mode.begin( ), mode.end ( ), mode.begin( ), []( unsigned char c )
                    { return static_cast<unsigned char>( std::tolower( c ) ); } );

    if( mode == "ascii" )
    {
        detail::writeVtu( filename, mesh, dataSetInfo, dataSetData, AsciiWriter { } );
    }
    else if( mode == "base64inline" )
    {
        detail::writeVtu( filename, mesh, dataSetInfo, dataSetData, Base64BinaryWriter { } );
    }
    else if( mode == "base64appended" )
    {
        detail::writeVtu( filename, mesh, dataSetInfo, dataSetData, Base64BinaryAppendedWriter { } );
    }
    else if( mode == "rawbinary" )
    {
        detail::writeVtu( filename, mesh, dataSetInfo, dataSetData, RawBinaryAppendedWriter { } );
    }
    else if( mode == "rawbinarycompressed" )
    {
        #ifdef VTU11_ENABLE_ZLIB
            detail::writeVtu( filename, mesh, dataSetInfo, dataSetData, CompressedRawBinaryAppendedWriter { } );
        #else
            detail::writeVtu( filename, mesh, dataSetInfo, dataSetData, RawBinaryAppendedWriter { } );
        #endif
    }
    else
    {
        VTU11_THROW( "Invalid write mode: \"" + writeMode + "\"." );
    }

} // writeVtu

namespace detail
{

struct PVtuDummyWriter
{
    void addHeaderAttributes( StringStringMap& ) { }
    void addDataAttributes( StringStringMap& ) { }
};

} // detail

inline void writePVtu( const std::string& path,
                       const std::string& baseName,
                       const std::vector<DataSetInfo>& dataSetInfo,
                       const size_t numberOfFiles )
{
    auto directory = vtu11fs::path { path } / baseName;
    auto pvtufile = vtu11fs::path { path } / ( baseName + ".pvtu" );

    if( !vtu11fs::exists( directory ) )
    {
        vtu11fs::create_directories( directory );
    }

    detail::PVtuDummyWriter writer;

    detail::writeVTUFile( pvtufile.string( ), "PUnstructuredGrid", writer,
                          [&]( std::ostream& output )
    {
        std::string ghostLevel = "0"; // Hardcoded to be 0

        ScopedXmlTag pUnstructuredGridFileTag( output,
            "PUnstructuredGrid", { { "GhostLevel", ghostLevel } } );

        {
            ScopedXmlTag pPointDataTag( output, "PPointData", { } );

            detail::writeDataSetPVtuHeaders( dataSetInfo, output, writer, DataSetType::PointData );

        } // PPointData

        {
            ScopedXmlTag pCellDataTag( output, "PCellData", { } );

            detail::writeDataSetPVtuHeaders( dataSetInfo, output, writer, DataSetType::CellData );

        } // PCellData

        {
            ScopedXmlTag pPointsTag( output, "PPoints", { } );
            StringStringMap attributes = { { "type", dataTypeString<double>( ) }, { "NumberOfComponents", "3" } };

            writer.addDataAttributes( attributes );

            writeEmptyTag( output, "PDataArray", attributes );

        } // PPoints

        for( size_t nFiles = 0; nFiles < numberOfFiles; ++nFiles )
        {
            std::string pieceName = baseName + "/" + baseName + "_" + std::to_string( nFiles ) + ".vtu";

            writeEmptyTag( output, "Piece", { { "Source", pieceName } } );

        } // Pieces

    } ); // writeVTUFile

} // writePVtu

template<typename MeshGenerator> inline
void writePartition( const std::string& path,
                     const std::string& baseName,
                     MeshGenerator& mesh,
                     const std::vector<DataSetInfo>& dataSetInfo,
                     const std::vector<DataSetData>& dataSetData,
                     size_t fileId,
                     const std::string& writeMode )
{
    auto vtuname = baseName + "_" + std::to_string( fileId ) + ".vtu";

    auto fullname = vtu11fs::path { path } /
                    vtu11fs::path { baseName } /
                    vtu11fs::path { vtuname };

    writeVtu( fullname.string( ), mesh, dataSetInfo, dataSetData, writeMode );

} // writePartition

} // namespace vtu11

#endif // VTU11_VTU11_IMPL_HPP

namespace Kratos
{

namespace { // helpers namespace

void GetNodalCoordinates(
    const ModelPart& rModelPart,
    std::vector<double>& rCoordinates,
    const bool WriteDeformedConfiguration)
{
    // NOTE: also in MPI all nodes (local and ghost) have to be written, because
    // they might be needed by the elements/conditions due to the connectivity

    if (rCoordinates.size() != rModelPart.NumberOfNodes()*3) {
        rCoordinates.resize(rModelPart.NumberOfNodes()*3);
    }

    std::size_t index(0);

    // Write nodes TODO loops should be OMP
    if (WriteDeformedConfiguration) {
        for (const auto& r_node : rModelPart.Nodes()) {
            rCoordinates[index++] = r_node.X();
            rCoordinates[index++] = r_node.Y();
            rCoordinates[index++] = r_node.Z();
        }
    } else {
        for (const auto& r_node : rModelPart.Nodes()) {
            rCoordinates[index++] = r_node.X0();
            rCoordinates[index++] = r_node.Y0();
            rCoordinates[index++] = r_node.Z0();
        }
    }
}

template <typename TContainerType>
void GetCellInformation(
    const TContainerType& rContainer,
    const std::unordered_map<int, int>& rKratosIdToVtuId,
    std::vector<vtu11::VtkIndexType>& rConnectivities,
    std::vector<vtu11::VtkIndexType>& rOffsets,
    std::vector<vtu11::VtkCellType>& rTypes)
{
        // IMPORTANT: The map geo_type_vtk_cell_type_map is to be extended to support new geometries
    // NOTE: See https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    const std::map<GeometryData::KratosGeometryType, int> geo_type_vtk_cell_type_map = {
        { GeometryData::KratosGeometryType::Kratos_Point2D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Point3D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Line2D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Line3D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Triangle2D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Triangle3D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    10 },
        { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     12 },
        { GeometryData::KratosGeometryType::Kratos_Prism3D6,         13 },
        { GeometryData::KratosGeometryType::Kratos_Pyramid3D5,       14 },
        { GeometryData::KratosGeometryType::Kratos_Line2D3,          21 },
        { GeometryData::KratosGeometryType::Kratos_Line3D3,          21 },
        { GeometryData::KratosGeometryType::Kratos_Triangle2D6,      22 },
        { GeometryData::KratosGeometryType::Kratos_Triangle3D6,      22 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, 23 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, 23 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   24 }
//         { GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    25 } // NOTE: Quadratic hexahedra (20) requires a conversor, order does not coincide with VTK
    };

    const std::size_t constainer_size(rContainer.size());
    rConnectivities.resize(0);
    if (constainer_size > 0) { // this can happen when a domain in MPI has no local Elements/Consitions
        // reserving a guess for the necessary entries
        rConnectivities.reserve(constainer_size * rContainer.begin()->GetGeometry().PointsNumber());
    }

    if (rOffsets.size() != constainer_size) {
        rOffsets.resize(constainer_size);
    }

    if (rTypes.size() != constainer_size) {
        rTypes.resize(constainer_size);
    }

    std::size_t container_index(0);
    for (const auto& r_entity : rContainer) {
        const auto& r_geom = r_entity.GetGeometry();
        rOffsets[container_index] = r_geom.PointsNumber();
        const auto& r_kratos_cell = r_geom.GetGeometryType();
        if (geo_type_vtk_cell_type_map.count(r_kratos_cell) > 0) {
            rTypes[container_index++] = geo_type_vtk_cell_type_map.at(r_kratos_cell);
        } else {
            KRATOS_ERROR << "Modelpart contains elements or conditions with "
             << "geometries for which no VTK-output is implemented!" << std::endl
             << "Cell type: " << static_cast<int>(r_kratos_cell) << std::endl;
        }
        for (const auto& r_node : r_geom) {
            rConnectivities.push_back(rKratosIdToVtuId.at(r_node.Id()));
        }
    }

    std::partial_sum(rOffsets.begin(), rOffsets.end(), rOffsets.begin());
}

vtu11::Vtu11UnstructuredMesh GetMesh(
    const ModelPart& rModelPart,
    const bool WriteDeformedConfiguration,
    const std::unordered_map<int, int>& rKratosIdToVtuId,
    std::vector<double>& rNodalCoordinates,
    std::vector<vtu11::VtkIndexType>& rConnectivities,
    std::vector<vtu11::VtkIndexType>& rOffsets,
    std::vector<vtu11::VtkCellType>& rTypes)
{
    GetNodalCoordinates(rModelPart, rNodalCoordinates, WriteDeformedConfiguration);

    const int num_elements = rModelPart.GetCommunicator().GlobalNumberOfElements();
    const int num_conditions = rModelPart.GetCommunicator().GlobalNumberOfConditions();

    if (num_elements > 0) {
        GetCellInformation(rModelPart.GetCommunicator().LocalMesh().Elements(), rKratosIdToVtuId, rConnectivities, rOffsets, rTypes);
    } else if (num_conditions > 0) {
        GetCellInformation(rModelPart.GetCommunicator().LocalMesh().Conditions(), rKratosIdToVtuId, rConnectivities, rOffsets, rTypes);
    }

    return {rNodalCoordinates, rConnectivities, rOffsets, rTypes};
}

// TODO this should use the AuxiliarModelPartUtilities
std::vector<double> GetScalarData(
    ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    KRATOS_TRY;

    rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    std::vector<double> results(rModelPart.NumberOfNodes());
    std::size_t counter = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        results[counter++] = r_node.FastGetSolutionStepValue(rVariable);
    }
    return results;

    KRATOS_CATCH("VTU GetScalarData");
}

std::vector<double> GetVectorData(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable)
{
    KRATOS_TRY;

    rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    std::vector<double> results(rModelPart.NumberOfNodes()*3);
    std::size_t counter = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const auto& r_vals = r_node.FastGetSolutionStepValue(rVariable);
        results[counter++] = r_vals[0];
        results[counter++] = r_vals[1];
        results[counter++] = r_vals[2];
    }
    return results;

    KRATOS_CATCH("VTU GetVectorData");
}
std::vector<double> GetScalarDataNonHist(
    ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    KRATOS_TRY;

    rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
    std::vector<double> results(rModelPart.NumberOfNodes());
    std::size_t counter = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        results[counter++] = r_node.GetValue(rVariable);
    }
    return results;

    KRATOS_CATCH("VTU GetScalarDataNonHist");
}

std::vector<double> GetVectorDataNonHist(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable)
{
    KRATOS_TRY;

    rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
    std::vector<double> results(rModelPart.NumberOfNodes()*3);
    std::size_t counter = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const auto& r_vals = r_node.GetValue(rVariable);
        results[counter++] = r_vals[0];
        results[counter++] = r_vals[1];
        results[counter++] = r_vals[2];
    }
    return results;

    KRATOS_CATCH("VTU GetVectorDataNonHist");
}

} // helpers namespace

Parameters VtuOutput::GetDefaultParameters()
{
    // IMPORTANT: when "output_control_type" is "time", then paraview will not be able to group them
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "file_format"                                 : "binary_raw_compressed",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_interval"                             : 1.0,
        "output_sub_model_parts"                      : false,
        "output_path"                                 : "VTU_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "write_deformed_configuration"                : false,
        "write_ids"                                   : false,
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : []
    })" );

    return default_parameters;
}


VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mOutputSettings(ThisParameters)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mOutputSettings.ValidateAndAssignDefaults(default_parameters);

    // Initialize other variables
    const std::string file_format = mOutputSettings["file_format"].GetString();
    if (file_format == "ascii") {
        mFileFormat = "ascii";
    } else if (file_format == "binary_raw") {
        mFileFormat = "rawbinary";
    } else if (file_format == "binary_raw_compressed") {
        mFileFormat = "rawbinarycompressed";
    } else if (file_format == "binary_base64") {
        mFileFormat = "base64inline";
    } else if (file_format == "binary_base64_appended") {
        mFileFormat = "base64appended";
    } else {
        KRATOS_ERROR << "Option for \"file_format\": " << file_format << " not recognised!\n Possible output formats options are: \"ascii\", \"binary_raw\", \"binary_raw_compressed\", \"binary_base64\", \"binary_base64_appended\"" << std::endl;
    }

    /*
    // Adding GP variables to nodal data variables list
    if(mOutputSettings["gauss_point_variables_extrapolated_to_nodes"].size() > 0) {
        Parameters gauss_intergration_param_non_hist = Parameters(R"(
        {
            "echo_level"                 : 0,
            "area_average"               : true,
            "average_variable"           : "NODAL_AREA",
            "list_of_variables"          : [],
            "extrapolate_non_historical" : true
        })");

        gauss_intergration_param_non_hist.SetValue("list_of_variables", mOutputSettings["gauss_point_variables_extrapolated_to_nodes"]);

        for(auto const& gauss_var : mOutputSettings["gauss_point_variables_extrapolated_to_nodes"])
            mOutputSettings["nodal_data_value_variables"].Append(gauss_var);

        // Making the gauss point to nodes process if any gauss point result is requested for
        mpGaussToNodesProcess = Kratos::make_unique<IntegrationValuesExtrapolationToNodesProcess>(rModelPart, gauss_intergration_param_non_hist);
    }*/

    const int num_elements = rModelPart.GetCommunicator().GlobalNumberOfElements();
    const int num_conditions = rModelPart.GetCommunicator().GlobalNumberOfConditions();

    KRATOS_WARNING_IF("VtuOutput", num_elements > 0 && num_conditions > 0) << "Modelpart \"" << rModelPart.Name() << "\" has both elements and conditions.\nGiving precedence to elements and writing only elements!" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::PrepareGaussPointResults()
{
    /*if(mOutputSettings["gauss_point_variables_extrapolated_to_nodes"].size() > 0){
        mpGaussToNodesProcess->Execute();
    }*/
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::PrintOutput(const std::string& rOutputFilename)
{
    // For Gauss point results
    PrepareGaussPointResults();

    // For whole model part
    WriteModelPartToFile(mrModelPart, false, rOutputFilename);

    // For sub model parts
    const bool print_sub_model_parts = mOutputSettings["output_sub_model_parts"].GetBool();
    if (print_sub_model_parts) {
        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {

            const auto& r_local_mesh = mrModelPart.GetCommunicator().LocalMesh();
            const auto& r_data_comm = mrModelPart.GetCommunicator().GetDataCommunicator();

            const int num_nodes = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfNodes()));
            const int num_elements = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
            const int num_conditions = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

            if (num_nodes == 0 && (num_elements != 0 || num_conditions != 0)) {
                WriteModelPartWithoutNodesToFile(r_sub_model_part, rOutputFilename);
            } else if (num_nodes != 0) {
                WriteModelPartToFile(r_sub_model_part, true, rOutputFilename);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::WriteModelPartToFile(ModelPart& rModelPart, const bool IsSubModelPart, const std::string& rOutputFilename)
{
    Initialize(rModelPart);

    std::string output_path, output_file_name;
    std::tie(output_path, output_file_name) = GetOutputFileName(rModelPart, IsSubModelPart, rOutputFilename);

    std::vector<double> nodal_coordinates;
    std::vector<vtu11::VtkIndexType> connectivity;
    std::vector<vtu11::VtkIndexType> offsets;
    std::vector<vtu11::VtkCellType> types;

    std::vector<vtu11::DataSetInfo> data_set_info;
    std::vector<vtu11::DataSetData> data_set_data;

    for (const auto& r_var_name : mOutputSettings["nodal_solution_step_data_variables"].GetStringArray()) {
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_var_name)) {
            const auto& r_var = KratosComponents<Variable<array_1d<double, 3>>>::Get(r_var_name);
            data_set_info.push_back({r_var_name, vtu11::DataSetType::PointData, 3});
            data_set_data.push_back(GetVectorData(rModelPart, r_var));

        } else if (KratosComponents<Variable<double>>::Has(r_var_name)) {
            const auto& r_var = KratosComponents<Variable<double>>::Get(r_var_name);
            data_set_info.push_back({r_var_name, vtu11::DataSetType::PointData, 1});
            data_set_data.push_back(GetScalarData(rModelPart, r_var));

        } else {
            KRATOS_WARNING("VtuOutput") << "Variable " << r_var_name << " is not suitable for vtu output and is skipped" << std::endl;
        }
    }

    for (const auto& r_var_name : mOutputSettings["nodal_data_value_variables"].GetStringArray()) {
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_var_name)) {
            const auto& r_var = KratosComponents<Variable<array_1d<double, 3>>>::Get(r_var_name);
            data_set_info.push_back({r_var_name, vtu11::DataSetType::PointData, 3});
            data_set_data.push_back(GetVectorDataNonHist(rModelPart, r_var));

        } else if (KratosComponents<Variable<double>>::Has(r_var_name)) {
            const auto& r_var = KratosComponents<Variable<double>>::Get(r_var_name);
            data_set_info.push_back({r_var_name, vtu11::DataSetType::PointData, 1});
            data_set_data.push_back(GetScalarDataNonHist(rModelPart, r_var));

        } else {
            KRATOS_WARNING("VtuOutput") << "Variable " << r_var_name << " is not suitable for vtu output and is skipped" << std::endl;
        }
    }

    vtu11::Vtu11UnstructuredMesh vtu11_mesh = GetMesh(
        rModelPart,
        mOutputSettings["write_deformed_configuration"].GetBool(),
        mKratosIdToVtkId,
        nodal_coordinates,
        connectivity,
        offsets,
        types);

    if (mrModelPart.IsDistributed()) {
        const auto& r_data_comm = mrModelPart.GetCommunicator().GetDataCommunicator();
        FilesystemExtensions::MPISafeCreateDirectories(Kratos::FilesystemExtensions::JoinPaths({output_path, output_file_name})); // to be sure that the dir exists (whould also be done in writePVtu)

        if (r_data_comm.Rank() == 0) {
            vtu11::writePVtu(
                output_path,
                output_file_name,
                data_set_info,
                r_data_comm.Size() );
        }
        r_data_comm.Barrier();

        vtu11::writePartition(
            output_path,
            output_file_name,
            vtu11_mesh,
            data_set_info,
            data_set_data,
            r_data_comm.Rank(),
            mFileFormat);

    } else {
        vtu11::writeVtu(
            Kratos::FilesystemExtensions::JoinPaths({output_path, output_file_name})+".vtu",
            vtu11_mesh,
            data_set_info,
            data_set_data,
            mFileFormat);
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::tuple<std::string, std::string> VtuOutput::GetOutputFileName(const ModelPart& rModelPart, const bool IsSubModelPart, const std::string& rOutputFilename) const
{
    std::string output_file_name = "";

    if (rOutputFilename != "") { // user specified file name externally
        output_file_name = rOutputFilename;
    } else {
        std::string model_part_name = rModelPart.Name();
        if (IsSubModelPart) {
            model_part_name = rModelPart.GetParentModelPart().Name() + "_" + model_part_name;
        }

        std::string label;
        std::stringstream ss;
        const std::string output_control = mOutputSettings["output_control_type"].GetString();
        if (output_control == "step") {
            ss << std::fixed << std::setprecision(mDefaultPrecision)<< std::setfill('0')
            << rModelPart.GetProcessInfo()[STEP];
            label = ss.str();
        } else if(output_control == "time") {
            ss << std::fixed << std::setprecision(mDefaultPrecision) << std::setfill('0')
            << rModelPart.GetProcessInfo()[TIME];
            label = ss.str();
        } else {
            KRATOS_ERROR << "Option for \"output_control_type\": " << output_control
                <<" not recognised!\nPossible output_control_type options "
                << "are: \"step\", \"time\"" << std::endl;
        }

        const std::string& r_custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
        const std::string& r_custom_name_postfix = mOutputSettings["custom_name_postfix"].GetString();
        output_file_name += r_custom_name_prefix + model_part_name + r_custom_name_postfix + "_" + label;
    }

    std::string output_path = "";
    if (mOutputSettings["save_output_files_in_folder"].GetBool()) {
        output_path = mOutputSettings["output_path"].GetString();

        // Create folder if it doesn't exist before
        FilesystemExtensions::MPISafeCreateDirectories(output_path);
    }

    return std::make_tuple(output_path, output_file_name);
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::Initialize(const ModelPart& rModelPart)
{
    CreateMapFromKratosIdToVTKId(rModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::CreateMapFromKratosIdToVTKId(const ModelPart& rModelPart)
{
    int vtk_id = 0;
    for(const auto& r_node : rModelPart.Nodes()) {
        mKratosIdToVtkId[r_node.Id()] = vtk_id++;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool VtuOutput::IsCompatibleVariable(const std::string& rVariableName) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::WriteNodalResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    // NOTE: also in MPI all nodes (local and ghost) have to be written, because
    // they might be needed by the elements/conditions due to the connectivity
    // Paraview needs a result on every node, therefore all results are written
    // this is why the synchronization is necessary

    // TODO perform synchronization of nodal results at the same time to
    // improve performance in MPI

    // write nodal results header
    Parameters nodal_solution_step_results = mOutputSettings["nodal_solution_step_data_variables"];
    Parameters nodal_variable_data_results = mOutputSettings["nodal_data_value_variables"];
    Parameters nodal_flags = mOutputSettings["nodal_flags"];
    const bool write_ids = mOutputSettings["write_ids"].GetBool();

    // Checking nodal_solution_step_results
    std::size_t counter_nodal_solution_step_results = 0;
    for (IndexType entry = 0; entry < nodal_solution_step_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_solution_step_results[entry].GetString();
        if (IsCompatibleVariable(r_nodal_result_name)) ++counter_nodal_solution_step_results;
    }

    // Checking nodal_variable_data_results
    std::size_t counter_nodal_variable_data_results = 0;
    for (IndexType entry = 0; entry < nodal_variable_data_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_variable_data_results[entry].GetString();
        if (IsCompatibleVariable(r_nodal_result_name)) ++counter_nodal_variable_data_results;
    }

    rFileStream << "POINT_DATA " << rModelPart.NumberOfNodes() << "\n";
    rFileStream << "FIELD FieldData " << counter_nodal_solution_step_results + counter_nodal_variable_data_results + nodal_flags.size() + (write_ids ? 1 : 0)  << "\n";

    // Writing nodal_solution_step_results
    for (IndexType entry = 0; entry < nodal_solution_step_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_solution_step_results[entry].GetString();
        WriteNodalContainerResults(r_nodal_result_name, rModelPart.Nodes(), true, rFileStream);
    }

    // Writing nodal_variable_data_results
    for (IndexType entry = 0; entry < nodal_variable_data_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_variable_data_results[entry].GetString();
        WriteNodalContainerResults(r_nodal_result_name, rModelPart.Nodes(), false, rFileStream);
    }

    // Writing nodal_flags
    if (nodal_flags.size() > 0) {
        mrModelPart.GetCommunicator().SynchronizeNodalFlags();
    }
    for (IndexType entry = 0; entry < nodal_flags.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_flags[entry].GetString();
        const Flags flag = KratosComponents<Flags>::Get(r_nodal_result_name);
        WriteFlagContainerVariable(rModelPart.Nodes(), flag, r_nodal_result_name, rFileStream);
    }

    // If we write ids
    if (write_ids) {
        WriteIdsToFile(rModelPart.Nodes(), "KRATOS_NODE_ID", rFileStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::WriteElementResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters element_data_value_variables = mOutputSettings["element_data_value_variables"];
    Parameters element_flags = mOutputSettings["element_flags"];
    Parameters gauss_point_variables_in_elements = mOutputSettings["gauss_point_variables_in_elements"];

    // Checking element_data_value_variables
    std::size_t counter_element_data_value_variables = 0;
    for (IndexType entry = 0; entry < element_data_value_variables.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_element_result_name = element_data_value_variables[entry].GetString();
        if (IsCompatibleVariable(r_element_result_name)) ++counter_element_data_value_variables;
    }

    // Checking gauss_point_variables_in_elements
    std::size_t counter_gauss_point_variables_in_elements = 0;
    for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_element_result_name = gauss_point_variables_in_elements[entry].GetString();
        if (IsCompatibleVariable(r_element_result_name)) ++counter_gauss_point_variables_in_elements;
    }

    const int num_elements = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));

    if (num_elements > 0) {
        // write cells header
        rFileStream << "CELL_DATA " << r_local_mesh.NumberOfElements() << "\n";
        const bool write_ids = mOutputSettings["write_ids"].GetBool();
        rFileStream << "FIELD FieldData " << counter_element_data_value_variables + element_flags.size() + (write_ids ? 2 : 0) + counter_gauss_point_variables_in_elements << "\n";
        for (IndexType entry = 0; entry < element_data_value_variables.size(); ++entry) {
            const std::string& r_element_result_name = element_data_value_variables[entry].GetString();
            WriteGeometricalContainerResults(r_element_result_name,r_local_mesh.Elements(),rFileStream);
        }

        // Writing element_flags
        if (element_flags.size() > 0) {
            mrModelPart.GetCommunicator().SynchronizeElementalFlags();
        }
        for (IndexType entry = 0; entry < element_flags.size(); ++entry) {
            // Write elemental flags results variable header
            const std::string& r_element_result_name = element_flags[entry].GetString();
            const Flags flag = KratosComponents<Flags>::Get(r_element_result_name);
            WriteFlagContainerVariable(r_local_mesh.Elements(), flag, r_element_result_name, rFileStream);
        }

        // If we write ids
        if (write_ids) {
            WritePropertiesIdsToFile(r_local_mesh.Elements(), rFileStream);
            WriteIdsToFile(r_local_mesh.Elements(), "KRATOS_ELEMENT_ID", rFileStream);
        }

        // Direct write GP values
        for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
            const std::string& r_element_result_name = gauss_point_variables_in_elements[entry].GetString();
            WriteGeometricalContainerIntegrationResults(r_element_result_name,r_local_mesh.Elements(),rFileStream);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::WriteConditionResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters condition_results = mOutputSettings["condition_data_value_variables"];
    Parameters condition_flags = mOutputSettings["condition_flags"];
    Parameters gauss_point_variables_in_elements = mOutputSettings["gauss_point_variables_in_elements"];

    // Checking condition_results
    std::size_t counter_condition_results = 0;
    for (IndexType entry = 0; entry < condition_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_condition_result_name = condition_results[entry].GetString();
        if (IsCompatibleVariable(r_condition_result_name)) ++counter_condition_results;
    }

    // Checking gauss_point_variables_in_elements
    std::size_t counter_gauss_point_variables_in_elements = 0;
    for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_element_result_name = gauss_point_variables_in_elements[entry].GetString();
        if (IsCompatibleVariable(r_element_result_name)) ++counter_gauss_point_variables_in_elements;
    }

    const int num_elements = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
    const int num_conditions = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(static_cast<int>(r_local_mesh.NumberOfConditions())));

    if (num_elements == 0 && num_conditions > 0) { // TODO: Can we have conditions and elements at the same time?
        // Write cells header
        rFileStream << "CELL_DATA " << r_local_mesh.NumberOfConditions() << "\n";
        const bool write_ids = mOutputSettings["write_ids"].GetBool();
        rFileStream << "FIELD FieldData " << counter_condition_results + condition_flags.size() + (write_ids ? 2 : 0) + counter_gauss_point_variables_in_elements << "\n";
        for (IndexType entry = 0; entry < condition_results.size(); ++entry) {
            const std::string& r_condition_result_name = condition_results[entry].GetString();
            WriteGeometricalContainerResults(r_condition_result_name,r_local_mesh.Conditions(),rFileStream);
        }

        // Writing condition_flags
        if (condition_flags.size() > 0) {
            // mrModelPart.GetCommunicator().SynchronizeConditionFlags(); // TODO implement this if at some point ghost-conditions are used
        }
        for (IndexType entry = 0; entry < condition_flags.size(); ++entry) {
            // Write conditional flags results variable header
            const std::string& r_condition_result_name = condition_flags[entry].GetString();
            const Flags flag = KratosComponents<Flags>::Get(r_condition_result_name);
            WriteFlagContainerVariable(r_local_mesh.Conditions(), flag, r_condition_result_name, rFileStream);
        }

        // If we write properties_id
        if (write_ids) {
            WritePropertiesIdsToFile(r_local_mesh.Conditions(), rFileStream);
            WriteIdsToFile(r_local_mesh.Conditions(), "KRATOS_CONDITION_ID", rFileStream);
        }

        // Direct write GP values
        for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
            const std::string& r_condition_result_name = gauss_point_variables_in_elements[entry].GetString();
            WriteGeometricalContainerIntegrationResults(r_condition_result_name,r_local_mesh.Conditions(),rFileStream);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::WriteNodalContainerResults(
    const std::string& rVariableName,
    const ModelPart::NodesContainerType& rNodes,
    const bool IsHistoricalValue,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<bool>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else {
        KRATOS_WARNING_ONCE(rVariableName) << mrModelPart.GetCommunicator().GetDataCommunicator() << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtuOutput, skipping it" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtuOutput::WriteGeometricalContainerResults(
    const std::string& rVariableName,
    const TContainerType& rContainer,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<bool>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<Flags>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Flags>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else {
        KRATOS_WARNING_ONCE(rVariableName) << mrModelPart.GetCommunicator().GetDataCommunicator() << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtuOutput, skipping it" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtuOutput::WriteGeometricalContainerIntegrationResults(
    const std::string& rVariableName,
    const TContainerType& rContainer,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteIntegrationScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<bool>>::Get(rVariableName);
        WriteIntegrationScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteIntegrationScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteIntegrationVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteIntegrationVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteIntegrationVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else {
        KRATOS_WARNING_ONCE(rVariableName) << mrModelPart.GetCommunicator().GetDataCommunicator() << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtuOutput, skipping it" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVarType>
void VtuOutput::WriteNodalScalarValues(
    const ModelPart::NodesContainerType& rNodes,
    const TVarType& rVariable,
    const bool IsHistoricalValue,
    std::ofstream& rFileStream) const
{
    if (IsHistoricalValue) {
        mrModelPart.GetCommunicator().SynchronizeVariable(rVariable);
        WriteScalarSolutionStepVariable(rNodes, rVariable, rFileStream);
    } else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
        WriteScalarContainerVariable(rNodes, rVariable, rFileStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVarType>
void VtuOutput::WriteNodalVectorValues(
    const ModelPart::NodesContainerType& rNodes,
    const TVarType& rVariable,
    const bool IsHistoricalValue,
    std::ofstream& rFileStream) const
{
    if (IsHistoricalValue) {
        mrModelPart.GetCommunicator().SynchronizeVariable(rVariable);
        WriteVectorSolutionStepVariable(rNodes, rVariable, rFileStream);
    } else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
        WriteVectorContainerVariable(rNodes, rVariable, rFileStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtuOutput::WriteScalarSolutionStepVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    rFileStream << rVariable.Name() << " 1 "
                << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.FastGetSolutionStepValue(rVariable);
        WriteScalarDataToFile((float)r_result, rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtuOutput::WriteVectorSolutionStepVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    if (rContainer.size() == 0) {
        return;
    }

    const int res_size = static_cast<int>((rContainer.begin()->FastGetSolutionStepValue(rVariable)).size());

    rFileStream << rVariable.Name() << " " << res_size
                << " " << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.FastGetSolutionStepValue(rVariable);
        WriteVectorDataToFile(r_result, rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtuOutput::WriteFlagContainerVariable(
    const TContainerType& rContainer,
    const Flags Flag,
    const std::string& rFlagName,
    std::ofstream& rFileStream) const
{
    rFileStream << rFlagName << " 1 "
                << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const float result = r_entity.IsDefined(Flag) ? float(r_entity.Is(Flag)) : -1.0;
        WriteScalarDataToFile(result, rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtuOutput::WriteScalarContainerVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    rFileStream << rVariable.Name() << " 1 "
                << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const double result = r_entity.GetValue(rVariable);
        WriteScalarDataToFile((float)result, rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtuOutput::WriteIntegrationScalarContainerVariable(
    const TContainerType& rContainer,
    const Variable<TVarType>& rVariable,
    std::ofstream& rFileStream) const
{
    rFileStream << rVariable.Name() << " 1 "
                << rContainer.size() << "  float\n";

    // Auxiliar values
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    auto& r_this_geometry_begin = (rContainer.begin())->GetGeometry();
    const GeometryData::IntegrationMethod this_integration_method = (rContainer.begin())->GetIntegrationMethod();
    const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
    const SizeType integration_points_number = r_integration_points.size();

    double aux_value;
    for (auto& r_entity : rContainer) { // TODO: CalculateOnIntegrationPoints should be const methods
        aux_value = 0.0;
        std::vector<TVarType> aux_result(integration_points_number);
        r_entity.CalculateOnIntegrationPoints(rVariable, aux_result, r_process_info);
        for (const double value : aux_result) {
            aux_value += value;
        }
        aux_value /= static_cast<double>(integration_points_number);
        WriteScalarDataToFile((float)aux_value, rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtuOutput::WriteVectorContainerVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    if (rContainer.size() == 0) {
        return;
    }

    const int res_size = static_cast<int>((rContainer.begin()->GetValue(rVariable)).size());

    rFileStream << rVariable.Name() << " " << res_size << " " << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.GetValue(rVariable);
        WriteVectorDataToFile(r_result, rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtuOutput::WriteIntegrationVectorContainerVariable(
    const TContainerType& rContainer,
    const Variable<TVarType>& rVariable,
    std::ofstream& rFileStream) const
{
    if (rContainer.size() == 0) {
        return;
    }

    // determining size of results
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    std::vector<TVarType> tmp_result;
    rContainer.begin()->CalculateOnIntegrationPoints(rVariable, tmp_result, r_process_info);
    const int res_size = tmp_result[0].size();

    rFileStream << rVariable.Name() << " " << res_size << " " << rContainer.size() << "  float\n";

    // Auxiliar values
    auto& r_this_geometry_begin = (rContainer.begin())->GetGeometry();
    const GeometryData::IntegrationMethod this_integration_method = (rContainer.begin())->GetIntegrationMethod();
    const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
    const SizeType integration_points_number = r_integration_points.size();

    TVarType aux_value;
    for (auto& r_entity : rContainer) { // TODO: CalculateOnIntegrationPoints should be const methods
        aux_value = ZeroVector(res_size);
        std::vector<TVarType> aux_result(integration_points_number);
        r_entity.CalculateOnIntegrationPoints(rVariable, aux_result, r_process_info);
        for (const TVarType& r_value : aux_result) {
            noalias(aux_value) += r_value;
        }
        aux_value /= static_cast<double>(integration_points_number);
        WriteVectorDataToFile(aux_value, rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::ForceBigEndian(unsigned char* pBytes) const
{
    if (mShouldSwap) {
        unsigned char tmp = pBytes[0];
        pBytes[0] = pBytes[3];
        pBytes[3] = tmp;
        tmp = pBytes[1];
        pBytes[1] = pBytes[2];
        pBytes[2] = tmp;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtuOutput::WritePropertiesIdsToFile(
    const TContainerType& rContainer,
    std::ofstream& rFileStream) const
{
    rFileStream << "PROPERTIES_ID" << " 1 "
                << rContainer.size() << "  int\n";

    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.GetProperties().Id(), rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtuOutput::WriteIdsToFile(
    const TContainerType& rContainer,
    const std::string DataName,
    std::ofstream& rFileStream) const
{
    rFileStream << DataName << " 1 "
                << rContainer.size() << "  int\n";

    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.Id(), rFileStream);
        if (mFileFormatLegacy == VtuOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}


/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::WriteModelPartWithoutNodesToFile(ModelPart& rModelPart, const std::string& rOutputFilename)
{
    // Getting model and creating auxiliar model part
    auto& r_model = mrModelPart.GetModel();
    const std::string& r_name_model_part = rModelPart.Name();
    auto& r_auxiliar_model_part = r_model.CreateModelPart("AUXILIAR_" + r_name_model_part);

    // Tranfering entities of the submodelpart
    FastTransferBetweenModelPartsProcess(r_auxiliar_model_part, rModelPart).Execute();

    // Tranfering nodes from root model part
    FastTransferBetweenModelPartsProcess(r_auxiliar_model_part, mrModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::NODES).Execute();

    // Marking to remove the nodes
    for (auto& r_node : r_auxiliar_model_part.Nodes()) {
        r_node.Set(TO_ERASE, true);
    }

    // Checking nodes from conditions
    for (auto& r_cond : r_auxiliar_model_part.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        for (auto& r_node : r_geometry) {
            r_node.Set(TO_ERASE, false);
        }
    }

    // Checking nodes from elements
    for (auto& r_elem : r_auxiliar_model_part.Elements()) {
        auto& r_geometry = r_elem.GetGeometry();
        for (auto& r_node : r_geometry) {
            r_node.Set(TO_ERASE, false);
        }
    }

    // Removing unused nodes
    r_auxiliar_model_part.RemoveNodes(TO_ERASE);

    // Actually writing the
    WriteModelPartToFile(r_auxiliar_model_part, true, rOutputFilename);

    // Deletin auxiliar modek part
    r_model.DeleteModelPart("AUXILIAR_" + r_name_model_part);
}

} // namespace Kratos
