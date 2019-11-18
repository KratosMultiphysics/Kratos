//          __        ____ ____
// ___  ___/  |_ __ _/_   /_   |
// \  \/ /\   __\  |  \   ||   |
//  \   /  |  | |  |  /   ||   |
//   \_/   |__| |____/|___||___|
//
//  License: BSD License ; see LICENSE
//

#ifndef VTU11_ZLIBWRITER_IMPL_HPP

#include "zlib.h"
#include "inc/utilities.hpp"

#include <iostream> // remove

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
