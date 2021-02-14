//          __        ____ ____
// ___  ___/  |_ __ _/_   /_   |
// \  \/ /\   __\  |  \   ||   |
//  \   /  |  | |  |  /   ||   |
//   \_/   |__| |____/|___||___|
//
//  License: BSD License ; see LICENSE
//

#ifndef VTU11_BASE64_IMPL_HPP
#define VTU11_BASE64_IMPL_HPP

#include <iostream>
#include <array>

namespace vtu11
{

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
    char byte = *( reinterpret_cast<const char*>(&(*it)) + byteIndex++ );

    if( byteIndex == size )
    {
      it++;
      byteIndex = 0;
    }

    return byte;
  };

  auto encodeTriplet = [&]( std::array<char, 3> bytes, size_t padding )
  {
    char tmp[5] = { base64Map[ ( bytes[0] & 0xfc ) >> 2 ],
                    base64Map[ ( ( bytes[0] & 0x03 ) << 4) + ( ( bytes[1] & 0xf0 ) >> 4 ) ],
                    base64Map[ ( ( bytes[1] & 0x0f ) << 2) + ( ( bytes[2] & 0xc0 ) >> 6 ) ],
                    base64Map[ bytes[2] & 0x3f ],
                    '\0' };

    std::fill( tmp + 4 - padding, tmp + 4, '=' );

    result += tmp;
  };

  // in steps of 3
  for( size_t i = 0; i < rawBytes / 3; ++i )
  {
    encodeTriplet( { next( ), next( ), next( ) }, 0 );
  }

  // cleanup
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

// http://www.cplusplus.com/forum/beginner/51572/

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


}// namespace vtu11

#endif // VTU11_BASE64_IMPL_HPP
