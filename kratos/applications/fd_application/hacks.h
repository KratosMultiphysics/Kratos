// BitHacks and other shady things
#include <sys/types.h>

uint interleave64(const uint &IX, const uint &IY) {

    // Strip the first half of the numbers
    uint X = IX & ( 0x00000000FFFFFFFF );
    uint Y = IY & ( 0x00000000FFFFFFFF );

    X = ( X | ( X << 16 ) ) & 0x0000FFFF0000FFFF;
    X = ( X | ( X << 8  ) ) & 0x00FF00FF00FF00FF;
    X = ( X | ( X << 4  ) ) & 0x0F0F0F0F0F0F0F0F;
    X = ( X | ( X << 2  ) ) & 0x3333333333333333;
    X = ( X | ( X << 1  ) ) & 0x5555555555555555;

    Y = ( Y | ( Y << 16 ) ) & 0x0000FFFF0000FFFF;
    Y = ( Y | ( Y << 8  ) ) & 0x00FF00FF00FF00FF;
    Y = ( Y | ( Y << 4  ) ) & 0x0F0F0F0F0F0F0F0F;
    Y = ( Y | ( Y << 2  ) ) & 0x3333333333333333;
    Y = ( Y | ( Y << 1  ) ) & 0x5555555555555555;

    return X | ( Y << 1 );
}

uint interleave64(const uint &IX, const uint &IY, const uint &IZ) {

    // Strip the first half of the numbers
    uint X = IX & ( 0x000000000000FFFF );
    uint Y = IY & ( 0x000000000000FFFF );
    uint Z = IZ & ( 0x000000000000FFFF );

    X = ( X | ( X << 16 ) ) & 0x00FF0000FF0000FF;
    X = ( X | ( X << 8  ) ) & 0xF00F00F00F00F00F;
    X = ( X | ( X << 4  ) ) & 0x30C30C30C30C30C3;
    X = ( X | ( X << 2  ) ) & 0x0249249249249249;

    Y = ( Y | ( Y << 16 ) ) & 0x00FF0000FF0000FF;
    Y = ( Y | ( Y << 8  ) ) & 0xF00F00F00F00F00F;
    Y = ( Y | ( Y << 4  ) ) & 0x30C30C30C30C30C3;
    Y = ( Y | ( Y << 2  ) ) & 0x0249249249249249;

    Z = ( Z | ( Z << 16 ) ) & 0x00FF0000FF0000FF;
    Z = ( Z | ( Z << 8  ) ) & 0xF00F00F00F00F00F;
    Z = ( Z | ( Z << 4  ) ) & 0x30C30C30C30C30C3;
    Z = ( Z | ( Z << 2  ) ) & 0x0249249249249249;

    return ( X | ( Y << 1 ) | ( Z << 2 ) ) & 0x0000FFFFFFFFFFFF; 
}