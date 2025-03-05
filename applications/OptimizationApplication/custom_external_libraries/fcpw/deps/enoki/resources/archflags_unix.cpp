#include <iostream>

int main(int argc, char *argv[]) {
#if defined(__AVX512DQ__)
    std::cout << "-march=skx" << std::endl;
#elif defined(__AVX512ER__)
    std::cout << "-march=knl" << std::endl;
#elif defined(__AVX2__)
    std::cout << "-mavx2" << std::endl;
#elif defined(__AVX__)
    std::cout << "-mavx" << std::endl;
#elif defined(__SSE4_2__)
    std::cout << "-msse4.2" << std::endl;
#elif defined(__aarch64__)
    std::cout << "-march=armv8-a+simd -mtune=cortex-a53" << std::endl;
#elif defined(__arm__)
    std::cout << "-march=armv7-a -mtune=cortex-a7 -mfpu=neon-vfpv4 -mfloat-abi=hard -mfp16-format=ieee" << std::endl;
#endif
    return 0;
}
