#!/bin/bash

# Set compiler
export CC=${CC:-gcc}
export CXX=${CXX:-g++}

# Configure
cmake                                                        \
-DCMAKE_BUILD_TYPE=Release                                   \
-DCMAKE_C_COMPILER=${CC}                                     \
-DCMAKE_CXX_COMPILER=${CXX}                                  \
-DCMAKE_C_FLAGS="-O3 -march=native -mtune=native -fopenmp"   \
-DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native -fopenmp" \
-DUSE_SCOTCH=OFF                                             \
-DLIBMMG_SHARED=ON                                           \
-DLIBMMG_STATIC=OFF                                          \
-DLIBMMGS_SHARED=ON                                          \
-DLIBMMGS_STATIC=OFF                                         \
-DLIBMMG2D_SHARED=ON                                         \
-DLIBMMG2D_STATIC=OFF                                        \
-DLIBMMG3D_SHARED=ON                                         \
-DLIBMMG3D_STATIC=OFF                                        \

# Build
cmake --build $(pwd) -- -j$(nproc)
