#!/bin/bash

# Set compiler
export CC=${CC:-gcc}
export CXX=${CXX:-g++}

# Configure
cmake ..                                                     \
-DCMAKE_BUILD_TYPE=Release                                   \
-DCMAKE_C_COMPILER=${CC}                                     \
-DCMAKE_CXX_COMPILER=${CXX}                                  \
-DCMAKE_C_FLAGS="-O3 -march=native -mtune=native -fopenmp"   \
-DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native -fopenmp" \
-DUSE_SCOTCH=OFF                                             \
-DUSE_ELAS=OFF                                               \
-DUSE_VTK=OFF                                                \
-DBUILD_SHARED_LIBS=ON                                       \

# Build
cmake --build $(pwd) -- -j$(nproc)
