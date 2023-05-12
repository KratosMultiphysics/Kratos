del CMakeCache.txt

cls

cmake .. -G"Visual Studio 16 2019" -A x64   ^
-DCMAKE_BUILD_TYPE=Release                  ^
-DCMAKE_CXX_FLAGS="-O3"                     ^
-DCMAKE_C_FLAGS="-O3"                       ^
-DUSE_SCOTCH=OFF                            ^
-DUSE_ELAS=OFF                              ^
-DUSE_VTK=OFF                               ^
-DBUILD_SHARED_LIBS=ON                                                          

cmake --build . --target install -- /property:configuration=Release /p:Platform=x64