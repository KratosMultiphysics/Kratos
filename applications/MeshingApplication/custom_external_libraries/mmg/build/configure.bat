del CMakeCache.txt

cls

cmake .. -G"Visual Studio 16 2019" -A x64   ^
-DCMAKE_BUILD_TYPE=Release                  ^
-DCMAKE_CXX_FLAGS="-O3"                     ^
-DCMAKE_C_FLAGS="-O3"                       ^
-DUSE_SCOTCH=OFF                            ^
-DLIBMMG_SHARED=OFF                         ^
-DLIBMMG_STATIC=ON                          ^
-DLIBMMGS_SHARED=OFF                        ^
-DLIBMMGS_STATIC=ON                         ^
-DLIBMMG2D_SHARED=OFF                       ^
-DLIBMMG2D_STATIC=ON                        ^
-DLIBMMG3D_SHARED=OFF                       ^
-DLIBMMG3D_STATIC=ON                                   

cmake --build . --target install -- /property:configuration=Release /p:Platform=x64
