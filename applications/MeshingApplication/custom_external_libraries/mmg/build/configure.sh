cmake .. \
-DCMAKE_BUILD_TYPE=Release                                  \
-DUSE_SCOTCH=OFF                                            \
-DSCOTCH_INCLUDE_DIR="/usr/include/scotch/"	            \
-DSCOTCH_LIBRARY="libscotchmetis-5.1.so"	            \
-DSCOTCHERR_LIBRARY="libscotcherr-5.1.so"	            \
-DCMAKE_CXX_FLAGS="-O3 -msse3 -fPIC -fopenmp"               \
-DCMAKE_C_FLAGS="-O3 -msse3 -fPIC -fopenmp"                 \
-DLIBMMG_SHARED=ON                                       \
-DLIBMMG_STATIC=ON                                       \
-DLIBMMGS_SHARED=ON                                      \
-DLIBMMGS_STATIC=ON                                      \
-DLIBMMG2D_SHARED=ON                                     \
-DLIBMMG2D_STATIC=ON                                     \
-DLIBMMG3D_SHARED=ON                                     \
-DLIBMMG3D_STATIC=ON                                     \
# If you ha modern processor use this instructions instead (look the whole list here https://software.intel.com/sites/landingpage/IntrinsicsGuide/)
# With relatively new (Intel i7 2013 works with -maxvx)
# -DCMAKE_CXX_FLAGS="-O3 -mavx2 -fPIC -fopenmp"               \
# -DCMAKE_C_FLAGS="-O3 -mavx2 -fPIC -fopenmp"                 \

#decomment this to have it verbose
# make VERBOSE=1 -j4
make -j2
# sudo make install
