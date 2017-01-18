cmake .. \
-DCMAKE_BUILD_TYPE=Release      \
-DUSE_SCOTCH=ON      \
-DSCOTCH_INCLUDE_DIR="/usr/include/scotch/"	            \
-DSCOTCH_LIBRARY="libscotchmetis-5.1.so"	            \
-DSCOTCHERR_LIBRARY="libscotcherr-5.1.so"	            \
-DCMAKE_CXX_FLAGS="-O3 -mavx2 -fPIC -fopenmp" \
-DCMAKE_C_FLAGS="-O3 -mavx2 -fPIC -fopenmp"
