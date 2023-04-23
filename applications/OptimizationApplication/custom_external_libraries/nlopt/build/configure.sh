# Set compiler
	export CC=${CC:-gcc}
	export CXX=${CXX:-g++}
	

	# Configure
	cmake ..                                                     \
	-DCMAKE_BUILD_TYPE=Release                                   \
	-DCMAKE_C_COMPILER=${CC}                                     \
	-DCMAKE_CXX_COMPILER=${CXX} 
-DCMAKE_C_FLAGS="-O3 -march=native -mtune=native -fopenmp -I~/nlopt-2.7.1/nlopt-2.7.1/build" \
-DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native -fopenmp -I~/nlopt-2.7.1/nlopt-2.7.1/build" \
-DCMAKE_EXE_LINKER_FLAGS="-L~/nlopt-2.7.1/nlopt-2.7.1/ -lnlopt"
# Build	
		cmake --build $(pwd) -- -j$(nproc)
