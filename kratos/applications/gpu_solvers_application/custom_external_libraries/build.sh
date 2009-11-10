# Compilation options:
#
#   -DUSE_UMUL24: Use __umul24 in GlobalIdx() instead of simple multiplication
#   -DUSE_TEXTURE_CACHING: Use texture caching in matrix vector multiplication kernels
#
# To run, simply issue: ./build.sh
#
nvcc gpu_sparse.cu -o gpu_sparse.o -c -g -Xcompiler -fPIC -DNDEBUG -DUSE_UMUL24 -DUSE_TEXTURE_CACHING -arch compute_13 -code compute_13 -lcublas -I/usr/local/lapackpp/include/lapackpp
ar cr libgpu_sparse.a gpu_sparse.o
