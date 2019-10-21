#include <thrust/sort.h>
#include <thrust/device_ptr.h>

//---------------------------------------------------------------------------
// NVCC is not yet able to compile C++11 code.
// Hence the need to keep Thrust and VexCL code in separate files.
//---------------------------------------------------------------------------
template <typename T>
void thrust_sort(T *begin, T *end) {
    thrust::sort(
            thrust::device_pointer_cast(begin),
            thrust::device_pointer_cast(end)
            );
}

//---------------------------------------------------------------------------
// Due to the code separation we also need to explicitly instantiate the
// necessary templates.
//---------------------------------------------------------------------------
#define VEXCL_INSTANTIATE_THRUST_SORT(T)                                       \
  template void thrust_sort<T>(T * begin, T * end)

VEXCL_INSTANTIATE_THRUST_SORT(int);

#undef VEXCL_INSTANTIATE_THRUST_SORT
