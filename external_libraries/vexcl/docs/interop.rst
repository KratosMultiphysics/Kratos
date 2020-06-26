Interoperability with other libraries
=====================================

VexCL does not try (too hard) to hide the implementation details from the user.
For example, in case of the OpenCL backend VexCL is based on the `Khronos C++
API`_, and the underlying OpenCL types are easily accessible. Hence, it should
be easy to interoperate with other OpenCL libraries. Similarly, in case of the
CUDA backend, VexCL backend types are thin wrappers around `CUDA Driver API`_.

When Boost.Compute_ backend is used, VexCL is based on the core classes of the
Boost.Compute_ library. It is very easy to apply Boost.Compute_ algorithms to
VexCL vectors and to use Boost.Compute_ buffers within VexCL expressions.

Here is an example:

.. code-block:: cpp

    #include <iostream>
    #include <boost/compute.hpp>

    #define VEXCL_BACKEND_COMPUTE
    #include <vexcl/vexcl.hpp>

    namespace compute = boost::compute;

    int main() {
	compute::command_queue bcq = compute::system::default_queue();

	const int n = 16;

	// Use boost.compute queue to allocate VexCL vectors:
	vex::vector<int> x({bcq}, n);
	x = 2 * vex::element_index();

	// Wrap boost.compute vectors into vexcl vectors (no data is copied):
	compute::vector<int> bcv(n, bcq.get_context());
	vex::vector<int> y({bcq}, bcv.get_buffer());
	y = x * 2;

	// Apply Boost.Compute algorithm to a vexcl vector:
	compute::sort(
	    compute::make_buffer_iterator<int>(x(0).raw_buffer(), 0),
	    compute::make_buffer_iterator<int>(x(0).raw_buffer(), n)
	    );
    }

.. _`Khronos C++ API`: https://www.khronos.org/registry/cl
.. _`CUDA Driver API`: http://docs.nvidia.com/cuda/cuda-driver-api/
.. _Boost.Compute: https://github.com/boostorg/compute

