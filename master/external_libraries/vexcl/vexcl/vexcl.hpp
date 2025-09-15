#ifndef VEXCL_VEXCL_HPP
#define VEXCL_VEXCL_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   vexcl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Vector expression template library for OpenCL.
 */

#include <vexcl/backend.hpp>

#include <vexcl/devlist.hpp>
#include <vexcl/constants.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/vector_view.hpp>
#include <vexcl/tensordot.hpp>
#include <vexcl/vector_pointer.hpp>
#include <vexcl/tagged_terminal.hpp>
#include <vexcl/temporary.hpp>
#include <vexcl/cast.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/spmat.hpp>
#include <vexcl/sparse/distributed.hpp>
#include <vexcl/sparse/matrix.hpp>
#include <vexcl/stencil.hpp>
#include <vexcl/gather.hpp>
#include <vexcl/random.hpp>
#include <vexcl/fft.hpp>
#include <vexcl/mba.hpp>
#include <vexcl/generator.hpp>
#include <vexcl/mba.hpp>
#include <vexcl/sort.hpp>
#include <vexcl/scan.hpp>
#include <vexcl/scan_by_key.hpp>
#include <vexcl/reduce_by_key.hpp>
#include <vexcl/profiler.hpp>
#include <vexcl/function.hpp>
#include <vexcl/logical.hpp>
#include <vexcl/enqueue.hpp>
#include <vexcl/image.hpp>
#include <vexcl/eval.hpp>

#ifndef VEXCL_BACKEND_CUDA
#include <vexcl/constant_address_space.hpp>
#endif

#endif
