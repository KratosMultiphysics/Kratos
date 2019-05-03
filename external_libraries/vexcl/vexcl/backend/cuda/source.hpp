#ifndef VEXCL_BACKEND_CUDA_SOURCE_HPP
#define VEXCL_BACKEND_CUDA_SOURCE_HPP

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
 * \file   vexcl/backend/cuda/source.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Helper class for CUDA source code generation.
 */

#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <cassert>

#include <vexcl/backend/common.hpp>
#include <vexcl/types.hpp>

namespace vex {

template <class T> struct global_ptr {};
template <class T> struct shared_ptr {};
template <class T> struct regstr_ptr {};

template <class T>
struct type_name_impl <global_ptr<T> > {
    static std::string get() {
        std::ostringstream s;
        s << type_name<T>() << " *";
        return s.str();
    }
};

template <class T>
struct type_name_impl < global_ptr<const T> > {
    static std::string get() {
        std::ostringstream s;
        s << "const " << type_name<T>() << " *";
        return s.str();
    }
};

template <class T>
struct type_name_impl <shared_ptr<T> > : type_name_impl< global_ptr<T> > { };

template <class T>
struct type_name_impl <shared_ptr<const T> > : type_name_impl< global_ptr<const T> > { };

template <class T>
struct type_name_impl <regstr_ptr<T> > : type_name_impl< global_ptr<T> > { };

template <class T>
struct type_name_impl <regstr_ptr<const T> > : type_name_impl< global_ptr<const T> > { };

template<typename T>
struct type_name_impl<T*>
{
    static std::string get() {
        return type_name_impl< global_ptr<T> >::get();
    }
};

namespace backend {
namespace cuda {

/// Returns standard CUDA program header.
/**
 * Defines pragmas necessary to work with double precision and anything
 * provided by the user with help of push_program_header().
 */
inline std::string standard_kernel_header(const command_queue &q) {
    std::ostringstream s;
    s
#if defined(_MSC_VER) || defined(__APPLE__)
      << "typedef unsigned char       uchar;\n"
         "typedef unsigned int        uint;\n"
         "typedef unsigned short      ushort;\n"
         "typedef unsigned long long  ulong;\n"
#endif
      << get_program_header(q);
    return s.str();
}

class source_generator {
    private:
        unsigned           indent;
        bool               first_prm;
        std::ostringstream src;

    public:
        source_generator() : indent(0), first_prm(true) { }

        source_generator(const command_queue &queue, bool include_standard_header = true)
            : indent(0), first_prm(true)
        {
            if (include_standard_header)
                src << standard_kernel_header(queue);
        }

        source_generator& new_line() {
            src << "\n" << std::string(2 * indent, ' ');
            return *this;
        }

        source_generator& open(const char *bracket) {
            new_line() << bracket;
            ++indent;
            return *this;
        }

        source_generator& close(const char *bracket) {
            assert(indent > 0);
            --indent;
            new_line() << bracket;
            return *this;
        }

        source_generator& begin_function(const std::string &return_type, const std::string &name) {
            new_line() << "__device__ " << return_type << " " << name;
            return *this;
        }

        template <class Return>
        source_generator& begin_function(const std::string &name) {
            return begin_function(type_name<Return>(), name);
        }

        source_generator& begin_function_parameters() {
            first_prm = true;
            return open("(");
        }

        source_generator& end_function_parameters() {
            return close(")").open("{");
        }

        source_generator& end_function() {
            return close("}");
        }

        source_generator& begin_kernel(const std::string &name) {
            new_line() << "extern \"C\" __global__ void " << name;
            return *this;
        }

        source_generator& begin_kernel_parameters() {
            first_prm = true;
            return open("(");
        }

        source_generator& end_kernel_parameters() {
            return close(")").open("{");
        }

        source_generator& end_kernel() {
            return close("}");
        }

        source_generator& parameter(const std::string &prm_name, const std::string &name) {
            prm_separator().new_line() << prm_name << " " << name;
            return *this;
        }

        template <class Prm>
        source_generator& parameter(const std::string &name) {
            return parameter(type_name<typename std::decay<Prm>::type>(), name);
        }

        template <class Prm>
        source_generator& smem_parameter(const std::string& = "smem") {
            return *this;
        }

        template <class Prm>
        source_generator& smem_declaration(const std::string &name = "smem") {
            new_line() << "extern __shared__ " << type_name<Prm>() << " " << name << "[];";
            return *this;
        }

        source_generator& smem_static_var(const std::string &type, const std::string &name) {
            new_line() << "__shared__ " << type <<  " " << name << ";";
            return *this;
        }

        source_generator& grid_stride_loop(
                const std::string &idx = "idx", const std::string &bnd = "n"
                )
        {
            new_line() << "for";
            open("(");
            new_line() << type_name<size_t>() << " " << idx << " = blockDim.x * blockIdx.x + threadIdx.x, "
                "grid_size = blockDim.x * gridDim.x;";
            new_line() << idx << " < " << bnd << ";";
            new_line() << idx << " += grid_size";
            close(")");

            return *this;
        }

        source_generator& barrier(bool /*global*/ = false) {
            src << "__syncthreads();";
            return *this;
        }

        std::string global_id(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "(threadIdx." << dim[d] << " + blockIdx." << dim[d]
              << " * blockDim." << dim[d] << ")";
            return s.str();
        }

        std::string global_size(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "(blockDim." << dim[d] << " * gridDim." << dim[d] << ")";
            return s.str();
        }

        std::string local_id(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "threadIdx." << dim[d];
            return s.str();
        }

        std::string local_size(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "blockDim." << dim[d];
            return s.str();
        }

        std::string group_id(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "blockIdx." << dim[d];
            return s.str();
        }

        std::string num_groups(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "gridDim." << dim[d];
            return s.str();
        }

        std::string str() const {
            return src.str();
        }
    private:
        template <class T>
        friend inline
        source_generator& operator<<(source_generator &src, const T &t) {
            src.src << t;
            return src;
        }

        source_generator& prm_separator() {
            if (first_prm)
                first_prm = false;
            else
                src << ",";

            return *this;
        }
};

} // namespace cuda
} // namespace backend

} // namespace vex

#endif
