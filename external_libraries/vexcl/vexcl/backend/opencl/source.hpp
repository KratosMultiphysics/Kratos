#ifndef VEXCL_BACKEND_OPENCL_SOURCE_HPP
#define VEXCL_BACKEND_OPENCL_SOURCE_HPP

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
 * \file   vexcl/backend/opencl/source.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Helper class for OpenCL source code generation.
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
template <class T> struct constant_ptr {};

template <class T>
struct type_name_impl <global_ptr<T> > {
    static std::string get() {
        std::ostringstream s;
        s << "global " << type_name<T>() << " *";
        return s.str();
    }
};

template <class T>
struct type_name_impl < global_ptr<const T> > {
    static std::string get() {
        std::ostringstream s;
        s << "global const " << type_name<T>() << " *";
        return s.str();
    }
};

template <class T>
struct type_name_impl <shared_ptr<T> > {
    static std::string get() {
        std::ostringstream s;
        s << "local " << type_name<T>() << " *";
        return s.str();
    }
};

template <class T>
struct type_name_impl <shared_ptr<const T> > {
    static std::string get() {
        std::ostringstream s;
        s << "local const " << type_name<T>() << " *";
        return s.str();
    }
};

template <class T>
struct type_name_impl <regstr_ptr<T> > {
    static std::string get() {
        std::ostringstream s;
        s << type_name<T>() << " *";
        return s.str();
    }
};

template <class T>
struct type_name_impl <regstr_ptr<const T> > {
    static std::string get() {
        std::ostringstream s;
        s << "const " << type_name<T>() << " *";
        return s.str();
    }
};

// Note that constant_ptr<T> and global_ptr<const T> are not the same.
// The former uses constant cache for read-only access to vector data
// while the latter simply informs the compiler that it is illegal to modify
// the vector data.
template <class T>
struct type_name_impl <constant_ptr<T> > {
    static std::string get() {
        std::ostringstream s;
        s << "constant "
          << type_name<typename std::decay<T>::type>()
          << " *";
        return s.str();
    }
};

template<typename T>
struct type_name_impl<T*>
{
    static std::string get() {
        return type_name_impl< global_ptr<T> >::get();
    }
};

namespace backend {
namespace opencl {

/// Returns standard OpenCL program header.
/**
 * Defines pragmas necessary to work with double precision and anything
 * provided by the user with help of push_program_header().
 */
inline std::string standard_kernel_header(const command_queue &q) {
    return std::string(
        "#if defined(cl_khr_fp64)\n"
        "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n"
        "#elif defined(cl_amd_fp64)\n"
        "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n"
        "#endif\n"
        )
        + get_program_header(q)
#ifdef VEXCL_AMD_SI_WORKAROUND
        + "kernel void __null_kernel() {}\n"
#endif
        ;
}

/// Helper class for OpenCL source code generation.
class source_generator {
    private:
        unsigned           indent;
        bool               first_prm, cpu;
        std::ostringstream src;

    public:
        source_generator() : indent(0), first_prm(true), cpu(false) { }

        source_generator(const command_queue &queue, bool include_standard_header = true)
            : indent(0), first_prm(true), cpu( is_cpu(queue) )
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
            new_line() << return_type << " " << name;
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
            new_line() << "kernel void " << name;
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

        source_generator& parameter(const std::string &prm_type, const std::string &name) {
            prm_separator().new_line() << prm_type << " " << name;
            return *this;
        }

        template <class Prm>
        source_generator& parameter(const std::string &name) {
            return parameter(type_name<typename std::decay<Prm>::type>(), name);
        }

        template <class Prm>
        source_generator& smem_parameter(const std::string &name = "smem") {
            return parameter< shared_ptr<Prm> >(name);
        }

        template <class Prm>
        source_generator& smem_declaration(const std::string& = "smem") {
            return *this;
        }

        source_generator& smem_static_var(const std::string &type, const std::string &name) {
            new_line() << "local " << type <<  " " << name << ";";
            return *this;
        }

        source_generator& grid_stride_loop(
                const std::string &idx = "idx", const std::string &bnd = "n"
                )
        {
            if (
                    cpu
#ifdef __APPLE__
                    && 0
#endif
               )
            {
                new_line() << type_name<size_t>() << " chunk_size  = (" << bnd
                           << " + get_global_size(0) - 1) / get_global_size(0);";
                new_line() << type_name<size_t>() << " chunk_start = get_global_id(0) * chunk_size;";
                new_line() << type_name<size_t>() << " chunk_end   = chunk_start + chunk_size;";
                new_line() << "if (" << bnd << " < chunk_end) chunk_end = " << bnd << ";";
                new_line() << "for(" << type_name<size_t>() << " "<< idx << " = chunk_start; "
                           << idx << " < chunk_end; ++" << idx << ")";
            } else {
                new_line() <<
                    "for(" << type_name<size_t>() << " " << idx << " = get_global_id(0);"
                    " " << idx << " < " << bnd << ";"
                    " " << idx << " += get_global_size(0))";
            }
            return *this;
        }

        source_generator& barrier(bool global = false) {
            if (global)
                src << "barrier(CLK_GLOBAL_MEM_FENCE);";
            else
                src << "barrier(CLK_LOCAL_MEM_FENCE);";
            return *this;
        }

        std::string global_id(int d) const {
            std::ostringstream s;
            s << "get_global_id(" << d << ")";
            return s.str();
        }

        std::string global_size(int d) const {
            std::ostringstream s;
            s << "get_global_size(" << d << ")";
            return s.str();
        }

        std::string local_id(int d) const {
            std::ostringstream s;
            s << "get_local_id(" << d << ")";
            return s.str();
        }

        std::string local_size(int d) const {
            std::ostringstream s;
            s << "get_local_size(" << d << ")";
            return s.str();
        }

        std::string group_id(int d) const {
            std::ostringstream s;
            s << "get_group_id(" << d << ")";
            return s.str();
        }

        std::string num_groups(int d) const {
            std::ostringstream s;
            s << "get_num_groups(" << d << ")";
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

} // namespace opencl
} // namespace backend

} // namespace vex

#endif
