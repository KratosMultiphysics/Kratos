#ifndef VEXCL_BACKEND_JIT_SOURCE_HPP
#define VEXCL_BACKEND_JIT_SOURCE_HPP

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
 * \file   vexcl/backend/jit/source.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Source code generation for the JIT backend.
 */

#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
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
struct type_name_impl <shared_ptr<T> >
  : type_name_impl<global_ptr<T>>
{};

template <class T>
struct type_name_impl <regstr_ptr<T> >
  : type_name_impl<global_ptr<T>>
{};

template <class T>
struct type_name_impl <constant_ptr<T> >
  : type_name_impl<global_ptr<const typename std::decay<T>::type> >
{};

template<typename T>
struct type_name_impl<T*>
  : type_name_impl<global_ptr<T>>
{};

namespace backend {
namespace jit {

inline std::string standard_kernel_header(const command_queue &q) {
    return std::string(R"(
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <boost/config.hpp>

using std::max;
using std::min;

#if defined(__APPLE__) || defined(__MACOSX)

#include <cstdint>
typedef uint64_t ulong;
typedef uint32_t uint;
typedef uint16_t ushort;
typedef uint8_t  uchar;

#define sincosf __sincosf
#define sincos __sincos

#endif

struct ndrange {
    size_t x,y,z;
};

template <typename T>
struct vector_type2 {
    union {
        T s[2];
        struct { T s0, s1; };
        struct { T x, y; };
    };

    template <class O>
    const vector_type2& operator=(O o) {
        s0 = o;
        s1 = o;
        return *this;
    }

    const vector_type2& operator*=(const vector_type2 &o) {
        s0 *= o.s0;
        s1 *= o.s1;
        return *this;
    }

    template <typename U>
    const vector_type2& operator*=(U o) {
        s0 *= o;
        s1 *= o;
        return *this;
    }

    const vector_type2& operator/=(const vector_type2 &o) {
        s0 /= o.s0;
        s1 /= o.s1;
        return *this;
    }

    template <typename U>
    const vector_type2& operator/=(U o) {
        s0 /= o;
        s1 /= o;
        return *this;
    }

    const vector_type2& operator+=(const vector_type2 &o) {
        s0 += o.s0;
        s1 += o.s1;
        return *this;
    }

    template <typename U>
    const vector_type2& operator+=(U o) {
        s0 += o;
        s1 += o;
        return *this;
    }

    const vector_type2& operator-=(const vector_type2 &o) {
        s0 -= o.s0;
        s1 -= o.s1;
        return *this;
    }

    template <typename U>
    const vector_type2& operator-=(U o) {
        s0 -= o;
        s1 -= o;
        return *this;
    }
};

template <typename T>
vector_type2<T> operator*(vector_type2<T> a, const vector_type2<T> &b) {
    a *= b;
    return a;
}

template <typename T, typename U>
vector_type2<T> operator*(vector_type2<T> a, U b) {
    a *= b;
    return a;
}

template <typename T>
vector_type2<T> operator/(vector_type2<T> a, const vector_type2<T> &b) {
    a /= b;
    return a;
}

template <typename T, typename U>
vector_type2<T> operator/(vector_type2<T> a, U b) {
    a /= b;
    return a;
}

template <typename T>
vector_type2<T> operator+(vector_type2<T> a, const vector_type2<T> &b) {
    a += b;
    return a;
}

template <typename T, typename U>
vector_type2<T> operator+(vector_type2<T> a, U b) {
    a += b;
    return a;
}

template <typename T>
vector_type2<T> operator-(vector_type2<T> a, const vector_type2<T> &b) {
    a -= b;
    return a;
}

template <typename T, typename U>
vector_type2<T> operator-(vector_type2<T> a, U b) {
    a -= b;
    return a;
}

template <typename T>
struct vector_type4 {
    union {
        T s[4];
        struct { T s0, s1, s2, s3; };
        struct { T x, y, z, w; };
    };

    template <class O>
    const vector_type4& operator=(O o) {
        s0 = o;
        s1 = o;
        s2 = o;
        s3 = o;
        return *this;
    }

    const vector_type4& operator*=(const vector_type4 &o) {
        s0 *= o.s0;
        s1 *= o.s1;
        s2 *= o.s2;
        s3 *= o.s3;
        return *this;
    }

    template <typename U>
    const vector_type4& operator*=(U o) {
        s0 *= o;
        s1 *= o;
        s2 *= o;
        s3 *= o;
        return *this;
    }

    const vector_type4& operator/=(const vector_type4 &o) {
        s0 /= o.s0;
        s1 /= o.s1;
        s2 /= o.s2;
        s3 /= o.s3;
        return *this;
    }

    template <typename U>
    const vector_type4& operator/=(U o) {
        s0 /= o;
        s1 /= o;
        s2 /= o;
        s3 /= o;
        return *this;
    }

    const vector_type4& operator+=(const vector_type4 &o) {
        s0 += o.s0;
        s1 += o.s1;
        s2 += o.s2;
        s3 += o.s3;
        return *this;
    }

    template <typename U>
    const vector_type4& operator+=(U o) {
        s0 += o;
        s1 += o;
        s2 += o;
        s3 += o;
        return *this;
    }

    const vector_type4& operator-=(const vector_type4 &o) {
        s0 -= o.s0;
        s1 -= o.s1;
        s2 -= o.s2;
        s3 -= o.s3;
        return *this;
    }

    template <typename U>
    const vector_type4& operator-=(U o) {
        s0 -= o;
        s1 -= o;
        s2 -= o;
        s3 -= o;
        return *this;
    }
};

template <typename T>
vector_type4<T> operator*(vector_type4<T> a, const vector_type4<T> &b) {
    a *= b;
    return a;
}

template <typename T, typename U>
vector_type4<T> operator*(vector_type4<T> a, U b) {
    a *= b;
    return a;
}

template <typename T>
vector_type4<T> operator/(vector_type4<T> a, const vector_type4<T> &b) {
    a /= b;
    return a;
}

template <typename T, typename U>
vector_type4<T> operator/(vector_type4<T> a, U b) {
    a /= b;
    return a;
}

template <typename T>
vector_type4<T> operator+(vector_type4<T> a, const vector_type4<T> &b) {
    a += b;
    return a;
}

template <typename T, typename U>
vector_type4<T> operator+(vector_type4<T> a, U b) {
    a += b;
    return a;
}

template <typename T>
vector_type4<T> operator-(vector_type4<T> a, const vector_type4<T> &b) {
    a -= b;
    return a;
}

template <typename T, typename U>
vector_type4<T> operator-(vector_type4<T> a, U b) {
    a -= b;
    return a;
}

#define VECTOR_TYPES(T) \
    typedef vector_type2<T> T ## 2; \
    typedef vector_type4<T> T ## 3; \
    typedef vector_type4<T> T ## 4;

typedef unsigned char uchar;

VECTOR_TYPES(float)
VECTOR_TYPES(double)
VECTOR_TYPES(char)
VECTOR_TYPES(uchar)
VECTOR_TYPES(short)
VECTOR_TYPES(ushort)
VECTOR_TYPES(int)
VECTOR_TYPES(uint)
VECTOR_TYPES(long)
VECTOR_TYPES(ulong)

template <typename Uint>
inline Uint mulhi(Uint a, Uint b) {
    const unsigned WHALF = std::numeric_limits<Uint>::digits/2;
    const Uint LOMASK = ((Uint)(~(Uint)0)) >> WHALF;
    Uint lo = a*b;
    Uint ahi = a>>WHALF;
    Uint alo = a& LOMASK;
    Uint bhi = b>>WHALF;
    Uint blo = b& LOMASK;
    Uint ahbl = ahi*blo;
    Uint albh = alo*bhi;
    Uint ahbl_albh = ((ahbl&LOMASK) + (albh&LOMASK));
    Uint hi = (ahi*bhi) + (ahbl>>WHALF) +  (albh>>WHALF);
    hi += ahbl_albh >> WHALF;
    hi += ((lo >> WHALF) < (ahbl_albh&LOMASK));
    return hi;
}

template <class T>
T atomic_add(T *p, T val) {
    T old;
#pragma omp atomic capture
    {
        old = *p; *p += val;
    }
    return old;
}

template <class T>
T atomic_sub(T *p, T val) {
    T old;
#pragma omp atomic capture
    {
        old = *p; *p -= val;
    }
    return old;
}

struct kernel_api {
    virtual void execute(const ndrange*, size_t, char*) const = 0;
};

#define KERNEL_PARAMETER(type, name) \
    type name = *reinterpret_cast<type*>(_p); _p+= sizeof(type)

)") + get_program_header(q);
}

class source_generator {
    private:
        unsigned indent;
        bool first_prm;

        enum {
            undefined,
            inside_function,
            inside_kernel
        } prm_state;

        std::ostringstream src;

    public:
        source_generator() : indent(0), first_prm(true), prm_state(undefined)
        { }

        source_generator(const command_queue &q, bool include_standard_header = true)
            : indent(0), first_prm(true), prm_state(undefined)
        {
            if (include_standard_header) src << standard_kernel_header(q);
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
            first_prm = true;
            new_line() << return_type << " " << name;
            return *this;
        }

        template <class Return>
        source_generator& begin_function(const std::string &name) {
            return begin_function(type_name<Return>(), name);
        }

        source_generator& begin_function_parameters() {
            prm_state = inside_function;
            first_prm = true;
            return open("(");
        }

        source_generator& end_function_parameters() {
            prm_state = undefined;
            return close(")").open("{");
        }

        source_generator& end_function() {
            return close("}");
        }

        source_generator& begin_kernel(const std::string &name) {
            new_line() << "struct " << name << "_t : public kernel_api"; open("{");
            new_line() << "void work(const ndrange*, const ndrange*, char*, char*) const;";
            new_line() << "void execute(const ndrange *dim, size_t smem_size, char *prm) const"; open("{");
            new_line() << "std::vector<char> smem(smem_size);";
            new_line() << "#pragma omp parallel for collapse(3) firstprivate(smem)";
            new_line() << "for(size_t id_z = 0; id_z < dim->z; ++id_z)"; open("{");
            new_line() << "for(size_t id_y = 0; id_y < dim->y; ++id_y)"; open("{");
            new_line() << "for(size_t id_x = 0; id_x < dim->x; ++id_x)"; open("{");
            new_line() << "ndrange id = {id_x, id_y, id_z};";
            new_line() << "work(dim, &id, smem.data(), prm);";
            close("}").close("}").close("}").close("}").close("};");
            new_line() << "extern \"C\" BOOST_SYMBOL_EXPORT " << name << "_t " << name << ";";
            new_line() << name << "_t " << name << ";";
            new_line() << "void " << name << "_t::work(const ndrange *_dim, const ndrange *_id, char *_smem, char *_p) const";
            open("{");
            return *this;
        }

        source_generator& begin_kernel_parameters() {
            prm_state = inside_kernel;
            return *this;
        }

        source_generator& end_kernel_parameters() {
            prm_state = undefined;
            return *this;
        }

        source_generator& end_kernel() {
            return close("}");
        }

        source_generator& parameter(const std::string &prm_type, const std::string &name) {
            switch(prm_state) {
                case inside_kernel:
                    return kernel_parameter(prm_type, name);
                case inside_function:
                    return function_parameter(prm_type, name);
                default:
                    throw std::logic_error("parameter definition outside of parameter block");
            }
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
            new_line() << type_name<shared_ptr<Prm>>() << " " << name
                << " = reinterpret_cast<" << type_name< shared_ptr<Prm> >()
                << ">(_smem);";
            return *this;
        }

        source_generator& grid_stride_loop(
                const std::string &idx = "idx", const std::string &bnd = "n"
                )
        {
            new_line() << "size_t chunk_size = (" << bnd << " + " << global_size(0) << " - 1) / " << global_size(0) << ";";
            new_line() << "size_t chunk_start = chunk_size * " << global_id(0) << ";";
            new_line() << "size_t chunk_end = chunk_start + chunk_size;";
            new_line() << "if (" << bnd << " < chunk_end) chunk_end = " << bnd << ";";
            new_line() << "for(size_t " << idx << " = chunk_start; " << idx << " < chunk_end; ++" << idx << ")";

            return *this;
        }

        std::string global_id(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "_id->" << dim[d];
            return s.str();
        }

        std::string global_size(int d) const {
            const char dim[] = {'x', 'y', 'z'};
            std::ostringstream s;
            s << "_dim->" << dim[d];
            return s.str();
        }

        std::string local_id(int d) const {
            return "0";
        }

        std::string local_size(int d) const {
            return "1";
        }

        std::string group_id(int d) const {
            return global_id(d);
        }

        std::string num_groups(int d) const {
            return global_size(d);
        }

        source_generator& barrier(bool /*global*/ = false) {
            return *this;
        }

        source_generator& smem_static_var(const std::string &type, const std::string &name) {
            new_line() << type <<  " " << name << ";";
            return *this;
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

        source_generator& function_parameter(const std::string &prm_type, const std::string &name) {
            prm_separator().new_line() << prm_type << " " << name;
            return *this;
        }

        source_generator& kernel_parameter(const std::string &prm_type, const std::string &name) {
            new_line() << "KERNEL_PARAMETER(" << prm_type << ", " << name << ");";
            return *this;
        }
};

} // namespace jit
} // namespace backend
} // namespace vex

#endif
