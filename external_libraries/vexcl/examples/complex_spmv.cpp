#include <iostream>
#include <vector>
#include <complex>

#include <vexcl/vexcl.hpp>

// The following enables use of std::complex<T> in vexcl expressions.
// It simply translates std::complex<T> on the host side into opencl vector
// types (float2/double2) on the device side.
//
// However, semantics of operations like multiplication, division, or
// math functions will be wrong. This is the main reason I hesitate doing
// this in the core library.
namespace vex {

template <typename T>
struct is_cl_native< std::complex<T> > : std::true_type {};

template <typename T>
struct type_name_impl< std::complex<T> >
{
    static std::string get() {
        std::ostringstream s;
        s << type_name<T>() << "2";
        return s.str();
    }
};

template <typename T>
struct cl_scalar_of< std::complex<T> > {
    typedef T type;
};

} // namespace vex

// Now we specialize a template from <vexcl/sparse/spmv_ops.hpp> that allows
// vexcl to generate semantically correct smpv kernels for complex values:
namespace vex {
namespace sparse {

template <typename T>
struct spmv_ops_impl<std::complex<T>, std::complex<T>>
{
    static void decl_accum_var(backend::source_generator &src, const std::string &name)
    {
        src.new_line() << type_name<T>() << "2 " << name << " = {0,0};";
    }

    static void append(backend::source_generator &src,
            const std::string &sum, const std::string &val)
    {
        src.new_line() << sum << " += " << val << ";";
    }

    static void append_product(backend::source_generator &src,
            const std::string &sum, const std::string &mat_val, const std::string &vec_val)
    {
        src.new_line() << sum << ".x += "
            << mat_val << ".x * " << vec_val << ".x - "
            << mat_val << ".y * " << vec_val << ".y;";
        src.new_line() << sum << ".y += "
            << mat_val << ".x * " << vec_val << ".y + "
            << mat_val << ".y * " << vec_val << ".x;";
    }
};

} // namespace sparse
} // namespace vex

int main() {
    vex::Context ctx(vex::Filter::Env && vex::Filter::Count(1));
    std::cout << ctx << std::endl;

    // 4x4 diagonal matrix in CSR format:
    std::vector<int> ptr = {0,1,2,3,4};
    std::vector<int> col = {0,1,2,3};
    std::vector<std::complex<double>> val = {
        {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0}};

    // complex vector:
    std::vector<std::complex<double>> x = {
        {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};

    // Device-side matrix and vectors:
    vex::sparse::matrix<std::complex<double>> A(ctx, 4, 4, ptr, col, val);
    vex::vector<std::complex<double>> X(ctx, x);
    vex::vector<std::complex<double>> Y(ctx, 4);


    Y = A * X;

    std::cout << Y << std::endl;
}
