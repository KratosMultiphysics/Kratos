#define BOOST_TEST_MODULE SparseMatrices
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/sparse/csr.hpp>
#include <vexcl/sparse/ell.hpp>
#include <vexcl/sparse/matrix.hpp>
#include <vexcl/sparse/distributed.hpp>

typedef std::array<std::array<double, 2>, 2> matrix_value;
typedef std::array<double, 2> vector_value;

matrix_value mconst(double c) {
    matrix_value a = {c, c, c, c};
    return a;
}

vector_value vconst(double c) {
    vector_value v = {c, c};
    return v;
}

namespace vex {

template <> struct is_cl_native< matrix_value > : std::true_type {};
template <> struct is_cl_native< vector_value > : std::true_type {};

template <> struct type_name_impl<matrix_value> {
    static std::string get() { return "double4"; }
};

template <> struct type_name_impl<vector_value> {
    static std::string get() { return "double2"; }
};

namespace sparse {

template <>
struct rhs_of< matrix_value > {
    typedef vector_value type;
};

template <>
struct spmv_ops_impl<matrix_value, vector_value> {
    static void decl_accum_var(backend::source_generator &src, const std::string &name)
    {
        src.new_line() << type_name<vector_value>() << " " << name << " = {0,0};";
    }

    static void append_product(backend::source_generator &src,
            const std::string &sum, const std::string &mat_val, const std::string &vec_val)
    {
        src.open("{");
        src.new_line() << type_name<vector_value>() << " v = " << vec_val << ";";
        src.new_line() << sum << ".x += " << mat_val << ".x * v.x + " << mat_val << ".y * v.y;";
        src.new_line() << sum << ".y += " << mat_val << ".z * v.x + " << mat_val << ".w * v.y;";
        src.close("}");
    }
};

} // namespace sparse
} // namespace vex

#include "context_setup.hpp"
#include "random_matrix.hpp"

BOOST_AUTO_TEST_CASE(csr)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q(1, ctx.queue(0));

    std::vector<int>    row;
    std::vector<int>    col;
    std::vector<double> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n);

    vex::sparse::csr<double> A(q, n, n, row, col, val);
    vex::vector<double> X(q, x);
    vex::vector<double> Y(q, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(int j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(ell)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q(1, ctx.queue(0));

    std::vector<int>    row;
    std::vector<int>    col;
    std::vector<double> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n);

    vex::sparse::ell<double> A(q, n, n, row, col, val);
    vex::vector<double> X(q, x);
    vex::vector<double> Y(q, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(int j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(matrix)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> q(1, ctx.queue(0));

    std::vector<int>    row;
    std::vector<int>    col;
    std::vector<double> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n);

    vex::sparse::matrix<double> A(q, n, n, row, col, val);
    vex::vector<double> X(q, x);
    vex::vector<double> Y(q, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(int j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(distributed)
{
    const int n = 1024;

    std::vector<int>    ptr;
    std::vector<int>    col;
    std::vector<double> val;

    ptr.push_back(0);
    for(int i = 0; i < n; ++i) {
        if (i > 0) {
            col.push_back(i-1);
            val.push_back(-1);
        }
        col.push_back(i);
        val.push_back(2);
        if (i + 1 < n) {
            col.push_back(i+1);
            val.push_back(-1);
        }

        ptr.push_back(static_cast<int>(col.size()));
    }

    vex::sparse::distributed<vex::sparse::ell<double>> A(ctx, n, n, ptr, col, val);

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = A * X;

    for(int i = 0; i < n; ++i) {
        double y = Y[i];
        double sum = 0;
        for(int j = ptr[i]; j < ptr[i + 1]; j++)
            sum += val[j] * x[col[j]];

        BOOST_CHECK_CLOSE(y, sum, 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(distributed_single)
{
    std::vector<vex::command_queue> q(1, ctx.queue(0));

    const int n = 1024;

    std::vector<int>    ptr;
    std::vector<int>    col;
    std::vector<double> val;

    ptr.push_back(0);
    for(int i = 0; i < n; ++i) {
        if (i > 0) {
            col.push_back(i-1);
            val.push_back(-1);
        }
        col.push_back(i);
        val.push_back(2);
        if (i + 1 < n) {
            col.push_back(i+1);
            val.push_back(-1);
        }

        ptr.push_back(static_cast<int>(col.size()));
    }

    vex::sparse::distributed<vex::sparse::ell<double>> A(q, n, n, ptr, col, val);

    std::vector<double> x = random_vector<double>(n);
    vex::vector<double> X(q, x);
    vex::vector<double> Y(q, n);

    Y = A * X;

    for(int i = 0; i < n; ++i) {
        double y = Y[i];
        double sum = 0;
        for(int j = ptr[i]; j < ptr[i + 1]; j++)
            sum += val[j] * x[col[j]];

        BOOST_CHECK_CLOSE(y, sum, 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(custom_values)
{
    const int n = 1024;
    std::vector<vex::command_queue> q(1, ctx.queue(0));

    std::vector<int> ptr, col;
    std::vector<matrix_value> val;

    ptr.push_back(0);
    for(int i = 0; i < n; ++i) {
        if (i > 0) {
            col.push_back(i-1);
            val.push_back(mconst(-1));
        }
        col.push_back(i);
        val.push_back(mconst(2));
        if (i + 1 < n) {
            col.push_back(i+1);
            val.push_back(mconst(-1));
        }

        ptr.push_back(static_cast<int>(col.size()));
    }

    vex::sparse::matrix<matrix_value> A(q, n, n, ptr, col, val);

    std::vector<vector_value> x(n, vconst(1));
    vex::vector<vector_value> X(q, x);
    vex::vector<vector_value> Y(q, n);

    Y = A * X;

    for(int i = 0; i < n; ++i) {
        vector_value y = Y[i];
        vector_value sum = {0,0};
        for(int j = ptr[i]; j < ptr[i + 1]; j++) {
            sum[0] += val[j][0][0] * x[col[j]][0] + val[j][0][1] * x[col[j]][1];
            sum[1] += val[j][1][0] * x[col[j]][0] + val[j][1][1] * x[col[j]][1];
        }

        BOOST_CHECK_CLOSE(y[0], sum[0], 1e-8);
        BOOST_CHECK_CLOSE(y[1], sum[1], 1e-8);
    }
}


BOOST_AUTO_TEST_SUITE_END()
