#define BOOST_TEST_MODULE SparseMatrixVectorProduct
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/spmat.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"
#include "random_matrix.hpp"

BOOST_AUTO_TEST_CASE(vector_product)
{
    const size_t n = 1024;

    std::vector<size_t> row;
    std::vector<size_t> col;
    std::vector<double> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n);

    vex::SpMat <double> A(ctx, n, n, row.data(), col.data(), val.data());
    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });

    Y -= A * X;

    check_sample(Y, [&](size_t, double a) { BOOST_CHECK_SMALL(a, 1e-8); });

    Y += 42 * (A * X);

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, 42 * sum, 1e-8);
            });

    Y = X + A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, x[idx] + sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(non_square_matrix)
{
    const size_t n = 1024;
    const size_t m = 2 * n;

    std::vector<size_t> row;
    std::vector<size_t> col;
    std::vector<double> val;

    random_matrix(n, m, 16, row, col, val);

    std::vector<double> x = random_vector<double>(m);

    vex::SpMat <double> A(ctx, n, m, row.data(), col.data(), val.data());
    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(non_default_types)
{
    const size_t n = 1024;

    std::vector<unsigned> row;
    std::vector<int>      col;
    std::vector<double>   val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n);

    vex::SpMat <double, int, unsigned> A(ctx, n, n, row.data(), col.data(), val.data());
    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(empty_rows)
{
    const size_t n = 1024;
    const size_t non_empty_part = 256;

    std::vector<size_t> row;
    std::vector<size_t> col;
    std::vector<double> val;

    row.reserve(n + 1);

    random_matrix(non_empty_part, n, 16, row, col, val);

    while(row.size() < n + 1) row.push_back(col.size());

    std::vector<double> x = random_vector<double>(n);

    vex::SpMat <double> A(ctx, n, n, row.data(), col.data(), val.data());
    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(ccsr_vector_product)
{
    const size_t n = 32;
    const double h2i = (n - 1) * (n - 1);

    std::vector<size_t> idx;
    std::vector<size_t> row(3);
    std::vector<int>    col(8);
    std::vector<double> val(8);

    idx.reserve(n * n * n);

    row[0] = 0;
    row[1] = 1;
    row[2] = 8;

    col[0] = 0;
    val[0] = 1;

    col[1] = -static_cast<int>(n * n);
    col[2] = -static_cast<int>(n);
    col[3] =    -1;
    col[4] =     0;
    col[5] =     1;
    col[6] =     n;
    col[7] =  (n * n);

    val[1] = -h2i;
    val[2] = -h2i;
    val[3] = -h2i;
    val[4] =  h2i * 6;
    val[5] = -h2i;
    val[6] = -h2i;
    val[7] = -h2i;

    for(size_t k = 0; k < n; k++) {
        for(size_t j = 0; j < n; j++) {
            for(size_t i = 0; i < n; i++) {
                if (
                        i == 0 || i == (n - 1) ||
                        j == 0 || j == (n - 1) ||
                        k == 0 || k == (n - 1)
                   )
                {
                    idx.push_back(0);
                } else {
                    idx.push_back(1);
                }
            }
        }
    }

    std::vector<double> x = random_vector<double>(n * n * n);

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    vex::SpMatCCSR<double,int> A(queue[0], x.size(), row.size() - 1,
            idx.data(), row.data(), col.data(), val.data());

    vex::vector<double> X(queue, x);
    vex::vector<double> Y(queue, x.size());

    Y = A * X;

    check_sample(Y, [&](size_t ii, double a) {
            double sum = 0;
            size_t i = idx[ii];
            for(size_t j = row[i]; j < row[i + 1]; j++)
                sum += val[j] * x[ii + col[j]];

            BOOST_CHECK_CLOSE(a, sum, 1e-8);
            });

    Y = X + A * X;

    check_sample(Y, [&](size_t ii, double a) {
            double sum = 0;
            size_t i = idx[ii];
            for(size_t j = row[i]; j < row[i + 1]; j++)
                sum += val[j] * x[ii + col[j]];

            BOOST_CHECK_CLOSE(a, x[ii] + sum, 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(inline_spmv)
{
    const size_t n = 1024;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    std::vector<size_t> row;
    std::vector<size_t> col;
    std::vector<double> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n);

    vex::SpMat <double> A(queue, n, n, row.data(), col.data(), val.data());
    vex::vector<double> X(queue, x);
    vex::vector<double> Y(queue, n);

    Y = sin(vex::make_inline(A * X));

    check_sample(Y, [&](size_t idx, double a) {
            double sum = 0;
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a, sin(sum), 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(multivector_product)
{
    const size_t n = 1024;
    const size_t m = 2;

    typedef std::array<double, m> elem_t;

    std::vector<size_t> row;
    std::vector<size_t> col;
    std::vector<double> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n * m);

    vex::SpMat <double> A(ctx, n, n, row.data(), col.data(), val.data());

    vex::multivector<double,m> X(ctx, x);
    vex::multivector<double,m> Y(ctx, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, elem_t a) {
            double sum[] = {0, 0};
            for(size_t j = row[idx]; j < row[idx + 1]; j++) {
                sum[0] += val[j] * x[0 + col[j]];
                sum[1] += val[j] * x[n + col[j]];
            }

            BOOST_CHECK_CLOSE(a[0], sum[0], 1e-8);
            BOOST_CHECK_CLOSE(a[1], sum[1], 1e-8);
            });

    Y = X + A * X;

    check_sample(Y, [&](size_t idx, elem_t a) {
            double sum[] = {0, 0};
            for(size_t j = row[idx]; j < row[idx + 1]; j++) {
                sum[0] += val[j] * x[0 + col[j]];
                sum[1] += val[j] * x[n + col[j]];
            }

            BOOST_CHECK_CLOSE(a[0], x[0 + idx] + sum[0], 1e-8);
            BOOST_CHECK_CLOSE(a[1], x[n + idx] + sum[1], 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(inline_multivector_product)
{
    const size_t n = 1024;
    const size_t m = 2;

    typedef std::array<double, m> elem_t;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    std::vector<size_t> row;
    std::vector<size_t> col;
    std::vector<double> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<double> x = random_vector<double>(n * m);

    vex::SpMat <double> A(queue, n, n, row.data(), col.data(), val.data());

    vex::multivector<double,m> X(queue, x);
    vex::multivector<double,m> Y(queue, n);

    Y = cos(vex::make_inline(A * X));

    check_sample(Y, [&](size_t idx, elem_t a) {
            double sum[] = {0, 0};
            for(size_t j = row[idx]; j < row[idx + 1]; j++) {
                sum[0] += val[j] * x[0 + col[j]];
                sum[1] += val[j] * x[n + col[j]];
            }

            BOOST_CHECK_CLOSE(a[0], cos(sum[0]), 1e-8);
            BOOST_CHECK_CLOSE(a[1], cos(sum[1]), 1e-8);
            });
}

BOOST_AUTO_TEST_CASE(ccsr_multivector_product)
{
    const size_t n = 32;
    const size_t N = n * n * n;
    const double h2i = (n - 1) * (n - 1);

    typedef std::array<double, 2> elem_t;

    std::vector<size_t> idx;
    std::vector<size_t> row(3);
    std::vector<int>    col(8);
    std::vector<double> val(8);

    idx.reserve(N);

    row[0] = 0;
    row[1] = 1;
    row[2] = 8;

    col[0] = 0;
    val[0] = 1;

    col[1] = -static_cast<int>(n * n);
    col[2] = -static_cast<int>(n);
    col[3] =    -1;
    col[4] =     0;
    col[5] =     1;
    col[6] =     n;
    col[7] =  (n * n);

    val[1] = -h2i;
    val[2] = -h2i;
    val[3] = -h2i;
    val[4] =  h2i * 6;
    val[5] = -h2i;
    val[6] = -h2i;
    val[7] = -h2i;

    for(size_t k = 0; k < n; k++) {
        for(size_t j = 0; j < n; j++) {
            for(size_t i = 0; i < n; i++) {
                if (
                        i == 0 || i + 1 == n ||
                        j == 0 || j + 1 == n ||
                        k == 0 || k + 1 == n
                   )
                {
                    idx.push_back(0);
                } else {
                    idx.push_back(1);
                }
            }
        }
    }

    std::vector<double> x = random_vector<double>(N * 2);

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    vex::SpMatCCSR<double,int> A(queue[0], N, row.size() - 1,
            idx.data(), row.data(), col.data(), val.data());

    vex::multivector<double,2> X(queue, x);
    vex::multivector<double,2> Y(queue, N);

    Y = A * X;

    check_sample(Y, [&](size_t ii, elem_t a) {
            double sum[] = {0, 0};
            size_t i = idx[ii];
            for(size_t j = row[i]; j < row[i + 1]; j++) {
                sum[0] += val[j] * x[0 + ii + col[j]];
                sum[1] += val[j] * x[N + ii + col[j]];
            }

            BOOST_CHECK_CLOSE(a[0], sum[0], 1e-8);
            BOOST_CHECK_CLOSE(a[1], sum[1], 1e-8);
            });

    Y = X + A * X;

    check_sample(Y, [&](size_t ii, elem_t a) {
            double sum[] = {0, 0};
            size_t i = idx[ii];
            for(size_t j = row[i]; j < row[i + 1]; j++) {
                sum[0] += val[j] * x[0 + ii + col[j]];
                sum[1] += val[j] * x[N + ii + col[j]];
            }

            BOOST_CHECK_CLOSE(a[0], x[0 + ii] + sum[0], 1e-8);
            BOOST_CHECK_CLOSE(a[1], x[N + ii] + sum[1], 1e-8);
            });
}

#if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
BOOST_AUTO_TEST_CASE(vector_valued_matrix)
{
    const size_t n = 1024;

    std::vector<size_t>     row;
    std::vector<size_t>     col;
    std::vector<cl_double2> val;

    random_matrix(n, n, 16, row, col, val);

    std::vector<cl_double2> x = random_vector<cl_double2>(n);

    vex::SpMat <cl_double2> A(ctx, n, n, row.data(), col.data(), val.data());
    vex::vector<cl_double2> X(ctx, x);
    vex::vector<cl_double2> Y(ctx, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, cl_double2 a) {
            cl_double2 sum = {{0, 0}};
            for(size_t j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(a.s[0], sum.s[0], 1e-8);
            BOOST_CHECK_CLOSE(a.s[1], sum.s[1], 1e-8);
            });
}

cl_double2 make_d2(double v) {
    cl_double2 d2 = {{v, v}};
    return d2;
}

BOOST_AUTO_TEST_CASE(vector_valued_ccsr_matrix)
{
    const size_t n = 32;
    const double h2i = (n - 1) * (n - 1);

    std::vector<size_t> idx;
    std::vector<size_t> row(3);
    std::vector<int>    col(8);
    std::vector<cl_double2> val(8);

    idx.reserve(n * n * n);

    row[0] = 0;
    row[1] = 1;
    row[2] = 8;

    col[0] = 0;
    val[0] = make_d2(1);

    col[1] = -static_cast<int>(n * n);
    col[2] = -static_cast<int>(n);
    col[3] =    -1;
    col[4] =     0;
    col[5] =     1;
    col[6] =     n;
    col[7] =  (n * n);

    val[1] = make_d2(-h2i);
    val[2] = make_d2(-h2i);
    val[3] = make_d2(-h2i);
    val[4] = make_d2( h2i * 6);
    val[5] = make_d2(-h2i);
    val[6] = make_d2(-h2i);
    val[7] = make_d2(-h2i);

    for(size_t k = 0; k < n; k++) {
        for(size_t j = 0; j < n; j++) {
            for(size_t i = 0; i < n; i++) {
                if (
                        i == 0 || i == (n - 1) ||
                        j == 0 || j == (n - 1) ||
                        k == 0 || k == (n - 1)
                   )
                {
                    idx.push_back(0);
                } else {
                    idx.push_back(1);
                }
            }
        }
    }

    std::vector<cl_double2> x = random_vector<cl_double2>(n * n * n);

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    vex::SpMatCCSR<cl_double2,int> A(queue[0], x.size(), row.size() - 1,
            idx.data(), row.data(), col.data(), val.data());

    vex::vector<cl_double2> X(queue, x);
    vex::vector<cl_double2> Y(queue, x.size());

    Y = A * X;

    check_sample(Y, [&](size_t ii, cl_double2 a) {
            cl_double2 sum = {{0, 0}};
            size_t i = idx[ii];
            for(size_t j = row[i]; j < row[i + 1]; j++)
                sum += val[j] * x[ii + col[j]];

            BOOST_CHECK_CLOSE(a.s[0], sum.s[0], 1e-8);
            BOOST_CHECK_CLOSE(a.s[1], sum.s[1], 1e-8);
            });
}
#endif

BOOST_AUTO_TEST_SUITE_END()
