#define BOOST_TEST_MODULE CUSPARSE
#include <cuda.h>
#if (CUDA_VERSION <= 10010)
#  define VEXCL_USE_CUSPARSE
#endif
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/spmat.hpp>
#include "random_matrix.hpp"
#include "context_setup.hpp"

#if (CUDA_VERSION <= 10010)
BOOST_AUTO_TEST_CASE(hyb_matrix)
{
    const size_t n = 1024;
    const size_t m = 2048;

    std::vector<int>    row;
    std::vector<int>    col;
    std::vector<double> val;
    std::vector<double> x = random_vector<double>(m);

    random_matrix(n, m, 16, row, col, val);

    vex::SpMat<double, int, int> A(ctx, n, m, row.data(), col.data(), val.data());

    vex::vector<double> X(ctx, x);
    vex::vector<double> Y(ctx, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double y) {
            double sum = 0;
            for(int j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(y, sum, 1e-8);
            });

    Y += 0.5 * (A * X);

    check_sample(Y, [&](size_t idx, double y) {
            double sum = 0;
            for(int j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(y, 1.5 * sum, 1e-8);
            });

}

BOOST_AUTO_TEST_CASE(crs_matrix)
{
    std::vector<vex::command_queue> single(1, ctx.queue(0));

    const size_t n = 1024;
    const size_t m = 2048;

    std::vector<int>    row;
    std::vector<int>    col;
    std::vector<double> val;
    std::vector<double> x = random_vector<double>(m);

    random_matrix(n, m, 16, row, col, val);

    vex::backend::cuda::spmat_crs<double>
        A(ctx.queue(0), n, m, row.data(), col.data(), val.data());

    vex::vector<double> X(single, x);
    vex::vector<double> Y(single, n);

    Y = A * X;

    check_sample(Y, [&](size_t idx, double y) {
            double sum = 0;
            for(int j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(y, sum, 1e-8);
            });

    Y += 0.5 * (A * X);

    check_sample(Y, [&](size_t idx, double y) {
            double sum = 0;
            for(int j = row[idx]; j < row[idx + 1]; j++)
                sum += val[j] * x[col[j]];

            BOOST_CHECK_CLOSE(y, 1.5 * sum, 1e-8);
            });

}
#endif

BOOST_AUTO_TEST_SUITE_END()
