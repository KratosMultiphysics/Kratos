#define BOOST_TEST_MODULE TensorDot
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/vector_view.hpp>
#include <vexcl/tensordot.hpp>
#include <vexcl/random.hpp>
#include <vexcl/reductor.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/function.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(axes_pairs) {
    auto ax = vex::axes_pairs(1, 2, 3, 4, 5, 6);

    BOOST_CHECK_EQUAL(ax[0][0], 1);
    BOOST_CHECK_EQUAL(ax[0][1], 2);
    BOOST_CHECK_EQUAL(ax[1][0], 3);
    BOOST_CHECK_EQUAL(ax[1][1], 4);
    BOOST_CHECK_EQUAL(ax[2][0], 5);
    BOOST_CHECK_EQUAL(ax[2][1], 6);
}

BOOST_AUTO_TEST_CASE(mat_mat) {
    using vex::_;
    using vex::extents;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    const size_t N = 32;

    vex::vector<double> a(queue, N * N);
    vex::vector<double> b(queue, N * N);

    vex::Random<double> rnd;
    vex::Reductor<double, vex::SUM> sum(queue);

    auto idx = vex::element_index();

    a = rnd(idx, 0);
    b = rnd(idx, 1);

    vex::slicer<2> dim(extents[N][N]);

    BOOST_CHECK_SMALL(
            sum(
                fabs(
                    vex::tensordot(
                        dim[_](a), dim[_](b), vex::axes_pairs(1, 0)
                        )
                    -
                    vex::reduce<vex::SUM>(
                        vex::extents[N][N][N],
                        vex::reshape(a, extents[N][N][N], extents[0][1]) *
                        vex::reshape(b, extents[N][N][N], extents[1][2]),
                        1)
                    )
               ),
            1e-8);
}

BOOST_AUTO_TEST_CASE(mat_vec) {
    using vex::_;
    using vex::extents;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    const size_t N = 32;

    vex::vector<double> a(queue, N * N);
    vex::vector<double> b(queue, N);

    vex::Random<double> rnd;
    vex::Reductor<double, vex::SUM> sum(queue);

    auto idx = vex::element_index();

    a = rnd(idx, 0);
    b = rnd(idx, 1);

    vex::slicer<2> mat(extents[N][N]);
    vex::slicer<1> vec(extents[N]);

    BOOST_CHECK_SMALL(
            sum(
                fabs(
                    vex::tensordot(
                        mat[_](a), vec[_](b), vex::axes_pairs(1, 0)
                        )
                    -
                    vex::reduce<vex::SUM>(
                        vex::extents[N][N],
                        a * vex::reshape(b, extents[N][N], extents[1]),
                        1)
                    )
               ),
            1e-8);
}

BOOST_AUTO_TEST_CASE(vec_mat) {
    using vex::_;
    using vex::extents;

    std::vector<vex::command_queue> queue(1, ctx.queue(0));

    const size_t N = 32;

    vex::vector<double> a(queue, N);
    vex::vector<double> b(queue, N * N);

    vex::Random<double> rnd;
    vex::Reductor<double, vex::SUM> sum(queue);

    auto idx = vex::element_index();

    a = rnd(idx, 0);
    b = rnd(idx, 1);

    vex::slicer<2> mat(extents[N][N]);
    vex::slicer<1> vec(extents[N]);

    BOOST_CHECK_SMALL(
            sum(
                fabs(
                    vex::tensordot(
                        vec[_](a), mat[_](b), vex::axes_pairs(0, 0)
                        )
                    -
                    vex::reduce<vex::SUM>(
                        vex::extents[N][N],
                        vex::reshape(a, extents[N][N], extents[0]) * b,
                        0)
                    )
               ),
            1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
