#define BOOST_TEST_MODULE FastFourierTransform
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/fft.hpp>
#include <vexcl/random.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/reductor.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(transform_expression)
{
    const size_t N = 1024;
    std::vector<vex::backend::command_queue> queue(1, ctx.queue(0));

    vex::vector<cl_float> data(queue, N);
    vex::FFT<cl_float> fft(queue, N);

    // should compile
    data += fft(data * data) * 5;
}

BOOST_AUTO_TEST_CASE(check_correctness)
{
#if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
    // Apple fails this test on CPUs
#   if defined(__APPLE__)
    if (vex::Filter::CPU(ctx.device(0)))
        return;
#   endif
#endif

    const size_t N = 1024;
    std::vector<vex::backend::command_queue> queue(1, ctx.queue(0));

    vex::vector<cl_double>  in  (queue, N);
    vex::vector<cl_double2> out (queue, N);
    vex::vector<cl_double>  back(queue, N);

    vex::Random<cl_double> rnd;

    in = rnd(vex::element_index(), std::rand());

    vex::FFT<cl_double,  cl_double2> fft (queue, N);
    vex::FFT<cl_double2, cl_double > ifft(queue, N, vex::fft::inverse);

    out  = fft (in );
    back = ifft(out);

    vex::Reductor<cl_double, vex::SUM> sum(queue);

    BOOST_CHECK(std::sqrt(sum(pow(in - back, 2.0)) / N) < 1e-3);
}

void test(const vex::Context &ctx, std::vector<size_t> ns, size_t batch) {
    std::cout << "FFT(C2C) size=" << ns[0];
    for(size_t i = 1; i < ns.size(); i++) std::cout << 'x' << ns[i];
    std::cout << " batch=" << batch << std::endl;

    std::vector<vex::backend::command_queue> queue(1, ctx.queue(0));

    size_t n1 = std::accumulate(ns.begin(), ns.end(), static_cast<size_t>(1), std::multiplies<size_t>());
    size_t n = n1 * batch;

    // random data.
    std::vector<cl_double2> inp_h = random_vector<cl_double2>(n);

    // test
    vex::vector<cl_double2> inp (queue, inp_h);
    vex::vector<cl_double2> out (queue, n);
    vex::vector<cl_double2> back(queue, n);

    std::vector<size_t> ns_(ns.begin(), ns.end());
    std::vector<vex::fft::direction> dirs (ns.size(), vex::fft::forward);
    std::vector<vex::fft::direction> idirs(ns.size(), vex::fft::inverse);
    if(batch != 1) {
        ns_.insert(ns_.begin(), batch);
        dirs.insert(dirs.begin(), vex::fft::none);
        idirs.insert(idirs.begin(), vex::fft::none);
    }
    vex::FFT<cl_double2> fft (queue, ns_, dirs);
    vex::FFT<cl_double2> ifft(queue, ns_, idirs);

    out  = fft (inp);
    back = ifft(out);

    auto rms = [&](const vex::vector<cl_double2> &a, const vex::vector<cl_double2> &b) -> double {
        // Required for CUDA compatibility:
        VEX_FUNCTION(cl_double2, minus, (cl_double2, a)(cl_double2, b),
                double2 r = {a.x - b.x, a.y - b.y};
                return r;
                );

        VEX_FUNCTION(double, dot2, (cl_double2, a)(cl_double2, b),
                return a.x * b.x + a.y * b.y;
                );

        static vex::Reductor<double, vex::SUM> sum(queue);
        return std::sqrt(sum(dot2(minus(a, b), minus(a, b)))) /
            std::sqrt(sum(dot2(b, b)));
    };

    BOOST_CHECK_SMALL(rms(back, inp), 1e-8);
}

// random dimension, mostly 1.
size_t random_dim(double p, double s) {
    static std::default_random_engine rng( std::rand() );
    static std::uniform_real_distribution<double> rnd(0.0, 1.0);

    return 1 + static_cast<size_t>( s * std::pow(rnd(rng), p) );
}

BOOST_AUTO_TEST_CASE(test_dimensions)
{
#if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
    // TODO: POCL fails this test.
    if (vex::Filter::Platform("Portable Computing Language")(ctx.device(0)))
        return;

    // Apple fails this test on CPUs
#if defined(__APPLE__)
    if (vex::Filter::CPU(ctx.device(0)))
        return;
#endif
#endif

    const size_t max = 1 << 12;

    vex::fft::planner p;

    for(size_t i = 0; i < 100; ++i) {
        // random number of dimensions, mostly 1.
        size_t dims = random_dim(3, 5);
        size_t batch = random_dim(5, 100);

        // random size.
        std::vector<size_t> n;
        size_t d_max = static_cast<size_t>(std::pow(static_cast<double>(max), 1.0 / dims));
        size_t total = batch;
        for(size_t d = 0 ; d < dims ; d++) {
            size_t sz = random_dim(dims == 1 ? 3 : 1, static_cast<double>(d_max));

            if(rand() % 3 != 0)
                sz = p.best_size(sz);

            n.push_back(sz);
            total *= sz;
        }

        // run
        if(total <= max) test(ctx, n, batch);
    }
}

BOOST_AUTO_TEST_SUITE_END()
