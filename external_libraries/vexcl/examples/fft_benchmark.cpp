#include <iostream>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>
#include <boost/io/ios_state.hpp>

#include <vexcl/devlist.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/fft.hpp>

using namespace vex;

typedef stopwatch<> watch;

#ifdef VEXCL_HAVE_CUDA
#  include <cufft.h>
#  include <cuda_runtime.h>

void check(cudaError_t status, const char *msg) {
    if (status != cudaSuccess)
        throw std::runtime_error(msg);
}

void check(cufftResult status, const char *msg) {
    if (status != CUFFT_SUCCESS)
        throw std::runtime_error(msg);
}

watch test_cufft(cl_float2 *data, size_t n, size_t m, size_t runs, bool) {
    size_t dataSize = sizeof(cufftComplex) * n * m;

    cufftHandle plan;
    if(m == 1)
        check(cufftPlan1d(&plan, (int)n, CUFFT_C2C, 1), "cufftPlan1d");
    else
        check(cufftPlan2d(&plan, (int)n, (int)m, CUFFT_C2C), "cufftPlan2d");

    cufftComplex *inData;
    check(cudaMalloc((void **)(&inData), dataSize), "cudaMalloc");

    cufftComplex *outData;
    check(cudaMalloc((void **)(&outData), dataSize), "cudaMalloc");

    // Send X to device
    check(cudaMemcpy(inData, data, dataSize, cudaMemcpyHostToDevice), "cudaMemcpy");

    watch w;
    for(size_t i = 0 ; i < runs ; i++) {
        cudaDeviceSynchronize();
        w.tic();
        cufftExecC2C(plan, inData, outData, CUFFT_FORWARD);
        cudaDeviceSynchronize();
        w.toc();
    }

    cufftDestroy(plan);
    cudaFree(inData);
    cudaFree(outData);
    return w;
}
#else
watch test_cufft(cl_float2 *, size_t, size_t, size_t, bool) {
    return watch();
}
#endif


#ifdef VEXCL_HAVE_FFTW
#  ifdef _OPENMP
#    include <omp.h>
#  endif
#  include <fftw3.h>

watch test_fftw(cl_float2 *data, size_t n, size_t m, size_t runs, bool dump_plan) {
    int sz[2] = {(int)n, (int)m};
    fftwf_complex *out = reinterpret_cast<fftwf_complex *>(
        fftwf_malloc(sizeof(fftwf_complex) * n * m));
    fftwf_plan p1 = fftwf_plan_dft(m == 1 ? 1 : 2, sz,
        reinterpret_cast<fftwf_complex *>(data),
        out, FFTW_FORWARD, FFTW_MEASURE);

    watch w;
    for(size_t i = 0 ; i < runs ; i++) {
        w.tic();
        fftwf_execute(p1);
        w.toc();
    }

    if(dump_plan) fftwf_fprint_plan(p1, stderr);
    fftwf_destroy_plan(p1);
    fftwf_free(out);
    return w;
}
#else
watch test_fftw(cl_float2 *, size_t, size_t, size_t, bool) {
    return watch();
}
#endif


watch test_clfft(Context &ctx, cl_float2 *data, size_t n, size_t m, size_t runs, bool dump_plan) {
    vector<cl_float2> a(ctx, n * m, data);
    vector<cl_float2> b(ctx, n * m);
    std::vector<size_t> sz; sz.push_back(n); if(m > 1) sz.push_back(m);
    FFT<cl_float2> fft(ctx, sz);

    watch w;
    for(size_t i = 0 ; i < runs ; i++) {
        ctx.finish();
        w.tic();
        b = fft(a);
        ctx.finish();
        w.toc();
    }

    if(dump_plan) std::cerr << fft.plan;
    return w;
}

void info(const watch &w, size_t size, size_t dim) {
    boost::io::ios_all_saver stream_state(std::cout);

    // FFT is O(n log n)
    double ops = dim == 1
        ? size * std::log(static_cast<double>(size)) // O(n log n)
        : 2.0 * size * size * std::log(static_cast<double>(size)); // O(n log n)[1D fft] * n[rows] * 2[transposed]
    std::cout << '\t';
    if(w.tics() == 0)
        std::cout << '-';
    else
        std::cout << std::scientific << (ops / w.average());
}

int main(int argc, char **argv) {
    using namespace boost::program_options;
    options_description desc("Options");
    bool add_prime = false, dump_plan = false;
    size_t composite = 0, min, max, runs;
    desc.add_options()
        ("help,h", "show help")
        ("runs,r", value(&runs)->default_value(1000),
            "run each FFT multiple times")
        ("min", value(&min)->default_value(1 << 4),
            "minimum size")
        ("max", value(&max)->default_value(1 << 20),
            "maximum size (all powers of two between)")
        ("prime,p", bool_switch(&add_prime),
            "add nearest prime")
        ("composite,c", value(&composite)->implicit_value(23),
            "add nearest multiple of this number (for each power of two)")
        ("dump-plan", bool_switch(&dump_plan),
            "show each CLFFT plan on stderr");
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if(vm.count("help") || min > max) {
        std::cerr << desc << std::endl;
        return EXIT_FAILURE;
    }

#if defined(_OPENMP) && defined(VEXCL_HAVE_FFTW)
    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
    Context ctx(Filter::Env && Filter::Count(1));
    std::cerr << ctx << std::endl;

    // power of two sizes
    std::vector<size_t> twos, ns;
    for(size_t n = min ; n <= max ; n *= 2) {
        twos.push_back(n);
        ns.push_back(n);
    }

    if(add_prime) {
        vex::fft::prime_generator prime;
        size_t prev = 0, next = 0;
        for(auto n = twos.begin() ; n != twos.end() ; n++) {
            while(next <= *n) { prev = next; next = prime(); }
            const size_t p = *n - prev < next - *n ? prev : next;
            ns.push_back(p);
        }
    }

    if(composite != 0) {
        for(auto n = twos.begin() ; n != twos.end() ; n++) {
            const size_t m = static_cast<size_t>(*n / static_cast<double>(composite) + 0.5);
            if(m > 1) ns.push_back(composite * m);
        }
    }

    const size_t max_len = *std::max_element(ns.begin(), ns.end());

    // random data
#ifdef VEXCL_HAVE_FFTW
    cl_float2 *data = reinterpret_cast<cl_float2 *>(
        fftwf_malloc(sizeof(cl_float2) * max_len));
#else
    cl_float2 *data = new cl_float2[max_len];
#endif
    std::minstd_rand gen;
    std::uniform_real_distribution<float> dist(-1000, 1000);
    for(size_t i = 0 ; i < 2 * max_len ; i++)
        reinterpret_cast<float*>(data)[i] = dist(gen);

    // 1D
    std::cout << "# prints `n log n / time`\n";
    std::cout << "#n\tfftw^1\tclfft^1\tcufft^1" << std::endl;
    for(auto n = ns.begin() ; n != ns.end() ; n++) {
        std::cout << *n;
        info(test_fftw (     data, *n, 1, runs, dump_plan), *n, 1);
        info(test_clfft(ctx, data, *n, 1, runs, dump_plan), *n, 1);
        info(test_cufft(     data, *n, 1, runs, dump_plan), *n, 1);
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // 2D
    std::cout << "# prints `2 n^2 log n / time`\n";
    std::cout << "#n\tfftw^2\tclfft^2\tcufft^2" << std::endl;
    for(auto n = ns.begin() ; n != ns.end() ; n++)
        if(*n * *n <= max_len) {
            std::cout << *n;
            info(test_fftw (     data, *n, *n, runs, dump_plan), *n, 2);
            info(test_clfft(ctx, data, *n, *n, runs, dump_plan), *n, 2);
            info(test_cufft(     data, *n, *n, runs, dump_plan), *n, 2);
            std::cout << std::endl;
        }

    return 0;
}
