#include <iostream>
#include <vexcl/devlist.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/fft.hpp>

using namespace vex;

const size_t runs = 100;
const size_t repeats = 10;


void profile(Context &ctx, std::vector<size_t> size) {
    size_t n = std::accumulate(size.begin(), size.end(),
	static_cast<size_t>(1), std::multiplies<size_t>());
    vector<cl_float2> a(ctx, n);
    vector<cl_float2> b(ctx, n);

    profiler<> prof(ctx);
    for(size_t i = 0 ; i < repeats ; i++) {
        prof.tic_cl("init");
        FFT<cl_float2> fft(ctx, size);
        fft.plan.profile = &prof;
        prof.toc("init");
        for(size_t j = 0 ; j < runs ; j++)
            b = fft(a);
    }
    std::cout << prof << std::endl;
}


int main() {
    try {
        Context ctx(Filter::Env && Filter::Count(1));
        std::cerr << ctx << std::endl;

        profile(ctx, std::vector<size_t>(2, 512));
        profile(ctx, std::vector<size_t>(2, 521));
    } catch(vex::backend::error &e) {
        std::cerr << e << std::endl;
    }

    return 0;
}
