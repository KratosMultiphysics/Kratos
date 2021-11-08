#include <vexcl/vexcl.hpp>
#include "thrust-sort.hpp"

template <class T>
int check(const vex::vector<T> &x, const char *who) {
    std::vector<T> y(x.size());

    vex::copy(x, y);
    if (!std::is_sorted(y.begin(), y.end())) {
        std::cerr << who << " has failed to sort a vector" << std::endl;
        return 0;
    }
    return 1;
}

int main() {
    vex::Context ctx(vex::Filter::Env && vex::Filter::Count(1));
    std::cout << ctx << std::endl;

    vex::profiler<> prof(ctx);

    typedef int T;
    const size_t n = 16 * 1024 * 1024;
    vex::vector<T> x(ctx, n);
    vex::Random<T> rnd;

    // Get raw pointers to the device memory.
    T *x_begin = x(0).raw_ptr();
    T *x_end   = x_begin + x.size();

    for(int pass = 0; pass < 11; ++pass) {
        std::cout << "." << std::flush;
        x = rnd(vex::element_index(), pass);

        if (pass) prof.tic_cl("Thrust");
        // Apply thrust algorithm.
        thrust_sort(x_begin, x_end);
        if (pass) prof.toc("Thrust");

        if (!check(x, "Thrust")) return 1;

        x = rnd(vex::element_index(), pass);
        if (pass) prof.tic_cl("VexCL");
        // Apply thrust algorithm.
        vex::sort(x);
        if (pass) prof.toc("VexCL");

        if (!check(x, "VexCL")) return 1;
    }
    std::cout << std::endl;

    std::cout << prof << std::endl;
}
