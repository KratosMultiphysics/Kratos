#include <enoki/python.h>
#include <enoki/complex.h>
#include <enoki/dynamic.h>

using namespace enoki;

Array<float, 4> a1() {
    auto result = Array<float, 4>(1, 2, 3, 4);
    std::cout << result << std::endl;
    return result;
}

Array<Packet<float>, 4> a2() {
    auto result = Array<Packet<float>, 4>(1, 2, 3, 4);
    std::cout << result << std::endl;
    return result;
}

void a3(Array<float, 4> value) {
    std::cout << value << std::endl;
}

void a4(Array<Packet<float>, 4> value) {
    std::cout << value << std::endl;
}

Complex<float> c1() {
    auto result = Complex<float>(1, 2);
    std::cout << result << std::endl;
    return result;
}

Complex<Packet<float>> c2() {
    auto result = Complex<Packet<float>>(1, 2);
    std::cout << result << std::endl;
    return result;
}

Complex<Array<Packet<float>, 4>> c2_b() {
    auto result = Complex<Array<Packet<float>, 4>>(1.f, 2.f);
    std::cout << result << std::endl;
    return result;
}


void c3(Complex<float> value) {
    std::cout << value << std::endl;
}

void c4(Complex<Packet<float>> value) {
    std::cout << value << std::endl;
}

template <typename Float> Float atan(Float x) {
    return enoki::atan(x);
}

PYBIND11_MODULE(test, m) {
    /* Real */
    m.def("a1", &a1);
    m.def("a2", &a2);
    m.def("a3", &a3);
    m.def("a4", &a4);

    /* Complex */
    m.def("c1", &c1);
    m.def("c2", &c2);
    m.def("c2_b", &c2_b);
    m.def("c3", &c3);
    m.def("c4", &c4);

    using FloatP = Packet<float>;
    m.def("atan", enoki::vectorize_wrapper(&atan<FloatP>));
}

