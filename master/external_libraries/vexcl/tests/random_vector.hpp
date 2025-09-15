#ifndef TESTS_RANDOM_VECTOR_HPP
#define TESTS_RANDOM_VECTOR_HPP

#include <random>
#include <vector>
#include <vexcl/types.hpp>

template <typename T, class Enable = void>
struct generator {};

template<typename T>
struct generator<T, typename std::enable_if<std::is_floating_point<T>::value>::type>
{
    static T get() {
        static std::default_random_engine rng( std::rand() );
        static std::uniform_real_distribution<T> rnd((T)0, (T)1);
        return rnd(rng);
    }
};

template<typename T>
struct generator<T, typename std::enable_if<std::is_integral<T>::value>::type>
{
    static T get() {
        static std::default_random_engine rng( std::rand() );
        static std::uniform_int_distribution<T> rnd(0, 100);
        return rnd(rng);
    }
};

template<typename T>
struct generator<T, typename std::enable_if<vex::is_cl_vector<T>::value>::type>
{
    static T get() {
        T r;
        for (unsigned i = 0; i < vex::cl_vector_length<T>::value; i++)
            r.s[i] = ::generator<typename vex::cl_scalar_of<T>::type>::get();
        return r;
    }
};

template<class T>
std::vector<T> random_vector(size_t n) {
    std::vector<T> x(n);
    std::generate(x.begin(), x.end(), ::generator<T>::get);
    return x;
}



#endif
