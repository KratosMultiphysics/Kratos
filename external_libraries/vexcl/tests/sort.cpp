#define BOOST_TEST_MODULE Sort
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/sort.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(sort_keys)
{
    const size_t n = 1000 * 1000;

    std::vector<float> k = random_vector<float>(n);
    vex::vector<float> keys(ctx, k);

    vex::sort(keys);
    vex::copy(keys, k);

    BOOST_CHECK( std::is_sorted(k.begin(), k.end()) );
}

BOOST_AUTO_TEST_CASE(sort_keys_vals_default)
{
    const size_t n = 1000 * 1000;

    std::vector<int  > k = random_vector<int  >(n);
    std::vector<float> v = random_vector<float>(n);
    std::vector<int>   p(n);

    vex::vector<int  > keys(ctx, k);
    vex::vector<float> vals(ctx, v);

    for(size_t i = 0; i < p.size(); ++i) p[i] = static_cast<int>(i);
    std::stable_sort(p.begin(), p.end(), [&](int i, int j) { return k[i] < k[j]; });

    vex::sort_by_key(keys, vals);

    check_sample(keys, [&](size_t pos, int val) {
            BOOST_CHECK_EQUAL(val, k[p[pos]]);
            });

    check_sample(vals, [&](size_t pos, float val) {
            BOOST_CHECK_EQUAL(val, v[p[pos]]);
            });
}

BOOST_AUTO_TEST_CASE(sort_keys_vals_custom_op)
{
    const size_t n = 1000 * 1000;

    std::vector<int  > k = random_vector<int  >(n);
    std::vector<float> v = random_vector<float>(n);
    std::vector<int>   p(n);

    vex::vector<int  > keys(ctx, k);
    vex::vector<float> vals(ctx, v);

    for(size_t i = 0; i < p.size(); ++i) p[i] = static_cast<int>(i);

    struct even_first_t {
        typedef bool result_type;
        even_first_t() {}

        VEX_DUAL_FUNCTOR(result_type, (int, a)(int, b),
            char bit1 = 1 & a;
            char bit2 = 1 & b;
            if (bit1 == bit2) return a < b;
            return bit1 < bit2;
            )
    } even_first;

    std::stable_sort(p.begin(), p.end(), [&](int i, int j) -> bool {
            int a = k[i];
            int b = k[j];
            return even_first(a, b);
            });

    vex::sort_by_key(keys, vals, even_first);

    check_sample(keys, [&](size_t pos, int val) {
            BOOST_CHECK_EQUAL(val, k[p[pos]]);
            });

    check_sample(vals, [&](size_t pos, float val) {
            BOOST_CHECK_EQUAL(val, v[p[pos]]);
            });
}

BOOST_AUTO_TEST_CASE(sort_keys_tuple)
{
    const size_t n = 1000 * 1000;

    std::vector<int>   k1 = random_vector<int>  (n);
    std::vector<float> k2 = random_vector<float>(n);

    vex::vector<int>   keys1(ctx, k1);
    vex::vector<float> keys2(ctx, k2);

    struct less_t {
        typedef bool result_type;
        less_t() {}

        VEX_DUAL_FUNCTOR(result_type, (int, a1)(float, a2)(int, b1)(float, b2),
            return (a1 == b1) ? (a2 < b2) : (a1 < b1);
            )
    } less;

    vex::sort(boost::fusion::vector_tie(keys1, keys2), less );
    vex::copy(keys1, k1);
    vex::copy(keys2, k2);

    BOOST_CHECK( std::is_sorted(
                boost::counting_iterator<size_t>(0),
                boost::counting_iterator<size_t>(n),
                [&](size_t i, size_t j) {
                    return std::make_tuple(k1[i], k2[i]) < std::make_tuple(k1[j], k2[j]);
                } ) );
}

BOOST_AUTO_TEST_CASE(sort_keys_vals_tuple)
{
    const size_t n = 1000 * 1000;

    std::vector<int>     k1 = random_vector<int>    (n);
    std::vector<float>   k2 = random_vector<float>  (n);
    std::vector<cl_long> v1 = random_vector<cl_long>(n);
    std::vector<short>   v2 = random_vector<short>  (n);

    vex::vector<int>     keys1(ctx, k1);
    vex::vector<float>   keys2(ctx, k2);
    vex::vector<cl_long> vals1(ctx, v1);
    vex::vector<short>   vals2(ctx, v2);

    std::vector<int> p(n);
    for(size_t i = 0; i < p.size(); ++i) p[i] = static_cast<int>(i);

    struct less_t {
        typedef bool result_type;
        less_t() {}

        VEX_DUAL_FUNCTOR(result_type, (int, a1)(float, a2)(int, b1)(float, b2),
            return (a1 == b1) ? (a2 < b2) : (a1 < b1);
            )
    } less;

    std::stable_sort(p.begin(), p.end(), [&](int i, int j) {
            return less(k1[i], k2[i], k1[j], k2[j]);
            });

    vex::sort_by_key(
            boost::fusion::vector_tie(keys1, keys2),
            boost::fusion::vector_tie(vals1, vals2),
            less
            );

    check_sample(keys1, [&](size_t pos, int val) {
            BOOST_CHECK_EQUAL(val, k1[p[pos]]);
            });

    check_sample(keys2, [&](size_t pos, float val) {
            BOOST_CHECK_EQUAL(val, k2[p[pos]]);
            });

    check_sample(vals1, [&](size_t pos, cl_long val) {
            BOOST_CHECK_EQUAL(val, v1[p[pos]]);
            });

    check_sample(vals2, [&](size_t pos, short val) {
            BOOST_CHECK_EQUAL(val, v2[p[pos]]);
            });
}

BOOST_AUTO_TEST_SUITE_END()
