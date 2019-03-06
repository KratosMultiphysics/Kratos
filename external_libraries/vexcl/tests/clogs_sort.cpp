#define BOOST_TEST_MODULE ClogsSort
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/external/clogs.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(sort_keys)
{
    const size_t n = 1000 * 1000;

    std::vector<cl_uint> k = random_vector<cl_uint>(n);
    vex::vector<cl_uint> keys(ctx, k);

    vex::clogs::sort(keys);
    std::sort(k.begin(), k.end());

    check_sample(keys, [&](size_t idx, cl_uint a) {
        BOOST_CHECK_EQUAL(a, k[idx]);
    });
}

BOOST_AUTO_TEST_CASE(sort_keys_vals)
{
    struct kv_pair {
        cl_ulong key;
        cl_float2 value;

        bool operator<(const kv_pair &other) const {
            return key < other.key;
        }
    };

    const size_t n = 1000 * 1000;

    std::vector<cl_ulong> k = random_vector<ulong>(n);
    std::vector<cl_float2> v = random_vector<cl_float2>(n);
    vex::vector<cl_ulong> keys(ctx, k);
    vex::vector<cl_float2> vals(ctx, v);

    std::vector<kv_pair> kv(n);
    for (size_t i = 0; i < n; i++)
        kv[i] = kv_pair{k[i], v[i]};
    std::stable_sort(kv.begin(), kv.end());

    vex::clogs::stable_sort_by_key(keys, vals);

    check_sample(keys, [&](size_t idx, cl_ulong a) {
        BOOST_CHECK_EQUAL(a, kv[idx].key);
    });
    check_sample(vals, [&](size_t idx, cl_float2 a) {
        for (int i = 0; i < 2; i++)
            BOOST_CHECK_EQUAL(a.s[i], kv[idx].value.s[i]);
    });
}

BOOST_AUTO_TEST_SUITE_END()
