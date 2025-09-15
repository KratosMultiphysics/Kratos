#define BOOST_TEST_MODULE VexIO
#include <boost/test/unit_test.hpp>
#if BOOST_VERSION >= 107100
#  include <boost/test/tools/output_test_stream.hpp>
#else
#  include <boost/test/output_test_stream.hpp>
#endif
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(stream_vector) {
    boost::test_tools::output_test_stream output;

    vex::vector<int> x(ctx, 32);
    x = vex::element_index();

    output << x;

    BOOST_CHECK(output.is_equal(
        "{\n"
        "     0:      0      1      2      3      4      5      6      7      8      9\n"
        "    10:     10     11     12     13     14     15     16     17     18     19\n"
        "    20:     20     21     22     23     24     25     26     27     28     29\n"
        "    30:     30     31\n"
        "}\n"
        ));
}

BOOST_AUTO_TEST_SUITE_END()
