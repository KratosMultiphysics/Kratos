#define BOOST_TEST_MODULE Logical
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/logical.hpp>
#include "context_setup.hpp"

BOOST_AUTO_TEST_CASE(logical)
{
    const size_t N = 1024;

    vex::vector<int> x(ctx, N);
    x = vex::element_index();

    vex::any_of any_of(ctx);
    vex::all_of all_of(ctx);

    BOOST_CHECK( any_of(x)       );
    BOOST_CHECK(!any_of(0 * x)   );
    BOOST_CHECK( any_of(x > N/2) );
    BOOST_CHECK(!any_of(x < 0)   );

    BOOST_CHECK(!all_of(x)           );
    BOOST_CHECK( all_of((x + 1) > 0) );
    BOOST_CHECK(!all_of(x > N/2)     );
}

BOOST_AUTO_TEST_SUITE_END()
