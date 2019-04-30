#define BOOST_TEST_MODULE MultipleObjects
#include <boost/test/unit_test.hpp>
#include <vexcl/vexcl.hpp>

bool empty_context();

BOOST_AUTO_TEST_CASE(multiple_objects)
{
    vex::Context ctx( vex::Filter::Env );

    BOOST_REQUIRE(!empty_context());

    std::cout << vex::current_context() << std::endl;
}
