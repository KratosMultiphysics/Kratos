#define BOOST_TEST_MODULE TypeSystem
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>

#define CHECK_IS_NATIVE(T) \
    BOOST_CHECK(vex::is_cl_native<T>::value); \
    BOOST_CHECK(vex::is_cl_native<cl_##T    >::value); \
    BOOST_CHECK(vex::is_cl_native<cl_##T##2 >::value); \
    BOOST_CHECK(vex::is_cl_native<cl_##T##4 >::value); \
    BOOST_CHECK(vex::is_cl_native<cl_##T##8 >::value); \
    BOOST_CHECK(vex::is_cl_native<cl_##T##16>::value)

BOOST_AUTO_TEST_CASE(cl_native_types)
{
    CHECK_IS_NATIVE( char   );
    CHECK_IS_NATIVE( short  );
    CHECK_IS_NATIVE( int    );
    CHECK_IS_NATIVE( uchar  );
    CHECK_IS_NATIVE( uint   );
    CHECK_IS_NATIVE( float  );
    CHECK_IS_NATIVE( double );
}

#define CHECK_TERMINAL(name) \
    BOOST_CHECK( vex::traits::is_vector_expr_terminal<name>::value )

BOOST_AUTO_TEST_CASE(vector_expression_terminals)
{
    CHECK_TERMINAL(char );
    CHECK_TERMINAL(int  );
    CHECK_TERMINAL(float);

    CHECK_TERMINAL(vex::vector_terminal );
    CHECK_TERMINAL(vex::elem_index   );
}
