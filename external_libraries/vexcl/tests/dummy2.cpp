/* Dummy file to test for multiple definition problems */
#include <vexcl/vexcl.hpp>

bool empty_context() {
    return vex::current_context().empty();
}
