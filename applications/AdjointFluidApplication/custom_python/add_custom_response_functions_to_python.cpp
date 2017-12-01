// External includes
#include <boost/python.hpp>

// Project includes

// Application includes
#include "custom_response_functions/drag_response_function.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

void AddCustomResponseFunctionsToPython()
{
  
    class_<DragResponseFunction<2>, bases<ResponseFunction>>
      ("DragResponseFunction2D", init<Parameters&>());

    class_<DragResponseFunction<3>, bases<ResponseFunction>>
      ("DragResponseFunction3D", init<Parameters&>());

}

} // namespace Python

} // namespace Kratos
