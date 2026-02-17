
//
//  Main authors:    Aron Noordam
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/uvec_utilities.h"

namespace Kratos
{

namespace Python
{
using namespace pybind11;

  void  AddCustomUtilitiesToPython(const pybind11::module& m)
  {
    class_<UvecUtilities, typename UvecUtilities::Pointer>(m, "UvecUtilities")
    .def_static("SetLoadOnCondition",&UvecUtilities::SetLoadOnCondition)
    .def_static("GetMovingConditionVariable", &UvecUtilities::GetMovingConditionVariable);

  }

}  // namespace Python.

} // Namespace Kratos