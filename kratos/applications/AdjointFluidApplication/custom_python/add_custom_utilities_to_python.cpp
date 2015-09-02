// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/time_averaged_primal_utility.h"

namespace Kratos
{

namespace Python
{
using namespace boost::python;

void AddCustomUtilitiesToPython()
{
  class_<TimeAveragedPrimalUtility>("TimeAveragedPrimalUtility", init<ModelPart&>())
    .def("Reset",&TimeAveragedPrimalUtility::Reset)
    .def("AddStep",&TimeAveragedPrimalUtility::AddStep)
    ;
}

}  // namespace Python

} // namespace Kratos
