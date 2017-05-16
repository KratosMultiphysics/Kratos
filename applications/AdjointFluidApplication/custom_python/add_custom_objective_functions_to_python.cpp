// External includes
#include <boost/python.hpp>

// Project includes

// Application includes
#include "custom_utilities/objective_function.h"
#include "custom_utilities/drag_objective_function.h"
#include "custom_python/add_custom_objective_functions_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

void AddCustomObjectiveFunctionsToPython()
{
    class_<ObjectiveFunction, boost::noncopyable>("ObjectiveFunction", init<>())
        .def("Initialize", &ObjectiveFunction::Initialize)
        .def("InitializeSolutionStep", &ObjectiveFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ObjectiveFunction::FinalizeSolutionStep)
        .def("Check", &ObjectiveFunction::Check)
        .def("CalculateAdjointVelocityContribution",
             &ObjectiveFunction::CalculateAdjointVelocityContribution)
        .def("CalculateAdjointAccelerationContribution",
             &ObjectiveFunction::CalculateAdjointAccelerationContribution)
        .def("CalculateSensitivityContribution", &ObjectiveFunction::CalculateSensitivityContribution)
        .def("Calculate", &ObjectiveFunction::Calculate);

    class_<DragObjectiveFunction<2>, bases<ObjectiveFunction>, boost::noncopyable>
    ("DragObjectiveFunction2D", init<Parameters&>());

    class_<DragObjectiveFunction<3>, bases<ObjectiveFunction>, boost::noncopyable>
    ("DragObjectiveFunction3D", init<Parameters&>());

}

} // namespace Python

} // namespace Kratos
