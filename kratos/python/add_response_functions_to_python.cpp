// External includes
#include <boost/python.hpp>

// Project includes
#include "python/add_response_functions_to_python.h"
#include "solving_strategies/response_functions/response_function.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

void AddResponseFunctionsToPython()
{
    class_<ResponseFunction, boost::noncopyable>("ResponseFunction", no_init)
        .def("Initialize", &ResponseFunction::Initialize)
        .def("InitializeSolutionStep", &ResponseFunction::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ResponseFunction::FinalizeSolutionStep)
        .def("Check", &ResponseFunction::Check)
        .def("Clear", &ResponseFunction::Clear)
        .def("CalculateValue", &ResponseFunction::CalculateValue)
        .def("UpdateSensitivities", &ResponseFunction::UpdateSensitivities);
}

} // namespace Python

} // namespace Kratos
