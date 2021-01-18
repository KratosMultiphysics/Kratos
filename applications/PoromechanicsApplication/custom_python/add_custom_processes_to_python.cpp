//
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_processes/apply_component_table_process.hpp"
#include "custom_processes/apply_double_table_process.hpp"
#include "custom_processes/apply_constant_hydrostatic_pressure_process.hpp"
#include "custom_processes/apply_hydrostatic_pressure_table_process.hpp"
#include "custom_processes/periodic_interface_process.hpp"


namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ApplyComponentTableProcess, ApplyComponentTableProcess::Pointer, Process>
    (m, "ApplyComponentTableProcess")
    .def( py::init< ModelPart&, Parameters>());
    py::class_<ApplyDoubleTableProcess, ApplyDoubleTableProcess::Pointer, Process>
    (m, "ApplyDoubleTableProcess")
    .def( py::init< ModelPart&, Parameters>());
    py::class_<ApplyConstantHydrostaticPressureProcess, ApplyConstantHydrostaticPressureProcess::Pointer, Process>
    (m, "ApplyConstantHydrostaticPressureProcess")
    .def( py::init< ModelPart&, Parameters>());
    py::class_<ApplyHydrostaticPressureTableProcess, ApplyHydrostaticPressureTableProcess::Pointer, Process>
    (m, "ApplyHydrostaticPressureTableProcess")
    .def( py::init< ModelPart&, Parameters>());
    py::class_<PeriodicInterfaceProcess, PeriodicInterfaceProcess::Pointer, Process>
    (m, "PeriodicInterfaceProcess")
    .def( py::init< ModelPart&, Parameters>());
}

}  // namespace Python.
} // Namespace Kratos
