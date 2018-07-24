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
    using namespace pybind11;

    class_<ApplyComponentTableProcess, ApplyComponentTableProcess::Pointer, Process>
    (m, "ApplyComponentTableProcess")
    .def(init < ModelPart&, Parameters>());
    class_<ApplyDoubleTableProcess, ApplyDoubleTableProcess::Pointer, Process>
    (m, "ApplyDoubleTableProcess")
    .def(init < ModelPart&, Parameters>());
    class_<ApplyConstantHydrostaticPressureProcess, ApplyConstantHydrostaticPressureProcess::Pointer, Process>
    (m, "ApplyConstantHydrostaticPressureProcess")
    .def(init < ModelPart&, Parameters>());
    class_<ApplyHydrostaticPressureTableProcess, ApplyHydrostaticPressureTableProcess::Pointer, Process>
    (m, "ApplyHydrostaticPressureTableProcess")
    .def(init < ModelPart&, Parameters>());
    class_<PeriodicInterfaceProcess, PeriodicInterfaceProcess::Pointer, Process>
    (m, "PeriodicInterfaceProcess")
    .def(init < ModelPart&, Parameters>());
}

}  // namespace Python.
} // Namespace Kratos
