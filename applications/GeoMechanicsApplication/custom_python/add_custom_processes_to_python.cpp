// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
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
#include "custom_processes/apply_constant_boundary_hydrostatic_pressure_process.hpp"
#include "custom_processes/apply_boundary_hydrostatic_pressure_table_process.hpp"
#include "custom_processes/apply_constant_phreatic_line_pressure_process.hpp"
#include "custom_processes/apply_phreatic_line_pressure_table_process.hpp"
#include "custom_processes/apply_constant_boundary_phreatic_line_pressure_process.hpp"
#include "custom_processes/apply_boundary_phreatic_line_pressure_table_process.hpp"
#include "custom_processes/apply_constant_phreatic_surface_pressure_process.hpp"
#include "custom_processes/apply_constant_interpolate_line_pressure_process.hpp"
#include "custom_processes/apply_time_dependent_interpolate_line_pressure_process.hpp"
#include "custom_processes/apply_phreatic_surface_pressure_table_process.hpp"
#include "custom_processes/apply_constant_boundary_phreatic_surface_pressure_process.hpp"
#include "custom_processes/apply_boundary_phreatic_surface_pressure_table_process.hpp"
#include "custom_processes/apply_excavation_process.hpp"
#include "custom_processes/apply_gradual_excavation_process.hpp"
#include "custom_processes/apply_write_result_scalar_process.hpp"
#include "custom_processes/gap_closure_interface_process.hpp"

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

        class_<ApplyBoundaryHydrostaticPressureTableProcess, ApplyBoundaryHydrostaticPressureTableProcess::Pointer, Process>
            (m, "ApplyBoundaryHydrostaticPressureTableProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyConstantBoundaryHydrostaticPressureProcess, ApplyConstantBoundaryHydrostaticPressureProcess::Pointer, Process>
            (m, "ApplyConstantBoundaryHydrostaticPressureProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyConstantPhreaticLinePressureProcess, ApplyConstantPhreaticLinePressureProcess::Pointer, Process>
            (m, "ApplyConstantPhreaticLinePressureProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyConstantInterpolateLinePressureProcess, ApplyConstantInterpolateLinePressureProcess::Pointer, Process>
            (m, "ApplyConstantInterpolateLinePressureProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyTimeDependentInterpolateLinePressureProcess, ApplyTimeDependentInterpolateLinePressureProcess::Pointer, Process>
            (m, "ApplyTimeDependentInterpolateLinePressureProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyPhreaticLinePressureTableProcess, ApplyPhreaticLinePressureTableProcess::Pointer, Process>
            (m, "ApplyPhreaticLinePressureTableProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyBoundaryPhreaticLinePressureTableProcess, ApplyBoundaryPhreaticLinePressureTableProcess::Pointer, Process>
            (m, "ApplyBoundaryPhreaticLinePressureTableProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyConstantBoundaryPhreaticLinePressureProcess, ApplyConstantBoundaryPhreaticLinePressureProcess::Pointer, Process>
            (m, "ApplyConstantBoundaryPhreaticLinePressureProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyConstantPhreaticSurfacePressureProcess, ApplyConstantPhreaticSurfacePressureProcess::Pointer, Process>
            (m, "ApplyConstantPhreaticSurfacePressureProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyPhreaticSurfacePressureTableProcess, ApplyPhreaticSurfacePressureTableProcess::Pointer, Process>
            (m, "ApplyPhreaticSurfacePressureTableProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyBoundaryPhreaticSurfacePressureTableProcess, ApplyBoundaryPhreaticSurfacePressureTableProcess::Pointer, Process>
            (m, "ApplyBoundaryPhreaticSurfacePressureTableProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyConstantBoundaryPhreaticSurfacePressureProcess, ApplyConstantBoundaryPhreaticSurfacePressureProcess::Pointer, Process>
            (m, "ApplyConstantBoundaryPhreaticSurfacePressureProcess")
            .def(init < ModelPart&, Parameters>());

        class_<ApplyExcavationProcess, ApplyExcavationProcess::Pointer, Process>
            (m, "ApplyExcavationProcess")
            .def(init < ModelPart&, Parameters&>());

        class_<ApplyGradualExcavationProcess, ApplyGradualExcavationProcess::Pointer, Process>
            (m, "ApplyGradualExcavationProcess")
            .def(init < ModelPart&, Parameters&>());

        class_<ApplyWriteScalarProcess, ApplyWriteScalarProcess::Pointer, Process>
            (m, "ApplyWriteScalarProcess")
            .def(init < ModelPart&, Parameters&>());

        class_<GapClosureInterfaceProcess, GapClosureInterfaceProcess::Pointer, Process>
            (m, "GapClosureInterfaceProcess")
            .def(init < ModelPart&, Parameters&>());

    }
}  // namespace Python.
} // Namespace Kratos
