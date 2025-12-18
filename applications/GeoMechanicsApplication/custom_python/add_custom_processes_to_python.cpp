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
#include "custom_python/add_custom_processes_to_python.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

#include "custom_processes/apply_boundary_hydrostatic_pressure_table_process.h"
#include "custom_processes/apply_boundary_phreatic_line_pressure_table_process.h"
#include "custom_processes/apply_boundary_phreatic_surface_pressure_table_process.h"
#include "custom_processes/apply_c_phi_reduction_process.h"
#include "custom_processes/apply_component_table_process.h"
#include "custom_processes/apply_constant_boundary_hydrostatic_pressure_process.h"
#include "custom_processes/apply_constant_boundary_phreatic_line_pressure_process.h"
#include "custom_processes/apply_constant_boundary_phreatic_surface_pressure_process.h"
#include "custom_processes/apply_constant_hydrostatic_pressure_process.h"
#include "custom_processes/apply_constant_interpolate_line_pressure_process.h"
#include "custom_processes/apply_constant_phreatic_line_pressure_process.h"
#include "custom_processes/apply_constant_phreatic_multi_line_pressure_process.h"
#include "custom_processes/apply_constant_phreatic_surface_pressure_process.h"
#include "custom_processes/apply_excavation_process.h"
#include "custom_processes/apply_final_stresses_of_previous_stage_to_initial_state.h"
#include "custom_processes/apply_hydrostatic_pressure_table_process.h"
#include "custom_processes/apply_initial_uniform_stress_field.h"
#include "custom_processes/apply_k0_procedure_process.h"
#include "custom_processes/apply_normal_load_table_process.h"
#include "custom_processes/apply_phreatic_line_pressure_table_process.h"
#include "custom_processes/apply_phreatic_multi_line_pressure_table_process.h"
#include "custom_processes/apply_phreatic_surface_pressure_table_process.h"
#include "custom_processes/apply_scalar_constraint_table_process.h"
#include "custom_processes/apply_vector_constraint_table_process.h"
#include "custom_processes/apply_write_result_scalar_process.h"
#include "custom_processes/calculate_incremental_motion_process.h"
#include "custom_processes/calculate_total_motion_process.h"
#include "custom_processes/deactivate_conditions_on_inactive_elements_process.hpp"
#include "custom_processes/find_neighbour_elements_of_conditions_process.h"
#include "custom_processes/find_neighbours_of_interfaces_process.h"
#include "custom_processes/geo_extrapolate_integration_point_values_to_nodes_process.h"
#include "custom_processes/set_absorbing_boundary_parameters_process.hpp"
#include "custom_processes/set_multiple_moving_loads.h"
#include "custom_processes/set_parameter_field_process.hpp"

namespace Kratos::Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ApplyComponentTableProcess, ApplyComponentTableProcess::Pointer, Process>(
        m, "ApplyComponentTableProcess", py::module_local())
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantHydrostaticPressureProcess, ApplyConstantHydrostaticPressureProcess::Pointer, Process>(
        m, "ApplyConstantHydrostaticPressureProcess", py::module_local())
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyHydrostaticPressureTableProcess, ApplyHydrostaticPressureTableProcess::Pointer, Process>(
        m, "ApplyHydrostaticPressureTableProcess", py::module_local())
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyBoundaryHydrostaticPressureTableProcess, ApplyBoundaryHydrostaticPressureTableProcess::Pointer, Process>(
        m, "ApplyBoundaryHydrostaticPressureTableProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantBoundaryHydrostaticPressureProcess, ApplyConstantBoundaryHydrostaticPressureProcess::Pointer, Process>(
        m, "ApplyConstantBoundaryHydrostaticPressureProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantPhreaticLinePressureProcess, ApplyConstantPhreaticLinePressureProcess::Pointer, Process>(
        m, "ApplyConstantPhreaticLinePressureProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantInterpolateLinePressureProcess, ApplyConstantInterpolateLinePressureProcess::Pointer, Process>(
        m, "ApplyConstantInterpolateLinePressureProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyPhreaticLinePressureTableProcess, ApplyPhreaticLinePressureTableProcess::Pointer, Process>(
        m, "ApplyPhreaticLinePressureTableProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantPhreaticMultiLinePressureProcess, ApplyConstantPhreaticMultiLinePressureProcess::Pointer, Process>(
        m, "ApplyConstantPhreaticMultiLinePressureProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyBoundaryPhreaticLinePressureTableProcess, ApplyBoundaryPhreaticLinePressureTableProcess::Pointer, Process>(
        m, "ApplyBoundaryPhreaticLinePressureTableProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantBoundaryPhreaticLinePressureProcess, ApplyConstantBoundaryPhreaticLinePressureProcess::Pointer, Process>(
        m, "ApplyConstantBoundaryPhreaticLinePressureProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantPhreaticSurfacePressureProcess, ApplyConstantPhreaticSurfacePressureProcess::Pointer, Process>(
        m, "ApplyConstantPhreaticSurfacePressureProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyPhreaticMultiLinePressureTableProcess, ApplyPhreaticMultiLinePressureTableProcess::Pointer, Process>(
        m, "ApplyPhreaticMultiLinePressureTableProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyPhreaticSurfacePressureTableProcess, ApplyPhreaticSurfacePressureTableProcess::Pointer, Process>(
        m, "ApplyPhreaticSurfacePressureTableProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyBoundaryPhreaticSurfacePressureTableProcess, ApplyBoundaryPhreaticSurfacePressureTableProcess::Pointer, Process>(
        m, "ApplyBoundaryPhreaticSurfacePressureTableProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyConstantBoundaryPhreaticSurfacePressureProcess, ApplyConstantBoundaryPhreaticSurfacePressureProcess::Pointer, Process>(
        m, "ApplyConstantBoundaryPhreaticSurfacePressureProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyExcavationProcess, ApplyExcavationProcess::Pointer, Process>(
        m, "ApplyExcavationProcess")
        .def(py::init<Model&, const Parameters&>());

    py::class_<ApplyWriteScalarProcess, ApplyWriteScalarProcess::Pointer, Process>(
        m, "ApplyWriteScalarProcess")
        .def(py::init<ModelPart&, Parameters&>());

    py::class_<ApplyK0ProcedureProcess, ApplyK0ProcedureProcess::Pointer, Process>(
        m, "ApplyK0ProcedureProcess")
        .def(py::init<Model&, Parameters>());

    py::class_<FindNeighbourElementsOfConditionsProcess, FindNeighbourElementsOfConditionsProcess::Pointer, Process>(
        m, "FindNeighbourElementsOfConditionsProcess")
        .def(py::init<ModelPart&>());

    py::class_<DeactivateConditionsOnInactiveElements, DeactivateConditionsOnInactiveElements::Pointer, Process>(
        m, "DeactivateConditionsOnInactiveElements")
        .def(py::init<ModelPart&>());

    py::class_<SetAbsorbingBoundaryParametersProcess, SetAbsorbingBoundaryParametersProcess::Pointer, Process>(
        m, "SetAbsorbingBoundaryParametersProcess")
        .def(py::init<ModelPart&, Parameters&>());

    py::class_<SetParameterFieldProcess, SetParameterFieldProcess::Pointer, Process>(
        m, "SetParameterFieldProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<SetMultipleMovingLoadsProcess, SetMultipleMovingLoadsProcess::Pointer, Process>(
        m, "SetMultipleMovingLoadsProcess")
        .def(py::init<ModelPart&, Parameters>());

    py::class_<ApplyVectorConstraintTableProcess, ApplyVectorConstraintTableProcess::Pointer, Process>(
        m, "ApplyVectorConstraintTableProcess")
        .def(py::init<Model&, const Parameters&>());

    py::class_<ApplyScalarConstraintTableProcess, ApplyScalarConstraintTableProcess::Pointer, Process>(
        m, "ApplyScalarConstraintTableProcess")
        .def(py::init<ModelPart&, const Parameters&>());

    py::class_<ApplyNormalLoadTableProcess, ApplyNormalLoadTableProcess::Pointer, Process>(
        m, "ApplyNormalLoadTableProcess")
        .def(py::init<ModelPart&, const Parameters&>());

    py::class_<GeoExtrapolateIntegrationPointValuesToNodesProcess, GeoExtrapolateIntegrationPointValuesToNodesProcess::Pointer, Process>(
        m, "GeoExtrapolateIntegrationPointValuesToNodesProcess")
        .def(py::init<Model&, const Parameters&>());

    py::class_<ApplyCPhiReductionProcess, ApplyCPhiReductionProcess::Pointer, Process>(
        m, "ApplyCPhiReductionProcess")
        .def(py::init<Model&, const Parameters&>());

    py::class_<CalculateIncrementalMotionProcess, CalculateIncrementalMotionProcess::Pointer, Process>(
        m, "CalculateIncrementalMotionProcess")
        .def(py::init<ModelPart&, const Parameters&>());

    py::class_<CalculateTotalMotionProcess, CalculateTotalMotionProcess::Pointer, Process>(
        m, "CalculateTotalMotionProcess")
        .def(py::init<ModelPart&, const Parameters&>());

    py::class_<ApplyFinalStressesOfPreviousStageToInitialState, ApplyFinalStressesOfPreviousStageToInitialState::Pointer, Process>(
        m, "ApplyFinalStressesOfPreviousStageToInitialState")
        .def(py::init<Model&, const Parameters&>());

    py::class_<ApplyInitialUniformStressField, ApplyInitialUniformStressField::Pointer, Process>(
        m, "ApplyInitialUniformStressField")
        .def(py::init<ModelPart&, const Parameters&>());

    py::class_<FindNeighboursOfInterfacesProcess, FindNeighboursOfInterfacesProcess::Pointer, Process>(
        m, "FindNeighboursOfInterfacesProcess")
        .def(py::init<Model&, const Parameters&>());
}

} // Namespace Kratos::Python.
