//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/kratos_parameters.h"

#include "processes/process.h"
#include "processes/output_process.h"
#include "python/add_processes_to_python.h"
#include "processes/calculate_embedded_nodal_variable_from_skin_process.h"
#include "processes/edge_based_gradient_recovery_process.h"
#include "processes/fast_transfer_between_model_parts_process.h"
#include "processes/find_nodal_h_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_conditions_neighbours_process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/flux_corrected_transport_convection_process.h"
#include "processes/calculate_nodal_area_process.h"
#include "processes/entity_erase_process.h"
#include "processes/eliminate_isolated_nodes_process.h"
#include "processes/calculate_distance_to_path_process.h"
#include "processes/calculate_signed_distance_to_3d_skin_process.h"
#include "processes/calculate_embedded_signed_distance_to_3d_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/tetrahedral_mesh_orientation_check.h"
#include "processes/variational_distance_calculation_process.h"
#include "processes/levelset_convection_process.h"
#include "processes/flux_corrected_transport_convection_process.h"
#include "processes/apply_constant_scalarvalue_process.h"
#include "processes/apply_constant_vectorvalue_process.h"
#include "processes/check_skin_process.h"
#include "processes/replace_elements_and_condition_process.h"
#include "processes/compute_nodal_gradient_process.h"
#include "processes/compute_nodal_normal_divergence_process.h"
#include "processes/assign_scalar_variable_to_entities_process.h"
#include "processes/assign_scalar_input_to_entities_process.h"
#include "processes/assign_scalar_field_to_entities_process.h"
#include "processes/reorder_and_optimize_modelpart_process.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "processes/apply_ray_casting_process.h"
#include "processes/apply_ray_casting_interface_recognition_process.h"
#include "processes/simple_mortar_mapper_process.h"
#include "processes/simple_mortar_mapper_wrapper_process.h"
#include "processes/skin_detection_process.h"
#include "processes/sub_model_part_skin_detection_process.h"
#include "processes/apply_periodic_boundary_condition_process.h"
#include "processes/integration_values_extrapolation_to_nodes_process.h"
#include "processes/time_averaging_process.h"
#include "processes/from_json_check_result_process.h"
#include "processes/set_initial_state_process.h"
#include "processes/split_internal_interfaces_process.h"
#include "processes/parallel_distance_calculation_process.h"
#include "processes/generic_find_elements_neighbours_process.h"
#include "processes/check_same_modelpart_using_skin_distance_process.h"
#include "processes/calculate_nodal_distance_to_skin_process.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos::Python
{

// Discontinuous distance computation auxiliar functions
template<std::size_t TDim>
void CalculateDiscontinuousEmbeddedVariableFromSkinDouble(
    CalculateDiscontinuousDistanceToSkinProcess<TDim> &rDiscDistProcess,
    const Variable<double> &rVariable,
    const Variable<double> &rEmbeddedVariable)
{
    rDiscDistProcess.CalculateEmbeddedVariableFromSkin(rVariable, rEmbeddedVariable);
}

template<std::size_t TDim>
void CalculateDiscontinuousEmbeddedVariableFromSkinArray(
    CalculateDiscontinuousDistanceToSkinProcess<TDim> &rDiscDistProcess,
    const Variable<array_1d<double,3>> &rVariable,
    const Variable<array_1d<double,3>> &rEmbeddedVariable)
{
    rDiscDistProcess.CalculateEmbeddedVariableFromSkin(rVariable, rEmbeddedVariable);
}

// Continuous distance computation auxiliar functions
template<std::size_t TDim>
void CalculateEmbeddedVariableFromSkinDouble(
    CalculateDistanceToSkinProcess<TDim> &rDistProcess,
    const Variable<double> &rVariable,
    const Variable<double> &rEmbeddedVariable)
{
    rDistProcess.CalculateEmbeddedVariableFromSkin(rVariable, rEmbeddedVariable);
}

template<std::size_t TDim>
void CalculateEmbeddedVariableFromSkinArray(
    CalculateDistanceToSkinProcess<TDim> &rDistProcess,
    const Variable<array_1d<double,3>> &rVariable,
    const Variable<array_1d<double,3>> &rEmbeddedVariable)
{
    rDistProcess.CalculateEmbeddedVariableFromSkin(rVariable, rEmbeddedVariable);
}

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

template <std::size_t TDim, std::size_t TNumNodes, class TVariableType, const SizeType TNumNodesMaster = TNumNodes>
void DefineSimpleMortarMapperProcess(pybind11::module& m, const std::string& className)
{
    using ProcessType = SimpleMortarMapperProcess<TDim, TNumNodes, TVariableType, TNumNodesMaster>;

    pybind11::class_<ProcessType, typename ProcessType::Pointer, Process>(m, className.c_str())
    .def(pybind11::init<ModelPart&, ModelPart&, TVariableType&>())
    .def(pybind11::init<ModelPart&, ModelPart&, Parameters>())
    .def(pybind11::init<ModelPart&, ModelPart&, TVariableType&, Parameters>())
    .def(pybind11::init<ModelPart&, ModelPart&, TVariableType&, Parameters, LinearSolverType::Pointer>())
    .def(pybind11::init<ModelPart&, ModelPart&, TVariableType&, TVariableType&>())
    .def(pybind11::init<ModelPart&, ModelPart&, TVariableType&, TVariableType&, Parameters>())
    .def(pybind11::init<ModelPart&, ModelPart&, TVariableType&, TVariableType&, Parameters, LinearSolverType::Pointer>())
    .def("Map", &ProcessType::Map)
    ;
}

void  AddProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Process, Process::Pointer, Flags>(m,"Process")
    .def(py::init<>())
    .def("Create",&Process::Create)
    .def("Execute",&Process::Execute)
    .def("ExecuteInitialize",&Process::ExecuteInitialize)
    .def("ExecuteBeforeSolutionLoop",&Process::ExecuteBeforeSolutionLoop)
    .def("ExecuteInitializeSolutionStep",&Process::ExecuteInitializeSolutionStep)
    .def("ExecuteFinalizeSolutionStep",&Process::ExecuteFinalizeSolutionStep)
    .def("ExecuteBeforeOutputStep",&Process::ExecuteBeforeOutputStep)
    .def("ExecuteAfterOutputStep",&Process::ExecuteAfterOutputStep)
    .def("ExecuteFinalize",&Process::ExecuteFinalize)
    .def("Check",&Process::Check)
    .def("Clear",&Process::Clear)
    .def("GetDefaultParameters",&Process::GetDefaultParameters)
    .def("Info",&Process::Info)
    .def("__str__", PrintObject<Process>)
    ;

    py::class_<OutputProcess, OutputProcess::Pointer, Process>(m,"OutputProcess")
    .def(py::init<>())
    .def("IsOutputStep",&OutputProcess::IsOutputStep)
    .def("PrintOutput",&OutputProcess::PrintOutput)
    ;

    py::class_<FindGlobalNodalNeighboursProcess, FindGlobalNodalNeighboursProcess::Pointer, Process>
        (m,"FindGlobalNodalNeighboursProcess")
    .def(py::init([](const DataCommunicator& rDataComm, ModelPart& rModelPart) {
        KRATOS_WARNING("FindGlobalNodalNeighboursProcess") << "Using deprecated constructor. Please use constructor without DataCommunicator.";
        return Kratos::make_shared<FindGlobalNodalNeighboursProcess>(rModelPart);
    }))
    .def(py::init([](ModelPart& rModelPart) {
        return Kratos::make_shared<FindGlobalNodalNeighboursProcess>(rModelPart);
    }))
    .def("ClearNeighbours",&FindGlobalNodalNeighboursProcess::ClearNeighbours)
    .def("GetNeighbourIds",&FindGlobalNodalNeighboursProcess::GetNeighbourIds)
    ;

    typedef FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> FindGlobalNodalNeighboursForConditionsProcess;
    py::class_<FindGlobalNodalNeighboursForConditionsProcess, FindGlobalNodalNeighboursForConditionsProcess::Pointer, Process>
        (m,"FindGlobalNodalNeighboursForConditionsProcess")
    .def(py::init([](const DataCommunicator& rDataComm, ModelPart& rModelPart) {
        KRATOS_WARNING("FindGlobalNodalNeighboursForConditionsProcess") << "Using deprecated constructor. Please use constructor without DataCommunicator.";
        return Kratos::make_shared<FindGlobalNodalNeighboursForConditionsProcess>(rModelPart, NEIGHBOUR_CONDITION_NODES);
    }))
    .def(py::init([](ModelPart& rModelPart) {
        return Kratos::make_shared<FindGlobalNodalNeighboursForConditionsProcess>(rModelPart, NEIGHBOUR_CONDITION_NODES);
    }))
    .def("ClearNeighbours",&FindGlobalNodalNeighboursForConditionsProcess::ClearNeighbours)
    .def("GetNeighbourIds",&FindGlobalNodalNeighboursForConditionsProcess::GetNeighbourIds)
    ;

    using FindGlobalNodalElementalNeighboursProcessType = FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>;
    py::class_<FindGlobalNodalElementalNeighboursProcessType, typename FindGlobalNodalElementalNeighboursProcessType::Pointer, Process>(m,"FindGlobalNodalElementalNeighboursProcess")
        .def(py::init([](const DataCommunicator& rDataComm, ModelPart& rModelPart) {
            KRATOS_WARNING("FindGlobalNodalElementalNeighboursProcess") << "Using deprecated constructor. Please use constructor without DataCommunicator.";
            return Kratos::make_shared<FindGlobalNodalElementalNeighboursProcessType>(rModelPart);
        }))
        .def(py::init<ModelPart&>())
        .def(py::init<Model&, Parameters>())
        .def("ClearNeighbours", [](FindGlobalNodalElementalNeighboursProcessType& rSelf){
            KRATOS_WARNING("FindGlobalNodalElementalNeighboursProcess") << "Using deprecated ClearNeighbours method. please use Clear().";
            rSelf.Clear();})
        .def("GetNeighbourIds",&FindGlobalNodalElementalNeighboursProcessType::GetNeighbourIds)
        ;

    using FindGlobalNodalConditionalNeighboursProcessType = FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType>;
    py::class_<FindGlobalNodalConditionalNeighboursProcessType, typename FindGlobalNodalConditionalNeighboursProcessType::Pointer, Process>(m,"FindGlobalNodalConditionNeighboursProcess")
        .def(py::init<Model&, Parameters>())
        .def(py::init<ModelPart&>())
        .def("ClearNeighbours", [](FindGlobalNodalConditionalNeighboursProcessType& rSelf){
            KRATOS_WARNING("FindGlobalNodalConditionNeighboursProcess") << "Using deprecated ClearNeighbours method. please use Clear().";
            rSelf.Clear();})
        .def("GetNeighbourIds",&FindGlobalNodalConditionalNeighboursProcessType::GetNeighbourIds)
        ;

    py::class_<FindIntersectedGeometricalObjectsProcess, FindIntersectedGeometricalObjectsProcess::Pointer, Process>
        (m, "FindIntersectedGeometricalObjectsProcess")
    .def(py::init<Model&, Parameters>())
    ;

    // Find NODAL_H (Historical variables stored)
    py::class_<FindNodalHProcess<FindNodalHSettings::SaveAsHistoricalVariable>, FindNodalHProcess<FindNodalHSettings::SaveAsHistoricalVariable>::Pointer, Process>(m,"FindNodalHProcess")
    .def(py::init<ModelPart&>())
    ;

    // Find NODAL_H (Non-historical variables stored)
    py::class_<FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>, FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"FindNodalHNonHistoricalProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<FindNodalNeighboursProcess, FindNodalNeighboursProcess::Pointer, Process>(m,"FindNodalNeighboursProcess")
            .def(py::init<ModelPart& >())
    .def("ClearNeighbours",&FindNodalNeighboursProcess::ClearNeighbours)
    .def(py::init<ModelPart&, unsigned int, unsigned int>())
    ;

    py::class_<FindConditionsNeighboursProcess, FindConditionsNeighboursProcess::Pointer, Process>(m,"FindConditionsNeighboursProcess")
        .def(py::init<Model&, Parameters>())
        .def(py::init<ModelPart&, const int, const unsigned int>())
        .def("ClearNeighbours",&FindConditionsNeighboursProcess::ClearNeighbours)
    ;

    py::class_<CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsHistoricalVariable>, CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsHistoricalVariable>::Pointer, Process>(m,"CalculateNodalAreaProcess")
    .def(py::init<Model&, Parameters>())
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, std::size_t>())
    ;

    py::class_<CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>, CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"CalculateNonHistoricalNodalAreaProcess")
    .def(py::init<Model&, Parameters>())
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, std::size_t>())
    ;

    py::class_<EntitiesEraseProcess<Node>, EntitiesEraseProcess<Node>::Pointer, Process>(m,"NodeEraseProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<EntitiesEraseProcess<Element>, EntitiesEraseProcess<Element>::Pointer, Process>(m,"ElementEraseProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<EntitiesEraseProcess<Condition>, EntitiesEraseProcess<Condition>::Pointer, Process>(m,"ConditionEraseProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<EntitiesEraseProcess<MasterSlaveConstraint>, EntitiesEraseProcess<MasterSlaveConstraint>::Pointer, Process>(m,"MasterSlaveConstraintEraseProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<EliminateIsolatedNodesProcess, EliminateIsolatedNodesProcess::Pointer, Process>(m,"EliminateIsolatedNodesProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<CalculateDistanceToPathProcess<CalculateDistanceToPathSettings::SaveAsHistoricalVariable>, CalculateDistanceToPathProcess<CalculateDistanceToPathSettings::SaveAsHistoricalVariable>::Pointer, Process>(m, "CalculateDistanceToPathProcess")
    .def(py::init<Model&, Parameters>())
    ;

    py::class_<CalculateDistanceToPathProcess<CalculateDistanceToPathSettings::SaveAsNonHistoricalVariable>, CalculateDistanceToPathProcess<CalculateDistanceToPathSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m, "CalculateDistanceToPathNonHistoricalProcess")
    .def(py::init<Model&, Parameters>())
    ;

    py::class_<CalculateSignedDistanceTo3DSkinProcess, CalculateSignedDistanceTo3DSkinProcess::Pointer, Process>(m,"CalculateSignedDistanceTo3DSkinProcess")
            .def(py::init<ModelPart&, ModelPart&>())
    .def("GenerateSkinModelPart",&CalculateSignedDistanceTo3DSkinProcess::GenerateSkinModelPart)
    .def("MappingPressureToStructure",&CalculateSignedDistanceTo3DSkinProcess::MappingPressureToStructure)
    ;

    py::class_<CalculateEmbeddedSignedDistanceTo3DSkinProcess, CalculateEmbeddedSignedDistanceTo3DSkinProcess::Pointer, Process>(m,"CalculateEmbeddedSignedDistanceTo3DSkinProcess")
            .def(py::init< ModelPart&, ModelPart& >())
    .def(py::init< ModelPart&, ModelPart&, bool>())
    ;

   py::class_<CalculateSignedDistanceTo3DConditionSkinProcess, CalculateSignedDistanceTo3DConditionSkinProcess::Pointer, Process>(m,"CalculateSignedDistanceTo3DConditionSkinProcess")
            .def(py::init<ModelPart&, ModelPart&>())
    ;

    py::class_<StructuredMeshGeneratorProcess, StructuredMeshGeneratorProcess::Pointer, Process>(m,"StructuredMeshGeneratorProcess")
            .def(py::init<const Geometry< Node >&, ModelPart&, Parameters>()) //TODO: VERIFY IF THE NEXT IS NEEDED: [with_custodian_and_ward<1, 2>()])
    ;

    auto orientation_check_interface = py::class_<TetrahedralMeshOrientationCheck, TetrahedralMeshOrientationCheck::Pointer, Process>(m,"TetrahedralMeshOrientationCheck")
            .def(py::init<ModelPart&, bool>())
            .def(py::init<ModelPart&, bool, Kratos::Flags>())
    .def("SwapAll",&TetrahedralMeshOrientationCheck::SwapAll)
    .def("SwapNegativeElements",&TetrahedralMeshOrientationCheck::SwapNegativeElements)
    ;
    orientation_check_interface.attr("ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS") = &TetrahedralMeshOrientationCheck::ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS;
    orientation_check_interface.attr("COMPUTE_NODAL_NORMALS") = &TetrahedralMeshOrientationCheck::COMPUTE_NODAL_NORMALS;
    orientation_check_interface.attr("COMPUTE_CONDITION_NORMALS") = &TetrahedralMeshOrientationCheck::COMPUTE_CONDITION_NORMALS;

    py::class_<VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"VariationalDistanceCalculationProcess2D")
            .def(py::init<ModelPart&, LinearSolverType::Pointer>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags, std::string>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags, std::string, double>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags, std::string, double, double>())
            .def_readonly_static("CALCULATE_EXACT_DISTANCES_TO_PLANE", &VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::CALCULATE_EXACT_DISTANCES_TO_PLANE)
    ;
    py::class_<VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"VariationalDistanceCalculationProcess3D")
            .def(py::init<ModelPart&, LinearSolverType::Pointer>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags, std::string>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags, std::string, double>())
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, Flags, std::string, double, double>())
            .def_readonly_static("CALCULATE_EXACT_DISTANCES_TO_PLANE", &VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::CALCULATE_EXACT_DISTANCES_TO_PLANE)
    ;

    py::class_<ParallelDistanceCalculationProcess<2>, ParallelDistanceCalculationProcess<2>::Pointer, Process>(m,"ParallelDistanceCalculationProcess2D")
        .def(py::init<ModelPart&, Parameters>())
        .def(py::init<Model&, Parameters>())
        .def("FindMaximumEdgeSize", &ParallelDistanceCalculationProcess<2>::FindMaximumEdgeSize)
    ;

    py::class_<ParallelDistanceCalculationProcess<3>, ParallelDistanceCalculationProcess<3>::Pointer, Process>(m,"ParallelDistanceCalculationProcess3D")
        .def(py::init<ModelPart&, Parameters>())
        .def(py::init<Model&, Parameters>())
        .def("FindMaximumEdgeSize", &ParallelDistanceCalculationProcess<3>::FindMaximumEdgeSize)
    ;

    py::class_<LevelSetConvectionProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, LevelSetConvectionProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"LevelSetConvectionProcess2D")
        .def(py::init<Model&, LinearSolverType::Pointer, Parameters>())
        .def(py::init<ModelPart&, LinearSolverType::Pointer, Parameters>())
    ;

    py::class_<LevelSetConvectionProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, LevelSetConvectionProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"LevelSetConvectionProcess3D")
        .def(py::init<Model&, LinearSolverType::Pointer, Parameters>())
        .def(py::init<ModelPart&, LinearSolverType::Pointer, Parameters>())
    ;

    py::class_<FluxCorrectedTransportConvectionProcess<2>, FluxCorrectedTransportConvectionProcess<2>::Pointer, Process>(m,"FluxCorrectedTransportConvectionProcess2D")
        .def(py::init<Model&, Parameters>())
    ;

    py::class_<FluxCorrectedTransportConvectionProcess<3>, FluxCorrectedTransportConvectionProcess<3>::Pointer, Process>(m,"FluxCorrectedTransportConvectionProcess3D")
        .def(py::init<Model&, Parameters>())
    ;

    py::class_<ApplyConstantScalarValueProcess, ApplyConstantScalarValueProcess::Pointer, Process>(m,"ApplyConstantScalarValueProcess")
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, const Variable<double>&, double, std::size_t, Flags>())
            .def(py::init< ModelPart&, Parameters >())
            .def(py::init<ModelPart&, const Variable<int>&, int, std::size_t, Flags>())
            .def(py::init<ModelPart&, const Variable<bool>&, bool, std::size_t, Flags>())
            .def("ExecuteInitialize", &ApplyConstantScalarValueProcess::ExecuteInitialize)
            .def_readonly_static("VARIABLE_IS_FIXED", &ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)
    ;

    py::class_<ApplyConstantVectorValueProcess, ApplyConstantVectorValueProcess::Pointer, Process>(m,"ApplyConstantVectorValueProcess")
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, const Variable<array_1d<double, 3 > >& , const double, const Vector , std::size_t, Flags>())
            .def(py::init< ModelPart&, Parameters >())
            .def_readonly_static("X_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::X_COMPONENT_FIXED)
            .def_readonly_static("Y_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Y_COMPONENT_FIXED)
            .def_readonly_static("Z_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Z_COMPONENT_FIXED)
    ;

    py::class_<CheckSkinProcess, CheckSkinProcess::Pointer, Process>(m,"CheckSkinProcess")
            .def(py::init<ModelPart&, Flags>())
    ;

    py::class_<ReplaceElementsAndConditionsProcess, ReplaceElementsAndConditionsProcess::Pointer, Process>(m,"ReplaceElementsAndConditionsProcess")
            .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<SetInitialStateProcess<3>, SetInitialStateProcess<3>::Pointer, Process>(m,"SetInitialStateProcess3D")
            .def(py::init<ModelPart&>())
            .def(py::init<ModelPart&, const Vector&, const Vector&, const Matrix&>())
            .def(py::init<ModelPart&, const Vector&, const int>())
            .def(py::init<ModelPart&, const Matrix&>())
    ;
    py::class_<SetInitialStateProcess<2>, SetInitialStateProcess<2>::Pointer, Process>(m,"SetInitialStateProcess2D")
            .def(py::init<ModelPart&>())
            .def(py::init<ModelPart&, const Vector&, const Vector&, const Matrix&>())
            .def(py::init<ModelPart&, const Vector&, const int>())
            .def(py::init<ModelPart&, const Matrix&>())
    ;

    /* Historical */
    py::class_<ComputeNodalGradientProcess< ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>, ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::Pointer, Process>(m,"ComputeNodalGradientProcess")
    .def(py::init<Model&, Parameters>())
    .def(py::init<ModelPart&, Parameters>())
    .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>&, const bool >())
    ;

    m.attr("ComputeNodalGradientProcess2D") = m.attr("ComputeNodalGradientProcess");
    m.attr("ComputeNodalGradientProcess3D") = m.attr("ComputeNodalGradientProcess");
    m.attr("ComputeNodalGradientProcessComp2D") = m.attr("ComputeNodalGradientProcess");
    m.attr("ComputeNodalGradientProcessComp3D") = m.attr("ComputeNodalGradientProcess");

    /* Non-Historical */
    py::class_<ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>, ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcess")
    .def(py::init<Model&, Parameters>())
    .def(py::init<ModelPart&, Parameters>())
    .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>&, const bool >())
    ;

    m.attr("ComputeNonHistoricalNodalGradientProcess2D") = m.attr("ComputeNonHistoricalNodalGradientProcess");
    m.attr("ComputeNonHistoricalNodalGradientProcess3D") = m.attr("ComputeNonHistoricalNodalGradientProcess");
    m.attr("ComputeNonHistoricalNodalGradientProcessComp2D") = m.attr("ComputeNonHistoricalNodalGradientProcess");
    m.attr("ComputeNonHistoricalNodalGradientProcessComp3D") = m.attr("ComputeNonHistoricalNodalGradientProcess");

    /* Historical */
    py::class_<ComputeNodalNormalDivergenceProcess< ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>, ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::Pointer, Process>(m,"ComputeNodalNormalDivergenceProcess")
    .def(py::init<ModelPart&, Variable<array_1d<double,3> >&, Variable<double>& , Variable<double>& >())
    .def(py::init<ModelPart&, Variable<array_1d<double,3> >&, Variable<double>& , Variable<double>&, const bool >())
    .def(py::init<ModelPart&, Variable<array_1d<double,3> >&, Variable<double>& , Variable<double>&, const bool, const bool >())
    ;

    /* Non-Historical */
    py::class_<ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>, ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"ComputeNonHistoricalNodalNormalDivergenceProcess")
    .def(py::init<ModelPart&, Variable<array_1d<double,3> >&, Variable<double>& , Variable<double>& >())
    .def(py::init<ModelPart&, Variable<array_1d<double,3> >&, Variable<double>& , Variable<double>&, const bool >())
    .def(py::init<ModelPart&, Variable<array_1d<double,3> >&, Variable<double>& , Variable<double>&, const bool, const bool >())
    ;

    // Check the same model part using skin distance
    py::class_<CheckSameModelPartUsingSkinDistanceProcess<2>, CheckSameModelPartUsingSkinDistanceProcess<2>::Pointer, Process>(m,"CheckSameModelPartUsingSkinDistanceProcess2D")
        .def(py::init<Model&, Parameters>())
        ;
    py::class_<CheckSameModelPartUsingSkinDistanceProcess<3>, CheckSameModelPartUsingSkinDistanceProcess<3>::Pointer, Process>(m,"CheckSameModelPartUsingSkinDistanceProcess3D")
        .def(py::init<Model&, Parameters>())
        ;

    // Discontinuous distance computation methods
    py::class_<CalculateDiscontinuousDistanceToSkinProcess<2>, CalculateDiscontinuousDistanceToSkinProcess<2>::Pointer, Process>(m,"CalculateDiscontinuousDistanceToSkinProcess2D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def(py::init<ModelPart&, ModelPart&, const Flags>())
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateDiscontinuousEmbeddedVariableFromSkinArray<2>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateDiscontinuousEmbeddedVariableFromSkinDouble<2>)
        .def_readonly_static("CALCULATE_ELEMENTAL_EDGE_DISTANCES", &CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES)
        .def_readonly_static("CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED", &CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED)
        .def_readonly_static("USE_POSITIVE_EPSILON_FOR_ZERO_VALUES", &CalculateDiscontinuousDistanceToSkinProcessFlags::USE_POSITIVE_EPSILON_FOR_ZERO_VALUES)
        ;

    py::class_<CalculateDiscontinuousDistanceToSkinProcess<3>, CalculateDiscontinuousDistanceToSkinProcess<3>::Pointer, Process>(m,"CalculateDiscontinuousDistanceToSkinProcess3D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def(py::init<ModelPart&, ModelPart&, const Flags>())
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateDiscontinuousEmbeddedVariableFromSkinArray<3>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateDiscontinuousEmbeddedVariableFromSkinDouble<3>)
        .def_readonly_static("CALCULATE_ELEMENTAL_EDGE_DISTANCES", &CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES)
        .def_readonly_static("CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED", &CalculateDiscontinuousDistanceToSkinProcessFlags::CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED)
        .def_readonly_static("USE_POSITIVE_EPSILON_FOR_ZERO_VALUES", &CalculateDiscontinuousDistanceToSkinProcessFlags::USE_POSITIVE_EPSILON_FOR_ZERO_VALUES)
        ;

    // Continuous distance computation methods
    py::class_<CalculateDistanceToSkinProcess<2>, CalculateDistanceToSkinProcess<2>::Pointer, Process>(m,"CalculateDistanceToSkinProcess2D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def(py::init<ModelPart&, ModelPart&, double>())
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinArray<2>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinDouble<2>)
    ;

    py::class_<CalculateDistanceToSkinProcess<3>, CalculateDistanceToSkinProcess<3>::Pointer, Process>(m,"CalculateDistanceToSkinProcess3D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def(py::init<ModelPart&, ModelPart&, double>())
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinArray<2>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinDouble<2>)
    ;

    // Continuous distance computation methods
    py::class_<ApplyRayCastingProcess<2>, ApplyRayCastingProcess<2>::Pointer, Process>(m,"ApplyRayCastingProcess2D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def(py::init<ModelPart&, ModelPart&, double>())
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
    ;

    py::class_<ApplyRayCastingProcess<3>, ApplyRayCastingProcess<3>::Pointer, Process>(m,"ApplyRayCastingProcess3D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def(py::init<ModelPart&, ModelPart&, double>())
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
    ;

    py::class_<ApplyRayCastingInterfaceRecognitionProcess<2>, 
      ApplyRayCastingInterfaceRecognitionProcess<2>::Pointer, Process>
      (m,"ApplyRayCastingInterfaceRecognitionProcess2D")
        .def(py::init<Model&, Parameters>())
    ;

    py::class_<ApplyRayCastingInterfaceRecognitionProcess<3>,
       ApplyRayCastingInterfaceRecognitionProcess<3>::Pointer, Process>
       (m,"ApplyRayCastingInterfaceRecognitionProcess3D")
        .def(py::init<Model&, Parameters>())
    ;

//     // Calculate embedded variable from skin processes
    py::class_<CalculateEmbeddedNodalVariableFromSkinProcess<double, SparseSpaceType, LocalSpaceType, LinearSolverType>, CalculateEmbeddedNodalVariableFromSkinProcess<double, SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer, Process>(
        m, "CalculateEmbeddedNodalVariableFromSkinProcessDouble")
        .def(py::init<Model &, Parameters>())
        ;

    py::class_<CalculateEmbeddedNodalVariableFromSkinProcess<array_1d<double, 3>, SparseSpaceType, LocalSpaceType, LinearSolverType>, CalculateEmbeddedNodalVariableFromSkinProcess<array_1d<double, 3>, SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer, Process>(
        m, "CalculateEmbeddedNodalVariableFromSkinProcessArray")
        .def(py::init<Model &, Parameters>())
        ;

    py::class_<ReorderAndOptimizeModelPartProcess, ReorderAndOptimizeModelPartProcess::Pointer, Process>(m,"ReorderAndOptimizeModelPartProcess")
            .def(py::init<ModelPart&, Parameters>())
            ;

    py::class_<AssignScalarVariableToEntitiesProcess<Node>, AssignScalarVariableToEntitiesProcess<Node>::Pointer, Process>(m,"AssignScalarVariableToNodesProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarVariableToEntitiesProcess<Condition>, AssignScalarVariableToEntitiesProcess<Condition>::Pointer, Process>(m,"AssignScalarVariableToConditionsProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarVariableToEntitiesProcess<Element>, AssignScalarVariableToEntitiesProcess<Element>::Pointer, Process>(m,"AssignScalarVariableToElementsProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>, AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>::Pointer, Process>(m,"AssignScalarVariableToMasterSlaveConstraintsProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarInputToEntitiesProcess<Node, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>, AssignScalarInputToEntitiesProcess<Node, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"AssignScalarInputToNodesProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarInputToEntitiesProcess<Node, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>, AssignScalarInputToEntitiesProcess<Node, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>::Pointer, Process>(m,"AssignScalarInputHistoricalToNodesProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>, AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"AssignScalarInputToConditionsProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>, AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"AssignScalarInputToElementsProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarFieldToEntitiesProcess<Node>, AssignScalarFieldToEntitiesProcess<Node>::Pointer, Process>(m,"AssignScalarFieldToNodesProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarFieldToEntitiesProcess<Condition>, AssignScalarFieldToEntitiesProcess<Condition>::Pointer, Process>(m,"AssignScalarFieldToConditionsProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarFieldToEntitiesProcess<Element>, AssignScalarFieldToEntitiesProcess<Element>::Pointer, Process>(m,"AssignScalarFieldToElementsProcess")
    .def(py::init<ModelPart&, Parameters >())
    ;

    /* Simple Mortar mapper */
    // Wrapper
    py::class_<SimpleMortarMapperProcessWrapper, SimpleMortarMapperProcessWrapper::Pointer, Process>(m, "SimpleMortarMapperProcess")
    .def(py::init<ModelPart&, ModelPart&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Parameters, LinearSolverType::Pointer>())
    ;

    // 2D
    DefineSimpleMortarMapperProcess<2, 2, Variable<double>>(m, "SimpleMortarMapperProcess2D2NDouble");
    DefineSimpleMortarMapperProcess<2, 2, Variable<array_1d<double, 3>>>(m, "SimpleMortarMapperProcess2D2NVector");

    // 3D - Triangle
    DefineSimpleMortarMapperProcess<3, 3, Variable<double>>(m, "SimpleMortarMapperProcess3D3NDouble");
    DefineSimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>>(m, "SimpleMortarMapperProcess3D3NVector");

    // 3D - Quadrilateral
    DefineSimpleMortarMapperProcess<3, 4, Variable<double>>(m, "SimpleMortarMapperProcess3D4NDouble");
    DefineSimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>>(m, "SimpleMortarMapperProcess3D4NVector");

    // 3D - Triangle - Quadrilateral
    DefineSimpleMortarMapperProcess<3, 3, Variable<double>, 4>(m, "SimpleMortarMapperProcess3D3N4NDouble");
    DefineSimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>, 4>(m, "SimpleMortarMapperProcess3D3N4NVector");

    // 3D - Quadrilateral - Triangle
    DefineSimpleMortarMapperProcess<3, 4, Variable<double>, 3>(m, "SimpleMortarMapperProcess3D4N3NDouble");
    DefineSimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>, 3>(m, "SimpleMortarMapperProcess3D4N3NVector");

    // Transfer between model parts
    py::class_<FastTransferBetweenModelPartsProcess, FastTransferBetweenModelPartsProcess::Pointer, Process> FastTransferBetweenModelPartsProcess_Scope(m, "FastTransferBetweenModelPartsProcess");

    FastTransferBetweenModelPartsProcess_Scope.def(py::init<ModelPart&, ModelPart&>());
    FastTransferBetweenModelPartsProcess_Scope.def(py::init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered>());
    FastTransferBetweenModelPartsProcess_Scope.def(py::init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered, const Flags >());
    FastTransferBetweenModelPartsProcess_Scope.def(py::init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered, const Flags, const bool >());

    // Adding FastTransferBetweenModelPartsProcess related enums
    py::enum_<FastTransferBetweenModelPartsProcess::EntityTransfered>(FastTransferBetweenModelPartsProcess_Scope, "EntityTransfered")
    .value("NODES", FastTransferBetweenModelPartsProcess::EntityTransfered::NODES)
    .value("ELEMENTS", FastTransferBetweenModelPartsProcess::EntityTransfered::ELEMENTS)
    .value("NODESANDELEMENTS", FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS)
    .value("CONDITIONS", FastTransferBetweenModelPartsProcess::EntityTransfered::CONDITIONS)
    .value("NODESANDCONDITIONS", FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDCONDITIONS)
    .value("CONSTRAINTS", FastTransferBetweenModelPartsProcess::EntityTransfered::CONSTRAINTS)
    .value("NODESANDCONSTRAINTS", FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDCONSTRAINTS)
    .value("ALL", FastTransferBetweenModelPartsProcess::EntityTransfered::ALL)
    ;

    py::class_<SkinDetectionProcess<2>, SkinDetectionProcess<2>::Pointer, Process>(m, "SkinDetectionProcess2D")
    .def(py::init<ModelPart&>())
    .def(py::init< ModelPart&, Parameters >())
    ;

    py::class_<SkinDetectionProcess<3>, SkinDetectionProcess<3>::Pointer, Process>(m, "SkinDetectionProcess3D")
    .def(py::init<ModelPart&>())
    .def(py::init< ModelPart&, Parameters >())
    ;

    py::class_<SubModelPartSkinDetectionProcess<2>, SubModelPartSkinDetectionProcess<2>::Pointer, SkinDetectionProcess<2>>(m, "SubModelPartSkinDetectionProcess2D")
    .def(py::init< ModelPart&, Parameters >())
    ;

    py::class_<SubModelPartSkinDetectionProcess<3>, SubModelPartSkinDetectionProcess<3>::Pointer, SkinDetectionProcess<3>>(m, "SubModelPartSkinDetectionProcess3D")
    .def(py::init< ModelPart&, Parameters >())
    ;

    py::class_<ApplyPeriodicConditionProcess, ApplyPeriodicConditionProcess::Pointer, Process>(m,"ApplyPeriodicConditionProcess")
            .def(py::init<ModelPart&,ModelPart&, Parameters>())
    ;

    // The process to recover internal variables
    py::class_<IntegrationValuesExtrapolationToNodesProcess, IntegrationValuesExtrapolationToNodesProcess::Pointer, Process>(m, "IntegrationValuesExtrapolationToNodesProcess")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<TimeAveragingProcess, TimeAveragingProcess::Pointer, Process>(m, "TimeAveragingProcess")
    .def(py::init<Model&, Parameters>())
    ;

    py::class_<SplitInternalInterfacesProcess, SplitInternalInterfacesProcess::Pointer, Process>(m, "SplitInternalInterfacesProcess")
    .def(py::init<Model&, Parameters>())
    ;

    auto from_json_check_result_process_interface =
    py::class_<FromJSONCheckResultProcess, FromJSONCheckResultProcess::Pointer, Process>(m, "FromJSONCheckResultProcess")
    .def(py::init<Model&>())
    .def(py::init<Model&, Parameters>())
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("IsCorrectResult", &FromJSONCheckResultProcess::IsCorrectResult)
    .def("GetErrorMessage", &FromJSONCheckResultProcess::GetErrorMessage)
    ;

    from_json_check_result_process_interface.attr("CORRECT_RESULT")                 = &FromJSONCheckResultProcess::CORRECT_RESULT;
    from_json_check_result_process_interface.attr("HISTORICAL_VALUE")               = &FromJSONCheckResultProcess::HISTORICAL_VALUE;
    from_json_check_result_process_interface.attr("CHECK_ONLY_LOCAL_ENTITIES")      = &FromJSONCheckResultProcess::CHECK_ONLY_LOCAL_ENTITIES;

    using ScalarEdgeBasedGradientRecoveryProcessType = EdgeBasedGradientRecoveryProcess<double,SparseSpaceType,LocalSpaceType,LinearSolverType>;
    py::class_<ScalarEdgeBasedGradientRecoveryProcessType, ScalarEdgeBasedGradientRecoveryProcessType::Pointer, Process>(m, "EdgeBasedGradientRecoveryProcessScalar")
        .def(py::init<Model&, LinearSolverType::Pointer, Parameters>())
    ;

    using ArrayEdgeBasedGradientRecoveryProcessType = EdgeBasedGradientRecoveryProcess<array_1d<double,3>,SparseSpaceType,LocalSpaceType,LinearSolverType>;
    py::class_<ArrayEdgeBasedGradientRecoveryProcessType, ArrayEdgeBasedGradientRecoveryProcessType::Pointer, Process>(m, "EdgeBasedGradientRecoveryProcessArray")
        .def(py::init<Model&, LinearSolverType::Pointer, Parameters>())
    ;

    py::class_<GenericFindElementalNeighboursProcess, GenericFindElementalNeighboursProcess::Pointer, Process> (m, "GenericFindElementalNeighboursProcess")
    .def(py::init<ModelPart&>())
    .def("HasNeighboursInFaces", &GenericFindElementalNeighboursProcess::HasNeighboursInFaces)
    ;

    py::class_<CalculateNodalDistanceToSkinProcess, CalculateNodalDistanceToSkinProcess::Pointer, Process> (m, "CalculateNodalDistanceToSkinProcess")
    .def(py::init<Model&, Parameters>())
    ;
}

}  // namespace Kratos::Python.