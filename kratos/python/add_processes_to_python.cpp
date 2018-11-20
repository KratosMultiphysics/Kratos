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
#include "python/add_processes_to_python.h"
#include "processes/fast_transfer_between_model_parts_process.h"
#include "processes/find_nodal_h_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_conditions_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/calculate_nodal_area_process.h"
#include "processes/node_erase_process.h"
#include "processes/element_erase_process.h"
#include "processes/condition_erase_process.h"
#include "processes/eliminate_isolated_nodes_process.h"
#include "processes/calculate_signed_distance_to_3d_skin_process.h"
#include "processes/calculate_embedded_signed_distance_to_3d_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "processes/translation_operation.h"
#include "processes/rotation_operation.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/tetrahedral_mesh_orientation_check.h"
#include "processes/compute_bdfcoefficients_process.h"
#include "processes/variational_distance_calculation_process.h"
#include "processes/levelset_convection_process.h"
#include "processes/apply_constant_scalarvalue_process.h"
#include "processes/apply_constant_vectorvalue_process.h"
#include "processes/check_skin_process.h"
#include "processes/replace_elements_and_condition_process.h"
#include "processes/compute_nodal_gradient_process.h"
#include "processes/assign_scalar_variable_to_conditions_process.h"
#include "processes/assign_scalar_field_to_conditions_process.h"
#include "processes/reorder_and_optimize_modelpart_process.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "processes/simple_mortar_mapper_process.h"
#include "processes/skin_detection_process.h"
#include "includes/node.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "utilities/python_function_callback_utility.h"


namespace Kratos
{

namespace Python
{
typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;

// Discontinuous distance computation auxiliar functions
template<std::size_t TDim>
void CalculateEmbeddedVariableFromSkinDouble(
    CalculateDiscontinuousDistanceToSkinProcess<TDim> &rDiscDistProcess,
    const Variable<double> &rVariable,
    const Variable<double> &rEmbeddedVariable)
{
    rDiscDistProcess.CalculateEmbeddedVariableFromSkin(rVariable, rEmbeddedVariable);
}

template<std::size_t TDim>
void CalculateEmbeddedVariableFromSkinArray(
    CalculateDiscontinuousDistanceToSkinProcess<TDim> &rDiscDistProcess,
    const Variable<array_1d<double,3>> &rVariable,
    const Variable<array_1d<double,3>> &rEmbeddedVariable)
{
    rDiscDistProcess.CalculateEmbeddedVariableFromSkin(rVariable, rEmbeddedVariable);
}

void  AddProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Process, Process::Pointer>(m,"Process")
    .def(py::init<>())
    .def("Execute",&Process::Execute)
    .def("ExecuteInitialize",&Process::ExecuteInitialize)
    .def("ExecuteBeforeSolutionLoop",&Process::ExecuteBeforeSolutionLoop)
    .def("ExecuteInitializeSolutionStep",&Process::ExecuteInitializeSolutionStep)
    .def("ExecuteFinalizeSolutionStep",&Process::ExecuteFinalizeSolutionStep)
    .def("ExecuteBeforeOutputStep",&Process::ExecuteBeforeOutputStep)
    .def("ExecuteAfterOutputStep",&Process::ExecuteAfterOutputStep)
    .def("ExecuteFinalize",&Process::ExecuteFinalize)
    .def("__str__", PrintObject<Process>)
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
            .def(py::init<ModelPart&, unsigned int, unsigned int>())
    .def("ClearNeighbours",&FindNodalNeighboursProcess::ClearNeighbours)
    ;

    py::class_<FindConditionsNeighboursProcess, FindConditionsNeighboursProcess::Pointer, Process>(m,"FindConditionsNeighboursProcess")
            .def(py::init<ModelPart&, int, unsigned int>())
    .def("ClearNeighbours",&FindConditionsNeighboursProcess::ClearNeighbours)
    ;

    py::class_<FindElementalNeighboursProcess, FindElementalNeighboursProcess::Pointer, Process>(m,"FindElementalNeighboursProcess")
            .def(py::init<ModelPart&, int, unsigned int>())
    .def("ClearNeighbours",&FindElementalNeighboursProcess::ClearNeighbours)
    ;

    py::class_<CalculateNodalAreaProcess, CalculateNodalAreaProcess::Pointer, Process>(m,"CalculateNodalAreaProcess")
            .def(py::init<ModelPart&, unsigned int>())
    ;

    py::class_<NodeEraseProcess, NodeEraseProcess::Pointer, Process>(m,"NodeEraseProcess")
            .def(py::init<ModelPart&>())
    ;

    py::class_<ElementEraseProcess, ElementEraseProcess::Pointer, Process>(m,"ElementEraseProcess")
            .def(py::init<ModelPart&>())
    ;

    py::class_<ConditionEraseProcess, ConditionEraseProcess::Pointer, Process>(m,"ConditionEraseProcess")
            .def(py::init<ModelPart&>())
    ;

    py::class_<EliminateIsolatedNodesProcess, EliminateIsolatedNodesProcess::Pointer, Process>(m,"EliminateIsolatedNodesProcess")
            .def(py::init<ModelPart&>())
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

    py::class_<TranslationOperation, TranslationOperation::Pointer, Process>(m,"TranslationOperation")
            .def(py::init<ModelPart&, DenseVector<int> ,DenseVector<int> ,unsigned int>())
    ;

    py::class_<RotationOperation, RotationOperation::Pointer, Process>(m,"RotationOperation")
            .def(py::init<ModelPart&, DenseVector<int> ,DenseVector<int> ,unsigned int>())
    ;

    py::class_<StructuredMeshGeneratorProcess, StructuredMeshGeneratorProcess::Pointer, Process>(m,"StructuredMeshGeneratorProcess")
            .def(py::init<const Geometry< Node<3> >&, ModelPart&, Parameters&>()) //TODO: VERIFY IF THE NEXT IS NEEDED: [with_custodian_and_ward<1, 2>()])
    ;

    auto orientation_check_interface = py::class_<TetrahedralMeshOrientationCheck, TetrahedralMeshOrientationCheck::Pointer, Process>(m,"TetrahedralMeshOrientationCheck")
            .def(py::init<ModelPart&, bool>())
            .def(py::init<ModelPart&, bool, Kratos::Flags>())
    .def("SwapAll",&TetrahedralMeshOrientationCheck::SwapAll)
    .def("SwapNegativeElements",&TetrahedralMeshOrientationCheck::SwapNegativeElements)
    ;
    orientation_check_interface.attr("ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS") = &TetrahedralMeshOrientationCheck::ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS;
    orientation_check_interface.attr("NOT_ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS") = &TetrahedralMeshOrientationCheck::ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS;
    orientation_check_interface.attr("COMPUTE_NODAL_NORMALS") = &TetrahedralMeshOrientationCheck::COMPUTE_NODAL_NORMALS;
    orientation_check_interface.attr("NOT_COMPUTE_NODAL_NORMALS") = &TetrahedralMeshOrientationCheck::NOT_COMPUTE_NODAL_NORMALS;
    orientation_check_interface.attr("COMPUTE_CONDITION_NORMALS") = &TetrahedralMeshOrientationCheck::COMPUTE_CONDITION_NORMALS;
    orientation_check_interface.attr("NOT_COMPUTE_CONDITION_NORMALS") = &TetrahedralMeshOrientationCheck::NOT_COMPUTE_CONDITION_NORMALS;

    py::class_<ComputeBDFCoefficientsProcess, ComputeBDFCoefficientsProcess::Pointer, Process>(m,"ComputeBDFCoefficientsProcess")
            .def(py::init<ModelPart&, const unsigned int>())
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    py::class_<VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"VariationalDistanceCalculationProcess2D")
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;
    py::class_<VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"VariationalDistanceCalculationProcess3D")
            .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;

    py::class_<LevelSetConvectionProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, LevelSetConvectionProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"LevelSetConvectionProcess2D")
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer>())
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double>())
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double>())
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double, const unsigned int>())
    ;
    py::class_<LevelSetConvectionProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, LevelSetConvectionProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"LevelSetConvectionProcess3D")
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer>())
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double>())
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double>())
        .def(py::init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double, const unsigned int>())
    ;

    py::class_<ApplyConstantScalarValueProcess, ApplyConstantScalarValueProcess::Pointer, Process>(m,"ApplyConstantScalarValueProcess")
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, const Variable<double>&, double, std::size_t, Flags>())
            .def(py::init< ModelPart&, Parameters& >())
            .def(py::init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, double, std::size_t, Flags>())
            .def(py::init<ModelPart&, const Variable<int>&, int, std::size_t, Flags>())
            .def(py::init<ModelPart&, const Variable<bool>&, bool, std::size_t, Flags>())
            .def("ExecuteInitialize", &ApplyConstantScalarValueProcess::ExecuteInitialize)
            .def_readonly_static("VARIABLE_IS_FIXED", &ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)
    ;

    py::class_<ApplyConstantVectorValueProcess, ApplyConstantVectorValueProcess::Pointer, Process>(m,"ApplyConstantVectorValueProcess")
            .def(py::init<ModelPart&, Parameters>())
            .def(py::init<ModelPart&, const Variable<array_1d<double, 3 > >& , const double, const Vector , std::size_t, Flags>())
            .def(py::init< ModelPart&, Parameters& >())
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

    /* Historical */
    // DOUBLE
    py::class_<ComputeNodalGradientProcess<2, Variable<double>, Historical>, ComputeNodalGradientProcess<2, Variable<double>, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcess2D")
            .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    py::class_<ComputeNodalGradientProcess<3, Variable<double>, Historical>, ComputeNodalGradientProcess<3, Variable<double>, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcess3D")
            .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    // COMPONENT
    py::class_<ComputeNodalGradientProcess<2, component_type, Historical>, ComputeNodalGradientProcess<2, component_type, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcessComp2D")
            .def(py::init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    py::class_<ComputeNodalGradientProcess<3, component_type, Historical>, ComputeNodalGradientProcess<3, component_type, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcessComp3D")
            .def(py::init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    /* Non-Historical */
    // DOUBLE
    py::class_<ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>, ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcess2D")
            .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    py::class_<ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>, ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcess3D")
            .def(py::init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    // COMPONENT
    py::class_<ComputeNodalGradientProcess<2, component_type, NonHistorical>, ComputeNodalGradientProcess<2, component_type, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcessComp2D")
            .def(py::init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    py::class_<ComputeNodalGradientProcess<3, component_type, NonHistorical>, ComputeNodalGradientProcess<3, component_type, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcessComp3D")
            .def(py::init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    // Discontinuous distance computation methods
    py::class_<CalculateDiscontinuousDistanceToSkinProcess<2>, CalculateDiscontinuousDistanceToSkinProcess<2>::Pointer, Process>(m,"CalculateDiscontinuousDistanceToSkinProcess2D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinArray<2>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinDouble<2>)
        ;

    py::class_<CalculateDiscontinuousDistanceToSkinProcess<3>, CalculateDiscontinuousDistanceToSkinProcess<3>::Pointer, Process>(m,"CalculateDiscontinuousDistanceToSkinProcess3D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinArray<3>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinDouble<3>)
        ;

    // Continuous distance computation methods
    py::class_<CalculateDistanceToSkinProcess<2>, CalculateDistanceToSkinProcess<2>::Pointer, Process>(m,"CalculateDistanceToSkinProcess2D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinArray<2>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinDouble<2>)
    ;

    py::class_<CalculateDistanceToSkinProcess<3>, CalculateDistanceToSkinProcess<3>::Pointer, Process>(m,"CalculateDistanceToSkinProcess3D")
        .def(py::init<ModelPart&, ModelPart&>())
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinArray<2>)
        .def("CalculateEmbeddedVariableFromSkin", CalculateEmbeddedVariableFromSkinDouble<2>)
    ;

    py::class_<ReorderAndOptimizeModelPartProcess, ReorderAndOptimizeModelPartProcess::Pointer, Process>(m,"ReorderAndOptimizeModelPartProcess")
            .def(py::init<ModelPart&, Parameters>())
            ;


    py::class_<AssignScalarVariableToConditionsProcess, AssignScalarVariableToConditionsProcess::Pointer, Process>(m,"AssignScalarVariableToConditionsProcess")
            .def(py::init<ModelPart&, Parameters >())
    ;

    py::class_<AssignScalarFieldToConditionsProcess, AssignScalarFieldToConditionsProcess::Pointer, Process>(m,"AssignScalarFieldToConditionsProcess")
            .def(py::init<ModelPart&, Parameters >())
    ;


    //typedef PointerVectorSet<Node<3>, IndexedObject> NodesContainerType;
    //typedef PointerVectorSet<Dof<double>, IndexedObject> DofsContainerType;

    //py::class_<AddDofsNodalProcess<Variable<double> >, AddDofsNodalProcess<Variable<double> >::Pointer, Process>(m,"AddDoubleDofsNodalProcess")
    // .def(py::init<Variable<double>, NodesContainerType&, DofsContainerType&>())
    // ;
    //py::class_<AddDofsNodalProcess<VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >, AddDofsNodalProcess<VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >::Pointer, Process>(m,"AddArrayComponentDofsNodalProcess")
    // ;

    /* Simple Mortar mapper */
    // 2D
    py::class_<SimpleMortarMapperProcess<2, 2, Variable<double>>, SimpleMortarMapperProcess<2, 2, Variable<double>>::Pointer, Process>(m, "SimpleMortarMapperProcess2D2NDouble")
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    py::class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >>, SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >>::Pointer, Process>(m, "SimpleMortarMapperProcess2D2NVector")
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Triangle
    py::class_<SimpleMortarMapperProcess<3, 3, Variable<double>>, SimpleMortarMapperProcess<3, 3, Variable<double>>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3NDouble")
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    py::class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >>, SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3NVector")
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Quadrilateral
    py::class_<SimpleMortarMapperProcess<3, 4, Variable<double>>, SimpleMortarMapperProcess<3, 4, Variable<double>>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4NDouble")
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    py::class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >>, SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4NVector")
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Triangle - Quadrilateral
    py::class_<SimpleMortarMapperProcess<3, 3, Variable<double>, 4>, SimpleMortarMapperProcess<3, 3, Variable<double>, 4>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3N4NDouble")
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    py::class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, 4>, SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, 4>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3N4NVector")
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Quadrilateral - Triangle
    py::class_<SimpleMortarMapperProcess<3, 4, Variable<double>, 3>, SimpleMortarMapperProcess<3, 4, Variable<double>, 3>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4N3NDouble")
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    py::class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, 3>, SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, 3>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4N3NVector")
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(py::init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

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
}

}  // namespace Python.

} // Namespace Kratos
