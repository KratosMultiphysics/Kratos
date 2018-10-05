//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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

void  AddProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<Process, Process::Pointer>(m,"Process")
    .def(init<>())
    .def("Execute",&Process::Execute)
    .def("ExecuteInitialize",&Process::ExecuteInitialize)
    .def("ExecuteBeforeSolutionLoop",&Process::ExecuteBeforeSolutionLoop)
    .def("ExecuteInitializeSolutionStep",&Process::ExecuteInitializeSolutionStep)
    .def("ExecuteFinalizeSolutionStep",&Process::ExecuteFinalizeSolutionStep)
    .def("ExecuteBeforeOutputStep",&Process::ExecuteBeforeOutputStep)
    .def("ExecuteAfterOutputStep",&Process::ExecuteAfterOutputStep)
    .def("ExecuteFinalize",&Process::ExecuteFinalize)
    .def("Check",&Process::Check)
    .def("__repr__", &Process::Info)
    ;

    // Find NODAL_H (Historical variables stored)
    class_<FindNodalHProcess<true>, FindNodalHProcess<true>::Pointer, Process>(m,"FindNodalHProcess")
    .def(init<ModelPart&>())
    ;

    // Find NODAL_H (Non-historical variables stored)
    class_<FindNodalHProcess<false>, FindNodalHProcess<false>::Pointer, Process>(m,"FindNodalHNonHistoricalProcess")
    .def(init<ModelPart&>())
    ;

    class_<FindNodalNeighboursProcess, FindNodalNeighboursProcess::Pointer, Process>(m,"FindNodalNeighboursProcess")
            .def(init<ModelPart&, unsigned int, unsigned int>())
    .def("ClearNeighbours",&FindNodalNeighboursProcess::ClearNeighbours)
    ;

    class_<FindConditionsNeighboursProcess, FindConditionsNeighboursProcess::Pointer, Process>(m,"FindConditionsNeighboursProcess")
            .def(init<ModelPart&, int, unsigned int>())
    .def("ClearNeighbours",&FindConditionsNeighboursProcess::ClearNeighbours)
    ;

    class_<FindElementalNeighboursProcess, FindElementalNeighboursProcess::Pointer, Process>(m,"FindElementalNeighboursProcess")
            .def(init<ModelPart&, int, unsigned int>())
    .def("ClearNeighbours",&FindElementalNeighboursProcess::ClearNeighbours)
    ;

    class_<CalculateNodalAreaProcess, CalculateNodalAreaProcess::Pointer, Process>(m,"CalculateNodalAreaProcess")
            .def(init<ModelPart&, unsigned int>())
    ;

    class_<NodeEraseProcess, NodeEraseProcess::Pointer, Process>(m,"NodeEraseProcess")
            .def(init<ModelPart&>())
    ;

    class_<ElementEraseProcess, ElementEraseProcess::Pointer, Process>(m,"ElementEraseProcess")
            .def(init<ModelPart&>())
    ;

    class_<ConditionEraseProcess, ConditionEraseProcess::Pointer, Process>(m,"ConditionEraseProcess")
            .def(init<ModelPart&>())
    ;

    class_<EliminateIsolatedNodesProcess, EliminateIsolatedNodesProcess::Pointer, Process>(m,"EliminateIsolatedNodesProcess")
            .def(init<ModelPart&>())
    ;

    class_<CalculateSignedDistanceTo3DSkinProcess, CalculateSignedDistanceTo3DSkinProcess::Pointer, Process>(m,"CalculateSignedDistanceTo3DSkinProcess")
            .def(init<ModelPart&, ModelPart&>())
    .def("GenerateSkinModelPart",&CalculateSignedDistanceTo3DSkinProcess::GenerateSkinModelPart)
    .def("MappingPressureToStructure",&CalculateSignedDistanceTo3DSkinProcess::MappingPressureToStructure)
    ;

    class_<CalculateEmbeddedSignedDistanceTo3DSkinProcess, CalculateEmbeddedSignedDistanceTo3DSkinProcess::Pointer, Process>(m,"CalculateEmbeddedSignedDistanceTo3DSkinProcess")
            .def(init< ModelPart&, ModelPart& >())
    .def(init< ModelPart&, ModelPart&, bool>())
    ;

   class_<CalculateSignedDistanceTo3DConditionSkinProcess, CalculateSignedDistanceTo3DConditionSkinProcess::Pointer, Process>(m,"CalculateSignedDistanceTo3DConditionSkinProcess")
            .def(init<ModelPart&, ModelPart&>())
    ;

    class_<TranslationOperation, TranslationOperation::Pointer, Process>(m,"TranslationOperation")
            .def(init<ModelPart&, DenseVector<int> ,DenseVector<int> ,unsigned int>())
    ;

    class_<RotationOperation, RotationOperation::Pointer, Process>(m,"RotationOperation")
            .def(init<ModelPart&, DenseVector<int> ,DenseVector<int> ,unsigned int>())
    ;

    class_<StructuredMeshGeneratorProcess, StructuredMeshGeneratorProcess::Pointer, Process>(m,"StructuredMeshGeneratorProcess")
            .def(init<const Geometry< Node<3> >&, ModelPart&, Parameters&>()) //TODO: VERIFY IF THE NEXT IS NEEDED: [with_custodian_and_ward<1, 2>()])
    ;

    class_<TetrahedralMeshOrientationCheck, TetrahedralMeshOrientationCheck::Pointer, Process>(m,"TetrahedralMeshOrientationCheck")
            .def(init<ModelPart&, bool>())
    .def("SwapAll",&TetrahedralMeshOrientationCheck::SwapAll)
    .def("SwapNegativeElements",&TetrahedralMeshOrientationCheck::SwapNegativeElements)
    ;

    class_<ComputeBDFCoefficientsProcess, ComputeBDFCoefficientsProcess::Pointer, Process>(m,"ComputeBDFCoefficientsProcess")
            .def(init<ModelPart&, const unsigned int>())
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_<VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"VariationalDistanceCalculationProcess2D")
            .def(init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;
    class_<VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"VariationalDistanceCalculationProcess3D")
            .def(init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;

    class_<LevelSetConvectionProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, LevelSetConvectionProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"LevelSetConvectionProcess2D")
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer>())
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double>())
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double>())
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double, const unsigned int>())
    ;
    class_<LevelSetConvectionProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, LevelSetConvectionProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"LevelSetConvectionProcess3D")
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer>())
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double>())
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double>())
        .def(init<Variable<double>&, ModelPart&, LinearSolverType::Pointer, const double, const double, const unsigned int>())
    ;

    class_<ApplyConstantScalarValueProcess, ApplyConstantScalarValueProcess::Pointer, Process>(m,"ApplyConstantScalarValueProcess")
            .def(init<ModelPart&, Parameters>())
            .def(init<ModelPart&, const Variable<double>&, double, std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, double, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<int>&, int, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<bool>&, bool, std::size_t, Flags>())
            .def("ExecuteInitialize", &ApplyConstantScalarValueProcess::ExecuteInitialize)
            .def_readonly_static("VARIABLE_IS_FIXED", &ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)
    ;

    class_<ApplyConstantVectorValueProcess, ApplyConstantVectorValueProcess::Pointer, Process>(m,"ApplyConstantVectorValueProcess")
            .def(init<ModelPart&, Parameters>())
            .def(init<ModelPart&, const Variable<array_1d<double, 3 > >& , const double, const Vector , std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def_readonly_static("X_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::X_COMPONENT_FIXED)
            .def_readonly_static("Y_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Y_COMPONENT_FIXED)
            .def_readonly_static("Z_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Z_COMPONENT_FIXED)
    ;

    class_<CheckSkinProcess, CheckSkinProcess::Pointer, Process>(m,"CheckSkinProcess")
            .def(init<ModelPart&, Flags>())
    ;

    class_<ReplaceElementsAndConditionsProcess, ReplaceElementsAndConditionsProcess::Pointer, Process>(m,"ReplaceElementsAndConditionsProcess")
            .def(init<ModelPart&, Parameters>())
    ;

    /* Historical */
    // DOUBLE
    class_<ComputeNodalGradientProcess<2, Variable<double>, Historical>, ComputeNodalGradientProcess<2, Variable<double>, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcess2D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, Variable<double>, Historical>, ComputeNodalGradientProcess<3, Variable<double>, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcess3D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    // COMPONENT
    class_<ComputeNodalGradientProcess<2, component_type, Historical>, ComputeNodalGradientProcess<2, component_type, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcessComp2D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, component_type, Historical>, ComputeNodalGradientProcess<3, component_type, Historical>::Pointer, Process>(m,"ComputeNodalGradientProcessComp3D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    /* Non-Historical */
    // DOUBLE
    class_<ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>, ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcess2D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    class_<ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>, ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcess3D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    // COMPONENT
    class_<ComputeNodalGradientProcess<2, component_type, NonHistorical>, ComputeNodalGradientProcess<2, component_type, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcessComp2D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, component_type, NonHistorical>, ComputeNodalGradientProcess<3, component_type, NonHistorical>::Pointer, Process>(m,"ComputeNonHistoricalNodalGradientProcessComp3D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<CalculateDiscontinuousDistanceToSkinProcess, CalculateDiscontinuousDistanceToSkinProcess::Pointer, Process>(m,"CalculateDiscontinuousDistanceToSkinProcess")
            .def(init<ModelPart&, ModelPart&>())
            ;

    class_<ReorderAndOptimizeModelPartProcess, ReorderAndOptimizeModelPartProcess::Pointer, Process>(m,"ReorderAndOptimizeModelPartProcess")
            .def(init<ModelPart&, Parameters>())
            ;


    class_<AssignScalarVariableToConditionsProcess, AssignScalarVariableToConditionsProcess::Pointer, Process>(m,"AssignScalarVariableToConditionsProcess")
            .def(init<ModelPart&, Parameters >())
    ;

    class_<AssignScalarFieldToConditionsProcess, AssignScalarFieldToConditionsProcess::Pointer, Process>(m,"AssignScalarFieldToConditionsProcess")
            .def(init<ModelPart&, Parameters >())
    ;


    //typedef PointerVectorSet<Node<3>, IndexedObject> NodesContainerType;
    //typedef PointerVectorSet<Dof<double>, IndexedObject> DofsContainerType;

    //class_<AddDofsNodalProcess<Variable<double> >, AddDofsNodalProcess<Variable<double> >::Pointer, Process>(m,"AddDoubleDofsNodalProcess")
    // .def(init<Variable<double>, NodesContainerType&, DofsContainerType&>())
    // ;
    //class_<AddDofsNodalProcess<VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >, AddDofsNodalProcess<VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >::Pointer, Process>(m,"AddArrayComponentDofsNodalProcess")
    // ;

    /* Simple Mortar mapper */
    // 2D
    class_<SimpleMortarMapperProcess<2, 2, Variable<double>>, SimpleMortarMapperProcess<2, 2, Variable<double>>::Pointer, Process>(m, "SimpleMortarMapperProcess2D2NDouble")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >>, SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >>::Pointer, Process>(m, "SimpleMortarMapperProcess2D2NVector")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Triangle
    class_<SimpleMortarMapperProcess<3, 3, Variable<double>>, SimpleMortarMapperProcess<3, 3, Variable<double>>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3NDouble")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >>, SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3NVector")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Quadrilateral
    class_<SimpleMortarMapperProcess<3, 4, Variable<double>>, SimpleMortarMapperProcess<3, 4, Variable<double>>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4NDouble")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >>, SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4NVector")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Triangle - Quadrilateral
    class_<SimpleMortarMapperProcess<3, 3, Variable<double>, 4>, SimpleMortarMapperProcess<3, 3, Variable<double>, 4>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3N4NDouble")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, 4>, SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, 4>::Pointer, Process>(m, "SimpleMortarMapperProcess3D3N4NVector")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // 3D - Quadrilateral - Triangle
    class_<SimpleMortarMapperProcess<3, 4, Variable<double>, 3>, SimpleMortarMapperProcess<3, 4, Variable<double>, 3>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4N3NDouble")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, 3>, SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, 3>::Pointer, Process>(m, "SimpleMortarMapperProcess3D4N3NVector")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    ;

    // Transfer between model parts
    class_<FastTransferBetweenModelPartsProcess, FastTransferBetweenModelPartsProcess::Pointer, Process> FastTransferBetweenModelPartsProcess_Scope(m, "FastTransferBetweenModelPartsProcess");

    FastTransferBetweenModelPartsProcess_Scope.def(init<ModelPart&, ModelPart&>());
    FastTransferBetweenModelPartsProcess_Scope.def(init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered>());
    FastTransferBetweenModelPartsProcess_Scope.def(init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered, const Flags >());
    FastTransferBetweenModelPartsProcess_Scope.def(init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered, const Flags, const bool >());

    // Adding FastTransferBetweenModelPartsProcess related enums
    enum_<FastTransferBetweenModelPartsProcess::EntityTransfered>(FastTransferBetweenModelPartsProcess_Scope, "EntityTransfered")
    .value("NODES", FastTransferBetweenModelPartsProcess::EntityTransfered::NODES)
    .value("ELEMENTS", FastTransferBetweenModelPartsProcess::EntityTransfered::ELEMENTS)
    .value("NODESANDELEMENTS", FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS)
    .value("CONDITIONS", FastTransferBetweenModelPartsProcess::EntityTransfered::CONDITIONS)
    .value("NODESANDCONDITIONS", FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDCONDITIONS)
    .value("CONSTRAINTS", FastTransferBetweenModelPartsProcess::EntityTransfered::CONSTRAINTS)
    .value("NODESANDCONSTRAINTS", FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDCONSTRAINTS)
    .value("ALL", FastTransferBetweenModelPartsProcess::EntityTransfered::ALL)
    ;

    class_<SkinDetectionProcess<2>, SkinDetectionProcess<2>::Pointer, Process>(m, "SkinDetectionProcess2D")
        .def(init<ModelPart&>())
        .def(init< ModelPart&, Parameters >())
        ;

    class_<SkinDetectionProcess<3>, SkinDetectionProcess<3>::Pointer, Process>(m, "SkinDetectionProcess3D")
        .def(init<ModelPart&>())
        .def(init< ModelPart&, Parameters >())
        ;
}

}  // namespace Python.

} // Namespace Kratos
