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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
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
#include "includes/node.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "utilities/python_function_callback_utility.h"


namespace Kratos
{

namespace Python
{
typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;

void  AddProcessesToPython()
{
    using namespace boost::python;

    class_<Process>("Process")
    .def("Execute",&Process::Execute)
    .def("ExecuteInitialize",&Process::ExecuteInitialize)
    .def("ExecuteBeforeSolutionLoop",&Process::ExecuteBeforeSolutionLoop)
    .def("ExecuteInitializeSolutionStep",&Process::ExecuteInitializeSolutionStep)
    .def("ExecuteFinalizeSolutionStep",&Process::ExecuteFinalizeSolutionStep)
    .def("ExecuteBeforeOutputStep",&Process::ExecuteBeforeOutputStep)
    .def("ExecuteAfterOutputStep",&Process::ExecuteAfterOutputStep)
    .def("ExecuteFinalize",&Process::ExecuteFinalize)
    .def(self_ns::str(self))
    ;

    class_<FastTransferBetweenModelPartsProcess, bases<Process> >("FastTransferBetweenModelPartsProcess",init<ModelPart&, ModelPart&, const std::string>())
    .def("Execute",&FastTransferBetweenModelPartsProcess::Execute)
    ;
    
    class_<FindNodalHProcess, bases<Process> >("FindNodalHProcess",init<ModelPart&>())
    .def("Execute",&FindNodalHProcess::Execute)
    ;
    
    class_<FindNodalNeighboursProcess, bases<Process> >("FindNodalNeighboursProcess",
            init<ModelPart&, int, int>())
    .def("ClearNeighbours",&FindNodalNeighboursProcess::ClearNeighbours)
    ;

    class_<FindConditionsNeighboursProcess, bases<Process> >("FindConditionsNeighboursProcess",
            init<ModelPart&, int, int>())
    .def("ClearNeighbours",&FindConditionsNeighboursProcess::ClearNeighbours)
    ;

    class_<FindElementalNeighboursProcess, bases<Process> >("FindElementalNeighboursProcess",
            init<ModelPart&, int, int>())
    .def("ClearNeighbours",&FindElementalNeighboursProcess::ClearNeighbours)
    ;

    class_<CalculateNodalAreaProcess, bases<Process> >("CalculateNodalAreaProcess",
            init<ModelPart&, unsigned int>())
    ;

    class_<NodeEraseProcess, bases<Process> >("NodeEraseProcess",
            init<ModelPart&>())
    ;

    class_<ElementEraseProcess, bases<Process> >("ElementEraseProcess",
            init<ModelPart&>())
    ;

    class_<ConditionEraseProcess, bases<Process> >("ConditionEraseProcess",
            init<ModelPart&>())
    ;

    class_<EliminateIsolatedNodesProcess, bases<Process> >("EliminateIsolatedNodesProcess",
            init<ModelPart&>())
    ;

    class_<CalculateSignedDistanceTo3DSkinProcess, bases<Process>, boost::noncopyable >("CalculateSignedDistanceTo3DSkinProcess",
            init<ModelPart&, ModelPart&>())
    .def("GenerateSkinModelPart",&CalculateSignedDistanceTo3DSkinProcess::GenerateSkinModelPart)
    .def("MappingPressureToStructure",&CalculateSignedDistanceTo3DSkinProcess::MappingPressureToStructure)
    ;

    class_<CalculateEmbeddedSignedDistanceTo3DSkinProcess, bases<Process>, boost::noncopyable >("CalculateEmbeddedSignedDistanceTo3DSkinProcess",
            init< ModelPart&, ModelPart& >())
    .def(init< ModelPart&, ModelPart&, bool>())
    ;

   class_<CalculateSignedDistanceTo3DConditionSkinProcess, bases<Process> >("CalculateSignedDistanceTo3DConditionSkinProcess",
            init<ModelPart&, ModelPart&>())
    ;

    class_<TranslationOperation, bases<Process> >("TranslationOperation",
            init<ModelPart&, boost::numeric::ublas::vector<int> ,boost::numeric::ublas::vector<int> ,unsigned int>())
    ;

    class_<RotationOperation, bases<Process> >("RotationOperation",
            init<ModelPart&, boost::numeric::ublas::vector<int> ,boost::numeric::ublas::vector<int> ,unsigned int>())
    ;

    class_<StructuredMeshGeneratorProcess, bases<Process>, boost::noncopyable >("StructuredMeshGeneratorProcess",
            init<const Geometry< Node<3> >&, ModelPart&, Parameters&>()[with_custodian_and_ward<1, 2>()])
    ;

    class_<TetrahedralMeshOrientationCheck, bases<Process>, boost::noncopyable >("TetrahedralMeshOrientationCheck",
            init<ModelPart&, bool>())
    .def("SwapAll",&TetrahedralMeshOrientationCheck::SwapAll)
    .def("SwapNegativeElements",&TetrahedralMeshOrientationCheck::SwapNegativeElements)
    ;

    class_<ComputeBDFCoefficientsProcess, bases<Process>, boost::noncopyable >("ComputeBDFCoefficientsProcess",
            init<ModelPart&, const unsigned int>())
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    class_<VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType > , bases<Process>, boost::noncopyable >("VariationalDistanceCalculationProcess2D",
            init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;
    class_<VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType > , bases<Process>, boost::noncopyable >("VariationalDistanceCalculationProcess3D",
            init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;

    class_<LevelSetConvectionProcess<2> , bases<Process>, boost::noncopyable >("LevelSetConvectionProcess2D",
            init<Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double >())
    .def(init< Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double, double>())
    .def(init< Variable<double>&, ModelPart&, LinearSolverType::Pointer, double, double,int>())
    ;
    class_<LevelSetConvectionProcess<3> , bases<Process>, boost::noncopyable >("LevelSetConvectionProcess3D",
            init<Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double>())
            .def(init< Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double, double>())
			.def(init< Variable<double>&, ModelPart&, LinearSolverType::Pointer, double, double,int>())
    ;

    class_<ApplyConstantScalarValueProcess , bases<Process>, boost::noncopyable >("ApplyConstantScalarValueProcess",
            init<ModelPart&, Parameters>())
            .def(init<ModelPart&, const Variable<double>&, double, std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, double, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<int>&, int, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<bool>&, bool, std::size_t, Flags>())
            .def("ExecuteInitialize", &ApplyConstantScalarValueProcess::ExecuteInitialize)
            .def_readonly("VARIABLE_IS_FIXED", &ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)
    ;

    class_<ApplyConstantVectorValueProcess , bases<Process>, boost::noncopyable >("ApplyConstantVectorValueProcess",
            init<ModelPart&, Parameters>())
            .def(init<ModelPart&, const Variable<array_1d<double, 3 > >& , const double, const Vector , std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def_readonly("X_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::X_COMPONENT_FIXED)
            .def_readonly("Y_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Y_COMPONENT_FIXED)
            .def_readonly("Z_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Z_COMPONENT_FIXED)
    ;

    class_<CheckSkinProcess , bases<Process>, boost::noncopyable >("CheckSkinProcess",
            init<ModelPart&, Flags>())
    ;

    class_<ReplaceElementsAndConditionsProcess , bases<Process>, boost::noncopyable >("ReplaceElementsAndConditionsProcess",
            init<ModelPart&, Parameters>())
    ;

    /* Historical */
    // DOUBLE
    class_<ComputeNodalGradientProcess<2, Variable<double>, Historical> , bases<Process>, boost::noncopyable >("ComputeNodalGradientProcess2D",
            init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, Variable<double>, Historical> , bases<Process>, boost::noncopyable >("ComputeNodalGradientProcess3D",
            init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    // COMPONENT
    class_<ComputeNodalGradientProcess<2, component_type, Historical> , bases<Process>, boost::noncopyable >("ComputeNodalGradientProcessComp2D",
            init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, component_type, Historical> , bases<Process>, boost::noncopyable >("ComputeNodalGradientProcessComp3D",
            init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;
    
    /* Non-Historical */
    // DOUBLE
    class_<ComputeNodalGradientProcess<2, Variable<double>, NonHistorical> , bases<Process>, boost::noncopyable >("ComputeNonHistoricalNodalGradientProcess2D",
            init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    class_<ComputeNodalGradientProcess<3, Variable<double>, NonHistorical> , bases<Process>, boost::noncopyable >("ComputeNonHistoricalNodalGradientProcess3D",
            init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    // COMPONENT
    class_<ComputeNodalGradientProcess<2, component_type, NonHistorical> , bases<Process>, boost::noncopyable >("ComputeNonHistoricalNodalGradientProcessComp2D",
            init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, component_type, NonHistorical> , bases<Process>, boost::noncopyable >("ComputeNonHistoricalNodalGradientProcessComp3D",
            init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<CalculateDiscontinuousDistanceToSkinProcess, bases<Process>, boost::noncopyable >("CalculateDiscontinuousDistanceToSkinProcess",
            init<ModelPart&, ModelPart&>())
            ;

    class_<ReorderAndOptimizeModelPartProcess, bases<Process>, boost::noncopyable >("ReorderAndOptimizeModelPartProcess",
            init<ModelPart&, Parameters>())
            ;


    class_<AssignScalarVariableToConditionsProcess, bases<Process>, boost::noncopyable >("AssignScalarVariableToConditionsProcess",
            init<ModelPart&, Parameters >())
    ;

    class_<AssignScalarFieldToConditionsProcess , bases<Process>, boost::noncopyable >("AssignScalarFieldToConditionsProcess",
            init<ModelPart&, Parameters >())
    ;


    //typedef PointerVectorSet<Node<3>, IndexedObject> NodesContainerType;
    //typedef PointerVectorSet<Dof<double>, IndexedObject> DofsContainerType;

    //class_<AddDofsNodalProcess<Variable<double> >, bases<Process> >("AddDoubleDofsNodalProcess")
    // .def(init<Variable<double>, NodesContainerType&, DofsContainerType&>())
    // ;
    //class_<AddDofsNodalProcess<VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >, bases<Process> >("AddArrayComponentDofsNodalProcess")
    // ;

    /* Simple Mortar mapper */
    // 2D 
    class_<SimpleMortarMapperProcess<2, 2, Variable<double>, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess2D2NDoubleHistorical", init<ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<double>, Historical>::Execute)
    ;

//     class_<SimpleMortarMapperProcess<2, 2, component_type, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess2D2NComponentHistorical", init<ModelPart&, component_type&>())
//     .def(init<ModelPart&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def(init<ModelPart&, component_type&, component_type&>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def("Execute",&SimpleMortarMapperProcess<2, 2, component_type, Historical>::Execute)
//     ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess2D2NVectorHistorical", init<ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess2D2NDoubleNonHistorical", init<ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical>::Execute)
    ;

//     class_<SimpleMortarMapperProcess<2, 2, component_type, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess2D2NComponentNonHistorical", init<ModelPart&, component_type&>())
//     .def(init<ModelPart&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def(init<ModelPart&, component_type&, component_type&>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def("Execute",&SimpleMortarMapperProcess<2, 2, component_type, NonHistorical>::Execute)
//     ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess2D2NVectorNonHistorical", init<ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, NonHistorical>::Execute)
    ;

    // 3D - Triangle
    class_<SimpleMortarMapperProcess<3, 3, Variable<double>, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D3NDoubleHistorical", init<ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<double>, Historical>::Execute)
    ;

//     class_<SimpleMortarMapperProcess<3, 3, component_type, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D3NComponentHistorical", init<ModelPart&, component_type&>())
//     .def(init<ModelPart&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def(init<ModelPart&, component_type&, component_type&>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def("Execute",&SimpleMortarMapperProcess<3, 3, component_type, Historical>::Execute)
//     ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D3NVectorHistorical", init<ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D3NDoubleNonHistorical", init<ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical>::Execute)
    ;

//     class_<SimpleMortarMapperProcess<3, 3, component_type, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D3NComponentNonHistorical", init<ModelPart&, component_type&>())
//     .def(init<ModelPart&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def(init<ModelPart&, component_type&, component_type&>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def("Execute",&SimpleMortarMapperProcess<3, 3, component_type, NonHistorical>::Execute)
//     ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D3NVectorNonHistorical", init<ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, NonHistorical>::Execute)
    ;

    // 3D - Quadrilateral
    class_<SimpleMortarMapperProcess<3, 4, Variable<double>, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D4NDoubleHistorical", init<ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<double>, Historical>::Execute)
    ;

//     class_<SimpleMortarMapperProcess<3, 4, component_type, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D4NComponentHistorical", init<ModelPart&, component_type&>())
//     .def(init<ModelPart&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def(init<ModelPart&, component_type&, component_type&>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def("Execute",&SimpleMortarMapperProcess<3, 4, component_type, Historical>::Execute)
//     ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, Historical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D4NVectorHistorical", init<ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D4NDoubleNonHistorical", init<ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical>::Execute)
    ;

//     class_<SimpleMortarMapperProcess<3, 4, component_type, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D4NComponentNonHistorical", init<ModelPart&, component_type&>())
//     .def(init<ModelPart&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def(init<ModelPart&, component_type&, component_type&>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters>())
//     .def(init<ModelPart&, component_type&, component_type&, Parameters, LinearSolverType::Pointer>())
//     .def("Execute",&SimpleMortarMapperProcess<3, 4, component_type, NonHistorical>::Execute)
//     ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, NonHistorical>, bases<Process>, boost::noncopyable >("SimpleMortarMapperProcess3D4NVectorNonHistorical", init<ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, NonHistorical>::Execute)
    ;

}

}  // namespace Python.

} // Namespace Kratos
