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

    class_<Process>(m,"Process")
    .def(init<>())
    .def("Execute",&Process::Execute)
    .def("ExecuteInitialize",&Process::ExecuteInitialize)
    .def("ExecuteBeforeSolutionLoop",&Process::ExecuteBeforeSolutionLoop)
    .def("ExecuteInitializeSolutionStep",&Process::ExecuteInitializeSolutionStep)
    .def("ExecuteFinalizeSolutionStep",&Process::ExecuteFinalizeSolutionStep)
    .def("ExecuteBeforeOutputStep",&Process::ExecuteBeforeOutputStep)
    .def("ExecuteAfterOutputStep",&Process::ExecuteAfterOutputStep)
    .def("ExecuteFinalize",&Process::ExecuteFinalize)
    .def("__repr__", &Process::Info)
    ;
    
    class_<FindNodalHProcess, Process >(m,"FindNodalHProcess")
    .def(init<ModelPart&>())
    .def("Execute",&FindNodalHProcess::Execute)
    ;
    
    class_<FindNodalNeighboursProcess, Process >(m,"FindNodalNeighboursProcess")
            .def(init<ModelPart&, unsigned int, unsigned int>())
    .def("ClearNeighbours",&FindNodalNeighboursProcess::ClearNeighbours)
    ;

    class_<FindConditionsNeighboursProcess, Process >(m,"FindConditionsNeighboursProcess")
            .def(init<ModelPart&, int, unsigned int>())
    .def("ClearNeighbours",&FindConditionsNeighboursProcess::ClearNeighbours)
    ;

    class_<FindElementalNeighboursProcess, Process >(m,"FindElementalNeighboursProcess")
            .def(init<ModelPart&, int, unsigned int>())
    .def("ClearNeighbours",&FindElementalNeighboursProcess::ClearNeighbours)
    ;

    class_<CalculateNodalAreaProcess, Process >(m,"CalculateNodalAreaProcess")
            .def(init<ModelPart&, unsigned int>())
    ;

    class_<NodeEraseProcess, Process >(m,"NodeEraseProcess")
            .def(init<ModelPart&>())
    ;

    class_<ElementEraseProcess, Process >(m,"ElementEraseProcess")
            .def(init<ModelPart&>())
    ;

    class_<ConditionEraseProcess, Process >(m,"ConditionEraseProcess")
            .def(init<ModelPart&>())
    ;

    class_<EliminateIsolatedNodesProcess, Process >(m,"EliminateIsolatedNodesProcess")
            .def(init<ModelPart&>())
    ;

    class_<CalculateSignedDistanceTo3DSkinProcess, Process>(m,"CalculateSignedDistanceTo3DSkinProcess")
            .def(init<ModelPart&, ModelPart&>())
    .def("GenerateSkinModelPart",&CalculateSignedDistanceTo3DSkinProcess::GenerateSkinModelPart)
    .def("MappingPressureToStructure",&CalculateSignedDistanceTo3DSkinProcess::MappingPressureToStructure)
    ;

    class_<CalculateEmbeddedSignedDistanceTo3DSkinProcess, Process>(m,"CalculateEmbeddedSignedDistanceTo3DSkinProcess")
            .def(init< ModelPart&, ModelPart& >())
    .def(init< ModelPart&, ModelPart&, bool>())
    ;

   class_<CalculateSignedDistanceTo3DConditionSkinProcess, Process >(m,"CalculateSignedDistanceTo3DConditionSkinProcess")
            .def(init<ModelPart&, ModelPart&>())
    ;

    class_<TranslationOperation, Process >(m,"TranslationOperation")
            .def(init<ModelPart&, boost::numeric::ublas::vector<int> ,boost::numeric::ublas::vector<int> ,unsigned int>())
    ;

    class_<RotationOperation, Process >(m,"RotationOperation")
            .def(init<ModelPart&, boost::numeric::ublas::vector<int> ,boost::numeric::ublas::vector<int> ,unsigned int>())
    ;

    class_<StructuredMeshGeneratorProcess, Process>(m,"StructuredMeshGeneratorProcess")
            .def(init<const Geometry< Node<3> >&, ModelPart&, Parameters&>()) //TODO: VERIFY IF THE NEXT IS NEEDED: [with_custodian_and_ward<1, 2>()])
    ;

    class_<TetrahedralMeshOrientationCheck, Process>(m,"TetrahedralMeshOrientationCheck")
            .def(init<ModelPart&, bool>())
    .def("SwapAll",&TetrahedralMeshOrientationCheck::SwapAll)
    .def("SwapNegativeElements",&TetrahedralMeshOrientationCheck::SwapNegativeElements)
    ;

    class_<ComputeBDFCoefficientsProcess, Process>(m,"ComputeBDFCoefficientsProcess")
            .def(init<ModelPart&, const unsigned int>())
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    class_<VariationalDistanceCalculationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType > , Process>(m,"VariationalDistanceCalculationProcess2D")
            .def(init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;
    class_<VariationalDistanceCalculationProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType > , Process>(m,"VariationalDistanceCalculationProcess3D")
            .def(init<ModelPart&, LinearSolverType::Pointer, unsigned int>())
    ;

    class_<LevelSetConvectionProcess<2> , Process>(m,"LevelSetConvectionProcess2D")
            .def(init<Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double >())
    .def(init< Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double, double>())
    .def(init< Variable<double>&, ModelPart&, LinearSolverType::Pointer, double, double,int>())
    ;
    class_<LevelSetConvectionProcess<3> , Process>(m,"LevelSetConvectionProcess3D")
            .def(init<Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double>())
            .def(init< Variable<double>& , ModelPart& , LinearSolverType::Pointer ,double, double>())
			.def(init< Variable<double>&, ModelPart&, LinearSolverType::Pointer, double, double,int>())
    ;

    class_<ApplyConstantScalarValueProcess , Process>(m,"ApplyConstantScalarValueProcess")
            .def(init<ModelPart&, Parameters>())
            .def(init<ModelPart&, const Variable<double>&, double, std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def(init<ModelPart&, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >&, double, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<int>&, int, std::size_t, Flags>())
            .def(init<ModelPart&, const Variable<bool>&, bool, std::size_t, Flags>())
            .def("ExecuteInitialize", &ApplyConstantScalarValueProcess::ExecuteInitialize)
            .def_readonly_static("VARIABLE_IS_FIXED", &ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)
    ;

    class_<ApplyConstantVectorValueProcess , Process>(m,"ApplyConstantVectorValueProcess")
            .def(init<ModelPart&, Parameters>())
            .def(init<ModelPart&, const Variable<array_1d<double, 3 > >& , const double, const Vector , std::size_t, Flags>())
            .def(init< ModelPart&, Parameters& >())
            .def_readonly_static("X_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::X_COMPONENT_FIXED)
            .def_readonly_static("Y_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Y_COMPONENT_FIXED)
            .def_readonly_static("Z_COMPONENT_FIXED", &ApplyConstantVectorValueProcess::Z_COMPONENT_FIXED)
    ;

    class_<CheckSkinProcess , Process>(m,"CheckSkinProcess")
            .def(init<ModelPart&, Flags>())
    ;

    class_<ReplaceElementsAndConditionsProcess , Process>(m,"ReplaceElementsAndConditionsProcess")
            .def(init<ModelPart&, Parameters>())
    ;

    /* Historical */
    // DOUBLE
    class_<ComputeNodalGradientProcess<2, Variable<double>, Historical> , Process>(m,"ComputeNodalGradientProcess2D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, Variable<double>, Historical> , Process>(m,"ComputeNodalGradientProcess3D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    // COMPONENT
    class_<ComputeNodalGradientProcess<2, component_type, Historical> , Process>(m,"ComputeNodalGradientProcessComp2D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, component_type, Historical> , Process>(m,"ComputeNodalGradientProcessComp3D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;
    
    /* Non-Historical */
    // DOUBLE
    class_<ComputeNodalGradientProcess<2, Variable<double>, NonHistorical> , Process>(m,"ComputeNonHistoricalNodalGradientProcess2D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    class_<ComputeNodalGradientProcess<3, Variable<double>, NonHistorical> , Process>(m,"ComputeNonHistoricalNodalGradientProcess3D")
            .def(init<ModelPart&, Variable<double>&, Variable<array_1d<double,3> >& , Variable<double>& >())
            ;

    // COMPONENT
    class_<ComputeNodalGradientProcess<2, component_type, NonHistorical> , Process>(m,"ComputeNonHistoricalNodalGradientProcessComp2D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<ComputeNodalGradientProcess<3, component_type, NonHistorical> , Process>(m,"ComputeNonHistoricalNodalGradientProcessComp3D")
            .def(init<ModelPart&, component_type&, Variable<array_1d<double,3> >& , Variable<double>& >())
    ;

    class_<CalculateDiscontinuousDistanceToSkinProcess, Process>(m,"CalculateDiscontinuousDistanceToSkinProcess")
            .def(init<ModelPart&, ModelPart&>())
            ;

    class_<ReorderAndOptimizeModelPartProcess, Process>(m,"ReorderAndOptimizeModelPartProcess")
            .def(init<ModelPart&, Parameters>())
            ;


    class_<AssignScalarVariableToConditionsProcess, Process>(m,"AssignScalarVariableToConditionsProcess")
            .def(init<ModelPart&, Parameters >())
    ;

    class_<AssignScalarFieldToConditionsProcess , Process>(m,"AssignScalarFieldToConditionsProcess")
            .def(init<ModelPart&, Parameters >())
    ;


    //typedef PointerVectorSet<Node<3>, IndexedObject> NodesContainerType;
    //typedef PointerVectorSet<Dof<double>, IndexedObject> DofsContainerType;

    //class_<AddDofsNodalProcess<Variable<double> >, Process >(m,"AddDoubleDofsNodalProcess")
    // .def(init<Variable<double>, NodesContainerType&, DofsContainerType&>())
    // ;
    //class_<AddDofsNodalProcess<VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >, Process >(m,"AddArrayComponentDofsNodalProcess")
    // ;

    /* Simple Mortar mapper */
    // 2D 
    class_<SimpleMortarMapperProcess<2, 2, Variable<double>, Historical>, Process>(m, "SimpleMortarMapperProcess2D2NDoubleHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<double>, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, Historical>, Process>(m, "SimpleMortarMapperProcess2D2NVectorHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical>, Process>(m, "SimpleMortarMapperProcess2D2NDoubleNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, NonHistorical>, Process>(m, "SimpleMortarMapperProcess2D2NVectorNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, NonHistorical>::Execute)
    ;

    // 3D - Triangle
    class_<SimpleMortarMapperProcess<3, 3, Variable<double>, Historical>, Process>(m, "SimpleMortarMapperProcess3D3NDoubleHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<double>, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, Historical>, Process>(m, "SimpleMortarMapperProcess3D3NVectorHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D3NDoubleNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D3NVectorNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, NonHistorical>::Execute)
    ;

    // 3D - Quadrilateral
    class_<SimpleMortarMapperProcess<3, 4, Variable<double>, Historical>, Process>(m, "SimpleMortarMapperProcess3D4NDoubleHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<double>, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, Historical>, Process>(m, "SimpleMortarMapperProcess3D4NVectorHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D4NDoubleNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D4NVectorNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, NonHistorical>::Execute)
    ;
    
    // 2D 
    class_<SimpleMortarMapperProcess<2, 2, Variable<double>, Historical, NonHistorical>, Process>(m, "SimpleMortarMapperProcess2D2NDoubleHistoricalToNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<double>, Historical, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, Historical, NonHistorical>, Process>(m, "SimpleMortarMapperProcess2D2NVectorHistoricalToNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, Historical, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical, Historical>, Process>(m, "SimpleMortarMapperProcess2D2NDoubleNonHistoricalToHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, NonHistorical, Historical>, Process>(m, "SimpleMortarMapperProcess2D2NVectorNonHistoricalToHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<2, 2, Variable<array_1d<double,3> >, NonHistorical, Historical>::Execute)
    ;

    // 3D - Triangle
    class_<SimpleMortarMapperProcess<3, 3, Variable<double>, Historical, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D3NDoubleHistoricalToNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<double>, Historical, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, Historical, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D3NVectorHistoricalToNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, Historical, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical, Historical>, Process>(m, "SimpleMortarMapperProcess3D3NDoubleNonHistoricalToHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, NonHistorical, Historical>, Process>(m, "SimpleMortarMapperProcess3D3NVectorNonHistoricalToHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 3, Variable<array_1d<double,3> >, NonHistorical, Historical>::Execute)
    ;

    // 3D - Quadrilateral
    class_<SimpleMortarMapperProcess<3, 4, Variable<double>, Historical, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D4NDoubleHistoricalToNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<double>, Historical, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, Historical, NonHistorical>, Process>(m, "SimpleMortarMapperProcess3D4NVectorHistoricalToNonHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, Historical, NonHistorical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical, Historical>, Process>(m, "SimpleMortarMapperProcess3D4NDoubleNonHistoricalToHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<double>&, Variable<double>&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical, Historical>::Execute)
    ;

    class_<SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, NonHistorical, Historical>, Process>(m, "SimpleMortarMapperProcess3D4NVectorNonHistoricalToHistorical")
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters>())
    .def(init<ModelPart&, ModelPart&, Variable<array_1d<double,3> >&, Variable<array_1d<double,3> >&, Parameters, LinearSolverType::Pointer>())
    .def("Execute",&SimpleMortarMapperProcess<3, 4, Variable<array_1d<double,3> >, NonHistorical, Historical>::Execute)
    ;

    class_<FastTransferBetweenModelPartsProcess, Process> FastTransferBetweenModelPartsProcess_Scope(m, "FastTransferBetweenModelPartsProcess");
    
    FastTransferBetweenModelPartsProcess_Scope.def(init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered>());
    FastTransferBetweenModelPartsProcess_Scope.def(init<ModelPart&, ModelPart&, const FastTransferBetweenModelPartsProcess::EntityTransfered, const Flags >());
    FastTransferBetweenModelPartsProcess_Scope.def("Execute",&FastTransferBetweenModelPartsProcess::Execute);

    // Adding FastTransferBetweenModelPartsProcess related enums
    enum_<FastTransferBetweenModelPartsProcess::EntityTransfered>(FastTransferBetweenModelPartsProcess_Scope, "EntityTransfered")
    .value("NODES", FastTransferBetweenModelPartsProcess::EntityTransfered::NODES)
    .value("ELEMENTS", FastTransferBetweenModelPartsProcess::EntityTransfered::ELEMENTS)
    .value("NODESANDELEMENTS", FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS)
    .value("CONDITIONS", FastTransferBetweenModelPartsProcess::EntityTransfered::CONDITIONS)
    .value("ALL", FastTransferBetweenModelPartsProcess::EntityTransfered::ALL)
    ;


}

}  // namespace Python.

} // Namespace Kratos
