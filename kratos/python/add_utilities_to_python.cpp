//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
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
#include "processes/process.h"
#include "python/add_utilities_to_python.h"
#include "utilities/variable_utils.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/body_normal_calculation_utils.h"
#include "utilities/body_distance_calculation_utils.h"
#include "utilities/signed_distance_calculation_utils.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "utilities/openmp_utils.h"
#include "utilities/brute_force_point_locator.h"
#include "utilities/deflation_utils.h"
#include "utilities/iso_printer.h"
#include "utilities/activation_utilities.h"
#include "utilities/convect_particles_utilities.h"
#include "utilities/condition_number_utility.h"
#include "utilities/mortar_utilities.h"
#include "utilities/read_materials_utility.h"


// #include "utilities/signed_distance_calculator_bin_based.h"
#include "utilities/divide_elem_utils.h"
#include "utilities/timer.h"

//#include "spatial_containers/bounding_box.h"
#include "utilities/bounding_box_utilities.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"
#include "utilities/geometry_tester.h"
#include "utilities/cutting_utility.h"

#include "utilities/python_function_callback_utility.h"
#include "utilities/interval_utility.h"
#include "utilities/table_stream_utility.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/sub_model_parts_list_utility.h"

namespace Kratos
{

namespace Python
{

/**
 * @brief Sets the current table utility on the process info
 * @param rCurrentProcessInfo The process info
 */
void SetOnProcessInfo(
    typename TableStreamUtility::Pointer pTable,
    ProcessInfo& rCurrentProcessInfo
    )
{
    rCurrentProcessInfo[TABLE_UTILITY] = pTable;
}


void AddUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    // NOTE: this function is special in that it accepts a "pyObject" - this is the reason for which it is defined in this same file
    class_<PythonGenericFunctionUtility,  PythonGenericFunctionUtility::Pointer >(m,"PythonGenericFunctionUtility")
    .def(init<const std::string&>() )
    .def(init<const std::string&, Parameters>())
    .def("UseLocalSystem", &PythonGenericFunctionUtility::UseLocalSystem)
    .def("DependsOnSpace", &PythonGenericFunctionUtility::DependsOnSpace)
    .def("RotateAndCallFunction", &PythonGenericFunctionUtility::RotateAndCallFunction)
    .def("CallFunction", &PythonGenericFunctionUtility::CallFunction)
    ;

    class_<ApplyFunctionToNodesUtility >(m,"ApplyFunctionToNodesUtility")
    .def(init<ModelPart::NodesContainerType&, PythonGenericFunctionUtility::Pointer >() )
    .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction< Variable<double> >)
    .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
    .def("ReturnFunction", &ApplyFunctionToNodesUtility::ReturnFunction)
    ;


    class_<DeflationUtils>(m,"DeflationUtils")
    .def(init<>())
    .def("VisualizeAggregates",&DeflationUtils::VisualizeAggregates)
    ;

    // This is required to recognize the different overloads of ConditionNumberUtility::GetConditionNumber
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef double (ConditionNumberUtility::*InputGetConditionNumber)(SparseSpaceType::MatrixType&, LinearSolverType::Pointer, LinearSolverType::Pointer);
    typedef double (ConditionNumberUtility::*DirectGetConditionNumber)(SparseSpaceType::MatrixType&);

    InputGetConditionNumber ThisGetConditionNumber = &ConditionNumberUtility::GetConditionNumber;
    DirectGetConditionNumber ThisDirectGetConditionNumber = &ConditionNumberUtility::GetConditionNumber;

    class_<ConditionNumberUtility>(m,"ConditionNumberUtility")
    .def(init<>())
    .def(init<LinearSolverType::Pointer, LinearSolverType::Pointer>())
    .def("GetConditionNumber", ThisGetConditionNumber)
    .def("GetConditionNumber", ThisDirectGetConditionNumber)
    ;

    class_<VariableUtils>(m, "VariableUtils")
        .def(init<>())
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<bool>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<double>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<array_1d<double,3>>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<Matrix>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<bool>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<double>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,3>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Vector>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Matrix>>)
        .def("SetVectorVar", &VariableUtils::SetVectorVar)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<Variable<double>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetNonHistoricalVectorVar", &VariableUtils::SetNonHistoricalVectorVar)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalScalarVar<Variable<double>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetVariable", &VariableUtils::SetVariable<bool>)
        .def("SetVariable", &VariableUtils::SetVariable<double>)
        .def("SetVariable", &VariableUtils::SetVariable<array_1d<double, 3>>)
        .def("SetVariable", &VariableUtils::SetVariable<Vector>)
        .def("SetVariable", &VariableUtils::SetVariable<Matrix>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::NodesContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ConditionsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ElementsContainerType>)
        .def("SaveVectorVar", &VariableUtils::SaveVectorVar)
        .def("SaveScalarVar", &VariableUtils::SaveScalarVar)
        .def("SaveVectorNonHistoricalVar", &VariableUtils::SaveVectorNonHistoricalVar)
        .def("SaveScalarNonHistoricalVar", &VariableUtils::SaveScalarNonHistoricalVar)
        .def("SelectNodeList", &VariableUtils::SelectNodeList)
        .def("CopyVectorVar", &VariableUtils::CopyVectorVar)
        .def("CopyScalarVar", &VariableUtils::CopyScalarVar)
        .def("SetToZero_VectorVar", &VariableUtils::SetToZero_VectorVar)
        .def("SetToZero_ScalarVar", &VariableUtils::SetToZero_ScalarVar)
        // .def("SetToZero_VelocityVectorVar", &VariableUtils::SetToZero_VelocityVectorVar)
        // .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<double> >)
        // .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > > )
        // .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 3> > > )
        .def("ApplyFixity", &VariableUtils::ApplyFixity<Variable<double>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<Variable<double>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumHistoricalNodeVectorVariable", &VariableUtils::SumHistoricalNodeVectorVariable)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumNonHistoricalNodeVectorVariable", &VariableUtils::SumNonHistoricalNodeVectorVariable)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<Variable<double>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumConditionVectorVariable", &VariableUtils::SumConditionVectorVariable)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<Variable<double>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumElementVectorVariable", &VariableUtils::SumElementVectorVariable)
        .def("AddDof", &VariableUtils::AddDof<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("CheckVariableKeys", &VariableUtils::CheckVariableKeys)
        .def("CheckDofs", &VariableUtils::CheckDofs);

    // This is required to recognize the different overloads of NormalCalculationUtils::CalculateOnSimplex
    typedef  void (NormalCalculationUtils::*CalcOnSimplexCondType)(NormalCalculationUtils::ConditionsArrayType&,int);
    typedef  void (NormalCalculationUtils::*CalcOnSimplexMPType)(ModelPart&,int);
    typedef  void (NormalCalculationUtils::*CalcOnSimplexWithDoubleVarType)(ModelPart&,int,Variable<double>&);
    typedef  void (NormalCalculationUtils::*CalcOnSimplexWithIntVarType)(ModelPart&,int,Variable<int>&);
//            typedef  void (NormalCalculationUtils::*CalcOnSimplexWithArrayVarType)(ModelPart&,int,Variable< array_1d<double,3> >&,const array_1d<double,3>&);
    typedef  void (NormalCalculationUtils::*CalcOnSimplexWithDoubleVarAlphaType)(ModelPart&,int,Variable<double>&,const double,const double);

    CalcOnSimplexCondType CalcOnSimplex_Cond = &NormalCalculationUtils::CalculateOnSimplex;
    CalcOnSimplexMPType CalcOnSimplex_ModelPart = &NormalCalculationUtils::CalculateOnSimplex;
    CalcOnSimplexWithDoubleVarType CalcOnSimplexWithDoubleVar = &NormalCalculationUtils::CalculateOnSimplex;
    CalcOnSimplexWithIntVarType CalcOnSimplexWithIntVar = &NormalCalculationUtils::CalculateOnSimplex;
    CalcOnSimplexWithDoubleVarAlphaType CalcOnSimplexWithDoubleVarAlpha = &NormalCalculationUtils::CalculateOnSimplex;

    class_<NormalCalculationUtils > (m,"NormalCalculationUtils")
    .def(init<>())
    .def("CalculateOnSimplex", CalcOnSimplex_Cond)
    .def("CalculateOnSimplex", CalcOnSimplex_ModelPart)
    .def("CalculateOnSimplex", CalcOnSimplexWithDoubleVar)
    .def("CalculateOnSimplex", CalcOnSimplexWithIntVar)
    .def("CalculateOnSimplex", CalcOnSimplexWithDoubleVarAlpha)
    .def("SwapNormals", &NormalCalculationUtils::SwapNormals)
//                    .def("CalculateOnSimplex", CalcOnSimplexWithArrayVar)
    ;

    class_<BodyNormalCalculationUtils > (m,"BodyNormalCalculationUtils")
    .def(init<>())
    .def("CalculateBodyNormals", &BodyNormalCalculationUtils::CalculateBodyNormals)
    ;

    class_<BodyDistanceCalculationUtils > (m,"BodyDistanceCalculationUtils")
    .def(init<>())
    .def("CalculateDistances2D", &BodyDistanceCalculationUtils::CalculateDistances < 2 >)
    .def("CalculateDistances3D", &BodyDistanceCalculationUtils::CalculateDistances < 3 >)
    ;

    class_<SignedDistanceCalculationUtils < 2 > >(m,"SignedDistanceCalculationUtils2D")
    .def(init<>())
    .def("CalculateDistances", &SignedDistanceCalculationUtils < 2 > ::CalculateDistances)
    .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 2 > ::FindMaximumEdgeSize)
    ;

    class_<SignedDistanceCalculationUtils < 3 > >(m,"SignedDistanceCalculationUtils3D")
    .def(init<>())
    .def("CalculateDistances", &SignedDistanceCalculationUtils < 3 > ::CalculateDistances)
    .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 3 > ::FindMaximumEdgeSize)
    ;

    class_<ParallelDistanceCalculator < 2 > >(m,"ParallelDistanceCalculator2D")
    .def(init<>())
    .def("CalculateDistances", &ParallelDistanceCalculator < 2 > ::CalculateDistances)
    .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 2 > ::CalculateInterfacePreservingDistances)
    .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 2 > ::CalculateDistancesLagrangianSurface)
    .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 2 > ::FindMaximumEdgeSize)
    ;

    class_<ParallelDistanceCalculator < 3 > >(m,"ParallelDistanceCalculator3D")
    .def(init<>())
    .def("CalculateDistances", &ParallelDistanceCalculator < 3 > ::CalculateDistances)
    .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 3 > ::CalculateInterfacePreservingDistances)
    .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 3 > ::CalculateDistancesLagrangianSurface)
    .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 3 > ::FindMaximumEdgeSize)
    ;

    class_<BruteForcePointLocator> (m, "BruteForcePointLocator")
    .def(init<ModelPart& >())
    .def("FindNode", &BruteForcePointLocator::FindNode)
    .def("FindElement", &BruteForcePointLocator::FindElement)
    .def("FindCondition", &BruteForcePointLocator::FindCondition)
    ;

    class_<ParticleConvectUtily<2> >(m,"ParticleConvectUtily2D")
    .def(init< BinBasedFastPointLocator < 2 >::Pointer >())
    .def("MoveParticles_Substepping", &ParticleConvectUtily<2>::MoveParticles_Substepping)
    .def("MoveParticles_RK4", &ParticleConvectUtily<2>::MoveParticles_RK4)
    ;

    class_<ParticleConvectUtily<3> >(m,"ParticleConvectUtily3D")
    .def(init< BinBasedFastPointLocator < 3 >::Pointer >())
    .def("MoveParticles_Substepping", &ParticleConvectUtily<3>::MoveParticles_Substepping)
    .def("MoveParticles_RK4", &ParticleConvectUtily<3>::MoveParticles_RK4)
    ;



    class_<IsosurfacePrinterApplication >(m,"IsosurfacePrinterApplication")
    .def(init<ModelPart& >() )
    .def("AddScalarVarIsosurface", &IsosurfacePrinterApplication::AddScalarVarIsosurface)
    .def("AddScalarVarIsosurfaceAndLower", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndLower)
    .def("AddScalarVarIsosurfaceAndHigher", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndHigher)
    .def("ClearData", &IsosurfacePrinterApplication::ClearData)
    .def("AddSkinConditions", &IsosurfacePrinterApplication::AddSkinConditions)
    .def("CreateNodesArray", &IsosurfacePrinterApplication::CreateNodesArray)
    ;


    // 	  class_<SignedDistanceCalculationBinBased<2> >(m,"SignedDistanceCalculationBinBased2D", init<>())
    // 			  .def("CalculateDistances",&SignedDistanceCalculationBinBased<2>::CalculateDistances )
    //                           .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<2>::FindMaximumEdgeSize )
    // 			  ;
    //
    // 	  class_<SignedDistanceCalculationBinBased<3> >(m,"SignedDistanceCalculationBinBased3D", init<>())
    // 			  .def("CalculateDistances",&SignedDistanceCalculationBinBased<3>::CalculateDistances )
    //                           .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<3>::FindMaximumEdgeSize )
    // 			  ;

    class_<DivideElemUtils >(m,"DivideElemUtils")
    .def(init<>())
    .def("DivideElement_2D", &DivideElemUtils::DivideElement_2D)
    ;

    class_<Timer >(m,"Timer")
    .def(init<>())
    .def_property("PrintOnScreen", &Timer::GetPrintOnScreen, &Timer::SetPrintOnScreen)
    .def_static("Start", &Timer::Start)
    .def_static("Stop", &Timer::Stop)
//     .staticmethod("Start")
//     .staticmethod("Stop")
    // 	    .def("PrintTimingInformation",Timer::PrintTimingInformation)
    .def("__repr__",&Timer::Info)
    ;




    class_<BoundingBoxUtilities >(m,"BoundingBoxUtilities")
    .def(init<ModelPart&, const unsigned int& >())
    .def("Test", &BoundingBoxUtilities::Test)
    ;


    //           class_<SplitElements, boost::noncopyable >
    //                     ("SplitElements", init<ModelPart&, int >() )
    //                     .def("Split", &SplitElements::Split)
    //                     ;


    // 	  def("PrintTimingInformation",Timer::PrintTimingInformation);

    class_<OpenMPUtils >(m,"OpenMPUtils")
    .def(init<>())
    .def_static("SetNumThreads", &OpenMPUtils::SetNumThreads)
//     .staticmethod("SetNumThreads")
    .def_static("GetNumThreads", &OpenMPUtils::GetNumThreads)
//     .staticmethod("GetNumThreads")
    .def_static("PrintOMPInfo", &OpenMPUtils::PrintOMPInfo)
//     .staticmethod("PrintOMPInfo")
    ;

    class_< BinBasedFastPointLocator < 2 > >(m,"BinBasedFastPointLocator2D")
    .def(init<ModelPart& >())
    .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabase)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
    .def("FindPointOnMesh", &BinBasedFastPointLocator < 2 > ::FindPointOnMeshSimplified)
    ;

    class_< BinBasedFastPointLocator < 3 > >(m,"BinBasedFastPointLocator3D")
    .def(init<ModelPart&  >())
    .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabase)
    .def("FindPointOnMesh", &BinBasedFastPointLocator < 3 > ::FindPointOnMeshSimplified)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
    ;

    class_< BinBasedNodesInElementLocator < 2 > >(m,"BinBasedNodesInElementLocator2D")
    .def(init<ModelPart& >())
    .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabase)
    .def("FindNodesInElement", &BinBasedNodesInElementLocator < 2 > ::FindNodesInElement)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
    ;

    class_< BinBasedNodesInElementLocator < 3 > >(m,"BinBasedNodesInElementLocator3D")
    .def(init<ModelPart&  >())
    .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabase)
    .def("FindNodesInElement", &BinBasedNodesInElementLocator < 3 > ::FindNodesInElement)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
    ;

    class_< ActivationUtilities >(m,"ActivationUtilities")
    .def(init< >())
    .def("ActivateElementsAndConditions", &ActivationUtilities::ActivateElementsAndConditions)
    ;

    class_< GeometryTesterUtility>(m,"GeometryTesterUtility")
    .def(init< >())
    .def("RunTest", &GeometryTesterUtility::RunTest)
    .def("TestTriangle2D3N", &GeometryTesterUtility::TestTriangle2D3N)
    .def("TestTriangle2D6N", &GeometryTesterUtility::TestTriangle2D6N)
    .def("TestTetrahedra3D4N", &GeometryTesterUtility::TestTetrahedra3D4N)
    .def("TestTetrahedra3D10N", &GeometryTesterUtility::TestTetrahedra3D10N)
    .def("TestHexahedra3D8N", &GeometryTesterUtility::TestHexahedra3D8N)
    .def("TestHexahedra3D27N", &GeometryTesterUtility::TestHexahedra3D27N)
    .def("TestHexahedra3D20N", &GeometryTesterUtility::TestHexahedra3D20N)
    ;

    class_<CuttingUtility >(m,"CuttingUtility")
    .def(init< >())
    .def("GenerateCut", &CuttingUtility::GenerateCut)
    .def("UpdateCutData", &CuttingUtility ::UpdateCutData)
    .def("AddSkinConditions", &CuttingUtility ::AddSkinConditions)
    .def("AddVariablesToCutModelPart", &CuttingUtility::AddVariablesToCutModelPart )
    .def("FindSmallestEdge", &CuttingUtility ::FindSmallestEdge)
    ;

    class_<IntervalUtility >(m,"IntervalUtility")
    .def(init<Parameters >())
    .def("GetIntervalBegin", &IntervalUtility::GetIntervalBegin)
    .def("GetIntervalEnd", &IntervalUtility::GetIntervalEnd)
    .def("IsInInterval", &IntervalUtility ::IsInInterval)
    ;

    // Adding table from table stream to python
    class_<TableStreamUtility, typename TableStreamUtility::Pointer>(m,"TableStreamUtility")
    .def(init<>())
    .def(init< bool >())
    .def("SetOnProcessInfo",SetOnProcessInfo)
    ;

    // Exact integration (for testing)
    class_<ExactMortarIntegrationUtility<2,2>>(m,"ExactMortarIntegrationUtility2D2N")
    .def(init<>())
    .def(init<const unsigned int>())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactIntegration)
    .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactAreaIntegration)
    ;
    class_<ExactMortarIntegrationUtility<3,3>>(m,"ExactMortarIntegrationUtility3D3N")
    .def(init<>())
    .def(init<const unsigned int>())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactIntegration)
    .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactAreaIntegration)
    .def("TestGiDDebug",&ExactMortarIntegrationUtility<3,3>::TestGiDDebug)
    ;

    class_<ExactMortarIntegrationUtility<3,4>>(m,"ExactMortarIntegrationUtility3D4N")
    .def(init<>())
    .def(init<const unsigned int>())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactIntegration)
    .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactAreaIntegration)
    .def("TestGiDDebug",&ExactMortarIntegrationUtility<3,4>::TestGiDDebug)
    ;

    // Sparse matrix multiplication utility
    class_<SparseMatrixMultiplicationUtility, typename SparseMatrixMultiplicationUtility::Pointer>(m, "SparseMatrixMultiplicationUtility")
    .def(init<>())
    .def("MatrixMultiplication",&SparseMatrixMultiplicationUtility::MatrixMultiplication<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
    .def("MatrixMultiplicationSaad",&SparseMatrixMultiplicationUtility::MatrixMultiplicationSaad<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
    .def("MatrixMultiplicationRMerge",&SparseMatrixMultiplicationUtility::MatrixMultiplicationRMerge<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
    .def("MatrixAdd",&SparseMatrixMultiplicationUtility::MatrixAdd<CompressedMatrix, CompressedMatrix>)
    ;

    // Mortar utilities
    class_<MortarUtilities, typename MortarUtilities::Pointer>(m, "MortarUtilities")
    .def(init<>())
    .def("ComputeNodesMeanNormalModelPart",&MortarUtilities::ComputeNodesMeanNormalModelPart)
    .def("InvertNormal",&MortarUtilities::InvertNormal<PointerVectorSet<Element, IndexedObject>>)
    .def("InvertNormal",&MortarUtilities::InvertNormal<PointerVectorSet<Condition, IndexedObject>>)
    ;

    // Read materials utility
    class_<ReadMaterialsUtility, typename ReadMaterialsUtility::Pointer>(m, "ReadMaterialsUtility")
    .def(init<Parameters, Model&>())
    ;

    // SubModelParts List Utility
    class_<SubModelPartsListUtility, typename SubModelPartsListUtility::Pointer>(m, "SubModelPartsListUtility")
    .def(init<ModelPart&>())
    .def("DebugComputeSubModelPartsList",&SubModelPartsListUtility::DebugComputeSubModelPartsList)
    .def("GetRecursiveSubModelPartNames",&SubModelPartsListUtility::GetRecursiveSubModelPartNames)
    .def("GetRecursiveSubModelPart",&SubModelPartsListUtility::GetRecursiveSubModelPart)
    ;
}

} // namespace Python.

} // Namespace Kratos
