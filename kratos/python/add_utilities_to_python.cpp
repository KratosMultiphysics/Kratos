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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "python/add_utilities_to_python.h"
#include "utilities/variable_utils.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/body_normal_calculation_utils.h"
#include "utilities/body_distance_calculation_utils.h"
#include "utilities/signed_distance_calculation_utils.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "utilities/openmp_utils.h"
#include "utilities/pointlocation.h"
#include "utilities/deflation_utils.h"
#include "utilities/iso_printer.h"
#include "utilities/activation_utilities.h"
#include "utilities/convect_particles_utilities.h"
#include "utilities/condition_number_utility.h"


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

namespace Kratos
{

namespace Python
{


void AddUtilitiesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    
    // NOTE: this function is special in that it accepts a "pyObject" - this is the reason for which it is defined in this same file
    class_<PythonGenericFunctionUtility,  PythonGenericFunctionUtility::Pointer >("PythonGenericFunctionUtility", init<const std::string&>() )
    .def(init<const std::string&, Parameters>())
    .def("UseLocalSystem", &PythonGenericFunctionUtility::UseLocalSystem)
    .def("DependsOnSpace", &PythonGenericFunctionUtility::DependsOnSpace)
    .def("RotateAndCallFunction", &PythonGenericFunctionUtility::RotateAndCallFunction)
    .def("CallFunction", &PythonGenericFunctionUtility::CallFunction)
    ;

    class_<ApplyFunctionToNodesUtility >("ApplyFunctionToNodesUtility", init<ModelPart::NodesContainerType&, PythonGenericFunctionUtility::Pointer >() )
    .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction< Variable<double> >)
    .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
    .def("ReturnFunction", &ApplyFunctionToNodesUtility::ReturnFunction)
    ;


    class_<DeflationUtils>("DeflationUtils", init<>())
    .def("VisualizeAggregates",&DeflationUtils::VisualizeAggregates)
    ;

    // This is required to recognize the different overloads of ConditionNumberUtility::GetConditionNumber
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef double (ConditionNumberUtility::*InputGetConditionNumber)(SparseSpaceType::MatrixType&, LinearSolverType::Pointer, LinearSolverType::Pointer);
    typedef double (ConditionNumberUtility::*DirectGetConditionNumber)(SparseSpaceType::MatrixType&);

    InputGetConditionNumber ThisGetConditionNumber = &ConditionNumberUtility::GetConditionNumber;
    DirectGetConditionNumber ThisDirectGetConditionNumber = &ConditionNumberUtility::GetConditionNumber;
    
    class_<ConditionNumberUtility>("ConditionNumberUtility", init<>())
    .def(init<LinearSolverType::Pointer, LinearSolverType::Pointer>())
    .def("GetConditionNumber", ThisGetConditionNumber)
    .def("GetConditionNumber", ThisDirectGetConditionNumber)
    ;

    class_<VariableUtils > ("VariableUtils", init<>())
    .def("SetVectorVar", &VariableUtils::SetVectorVar)
    .def("SetScalarVar", &VariableUtils::SetScalarVar< Variable<double> >)
    .def("SetScalarVar", &VariableUtils::SetScalarVar< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >)
    .def("SetFlag", &VariableUtils::SetFlag< ModelPart::NodesContainerType >)
    .def("SetFlag", &VariableUtils::SetFlag< ModelPart::ConditionsContainerType >)
    .def("SetFlag", &VariableUtils::SetFlag< ModelPart::ElementsContainerType >)
    .def("SaveVectorVar", &VariableUtils::SaveVectorVar)
    .def("SaveScalarVar", &VariableUtils::SaveScalarVar)
    .def("SelectNodeList", &VariableUtils::SelectNodeList)
    .def("CopyVectorVar", &VariableUtils::CopyVectorVar)
    .def("CopyScalarVar", &VariableUtils::CopyScalarVar)
    .def("SetToZero_VectorVar", &VariableUtils::SetToZero_VectorVar)
    .def("SetToZero_ScalarVar", &VariableUtils::SetToZero_ScalarVar)
    // .def("SetToZero_VelocityVectorVar", &VariableUtils::SetToZero_VelocityVectorVar)
    // .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<double> >)
    // .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > > )
    // .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 3> > > )
    .def("ApplyFixity", &VariableUtils::ApplyFixity< Variable<double> >)
    .def("ApplyFixity", &VariableUtils::ApplyFixity< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def("ApplyVector", &VariableUtils::ApplyVector< Variable<double> >)
    .def("ApplyVector", &VariableUtils::ApplyVector< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable< Variable<double> > )
    .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def("SumHistoricalNodeVectorVariable", &VariableUtils::SumHistoricalNodeVectorVariable)
    .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable< Variable<double> > )
    .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def("SumNonHistoricalNodeVectorVariable", &VariableUtils::SumNonHistoricalNodeVectorVariable)
    .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable< Variable<double> > )
    .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def("SumConditionVectorVariable", &VariableUtils::SumConditionVectorVariable)
    .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable< Variable<double> > )
    .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def("SumElementVectorVariable", &VariableUtils::SumElementVectorVariable)
    .def("AddDof", &VariableUtils::AddDof< Variable<double> > )
    .def("AddDof", &VariableUtils::AddDof< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
    .def("AddDof", &VariableUtils::AddDofWithReaction< Variable<double> > )
    .def("AddDof", &VariableUtils::AddDofWithReaction< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > > )
	.def("CheckVariableKeys", &VariableUtils::CheckVariableKeys)
	.def("CheckDofs", &VariableUtils::CheckDofs)
    ;

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

    class_<NormalCalculationUtils > ("NormalCalculationUtils", init<>())
    .def("CalculateOnSimplex", CalcOnSimplex_Cond)
    .def("CalculateOnSimplex", CalcOnSimplex_ModelPart)
    .def("CalculateOnSimplex", CalcOnSimplexWithDoubleVar)
    .def("CalculateOnSimplex", CalcOnSimplexWithIntVar)
    .def("CalculateOnSimplex", CalcOnSimplexWithDoubleVarAlpha)
    .def("SwapNormals", &NormalCalculationUtils::SwapNormals)
//                    .def("CalculateOnSimplex", CalcOnSimplexWithArrayVar)
    ;

    class_<BodyNormalCalculationUtils > ("BodyNormalCalculationUtils", init<>())
    .def("CalculateBodyNormals", &BodyNormalCalculationUtils::CalculateBodyNormals)
    ;

    class_<BodyDistanceCalculationUtils > ("BodyDistanceCalculationUtils", init<>())
    .def("CalculateDistances2D", &BodyDistanceCalculationUtils::CalculateDistances < 2 >)
    .def("CalculateDistances3D", &BodyDistanceCalculationUtils::CalculateDistances < 3 >)
    ;

    class_<SignedDistanceCalculationUtils < 2 > >("SignedDistanceCalculationUtils2D", init<>())
    .def("CalculateDistances", &SignedDistanceCalculationUtils < 2 > ::CalculateDistances)
    .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 2 > ::FindMaximumEdgeSize)
    ;

    class_<SignedDistanceCalculationUtils < 3 > >("SignedDistanceCalculationUtils3D", init<>())
    .def("CalculateDistances", &SignedDistanceCalculationUtils < 3 > ::CalculateDistances)
    .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 3 > ::FindMaximumEdgeSize)
    ;

    class_<ParallelDistanceCalculator < 2 >, boost::noncopyable > ("ParallelDistanceCalculator2D", init<>())
    .def("CalculateDistances", &ParallelDistanceCalculator < 2 > ::CalculateDistances)
    .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 2 > ::CalculateInterfacePreservingDistances)
    .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 2 > ::CalculateDistancesLagrangianSurface)
    .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 2 > ::FindMaximumEdgeSize)
    ;

    class_<ParallelDistanceCalculator < 3 >, boost::noncopyable > ("ParallelDistanceCalculator3D", init<>())
    .def("CalculateDistances", &ParallelDistanceCalculator < 3 > ::CalculateDistances)
    .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 3 > ::CalculateInterfacePreservingDistances)
    .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 3 > ::CalculateDistancesLagrangianSurface)
    .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 3 > ::FindMaximumEdgeSize)
    ;

    class_<PointLocation > ("PointLocation", init<ModelPart& >())
    .def("Find", &PointLocation::Find)
    .def("Find2D", &PointLocation::Find2D)
    .def("Find3D", &PointLocation::Find3D)
    .def("found", &PointLocation::found)
    .def("ReturnDefaultPointData_scalar", &PointLocation::ReturnDefaultPointData_scalar)
    .def("ReturnDefaultPointData_vector", &PointLocation::ReturnDefaultPointData_vector)
    .def("ReturnCustomPointData_scalar", &PointLocation::ReturnCustomPointData_scalar)
    .def("ReturnCustomPointData_vector", &PointLocation::ReturnCustomPointData_vector)
    ;

    class_<ParticleConvectUtily<2> > ("ParticleConvectUtily2D", init< BinBasedFastPointLocator < 2 >::Pointer >())
    .def("MoveParticles_Substepping", &ParticleConvectUtily<2>::MoveParticles_Substepping)
    .def("MoveParticles_RK4", &ParticleConvectUtily<2>::MoveParticles_RK4)
    ;

    class_<ParticleConvectUtily<3> > ("ParticleConvectUtily3D", init< BinBasedFastPointLocator < 3 >::Pointer >())
    .def("MoveParticles_Substepping", &ParticleConvectUtily<3>::MoveParticles_Substepping)
    .def("MoveParticles_RK4", &ParticleConvectUtily<3>::MoveParticles_RK4)
    ;



    class_<IsosurfacePrinterApplication, boost::noncopyable >
    ("IsosurfacePrinterApplication",
     init<ModelPart& >() )
    .def("AddScalarVarIsosurface", &IsosurfacePrinterApplication::AddScalarVarIsosurface)
    .def("AddScalarVarIsosurfaceAndLower", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndLower)
    .def("AddScalarVarIsosurfaceAndHigher", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndHigher)
    .def("ClearData", &IsosurfacePrinterApplication::ClearData)
    .def("AddSkinConditions", &IsosurfacePrinterApplication::AddSkinConditions)
    .def("CreateNodesArray", &IsosurfacePrinterApplication::CreateNodesArray)
    ;


    // 	  class_<SignedDistanceCalculationBinBased<2> >("SignedDistanceCalculationBinBased2D", init<>())
    // 			  .def("CalculateDistances",&SignedDistanceCalculationBinBased<2>::CalculateDistances )
    //                           .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<2>::FindMaximumEdgeSize )
    // 			  ;
    //
    // 	  class_<SignedDistanceCalculationBinBased<3> >("SignedDistanceCalculationBinBased3D", init<>())
    // 			  .def("CalculateDistances",&SignedDistanceCalculationBinBased<3>::CalculateDistances )
    //                           .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<3>::FindMaximumEdgeSize )
    // 			  ;

    class_<DivideElemUtils > ("DivideElemUtils", init<>())
    .def("DivideElement_2D", &DivideElemUtils::DivideElement_2D)
    ;

    class_<Timer > ("Timer", init<>())
    .add_property("PrintOnScreen", &Timer::GetPrintOnScreen, &Timer::SetPrintOnScreen)
    .def("Start", &Timer::Start)
    .def("Stop", &Timer::Stop)
    .staticmethod("Start")
    .staticmethod("Stop")
    // 	    .def("PrintTimingInformation",Timer::PrintTimingInformation)
    .def(self_ns::str(self))
    ;




    class_<BoundingBoxUtilities > ("BoundingBoxUtilities", init<ModelPart&, const unsigned int& >())
    .def("Test", &BoundingBoxUtilities::Test)
    ;


    //           class_<SplitElements, boost::noncopyable >
    //                     ("SplitElements", init<ModelPart&, int >() )
    //                     .def("Split", &SplitElements::Split)
    //                     ;


    // 	  def("PrintTimingInformation",Timer::PrintTimingInformation);

    class_<OpenMPUtils > ("OpenMPUtils", init<>())
    .def("SetNumThreads", &OpenMPUtils::SetNumThreads)
    .staticmethod("SetNumThreads")
    .def("GetNumThreads", &OpenMPUtils::GetNumThreads)
    .staticmethod("GetNumThreads")
    .def("PrintOMPInfo", &OpenMPUtils::PrintOMPInfo)
    .staticmethod("PrintOMPInfo")
    ;

    class_< BinBasedFastPointLocator < 2 > > ("BinBasedFastPointLocator2D", init<ModelPart& >())
    .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabase)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
    .def("FindPointOnMesh", &BinBasedFastPointLocator < 2 > ::FindPointOnMeshSimplified)
    ;

    class_< BinBasedFastPointLocator < 3 > > ("BinBasedFastPointLocator3D", init<ModelPart&  >())
    .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabase)
    .def("FindPointOnMesh", &BinBasedFastPointLocator < 3 > ::FindPointOnMeshSimplified)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
    ;

    class_< BinBasedNodesInElementLocator < 2 > > ("BinBasedNodesInElementLocator2D", init<ModelPart& >())
    .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabase)
    .def("FindNodesInElement", &BinBasedNodesInElementLocator < 2 > ::FindNodesInElement)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
    ;

    class_< BinBasedNodesInElementLocator < 3 > > ("BinBasedNodesInElementLocator3D", init<ModelPart&  >())
    .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabase)
    .def("FindNodesInElement", &BinBasedNodesInElementLocator < 3 > ::FindNodesInElement)
    .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
    ;

    class_< ActivationUtilities > ("ActivationUtilities", init< >())
    .def("ActivateElementsAndConditions", &ActivationUtilities::ActivateElementsAndConditions)
    ;

    class_< GeometryTesterUtility, boost::noncopyable> ("GeometryTesterUtility", init< >())
    .def("RunTest", &GeometryTesterUtility::RunTest)
    .def("TestTriangle2D3N", &GeometryTesterUtility::TestTriangle2D3N)
    .def("TestTriangle2D6N", &GeometryTesterUtility::TestTriangle2D6N)
    .def("TestTetrahedra3D4N", &GeometryTesterUtility::TestTetrahedra3D4N)
    .def("TestTetrahedra3D10N", &GeometryTesterUtility::TestTetrahedra3D10N)
    .def("TestHexahedra3D8N", &GeometryTesterUtility::TestHexahedra3D8N)
    .def("TestHexahedra3D27N", &GeometryTesterUtility::TestHexahedra3D27N)
    .def("TestHexahedra3D20N", &GeometryTesterUtility::TestHexahedra3D20N)
    ;

    class_<CuttingUtility >("CuttingUtility", init< >())
    .def("GenerateCut", &CuttingUtility::GenerateCut)
    .def("UpdateCutData", &CuttingUtility ::UpdateCutData)
    .def("AddSkinConditions", &CuttingUtility ::AddSkinConditions)
    .def("FindSmallestEdge", &CuttingUtility ::FindSmallestEdge)
    ;

    class_<IntervalUtility >("IntervalUtility", init<Parameters >())
    .def("GetIntervalBegin", &IntervalUtility::GetIntervalBegin)
    .def("GetIntervalEnd", &IntervalUtility::GetIntervalEnd)
    .def("IsInInterval", &IntervalUtility ::IsInInterval)
    ;
    
    // Adding table from table stream to python
    class_<TableStreamUtility>("TableStreamUtility", init<>())
    .def(init< bool >())
    ;
    
    // Exact integration (for testing)
    class_<ExactMortarIntegrationUtility<2,2>>("ExactMortarIntegrationUtility2D2N", init<>())
    .def(init<const unsigned int>())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactIntegration)
    .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactAreaIntegration)
    ;
    class_<ExactMortarIntegrationUtility<3,3>>("ExactMortarIntegrationUtility3D3N", init<>())
    .def(init<const unsigned int>())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactIntegration)
    .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactAreaIntegration)
    ;
    class_<ExactMortarIntegrationUtility<3,4>>("ExactMortarIntegrationUtility3D4N", init<>())
    .def(init<const unsigned int>())
    .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactIntegration)
    .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactAreaIntegration)
    ;
}

} // namespace Python.

} // Namespace Kratos
