//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
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
#include "python/add_other_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "includes/global_pointer_variables.h"
#include "utilities/python_function_callback_utility.h"

//Other utilities
#include "utilities/convect_particles_utilities.h"
#include "utilities/condition_number_utility.h"
#include "utilities/mortar_utilities.h"
// #include "utilities/signed_distance_calculator_bin_based.h"  // Why is this commented?
#include "utilities/timer.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/merge_variable_lists_utility.h"
#include "utilities/variable_redistribution_utility.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/time_discretization.h"

namespace Kratos {
namespace Python {

// Auxiliar ModelPart Utility
void ModelPartRemoveElementAndBelongings1(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongings(ElementId, IdentifierFlag);
}
void ModelPartRemoveElementAndBelongings2(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongings(ElementId, IdentifierFlag, ThisIndex);
}
void ModelPartRemoveElementAndBelongings3(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongings(pThisElement, IdentifierFlag);
}

void ModelPartRemoveElementAndBelongings4(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongings(pThisElement, IdentifierFlag, ThisIndex);
}

void ModelPartRemoveElementAndBelongingsFromAllLevels1(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(ElementId, IdentifierFlag);
}

void ModelPartRemoveElementAndBelongingsFromAllLevels2(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(ElementId, IdentifierFlag, ThisIndex);
}

void ModelPartRemoveElementAndBelongingsFromAllLevels3(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(pThisElement, IdentifierFlag);
}

void ModelPartRemoveElementAndBelongingsFromAllLevels4(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(pThisElement, IdentifierFlag, ThisIndex);
}

void ModelPartRemoveConditionAndBelongings1(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(ConditionId, IdentifierFlag);
}

void ModelPartRemoveConditionAndBelongings2(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(ConditionId, IdentifierFlag, ThisIndex);
}

void ModelPartRemoveConditionAndBelongings3(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(pThisCondition, IdentifierFlag);
}

void ModelPartRemoveConditionAndBelongings4(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(pThisCondition, IdentifierFlag, ThisIndex);
}

void ModelPartRemoveConditionAndBelongingsFromAllLevels1(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(ConditionId, IdentifierFlag);
}

void ModelPartRemoveConditionAndBelongingsFromAllLevels2(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(ConditionId, IdentifierFlag, ThisIndex);
}

void ModelPartRemoveConditionAndBelongingsFromAllLevels3(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(pThisCondition, IdentifierFlag);
}

void ModelPartRemoveConditionAndBelongingsFromAllLevels4(AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag, ModelPart::IndexType ThisIndex)
{
    rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(pThisCondition, IdentifierFlag, ThisIndex);
}


// Timer
void PrintTimingInformation(Timer& rTimer)
{
    rTimer.PrintTimingInformation();
}


// Mortar
void ComputeNodesTangentModelPartWithSlipVariable(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rSlipVariable,
    const double SlipCoefficient,
    const bool SlipAlways
    )
{
    MortarUtilities::ComputeNodesTangentModelPart(rModelPart, &rSlipVariable, SlipCoefficient, SlipAlways);
}

void ComputeNodesTangentModelPartWithOutSlipVariable(
    ModelPart& rModelPart,
    const double SlipCoefficient,
    const bool SlipAlways
    )
{
    MortarUtilities::ComputeNodesTangentModelPart(rModelPart, NULL, SlipCoefficient, SlipAlways);
}

void ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlip(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rSlipVariable,
    const double SlipCoefficient
    )
{
    MortarUtilities::ComputeNodesTangentModelPart(rModelPart, &rSlipVariable, SlipCoefficient, false);
}

void ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlip(
    ModelPart& rModelPart,
    const double SlipCoefficient
    )
{
    MortarUtilities::ComputeNodesTangentModelPart(rModelPart, NULL, SlipCoefficient, false);
}

void ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlipUnitary(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rSlipVariable
    )
{
    MortarUtilities::ComputeNodesTangentModelPart(rModelPart, &rSlipVariable, 1.0, false);
}

void ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlipUnitary(ModelPart& rModelPart)
{
    MortarUtilities::ComputeNodesTangentModelPart(rModelPart, NULL, 1.0, false);
}



void AddOtherUtilitiesToPython(pybind11::module &m)  
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    // NOTE: this function is special in that it accepts a "pyObject" - this is the reason for which it is defined in this same file
    py::class_<PythonGenericFunctionUtility,  PythonGenericFunctionUtility::Pointer >(m,"PythonGenericFunctionUtility")
        .def(py::init<const std::string&>() )
        .def(py::init<const std::string&, Parameters>())
        .def("UseLocalSystem", &PythonGenericFunctionUtility::UseLocalSystem)
        .def("DependsOnSpace", &PythonGenericFunctionUtility::DependsOnSpace)
        .def("RotateAndCallFunction", &PythonGenericFunctionUtility::RotateAndCallFunction)
        .def("CallFunction", &PythonGenericFunctionUtility::CallFunction)
        ;

    py::class_<ApplyFunctionToNodesUtility >(m,"ApplyFunctionToNodesUtility")
        .def(py::init<ModelPart::NodesContainerType&, PythonGenericFunctionUtility::Pointer >() )
        .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction< Variable<double> >)
        .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
        .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >)
        .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >)
        .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >)
        .def("ReturnFunction", &ApplyFunctionToNodesUtility::ReturnFunction)
        ;
    
    // This is required to recognize the different overloads of ConditionNumberUtility::GetConditionNumber
    typedef double (ConditionNumberUtility::*InputGetConditionNumber)(SparseSpaceType::MatrixType&, LinearSolverType::Pointer, LinearSolverType::Pointer);
    typedef double (ConditionNumberUtility::*DirectGetConditionNumber)(SparseSpaceType::MatrixType&);

    InputGetConditionNumber ThisGetConditionNumber = &ConditionNumberUtility::GetConditionNumber;
    DirectGetConditionNumber ThisDirectGetConditionNumber = &ConditionNumberUtility::GetConditionNumber;

    py::class_<ConditionNumberUtility,ConditionNumberUtility::Pointer>(m,"ConditionNumberUtility")
        .def(py::init<>())
        .def(py::init<LinearSolverType::Pointer, LinearSolverType::Pointer>())
        .def("GetConditionNumber", ThisGetConditionNumber)
        .def("GetConditionNumber", ThisDirectGetConditionNumber)
        ;

    py::class_<ParticleConvectUtily<2> >(m,"ParticleConvectUtily2D")
        .def(py::init< BinBasedFastPointLocator < 2 >::Pointer >())
        .def("MoveParticles_Substepping", &ParticleConvectUtily<2>::MoveParticles_Substepping)
        .def("MoveParticles_RK4", &ParticleConvectUtily<2>::MoveParticles_RK4)
        ;

    py::class_<ParticleConvectUtily<3> >(m,"ParticleConvectUtily3D")
        .def(py::init< BinBasedFastPointLocator < 3 >::Pointer >())
        .def("MoveParticles_Substepping", &ParticleConvectUtily<3>::MoveParticles_Substepping)
        .def("MoveParticles_RK4", &ParticleConvectUtily<3>::MoveParticles_RK4)
        ;

    py::class_<Timer >(m,"Timer")
        .def(py::init<>())
        .def_property("PrintOnScreen", &Timer::GetPrintOnScreen, &Timer::SetPrintOnScreen)
        .def_property("PrintIntervalInformation", &Timer::GetPrintIntervalInformation, &Timer::SetPrintIntervalInformation)
        .def_static("Start", &Timer::Start)
        .def_static("Stop", &Timer::Stop)
        .def_static("GetTime", &Timer::GetTime)
        .def_static("SetOuputFile", &Timer::SetOuputFile)
        .def_static("CloseOuputFile", &Timer::CloseOuputFile)
        .def_static("GetPrintOnScreen", &Timer::GetPrintOnScreen)
        .def_static("SetPrintOnScreen", &Timer::SetPrintOnScreen)
        .def_static("GetPrintIntervalInformation", &Timer::GetPrintIntervalInformation)
        .def_static("SetPrintIntervalInformation", &Timer::SetPrintIntervalInformation)
        .def_static("PrintTimingInformation", PrintTimingInformation)
        .def("__str__", PrintObject<Timer>)
        ;

    // Exact integration (for testing)
    py::class_<ExactMortarIntegrationUtility<2,2>>(m,"ExactMortarIntegrationUtility2D2N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactAreaIntegration)
        ;

    py::class_<ExactMortarIntegrationUtility<3,3>>(m,"ExactMortarIntegrationUtility3D3N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactAreaIntegration)
        .def("TestGiDDebug",&ExactMortarIntegrationUtility<3,3>::TestGiDDebug)
        ;

    py::class_<ExactMortarIntegrationUtility<3,4>>(m,"ExactMortarIntegrationUtility3D4N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactAreaIntegration)
        .def("TestGiDDebug",&ExactMortarIntegrationUtility<3,4>::TestGiDDebug)
        ;

    py::class_<ExactMortarIntegrationUtility<3,3,false,4>>(m,"ExactMortarIntegrationUtility3D3N4N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3,false,4>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,3,false,4>::TestGetExactAreaIntegration)
        .def("TestGiDDebug",&ExactMortarIntegrationUtility<3,3,false,4>::TestGiDDebug)
        ;

    py::class_<ExactMortarIntegrationUtility<3,4,false,3>>(m,"ExactMortarIntegrationUtility3D4N3N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4,false,3>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,4,false,3>::TestGetExactAreaIntegration)
        .def("TestGiDDebug",&ExactMortarIntegrationUtility<3,4,false,3>::TestGiDDebug)
        ;


    // Mortar utilities
    auto mortar_utilities = m.def_submodule("MortarUtilities");
    mortar_utilities.def("ComputeNodesMeanNormalModelPart",&MortarUtilities::ComputeNodesMeanNormalModelPart);
    mortar_utilities.def("ComputeNodesTangentFromNormalModelPart",&MortarUtilities::ComputeNodesTangentFromNormalModelPart);
    mortar_utilities.def("ComputeNodesTangentModelPart",ComputeNodesTangentModelPartWithSlipVariable);
    mortar_utilities.def("ComputeNodesTangentModelPart",ComputeNodesTangentModelPartWithOutSlipVariable);
    mortar_utilities.def("ComputeNodesTangentModelPart",ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlip);
    mortar_utilities.def("ComputeNodesTangentModelPart",ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlip);
    mortar_utilities.def("ComputeNodesTangentModelPart",ComputeNodesTangentModelPartWithSlipVariableNotAlwaysSlipUnitary);
    mortar_utilities.def("ComputeNodesTangentModelPart",ComputeNodesTangentModelPartWithOutSlipVariableNotAlwaysSlipUnitary);
    mortar_utilities.def("InvertNormal",&MortarUtilities::InvertNormal<PointerVectorSet<Element, IndexedObject>>);
    mortar_utilities.def("InvertNormal",&MortarUtilities::InvertNormal<PointerVectorSet<Condition, IndexedObject>>);
    mortar_utilities.def("InvertNormal",&MortarUtilities::InvertNormalForFlag<PointerVectorSet<Element, IndexedObject>>);
    mortar_utilities.def("InvertNormal",&MortarUtilities::InvertNormalForFlag<PointerVectorSet<Condition, IndexedObject>>);


    // AssignUniqueModelPartCollectionTagUtility
    py::class_<AssignUniqueModelPartCollectionTagUtility, typename AssignUniqueModelPartCollectionTagUtility::Pointer>(m, "AssignUniqueModelPartCollectionTagUtility")
        .def(py::init<ModelPart&>())
        .def("DebugAssignUniqueModelPartCollectionTag",&AssignUniqueModelPartCollectionTagUtility::DebugAssignUniqueModelPartCollectionTag)
        .def("GetRecursiveSubModelPartNames",&AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames)
        .def("GetRecursiveSubModelPart",&AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart)
        ;

    py::class_<MergeVariableListsUtility, typename MergeVariableListsUtility::Pointer>(m, "MergeVariableListsUtility")
        .def(py::init<>())
        .def("Merge",&MergeVariableListsUtility::Merge)
        ;
    


    // VariableRedistributionUtility
    typedef void (*DistributePointDoubleType)(ModelPart&, const Variable< double >&, const Variable< double >&, double, unsigned int);
    typedef void (*DistributePointArrayType)(ModelPart&, const Variable< array_1d<double,3> >&, const Variable< array_1d<double,3> >&,double, unsigned int);

    DistributePointDoubleType DistributePointDouble = &VariableRedistributionUtility::DistributePointValues;
    DistributePointArrayType  DistributePointArray  = &VariableRedistributionUtility::DistributePointValues;

    typedef void (*ConvertDistributedDoubleType)(ModelPart&, const Variable< double >&, const Variable< double >&);
    typedef void (*ConvertDistributedArrayType)(ModelPart&, const Variable< array_1d<double,3> >&, const Variable< array_1d<double,3> >&);

    ConvertDistributedDoubleType ConvertDistributedDouble = &VariableRedistributionUtility::ConvertDistributedValuesToPoint;
    ConvertDistributedArrayType  ConvertDistributedArray  = &VariableRedistributionUtility::ConvertDistributedValuesToPoint;

    // Note: The StaticMethod thing should be done only once for each set of overloads
    py::class_< VariableRedistributionUtility >(m,"VariableRedistributionUtility")
        .def_static("DistributePointValues",DistributePointDouble)
        .def_static("DistributePointValues",DistributePointArray)
        .def_static("ConvertDistributedValuesToPoint",ConvertDistributedDouble)
        .def_static("ConvertDistributedValuesToPoint",ConvertDistributedArray)
        ;

    // Auxiliar ModelPart Utility

    py::class_<AuxiliarModelPartUtilities, typename AuxiliarModelPartUtilities::Pointer>(m, "AuxiliarModelPartUtilities")
        .def(py::init<ModelPart&>())
        .def("RemoveElementAndBelongings", ModelPartRemoveElementAndBelongings1)
        .def("RemoveElementAndBelongings", ModelPartRemoveElementAndBelongings2)
        .def("RemoveElementAndBelongings", ModelPartRemoveElementAndBelongings3)
        .def("RemoveElementAndBelongings", ModelPartRemoveElementAndBelongings4)
        .def("RemoveElementsAndBelongings", &Kratos::AuxiliarModelPartUtilities::RemoveElementsAndBelongings)
        .def("RemoveElementAndBelongingsFromAllLevels", ModelPartRemoveElementAndBelongingsFromAllLevels1)
        .def("RemoveElementAndBelongingsFromAllLevels", ModelPartRemoveElementAndBelongingsFromAllLevels2)
        .def("RemoveElementAndBelongingsFromAllLevels", ModelPartRemoveElementAndBelongingsFromAllLevels3)
        .def("RemoveElementAndBelongingsFromAllLevels", ModelPartRemoveElementAndBelongingsFromAllLevels4)
        .def("RemoveElementsAndBelongingsFromAllLevels", &Kratos::AuxiliarModelPartUtilities::RemoveElementsAndBelongingsFromAllLevels)
        .def("RemoveConditionAndBelongings", ModelPartRemoveConditionAndBelongings1)
        .def("RemoveConditionAndBelongings", ModelPartRemoveConditionAndBelongings2)
        .def("RemoveConditionAndBelongings", ModelPartRemoveConditionAndBelongings3)
        .def("RemoveConditionAndBelongings", ModelPartRemoveConditionAndBelongings4)
        .def("RemoveConditionsAndBelongings", &Kratos::AuxiliarModelPartUtilities::RemoveConditionsAndBelongings)
        .def("RemoveConditionAndBelongingsFromAllLevels", ModelPartRemoveConditionAndBelongingsFromAllLevels1)
        .def("RemoveConditionAndBelongingsFromAllLevels", ModelPartRemoveConditionAndBelongingsFromAllLevels2)
        .def("RemoveConditionAndBelongingsFromAllLevels", ModelPartRemoveConditionAndBelongingsFromAllLevels3)
        .def("RemoveConditionAndBelongingsFromAllLevels", ModelPartRemoveConditionAndBelongingsFromAllLevels4)
        .def("RemoveConditionsAndBelongingsFromAllLevels", &Kratos::AuxiliarModelPartUtilities::RemoveConditionsAndBelongingsFromAllLevels)
        ;

    // Sparse matrix multiplication utility
    py::class_<SparseMatrixMultiplicationUtility, typename SparseMatrixMultiplicationUtility::Pointer>(m, "SparseMatrixMultiplicationUtility")
        .def(py::init<>())
        .def_static("MatrixMultiplication",&SparseMatrixMultiplicationUtility::MatrixMultiplication<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def_static("MatrixMultiplicationSaad",&SparseMatrixMultiplicationUtility::MatrixMultiplicationSaad<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def_static("MatrixMultiplicationRMerge",&SparseMatrixMultiplicationUtility::MatrixMultiplicationRMerge<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def_static("MatrixAdd",&SparseMatrixMultiplicationUtility::MatrixAdd<CompressedMatrix, CompressedMatrix>)
        .def_static("TransposeMatrix",&SparseMatrixMultiplicationUtility::TransposeMatrix<CompressedMatrix, CompressedMatrix>)
        ;


    // TimeDiscretization
    auto mod_time_discretization = m.def_submodule("TimeDiscretization");

    py::class_<TimeDiscretization::BDF>(mod_time_discretization, "BDF")
        .def(py::init<const unsigned int>())
        .def("GetTimeOrder", (unsigned int (TimeDiscretization::BDF::*)() const) & TimeDiscretization::BDF::GetTimeOrder)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF::*)(double) const) & TimeDiscretization::BDF::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF::*)(double, double) const) & TimeDiscretization::BDF::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF::*)(const ProcessInfo &) const) & TimeDiscretization::BDF::ComputeBDFCoefficients)
        .def("ComputeAndSaveBDFCoefficients", (void (TimeDiscretization::BDF::*)(ProcessInfo &) const) & TimeDiscretization::BDF::ComputeAndSaveBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF1>(mod_time_discretization, "BDF1")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF1::*)(double) const) &TimeDiscretization::BDF1::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF1::*)(const ProcessInfo&) const) &TimeDiscretization::BDF1::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF2>(mod_time_discretization, "BDF2")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF2::*)(double, double) const) &TimeDiscretization::BDF2::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF2::*)(const ProcessInfo&) const) &TimeDiscretization::BDF2::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF3>(mod_time_discretization, "BDF3")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF3::*)(double) const) &TimeDiscretization::BDF3::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF3::*)(const ProcessInfo&) const) &TimeDiscretization::BDF3::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF4>(mod_time_discretization, "BDF4")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF4::*)(double) const) &TimeDiscretization::BDF4::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF4::*)(const ProcessInfo&) const) &TimeDiscretization::BDF4::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF5>(mod_time_discretization, "BDF5")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF5::*)(double) const) &TimeDiscretization::BDF5::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF5::*)(const ProcessInfo&) const) &TimeDiscretization::BDF5::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF6>(mod_time_discretization, "BDF6")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF6::*)(double) const) &TimeDiscretization::BDF6::ComputeBDFCoefficients)
        .def("ComputeBDFCoefficients", (std::vector<double> (TimeDiscretization::BDF6::*)(const ProcessInfo&) const) &TimeDiscretization::BDF6::ComputeBDFCoefficients)
        ;

    py::class_<TimeDiscretization::Newmark>(mod_time_discretization, "Newmark")
        .def(py::init<>())
        .def(py::init<const double, const double>())
        .def("GetBeta", &TimeDiscretization::Newmark::GetBeta)
        .def("GetGamma", &TimeDiscretization::Newmark::GetGamma)
        ;
    py::class_<TimeDiscretization::Bossak>(mod_time_discretization, "Bossak")
        .def(py::init<>())
        .def(py::init<const double>())
        .def(py::init<const double, const double, const double>())
        .def("GetBeta", &TimeDiscretization::Bossak::GetBeta)
        .def("GetGamma", &TimeDiscretization::Bossak::GetGamma)
        .def("GetAlphaM", &TimeDiscretization::Bossak::GetAlphaM)
        ;
    py::class_<TimeDiscretization::GeneralizedAlpha>(mod_time_discretization, "GeneralizedAlpha")
        .def(py::init<>())
        .def(py::init<const double, const double>())
        .def(py::init<const double, const double, const double, const double>())
        .def("GetBeta", &TimeDiscretization::GeneralizedAlpha::GetBeta)
        .def("GetGamma", &TimeDiscretization::GeneralizedAlpha::GetGamma)
        .def("GetAlphaM", &TimeDiscretization::GeneralizedAlpha::GetAlphaM)
        .def("GetAlphaF", &TimeDiscretization::GeneralizedAlpha::GetAlphaF)
        ;

    std::size_t (*GetMinimumBufferSizeBDF)(const TimeDiscretization::BDF&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF1)(const TimeDiscretization::BDF1&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF2)(const TimeDiscretization::BDF2&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF3)(const TimeDiscretization::BDF3&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF4)(const TimeDiscretization::BDF4&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF5)(const TimeDiscretization::BDF5&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF6)(const TimeDiscretization::BDF6&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeNewmark)(const TimeDiscretization::Newmark&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBossak)(const TimeDiscretization::Bossak&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeGeneralizedAlpha)(const TimeDiscretization::GeneralizedAlpha&) = &TimeDiscretization::GetMinimumBufferSize;

    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF1 );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF2 );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF3 );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF4 );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF5 );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF6 );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeNewmark );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeBossak );
    mod_time_discretization.def("GetMinimumBufferSize", GetMinimumBufferSizeGeneralizedAlpha );



}

} // namespace Python.
} // Namespace Kratos
