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
#include <pybind11/stl.h>

// External includes

// Project includes
#include "python/add_other_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "includes/global_pointer_variables.h"

//Other utilities
#include "utilities/python_function_callback_utility.h"
#include "utilities/condition_number_utility.h"
#include "utilities/mortar_utilities.h"
#include "utilities/deflation_utils.h"
#include "utilities/timer.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/merge_variable_lists_utility.h"
#include "utilities/variable_redistribution_utility.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/time_discretization.h"
#include "utilities/table_stream_utility.h"
#include "utilities/read_materials_utility.h"
#include "utilities/activation_utilities.h"
#include "utilities/sensitivity_builder.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/constraint_utilities.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "utilities/properties_utilities.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/file_name_data_collector.h"
#include "utilities/sensitivity_utilities.h"

namespace Kratos {
namespace Python {

/**
 * @brief A thin wrapper for GetSortedListOfFileNameData. The reason for having the wrapper is to replace the original lambda implementation as it causes gcc 4.8 to generate bad code on Centos7 which leads to memory corruption.
 */
pybind11::list GetSortedListOfFileNameDataHelper(
    std::vector<FileNameDataCollector::FileNameData>& rFileNameDataList,
    const std::vector<std::string> & rSortingFlagsOrder
    )
{
    FileNameDataCollector::SortListOfFileNameData(rFileNameDataList, rSortingFlagsOrder);
    pybind11::list result;
    for (unsigned int j = 0; j < rFileNameDataList.size(); j++)
    {
        result.append(rFileNameDataList[j]);
    }
    return result;
}

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

//timer
void PrintTimingInformation(Timer& rTimer)
{
    rTimer.PrintTimingInformation();
}

//mortar
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


//compare elements and conditions utilities
std::string GetRegisteredNameElement(const Element& rElement)
{
    std::string name;
    CompareElementsAndConditionsUtility::GetRegisteredName(rElement, name);
    return name;
}

std::string GetRegisteredNameCondition(const Condition& rCondition)
{
    std::string name;
    CompareElementsAndConditionsUtility::GetRegisteredName(rCondition, name);
    return name;
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
        .def("FunctionBody", &PythonGenericFunctionUtility::FunctionBody)
        .def("RotateAndCallFunction", &PythonGenericFunctionUtility::RotateAndCallFunction)
        .def("CallFunction", &PythonGenericFunctionUtility::CallFunction)
        ;

    py::class_<ApplyFunctionToNodesUtility >(m,"ApplyFunctionToNodesUtility")
        .def(py::init<ModelPart::NodesContainerType&, PythonGenericFunctionUtility::Pointer >() )
        .def("ApplyFunction", &ApplyFunctionToNodesUtility::ApplyFunction< Variable<double> >)
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

    //deflation utilities
    py::class_<DeflationUtils>(m,"DeflationUtils")
        .def(py::init<>())
        .def("VisualizeAggregates",&DeflationUtils::VisualizeAggregates)
        ;

    //timer
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
        .def(py::init<const std::size_t, const double, const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double, const bool>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<2,2>::TestGetExactAreaIntegration)
        ;

    py::class_<ExactMortarIntegrationUtility<3,3>>(m,"ExactMortarIntegrationUtility3D3N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double, const bool>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,3>::TestGetExactAreaIntegration)
        .def("TestIODebug",&ExactMortarIntegrationUtility<3,3>::TestIODebug)
        ;

    py::class_<ExactMortarIntegrationUtility<3,4>>(m,"ExactMortarIntegrationUtility3D4N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double, const bool>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,4>::TestGetExactAreaIntegration)
        .def("TestIODebug",&ExactMortarIntegrationUtility<3,4>::TestIODebug)
        ;

    py::class_<ExactMortarIntegrationUtility<3,3,false,4>>(m,"ExactMortarIntegrationUtility3D3N4N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double, const bool>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,3,false,4>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,3,false,4>::TestGetExactAreaIntegration)
        .def("TestIODebug",&ExactMortarIntegrationUtility<3,3,false,4>::TestIODebug)
        ;

    py::class_<ExactMortarIntegrationUtility<3,4,false,3>>(m,"ExactMortarIntegrationUtility3D4N3N")
        .def(py::init<>())
        .def(py::init<const std::size_t>())
        .def(py::init<const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double>())
        .def(py::init<const std::size_t, const double, const std::size_t, const double, const bool>())
        .def("TestGetExactIntegration",&ExactMortarIntegrationUtility<3,4,false,3>::TestGetExactIntegration)
        .def("TestGetExactAreaIntegration",&ExactMortarIntegrationUtility<3,4,false,3>::TestGetExactAreaIntegration)
        .def("TestIODebug",&ExactMortarIntegrationUtility<3,4,false,3>::TestIODebug)
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

    // Merge variable list utility
    py::class_<MergeVariableListsUtility, typename MergeVariableListsUtility::Pointer>(m, "MergeVariableListsUtility")
        .def(py::init<>())
        .def("Merge",&MergeVariableListsUtility::Merge)
        ;

    // VariableRedistributionUtility
    py::class_< VariableRedistributionUtility >(m,"VariableRedistributionUtility")
        .def_static("ConvertDistributedValuesToPoint", [](ModelPart& rModelPart, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPoint(rModelPart, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPoint", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPoint(rModelPart, rElements, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPoint", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPoint(rModelPart, rConditions, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPoint", [](ModelPart& rModelPart, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPoint(rModelPart, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPoint", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPoint(rModelPart, rElements, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPoint", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPoint(rModelPart, rConditions, rPointVariable, rDistributedVariable);})
        .def_static("DistributePointValues", [](ModelPart& rModelPart, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValues(rModelPart, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValues", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValues(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValues", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValues(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValues", [](ModelPart& rModelPart, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValues(rModelPart, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValues", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValues(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValues", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValues(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("ConvertDistributedValuesToPointNonHistorical", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(rModelPart, rElements, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPointNonHistorical", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(rModelPart, rConditions, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPointNonHistorical", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(rModelPart, rElements, rPointVariable, rDistributedVariable);})
        .def_static("ConvertDistributedValuesToPointNonHistorical", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable){VariableRedistributionUtility::ConvertDistributedValuesToPointNonHistorical(rModelPart, rConditions, rPointVariable, rDistributedVariable);})
        .def_static("DistributePointValuesNonHistorical", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValuesNonHistorical(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValuesNonHistorical", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<double>& rPointVariable, const Variable<double>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValuesNonHistorical(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValuesNonHistorical", [](ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValuesNonHistorical(rModelPart, rElements, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        .def_static("DistributePointValuesNonHistorical", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType& rConditions, const Variable<array_1d<double,3>>& rPointVariable, const Variable<array_1d<double,3>>& rDistributedVariable, double Tolerance, double MaximumIterations){VariableRedistributionUtility::DistributePointValuesNonHistorical(rModelPart, rConditions, rPointVariable, rDistributedVariable, Tolerance, MaximumIterations);})
        ;

    // Auxiliar ModelPart Utility
    py::class_<AuxiliarModelPartUtilities, typename AuxiliarModelPartUtilities::Pointer>(m, "AuxiliarModelPartUtilities")
        .def(py::init<ModelPart&>())
        .def("RecursiveEnsureModelPartOwnsProperties", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities) { rAuxiliarModelPartUtilities.RecursiveEnsureModelPartOwnsProperties();})
        .def("RecursiveEnsureModelPartOwnsProperties", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, const bool RemovePreviousProperties) { rAuxiliarModelPartUtilities.RecursiveEnsureModelPartOwnsProperties(RemovePreviousProperties);})
        .def("EnsureModelPartOwnsProperties", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities) { rAuxiliarModelPartUtilities.EnsureModelPartOwnsProperties();})
        .def("EnsureModelPartOwnsProperties", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, const bool RemovePreviousProperties) { rAuxiliarModelPartUtilities.EnsureModelPartOwnsProperties(RemovePreviousProperties);})
        .def("RemoveElementAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveElementAndBelongings(ElementId, IdentifierFlag);})
        .def("RemoveElementAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveElementAndBelongings(ElementId, IdentifierFlag, ThisIndex);})
        .def("RemoveElementAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveElementAndBelongings(pThisElement, IdentifierFlag);})
        .def("RemoveElementAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveElementAndBelongings(pThisElement, IdentifierFlag, ThisIndex);})
        .def("RemoveElementsAndBelongings", &Kratos::AuxiliarModelPartUtilities::RemoveElementsAndBelongings)
        .def("RemoveElementAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(ElementId, IdentifierFlag);})
        .def("RemoveElementAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ElementId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(ElementId, IdentifierFlag, ThisIndex);})
        .def("RemoveElementAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(pThisElement, IdentifierFlag);})
        .def("RemoveElementAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ElementType::Pointer pThisElement, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveElementAndBelongingsFromAllLevels(pThisElement, IdentifierFlag, ThisIndex);})
        .def("RemoveElementsAndBelongingsFromAllLevels", &Kratos::AuxiliarModelPartUtilities::RemoveElementsAndBelongingsFromAllLevels)
        .def("RemoveConditionAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(ConditionId, IdentifierFlag);})
        .def("RemoveConditionAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(ConditionId, IdentifierFlag, ThisIndex);})
        .def("RemoveConditionAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(pThisCondition, IdentifierFlag);})
        .def("RemoveConditionAndBelongings", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongings(pThisCondition, IdentifierFlag, ThisIndex);})
        .def("RemoveConditionsAndBelongings", &Kratos::AuxiliarModelPartUtilities::RemoveConditionsAndBelongings)
        .def("RemoveConditionAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(ConditionId, IdentifierFlag);})
        .def("RemoveConditionAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::IndexType ConditionId, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(ConditionId, IdentifierFlag, ThisIndex);})
        .def("RemoveConditionAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(pThisCondition, IdentifierFlag);})
        .def("RemoveConditionAndBelongingsFromAllLevels", [](AuxiliarModelPartUtilities& rAuxiliarModelPartUtilities, ModelPart::ConditionType::Pointer pThisCondition, Flags IdentifierFlag, ModelPart::IndexType ThisIndex) { rAuxiliarModelPartUtilities.RemoveConditionAndBelongingsFromAllLevels(pThisCondition, IdentifierFlag, ThisIndex);})
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


    // Adding table from table stream to python
    py::class_<TableStreamUtility, typename TableStreamUtility::Pointer>(m,"TableStreamUtility")
        .def(py::init<>())
        .def(py::init< bool >())
        .def("SetOnProcessInfo",SetOnProcessInfo)
        ;

    // Read materials utility
    py::class_<ReadMaterialsUtility, typename ReadMaterialsUtility::Pointer>(m, "ReadMaterialsUtility")
    .def(py::init<Model&>())
    .def(py::init<Parameters, Model&>())
    .def("ReadMaterials",&ReadMaterialsUtility::ReadMaterials)
    ;

    //activation utilities
    py::class_< ActivationUtilities >(m,"ActivationUtilities")
        .def(py::init< >())
        .def("ActivateElementsAndConditions", &ActivationUtilities::ActivateElementsAndConditions)
        ;

    //sensitivity builder
    py::class_<SensitivityBuilder>(m, "SensitivityBuilder")
        .def(py::init<Parameters, ModelPart&, AdjointResponseFunction::Pointer>())
        .def(py::init<Parameters, ModelPart&, AdjointResponseFunction::Pointer, SensitivityBuilderScheme::Pointer>())
        .def("Initialize", &SensitivityBuilder::Initialize)
        .def("InitializeSolutionStep", &SensitivityBuilder::InitializeSolutionStep)
        .def("UpdateSensitivities", &SensitivityBuilder::UpdateSensitivities)
        .def("FinalizeSolutionStep", &SensitivityBuilder::FinalizeSolutionStep)
        .def("Finalize", &SensitivityBuilder::Finalize)
        ;

    //Sensitivity utilities
    py::class_<SensitivityUtilities>(m,"SensitivityUtilities")
        .def_static("AssignConditionDerivativesToNodes", &SensitivityUtilities::AssignEntityDerivativesToNodes<ModelPart::ConditionsContainerType>)
        .def_static("AssignElementDerivativesToNodes", &SensitivityUtilities::AssignEntityDerivativesToNodes<ModelPart::ElementsContainerType>)
        ;

    //OpenMP utilities
    py::class_<OpenMPUtils >(m,"OpenMPUtils")
        .def(py::init<>())
        .def_static("SetNumThreads", [](const int NumThreads){
            KRATOS_WARNING("OpenMPUtils") << "\"SetNumThreads\" is deprecated, please use ParallelUtilities.SetNumThreads instead" << std::endl;
            ParallelUtilities::SetNumThreads(NumThreads);})
        .def_static("GetNumThreads", [](){
            KRATOS_WARNING("OpenMPUtils") << "\"GetNumThreads\" is deprecated, please use ParallelUtilities.GetNumThreads instead" << std::endl;
            return ParallelUtilities::GetNumThreads();})
        .def_static("GetNumberOfProcessors", [](){
            KRATOS_WARNING("OpenMPUtils") << "\"GetNumberOfProcessors\" is deprecated, please use ParallelUtilities.GetNumProcs instead" << std::endl;
            return ParallelUtilities::GetNumProcs();})
        .def_static("PrintOMPInfo", &OpenMPUtils::PrintOMPInfo)
        ;

    // ParallelUtilities
    py::class_<ParallelUtilities >(m,"ParallelUtilities")
        .def_static("GetNumThreads", &ParallelUtilities::GetNumThreads)
        .def_static("SetNumThreads", &ParallelUtilities::SetNumThreads)
        .def_static("GetNumProcs",   &ParallelUtilities::GetNumProcs)
        ;

    // EntitiesUtilities
    auto entities_utilities = m.def_submodule("EntitiesUtilities");
    entities_utilities.def("InitializeAllEntities", &EntitiesUtilities::InitializeAllEntities );
    entities_utilities.def("InitializeConditions", &EntitiesUtilities::InitializeEntities<Condition> );
    entities_utilities.def("InitializeElements", &EntitiesUtilities::InitializeEntities<Element> );
    entities_utilities.def("InitializeMasterSlaveConstraints", &EntitiesUtilities::InitializeEntities<MasterSlaveConstraint> );

    // ConstraintUtilities
    auto constraint_utilities = m.def_submodule("ConstraintUtilities");
    constraint_utilities.def("ResetSlaveDofs", &ConstraintUtilities::ResetSlaveDofs );
    constraint_utilities.def("ApplyConstraints", &ConstraintUtilities::ApplyConstraints );

    // Compare elements and conditions utility
    auto mod_compare_elem_cond_utils = m.def_submodule("CompareElementsAndConditionsUtility");
    mod_compare_elem_cond_utils.def("GetRegisteredName", GetRegisteredNameElement );
    mod_compare_elem_cond_utils.def("GetRegisteredName", GetRegisteredNameCondition );

    // PropertiesUtilities
    auto mod_prop_utils = m.def_submodule("PropertiesUtilities");
    mod_prop_utils.def("CopyPropertiesValues", &PropertiesUtilities::CopyPropertiesValues);

    // coordinate transformation utilities
    typedef CoordinateTransformationUtils<LocalSpaceType::MatrixType, LocalSpaceType::VectorType, double> CoordinateTransformationUtilsType;
    py::class_<
        CoordinateTransformationUtilsType,
        CoordinateTransformationUtilsType::Pointer>
        (m,"CoordinateTransformationUtils")
        .def(py::init<const unsigned int, const unsigned int, const Kratos::Flags&>())
        .def("Rotate", (void(CoordinateTransformationUtilsType::*)(LocalSpaceType::MatrixType&, LocalSpaceType::VectorType&, ModelPart::GeometryType&)const)(&CoordinateTransformationUtilsType::Rotate))
        .def("Rotate", (void(CoordinateTransformationUtilsType::*)(LocalSpaceType::VectorType&, ModelPart::GeometryType&)const)(&CoordinateTransformationUtilsType::Rotate))
        .def("ApplySlipCondition", (void(CoordinateTransformationUtilsType::*)(LocalSpaceType::MatrixType&, LocalSpaceType::VectorType&, ModelPart::GeometryType&)const)(&CoordinateTransformationUtilsType::ApplySlipCondition))
        .def("ApplySlipCondition", (void(CoordinateTransformationUtilsType::*)(LocalSpaceType::VectorType&, ModelPart::GeometryType&)const)(&CoordinateTransformationUtilsType::ApplySlipCondition))
        .def("RotateVelocities", &CoordinateTransformationUtilsType::RotateVelocities)
        .def("RecoverVelocities", &CoordinateTransformationUtilsType::RecoverVelocities)
        .def("CalculateRotationOperatorPure", (void(CoordinateTransformationUtilsType::*)(LocalSpaceType::MatrixType&, const ModelPart::GeometryType::PointType&)const)(&CoordinateTransformationUtilsType::CalculateRotationOperatorPure))
        .def("CalculateRotationOperatorPureShapeSensitivities", (void(CoordinateTransformationUtilsType::*)(LocalSpaceType::MatrixType&, const std::size_t, const std::size_t, const ModelPart::GeometryType::PointType&)const)(&CoordinateTransformationUtilsType::CalculateRotationOperatorPureShapeSensitivities))
        ;

    // add FileNameDataCollector
    auto file_name_data_collector = py::class_<
        FileNameDataCollector,
        FileNameDataCollector::Pointer>
        (m, "FileNameDataCollector")
        .def(py::init<const ModelPart&, const std::string&, const std::unordered_map<std::string, std::string>&>())
        .def("GetFileName", &FileNameDataCollector::GetFileName)
        .def("GetPath", &FileNameDataCollector::GetPath)
        .def("GetSortedFileNamesList", &FileNameDataCollector::GetSortedFileNamesList)
        .def("RetrieveFileNameData", &FileNameDataCollector::RetrieveFileNameData)
        .def("GetFileNameDataList", &FileNameDataCollector::GetFileNameDataList)
        .def_static("ExtractFileNamePattern", &FileNameDataCollector::ExtractFileNamePattern)
        .def_static("GetSortedListOfFileNameData", &GetSortedListOfFileNameDataHelper)
        ;

    // add FileNameData holder
    py::class_<
        FileNameDataCollector::FileNameData,
        FileNameDataCollector::FileNameData::Pointer>
        (file_name_data_collector, "FileNameData")
        .def(py::init<>())
        .def(py::init<const std::string&, int, int, double>())
        .def("SetFileName", &FileNameDataCollector::FileNameData::SetFileName)
        .def("GetFileName", &FileNameDataCollector::FileNameData::GetFileName)
        .def("SetRank", &FileNameDataCollector::FileNameData::SetRank)
        .def("GetRank", &FileNameDataCollector::FileNameData::GetRank)
        .def("SetStep", &FileNameDataCollector::FileNameData::SetStep)
        .def("GetStep", &FileNameDataCollector::FileNameData::GetStep)
        .def("SetTime", &FileNameDataCollector::FileNameData::SetTime)
        .def("GetTime", &FileNameDataCollector::FileNameData::GetTime)
        .def("Clear", &FileNameDataCollector::FileNameData::Clear)
        .def("__eq__", &FileNameDataCollector::FileNameData::operator==)
        ;
}

} // namespace Python.
} // Namespace Kratos
