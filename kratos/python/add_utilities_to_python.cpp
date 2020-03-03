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
#include "includes/global_pointer_variables.h"

// #include "utilities/signed_distance_calculator_bin_based.h"
#include "utilities/timer.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_fast_point_locator_conditions.h"
#include "utilities/binbased_nodes_in_element_locator.h"
#include "utilities/embedded_skin_utility.h"
#include "utilities/geometry_tester.h"
#include "utilities/cutting_utility.h"
#include "utilities/python_function_callback_utility.h"
#include "utilities/interval_utility.h"
#include "utilities/table_stream_utility.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/merge_variable_lists_utility.h"
#include "utilities/variable_redistribution_utility.h"
#include "utilities/sensitivity_builder.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/time_discretization.h"
#include "utilities/geometrical_transformation_utilities.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos {
namespace Python {
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

// Embedded skin utility auxiliar functions
template<std::size_t TDim>
void InterpolateMeshVariableToSkinDouble(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<double> &rVariable,
    const Variable<double> &rEmbeddedVariable)
{
    rEmbeddedSkinUtility.InterpolateMeshVariableToSkin(rVariable, rEmbeddedVariable);
}

template<std::size_t TDim>
void InterpolateMeshVariableToSkinArray(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<array_1d<double,3>> &rVariable,
    const Variable<array_1d<double,3>> &rEmbeddedVariable)
{
    rEmbeddedSkinUtility.InterpolateMeshVariableToSkin(rVariable, rEmbeddedVariable);
}

template<std::size_t TDim>
void InterpolateDiscontinuousMeshVariableToSkinDouble(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<double> &rVariable,
    const Variable<double> &rEmbeddedVariable,
    const std::string &rInterfaceSide)
{
    rEmbeddedSkinUtility.InterpolateDiscontinuousMeshVariableToSkin(rVariable, rEmbeddedVariable, rInterfaceSide);
}

template<std::size_t TDim>
void InterpolateDiscontinuousMeshVariableToSkinArray(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<array_1d<double,3>> &rVariable,
    const Variable<array_1d<double,3>> &rEmbeddedVariable,
    const std::string &rInterfaceSide)
{
    rEmbeddedSkinUtility.InterpolateDiscontinuousMeshVariableToSkin(rVariable, rEmbeddedVariable, rInterfaceSide);
}

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

void CalculateDistancesDefault2D(ParallelDistanceCalculator<2>& rParallelDistanceCalculator,ModelPart& rModelPart, const Variable<double>& rDistanceVar, const Variable<double>& rAreaVar, const unsigned int max_levels, const double max_distance)
{
    rParallelDistanceCalculator.CalculateDistances(rModelPart, rDistanceVar, rAreaVar, max_levels, max_distance);
}

void CalculateDistancesFlag2D(ParallelDistanceCalculator<2>& rParallelDistanceCalculator, ModelPart& rModelPart, const Variable<double>& rDistanceVar, const Variable<double>& rAreaVar, const unsigned int max_levels, const double max_distance, Flags Options)
{
    rParallelDistanceCalculator.CalculateDistances(rModelPart, rDistanceVar, rAreaVar, max_levels, max_distance, Options);
}

void CalculateDistancesDefault3D(ParallelDistanceCalculator<3>& rParallelDistanceCalculator,ModelPart& rModelPart, const Variable<double>& rDistanceVar, const Variable<double>& rAreaVar, const unsigned int max_levels, const double max_distance)
{
    rParallelDistanceCalculator.CalculateDistances(rModelPart, rDistanceVar, rAreaVar, max_levels, max_distance);
}

void CalculateDistancesFlag3D(ParallelDistanceCalculator<3>& rParallelDistanceCalculator, ModelPart& rModelPart, const Variable<double>& rDistanceVar, const Variable<double>& rAreaVar, const unsigned int max_levels, const double max_distance, Flags Options)
{
    rParallelDistanceCalculator.CalculateDistances(rModelPart, rDistanceVar, rAreaVar, max_levels, max_distance, Options);
}

void VariableUtilsUpdateCurrentPosition(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes);
}

void VariableUtilsUpdateCurrentPositionWithVariable(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes,
    const VariableUtils::ArrayVarType &rUpdateVariable
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes, rUpdateVariable);
}

void VariableUtilsUpdateCurrentPositionWithVariableAndPosition(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes,
    const VariableUtils::ArrayVarType &rUpdateVariable,
    const IndexType BufferPosition
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes, rUpdateVariable, BufferPosition);
}

template<class TVarType>
void VariableUtilsCopyModelPartNodalVar(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void VariableUtilsCopyModelPartNodalVarWithDestination(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TVarType &rDestinationVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void CopyModelPartNodalVarToNonHistoricalVar(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void CopyModelPartNodalVarToNonHistoricalVarWithDestination(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TVarType &rDestinationVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

/**
 * @brief Auxiliary set variable export function
 * This function is required to export the SetVariable overloaded method with a unique name
 * @tparam TDataType The variable data type
 * @tparam Variable<TDataType> The variable type
 * @param rVariableUtils Reference to the self variable utils class
 * @param rVariable Reference to the variable to be set
 * @param rValue Reference to the value to set
 * @param rNodes Reference to the nodes container
 */
template<class TDataType, class TVarType = Variable<TDataType>>
void VariableUtilsSetVariable(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    NodesContainerType &rNodes)
{
    rVariableUtils.SetVariable(rVariable, rValue, rNodes);
}


/**
 * @brief Auxiliary set variable export function
 * This function is required to export the SetVariable with flag overloaded method with a unique name
 * @tparam TDataType The variable data type
 * @tparam Variable<TDataType> The variable type
 * @param rVariableUtils Reference to the self variable utils class
 * @param rVariable Reference to the variable to be set
 * @param rValue Reference to the value to set
 * @param rNodes Reference to the nodes container
 * @param Flag Flag to filter the nodes that are set
 * @param CheckValue Flag value to be checked
 */
template <class TDataType, class TVarType = Variable<TDataType>>
void VariableUtilsSetVariableForFlag(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    NodesContainerType &rNodes,
    const Flags Flag,
    const bool CheckValue = true)
{
    rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);
}

template <class TDataType, class TContainerType, class TVarType = Variable<TDataType>>
void VariableUtilsSetNonHistoricalVariable(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    TContainerType &rContainer)
{
    rVariableUtils.SetNonHistoricalVariable(rVariable, rValue, rContainer);
}

template <class TDataType, class TContainerType, class TVarType = Variable<TDataType>>
void VariableUtilsSetNonHistoricalVariableForFlag(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    TContainerType &rContainer,
    const Flags Flag,
    const bool CheckValue = true)
{
    rVariableUtils.SetNonHistoricalVariable(rVariable, rValue, rContainer, Flag, CheckValue);
}

void PrintTimingInformation(Timer& rTimer)
{
    rTimer.PrintTimingInformation();
}

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

void AddUtilitiesToPython(pybind11::module &m)
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


    py::class_<DeflationUtils>(m,"DeflationUtils")
        .def(py::init<>())
        .def("VisualizeAggregates",&DeflationUtils::VisualizeAggregates)
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

    py::class_<VariableUtils>(m, "VariableUtils")
        .def(py::init<>())
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<bool>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<double>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<Matrix>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<bool>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<double>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<Vector>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<Matrix>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<bool>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<double>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<Matrix>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<bool>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<double>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<Vector>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<Matrix>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<bool>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<double>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Vector>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Matrix>>)
        .def("SetVectorVar", VariableUtilsSetVariable<array_1d<double,3>>)
        .def("SetVectorVar", VariableUtilsSetVariableForFlag<array_1d<double, 3>>)
        .def("SetScalarVar", VariableUtilsSetVariable<double>)
        .def("SetScalarVar", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetScalarVar", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetScalarVar", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetScalarVar", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetScalarVar", VariableUtilsSetVariableForFlag<double>)
        .def("SetScalarVar", VariableUtilsSetVariableForFlag<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetScalarVar", VariableUtilsSetVariableForFlag<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetScalarVar", VariableUtilsSetVariableForFlag<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetScalarVar", VariableUtilsSetVariableForFlag<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetNonHistoricalVectorVar", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetVariable", VariableUtilsSetVariable<int>)
        .def("SetVariable", VariableUtilsSetVariable<bool>)
        .def("SetVariable", VariableUtilsSetVariable<double>)
        .def("SetVariable", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetVariable", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetVariable", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetVariable", VariableUtilsSetVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetVariable", VariableUtilsSetVariable<array_1d<double, 3>>)
        .def("SetVariable", VariableUtilsSetVariable<Vector>)
        .def("SetVariable", VariableUtilsSetVariable<Matrix>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<bool>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<double>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<array_1d<double, 3>>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<array_1d<double, 4>>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<array_1d<double, 6>>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<array_1d<double, 9>>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<Vector>)
        .def("SetVariable", VariableUtilsSetVariableForFlag<Matrix>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<int>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<double>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<array_1d<double, 3>>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Matrix, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Matrix, ModelPart::ElementsContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::NodesContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ConditionsContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::NodesContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ConditionsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::NodesContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::ConditionsContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::ElementsContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::NodesContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::ConditionsContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::ElementsContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("SaveScalarVar", &VariableUtils::SaveVariable<double>)              // To be removed
        .def("SaveVectorVar", &VariableUtils::SaveVariable<array_1d<double, 3>>) // To be removed
        .def("SaveVariable", &VariableUtils::SaveVariable<bool>)
        .def("SaveVariable", &VariableUtils::SaveVariable<double>)
        .def("SaveVariable", &VariableUtils::SaveVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 3>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 4>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 6>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 9>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<Vector>)
        .def("SaveVariable", &VariableUtils::SaveVariable<Matrix>)
        .def("SaveScalarNonHistoricalVar", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType>)              // To be removed
        .def("SaveVectorNonHistoricalVar", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>) // To be removed
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<bool, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Vector, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Matrix, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<bool, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ConditionsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ConditionsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ConditionsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ConditionsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Vector, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Matrix, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<bool, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ElementsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ElementsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ElementsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ElementsContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Vector, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Matrix, ModelPart::ElementsContainerType>)
        .def("SelectNodeList", &VariableUtils::SelectNodeList)
        .def("CopyScalarVar", &VariableUtils::CopyVariable<double>)                                                                    // To be removed
        .def("CopyVectorVar", &VariableUtils::CopyVariable<array_1d<double, 3>>)                                                       // To be removed
        .def("CopyComponentVar", &VariableUtils::CopyVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>) // To be removed
        .def("CopyVariable", &VariableUtils::CopyVariable<bool>)
        .def("CopyVariable", &VariableUtils::CopyVariable<double>)
        .def("CopyVariable", &VariableUtils::CopyVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 3>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 4>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 6>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 9>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<Vector>)
        .def("CopyVariable", &VariableUtils::CopyVariable<Matrix>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<Variable<double>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<Variable<double>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalVariable<double>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalVariable<double, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumHistoricalNodeVectorVariable", &VariableUtils::SumHistoricalVariable<array_1d<double, 3>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumNonHistoricalNodeVectorVariable", &VariableUtils::SumNonHistoricalNodeVectorVariable)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<Variable<double>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumConditionVectorVariable", &VariableUtils::SumConditionVectorVariable)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<Variable<double>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumElementVectorVariable", &VariableUtils::SumElementVectorVariable)
        .def("AddDof", &VariableUtils::AddDof<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("AddDof", &VariableUtils::AddDof<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("CheckVariableKeys", &VariableUtils::CheckVariableKeys)
        .def("CheckDofs", &VariableUtils::CheckDofs)
        .def("UpdateCurrentToInitialConfiguration", &VariableUtils::UpdateCurrentToInitialConfiguration)
        .def("UpdateInitialToCurrentConfiguration", &VariableUtils::UpdateInitialToCurrentConfiguration)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPosition)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariable)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariableAndPosition);

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

    py::class_<NormalCalculationUtils > (m,"NormalCalculationUtils")
        .def(py::init<>())
        .def("CalculateOnSimplex", CalcOnSimplex_Cond)
        .def("CalculateOnSimplex", CalcOnSimplex_ModelPart)
        .def("CalculateOnSimplex", CalcOnSimplexWithDoubleVar)
        .def("CalculateOnSimplex", CalcOnSimplexWithIntVar)
        .def("CalculateOnSimplex", CalcOnSimplexWithDoubleVarAlpha)
        .def("SwapNormals", &NormalCalculationUtils::SwapNormals)
    //                    .def("CalculateOnSimplex", CalcOnSimplexWithArrayVar)
        ;

    py::class_<BodyNormalCalculationUtils > (m,"BodyNormalCalculationUtils")
        .def(py::init<>())
        .def("CalculateBodyNormals", &BodyNormalCalculationUtils::CalculateBodyNormals)
        ;

    py::class_<BodyDistanceCalculationUtils > (m,"BodyDistanceCalculationUtils")
        .def(py::init<>())
        .def("CalculateDistances2D", &BodyDistanceCalculationUtils::CalculateDistances < 2 >)
        .def("CalculateDistances3D", &BodyDistanceCalculationUtils::CalculateDistances < 3 >)
        ;

    py::class_<SignedDistanceCalculationUtils < 2 > >(m,"SignedDistanceCalculationUtils2D")
        .def(py::init<>())
        .def("CalculateDistances", &SignedDistanceCalculationUtils < 2 > ::CalculateDistances)
        .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 2 > ::FindMaximumEdgeSize)
        ;

    py::class_<SignedDistanceCalculationUtils < 3 > >(m,"SignedDistanceCalculationUtils3D")
        .def(py::init<>())
        .def("CalculateDistances", &SignedDistanceCalculationUtils < 3 > ::CalculateDistances)
        .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 3 > ::FindMaximumEdgeSize)
        ;

    py::class_<ParallelDistanceCalculator < 2 > >(m,"ParallelDistanceCalculator2D")
        .def(py::init<>())
        .def("CalculateDistances", CalculateDistancesDefault2D)
        .def("CalculateDistances", CalculateDistancesFlag2D)
        .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 2 > ::CalculateInterfacePreservingDistances)
        .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 2 > ::CalculateDistancesLagrangianSurface)
        .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 2 > ::FindMaximumEdgeSize)
        .def_readonly_static("CALCULATE_EXACT_DISTANCES_TO_PLANE", &ParallelDistanceCalculator<2>::CALCULATE_EXACT_DISTANCES_TO_PLANE)
        ;

    py::class_<ParallelDistanceCalculator < 3 > >(m,"ParallelDistanceCalculator3D")
        .def(py::init<>())
        .def("CalculateDistances", CalculateDistancesDefault3D)
        .def("CalculateDistances", CalculateDistancesFlag3D)
        .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 3 > ::CalculateInterfacePreservingDistances)
        .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 3 > ::CalculateDistancesLagrangianSurface)
        .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 3 > ::FindMaximumEdgeSize)
        .def_readonly_static("CALCULATE_EXACT_DISTANCES_TO_PLANE", &ParallelDistanceCalculator<3>::CALCULATE_EXACT_DISTANCES_TO_PLANE)
        ;

    py::class_<BruteForcePointLocator> (m, "BruteForcePointLocator")
        .def(py::init<ModelPart& >())
        .def("FindNode", &BruteForcePointLocator::FindNode)
        .def("FindElement", &BruteForcePointLocator::FindElement)
        .def("FindCondition", &BruteForcePointLocator::FindCondition)
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

    py::class_<IsosurfacePrinterApplication >(m,"IsosurfacePrinterApplication")
        .def(py::init<ModelPart& >() )
        .def("AddScalarVarIsosurface", &IsosurfacePrinterApplication::AddScalarVarIsosurface)
        .def("AddScalarVarIsosurfaceAndLower", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndLower)
        .def("AddScalarVarIsosurfaceAndHigher", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndHigher)
        .def("ClearData", &IsosurfacePrinterApplication::ClearData)
        .def("AddSkinConditions", &IsosurfacePrinterApplication::AddSkinConditions)
        .def("CreateNodesArray", &IsosurfacePrinterApplication::CreateNodesArray)
        ;

//     py::class_<SignedDistanceCalculationBinBased<2> >(m,"SignedDistanceCalculationBinBased2D", init<>())
//             .def("CalculateDistances",&SignedDistanceCalculationBinBased<2>::CalculateDistances )
//                             .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<2>::FindMaximumEdgeSize )
//             ;
//
//     py::class_<SignedDistanceCalculationBinBased<3> >(m,"SignedDistanceCalculationBinBased3D", init<>())
//             .def("CalculateDistances",&SignedDistanceCalculationBinBased<3>::CalculateDistances )
//                             .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<3>::FindMaximumEdgeSize )
//             ;

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

    py::class_<OpenMPUtils >(m,"OpenMPUtils")
        .def(py::init<>())
        .def_static("SetNumThreads", &OpenMPUtils::SetNumThreads)
    //     .staticmethod("SetNumThreads")
        .def_static("GetNumThreads", &OpenMPUtils::GetNumThreads)
    //     .staticmethod("GetNumThreads")
        .def_static("PrintOMPInfo", &OpenMPUtils::PrintOMPInfo)
    //     .staticmethod("PrintOMPInfo")
        ;

    py::class_< BinBasedFastPointLocator < 2 >, BinBasedFastPointLocator < 2 >::Pointer >(m,"BinBasedFastPointLocator2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", &BinBasedFastPointLocator < 2 > ::FindPointOnMeshSimplified)
        ;

    py::class_< BinBasedFastPointLocator < 3 >, BinBasedFastPointLocator < 3 >::Pointer >(m,"BinBasedFastPointLocator3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabase)
        .def("FindPointOnMesh", &BinBasedFastPointLocator < 3 > ::FindPointOnMeshSimplified)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
        ;

    py::class_< BinBasedFastPointLocatorConditions < 2 > >(m,"BinBasedFastPointLocatorConditions2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocatorConditions < 2 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocatorConditions < 2 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", &BinBasedFastPointLocatorConditions < 2 > ::FindPointOnMeshSimplified)
        ;

    py::class_< BinBasedFastPointLocatorConditions < 3 > >(m,"BinBasedFastPointLocatorConditions3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocatorConditions < 3 > ::UpdateSearchDatabase)
        .def("FindPointOnMesh", &BinBasedFastPointLocatorConditions < 3 > ::FindPointOnMeshSimplified)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocatorConditions < 3 > ::UpdateSearchDatabaseAssignedSize)
        ;

    py::class_< BinBasedNodesInElementLocator < 2 > >(m,"BinBasedNodesInElementLocator2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabase)
        .def("FindNodesInElement", &BinBasedNodesInElementLocator < 2 > ::FindNodesInElement)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
        ;

    py::class_< BinBasedNodesInElementLocator < 3 > >(m,"BinBasedNodesInElementLocator3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabase)
        .def("FindNodesInElement", &BinBasedNodesInElementLocator < 3 > ::FindNodesInElement)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
        ;

    py::class_< ActivationUtilities >(m,"ActivationUtilities")
        .def(py::init< >())
        .def("ActivateElementsAndConditions", &ActivationUtilities::ActivateElementsAndConditions)
        ;

    py::class_< EmbeddedSkinUtility < 2 > >(m,"EmbeddedSkinUtility2D")
        .def(py::init< ModelPart&, ModelPart&, const std::string >())
        .def("GenerateSkin", &EmbeddedSkinUtility < 2 > ::GenerateSkin)
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinArray< 2 > )
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinDouble< 2 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinArray< 2 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinDouble< 2 > )
        ;

    py::class_< EmbeddedSkinUtility <3 > >(m,"EmbeddedSkinUtility3D")
        .def(py::init< ModelPart&, ModelPart&, const std::string >())
        .def("GenerateSkin", &EmbeddedSkinUtility < 3 > ::GenerateSkin)
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinArray< 3 > )
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinDouble< 3 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinArray< 3 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinDouble< 3 > )
        ;

    py::class_< GeometryTesterUtility>(m,"GeometryTesterUtility")
        .def(py::init< >())
        .def("RunTest", &GeometryTesterUtility::RunTest)
        .def("TestTriangle2D3N", &GeometryTesterUtility::TestTriangle2D3N)
        .def("TestTriangle2D6N", &GeometryTesterUtility::TestTriangle2D6N)
        .def("TestTetrahedra3D4N", &GeometryTesterUtility::TestTetrahedra3D4N)
        .def("TestTetrahedra3D10N", &GeometryTesterUtility::TestTetrahedra3D10N)
        .def("TestHexahedra3D8N", &GeometryTesterUtility::TestHexahedra3D8N)
        .def("TestHexahedra3D27N", &GeometryTesterUtility::TestHexahedra3D27N)
        .def("TestHexahedra3D20N", &GeometryTesterUtility::TestHexahedra3D20N)
        .def("TestQuadrilateralInterface2D4N", &GeometryTesterUtility::TestQuadrilateralInterface2D4N)
        .def("TestPrismInterface3D6N", &GeometryTesterUtility::TestPrismInterface3D6N)
        .def("TestHexahedraInterface3D8N", &GeometryTesterUtility::TestHexahedraInterface3D8N)
        ;

    py::class_<CuttingUtility >(m,"CuttingUtility")
        .def(py::init< >())
        .def("GenerateCut", &CuttingUtility::GenerateCut)
        .def("UpdateCutData", &CuttingUtility ::UpdateCutData)
        .def("AddSkinConditions", &CuttingUtility ::AddSkinConditions)
        .def("AddVariablesToCutModelPart", &CuttingUtility::AddVariablesToCutModelPart )
        .def("FindSmallestEdge", &CuttingUtility ::FindSmallestEdge)
        ;

    py::class_<IntervalUtility >(m,"IntervalUtility")
        .def(py::init<Parameters >())
        .def("GetIntervalBegin", &IntervalUtility::GetIntervalBegin)
        .def("GetIntervalEnd", &IntervalUtility::GetIntervalEnd)
        .def("IsInInterval", &IntervalUtility ::IsInInterval)
        ;

    // Adding table from table stream to python
    py::class_<TableStreamUtility, typename TableStreamUtility::Pointer>(m,"TableStreamUtility")
        .def(py::init<>())
        .def(py::init< bool >())
        .def("SetOnProcessInfo",SetOnProcessInfo)
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

    // Sparse matrix multiplication utility
    py::class_<SparseMatrixMultiplicationUtility, typename SparseMatrixMultiplicationUtility::Pointer>(m, "SparseMatrixMultiplicationUtility")
        .def(py::init<>())
        .def_static("MatrixMultiplication",&SparseMatrixMultiplicationUtility::MatrixMultiplication<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def_static("MatrixMultiplicationSaad",&SparseMatrixMultiplicationUtility::MatrixMultiplicationSaad<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def_static("MatrixMultiplicationRMerge",&SparseMatrixMultiplicationUtility::MatrixMultiplicationRMerge<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def_static("MatrixAdd",&SparseMatrixMultiplicationUtility::MatrixAdd<CompressedMatrix, CompressedMatrix>)
        .def_static("TransposeMatrix",&SparseMatrixMultiplicationUtility::TransposeMatrix<CompressedMatrix, CompressedMatrix>)
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

    // Read materials utility
    py::class_<ReadMaterialsUtility, typename ReadMaterialsUtility::Pointer>(m, "ReadMaterialsUtility")
    .def(py::init<Model&>())
    .def(py::init<Parameters, Model&>())
    .def("ReadMaterials",&ReadMaterialsUtility::ReadMaterials)
    ;

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

    py::class_<SensitivityBuilder>(m, "SensitivityBuilder")
        .def(py::init<Parameters, ModelPart&, AdjointResponseFunction::Pointer>())
        .def("Initialize", &SensitivityBuilder::Initialize)
        .def("UpdateSensitivities", &SensitivityBuilder::UpdateSensitivities);

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

    // GeometricalTransformationUtilities
    auto mod_geom_trans_utils = m.def_submodule("GeometricalTransformationUtilities");
    mod_geom_trans_utils.def("CalculateTranslationMatrix", &GeometricalTransformationUtilities::CalculateTranslationMatrix );
    mod_geom_trans_utils.def("CalculateRotationMatrix", &GeometricalTransformationUtilities::CalculateRotationMatrix );

    // GeometricalTransformationUtilities
    auto mod_compare_elem_cond_utils = m.def_submodule("CompareElementsAndConditionsUtility");
    mod_compare_elem_cond_utils.def("GetRegisteredName", GetRegisteredNameElement );
    mod_compare_elem_cond_utils.def("GetRegisteredName", GetRegisteredNameCondition );
}

} // namespace Python.
} // Namespace Kratos
