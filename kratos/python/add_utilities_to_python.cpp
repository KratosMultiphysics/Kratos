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

// #include "utilities/signed_distance_calculator_bin_based.h"
#include "utilities/divide_elem_utils.h"
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
#include "utilities/sub_model_parts_list_utility.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "utilities/merge_variable_lists_utility.h"
#include "utilities/variable_redistribution_utility.h"
#include "utilities/sensitivity_builder.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/time_discretization.h"

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

void AddUtilitiesToPython(pybind11::module& m)
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

    py::class_<ConditionNumberUtility>(m,"ConditionNumberUtility")
        .def(py::init<>())
        .def(py::init<LinearSolverType::Pointer, LinearSolverType::Pointer>())
        .def("GetConditionNumber", ThisGetConditionNumber)
        .def("GetConditionNumber", ThisDirectGetConditionNumber)
        ;

    py::class_<VariableUtils>(m, "VariableUtils")
        .def(py::init<>())
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<bool>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<double>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<array_1d<double,3>>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<array_1d<double,4>>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<array_1d<double,6>>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<array_1d<double,9>>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVar", &VariableUtils::CopyModelPartNodalVar<Variable<Matrix>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<bool>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<double>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,3>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,4>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,6>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double,9>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Vector>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Matrix>>)
        .def("SetVectorVar", &VariableUtils::SetVectorVar)
        .def("SetVectorVar", &VariableUtils::SetVectorVarForFlag)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<Variable<double>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVar<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<Variable<double>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetScalarVar", &VariableUtils::SetScalarVarForFlag<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SetNonHistoricalScalarVar", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SetNonHistoricalVectorVar", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetVariable", &VariableUtils::SetVariable<bool>)
        .def("SetVariable", &VariableUtils::SetVariable<double>)
        .def("SetVariable", &VariableUtils::SetVariable<array_1d<double, 3>>)
        .def("SetVariable", &VariableUtils::SetVariable<Vector>)
        .def("SetVariable", &VariableUtils::SetVariable<Matrix>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<bool>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<double>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 3>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 4>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 6>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<array_1d<double, 9>>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<Vector>)
        .def("SetVariable", &VariableUtils::SetVariableForFlag<Matrix>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariable<Matrix, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", &VariableUtils::SetNonHistoricalVariableForFlag<Matrix, ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::NodesContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ConditionsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("SaveVectorVar", &VariableUtils::SaveVectorVar)
        .def("SaveScalarVar", &VariableUtils::SaveScalarVar)
        .def("SaveVectorNonHistoricalVar", &VariableUtils::SaveVectorNonHistoricalVar)
        .def("SaveScalarNonHistoricalVar", &VariableUtils::SaveScalarNonHistoricalVar)
        .def("SelectNodeList", &VariableUtils::SelectNodeList)
        .def("CopyVectorVar", &VariableUtils::CopyVectorVar)
        .def("CopyComponentVar", &VariableUtils::CopyComponentVar)
        .def("CopyScalarVar", &VariableUtils::CopyScalarVar)
        .def("SetToZero_VectorVar", &VariableUtils::SetToZero_VectorVar)
        .def("SetToZero_ScalarVar", &VariableUtils::SetToZero_ScalarVar)
//         .def("SetToZero_VelocityVectorVar", &VariableUtils::SetToZero_VelocityVectorVar)
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<double> >)
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 4> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 6> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< VariableComponent< VectorComponentAdaptor<array_1d<double, 9> > > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 3> > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 4> > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 6> > > )
//         .def("CheckVariableExists", &VariableUtils::CheckVariableExists< Variable<array_1d<double, 9> > > )
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
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalNodeScalarVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>)
        .def("SumHistoricalNodeVectorVariable", &VariableUtils::SumHistoricalNodeVectorVariable)
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
        .def("CalculateDistances", &ParallelDistanceCalculator < 2 > ::CalculateDistances)
        .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 2 > ::CalculateInterfacePreservingDistances)
        .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 2 > ::CalculateDistancesLagrangianSurface)
        .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 2 > ::FindMaximumEdgeSize)
        ;

    py::class_<ParallelDistanceCalculator < 3 > >(m,"ParallelDistanceCalculator3D")
        .def(py::init<>())
        .def("CalculateDistances", &ParallelDistanceCalculator < 3 > ::CalculateDistances)
        .def("CalculateInterfacePreservingDistances", &ParallelDistanceCalculator < 3 > ::CalculateInterfacePreservingDistances)
        .def("CalculateDistancesLagrangianSurface", &ParallelDistanceCalculator < 3 > ::CalculateDistancesLagrangianSurface)
        .def("FindMaximumEdgeSize", &ParallelDistanceCalculator < 3 > ::FindMaximumEdgeSize)
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

    py::class_<DivideElemUtils >(m,"DivideElemUtils")
        .def(py::init<>())
        .def("DivideElement_2D", &DivideElemUtils::DivideElement_2D)
        ;

    py::class_<Timer >(m,"Timer")
        .def(py::init<>())
        .def_property("PrintOnScreen", &Timer::GetPrintOnScreen, &Timer::SetPrintOnScreen)
        .def_static("Start", &Timer::Start)
        .def_static("Stop", &Timer::Stop)
    //     .staticmethod("Start")
    //     .staticmethod("Stop")
        //      .def("PrintTimingInformation",Timer::PrintTimingInformation)
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

    py::class_< BinBasedFastPointLocator < 2 > >(m,"BinBasedFastPointLocator2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", &BinBasedFastPointLocator < 2 > ::FindPointOnMeshSimplified)
        ;

    py::class_< BinBasedFastPointLocator < 3 > >(m,"BinBasedFastPointLocator3D")
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
        ;

    py::class_< EmbeddedSkinUtility <3 > >(m,"EmbeddedSkinUtility3D")
        .def(py::init< ModelPart&, ModelPart&, const std::string >())
        .def("GenerateSkin", &EmbeddedSkinUtility < 3 > ::GenerateSkin)
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinArray< 3 > )
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinDouble< 3 > )
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
        .def("MatrixMultiplication",&SparseMatrixMultiplicationUtility::MatrixMultiplication<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def("MatrixMultiplicationSaad",&SparseMatrixMultiplicationUtility::MatrixMultiplicationSaad<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def("MatrixMultiplicationRMerge",&SparseMatrixMultiplicationUtility::MatrixMultiplicationRMerge<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
        .def("MatrixAdd",&SparseMatrixMultiplicationUtility::MatrixAdd<CompressedMatrix, CompressedMatrix>)
        .def("TransposeMatrix",&SparseMatrixMultiplicationUtility::TransposeMatrix<CompressedMatrix, CompressedMatrix>)
        ;

    // Mortar utilities
    py::class_<MortarUtilities, typename MortarUtilities::Pointer>(m, "MortarUtilities")
        .def(py::init<>())
        .def("ComputeNodesMeanNormalModelPart",&MortarUtilities::ComputeNodesMeanNormalModelPart)
        .def("InvertNormal",&MortarUtilities::InvertNormal<PointerVectorSet<Element, IndexedObject>>)
        .def("InvertNormal",&MortarUtilities::InvertNormal<PointerVectorSet<Condition, IndexedObject>>)
        ;

    // Read materials utility
    py::class_<ReadMaterialsUtility, typename ReadMaterialsUtility::Pointer>(m, "ReadMaterialsUtility")
    .def(py::init<Model&>())
    .def(py::init<Parameters, Model&>())
    .def("ReadMaterials",&ReadMaterialsUtility::ReadMaterials)
    ;

    // SubModelParts List Utility
    py::class_<SubModelPartsListUtility, typename SubModelPartsListUtility::Pointer>(m, "SubModelPartsListUtility")
        .def(py::init<ModelPart&>())
        .def("DebugComputeSubModelPartsList",&SubModelPartsListUtility::DebugComputeSubModelPartsList)
        .def("GetRecursiveSubModelPartNames",&SubModelPartsListUtility::GetRecursiveSubModelPartNames)
        .def("GetRecursiveSubModelPart",&SubModelPartsListUtility::GetRecursiveSubModelPart)
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

    py::class_<TimeDiscretization::BDF1>(m, "BDF1")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", &TimeDiscretization::BDF1::ComputeBDFCoefficients)
        ;

    void (TimeDiscretization::BDF2::*ComputeBDF2Coefficients)(double, std::array<double, 3>&) const
        = &TimeDiscretization::BDF2::ComputeBDFCoefficients;
    void (TimeDiscretization::BDF2::*ComputeBDF2CoefficientsVariable)(double, double, std::array<double, 3>&) const
        = &TimeDiscretization::BDF2::ComputeBDFCoefficients;

    py::class_<TimeDiscretization::BDF2>(m, "BDF2")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", ComputeBDF2Coefficients)
        .def("ComputeBDFCoefficients", ComputeBDF2CoefficientsVariable)
        ;
    py::class_<TimeDiscretization::BDF3>(m, "BDF3")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", &TimeDiscretization::BDF3::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF4>(m, "BDF4")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", &TimeDiscretization::BDF4::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF5>(m, "BDF5")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", &TimeDiscretization::BDF5::ComputeBDFCoefficients)
        ;
    py::class_<TimeDiscretization::BDF6>(m, "BDF6")
        .def(py::init<>())
        .def("ComputeBDFCoefficients", &TimeDiscretization::BDF6::ComputeBDFCoefficients)
        ;

    py::class_<TimeDiscretization::Newmark>(m, "Newmark")
        .def(py::init<>())
        .def("GetBeta", &TimeDiscretization::Newmark::GetBeta)
        .def("GetGamma", &TimeDiscretization::Newmark::GetGamma)
        ;
    py::class_<TimeDiscretization::Bossak>(m, "Bossak")
        .def(py::init<const double>())
        .def("GetBeta", &TimeDiscretization::Bossak::GetBeta)
        .def("GetGamma", &TimeDiscretization::Bossak::GetGamma)
        ;
    py::class_<TimeDiscretization::GeneralizedAlpha>(m, "GeneralizedAlpha")
        .def(py::init<const double, const double>())
        .def("GetBeta", &TimeDiscretization::GeneralizedAlpha::GetBeta)
        .def("GetGamma", &TimeDiscretization::GeneralizedAlpha::GetGamma)
        ;

    std::size_t (*GetMinimumBufferSizeBDF1)(const TimeDiscretization::BDF1&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF2)(const TimeDiscretization::BDF2&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF3)(const TimeDiscretization::BDF3&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF4)(const TimeDiscretization::BDF4&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF5)(const TimeDiscretization::BDF5&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBDF6)(const TimeDiscretization::BDF6&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeNewmark)(const TimeDiscretization::Newmark&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeBossak)(const TimeDiscretization::Bossak&) = &TimeDiscretization::GetMinimumBufferSize;
    std::size_t (*GetMinimumBufferSizeGneralizedAlpha)(const TimeDiscretization::GeneralizedAlpha&) = &TimeDiscretization::GetMinimumBufferSize;

    m.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF1 );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF2 );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF3 );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF4 );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF5 );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeBDF6 );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeNewmark );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeBossak );
    m.def("GetMinimumBufferSize", GetMinimumBufferSizeGneralizedAlpha );
}

} // namespace Python.
} // Namespace Kratos
