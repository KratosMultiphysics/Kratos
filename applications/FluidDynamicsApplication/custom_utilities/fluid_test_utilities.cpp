//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <functional>
#include <random>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/model_part.h"

// Application includes

// Include base h
#include "fluid_test_utilities.h"

namespace Kratos
{
//// Static Operations ///////////////////////////////////////////////////////////////////

template <>
void FluidTestUtilities::AssignRandomValues(
    double& rValue,
    const std::string& rSeed,
    const int DomainSize,
    const double MinValue,
    const double MaxValue)
{
    std::seed_seq seed(rSeed.begin(), rSeed.end());

    // using the following specific generator to have consistent random numbers
    // generator across different platforms. No c++ distributions should be used
    // because they does not guarantee the multi-platform consistency.
    std::mt19937 generator(seed);
    rValue = MinValue + (generator() - generator.min()) * (MaxValue - MinValue) / (generator.max() - generator.min());
}

template <>
void FluidTestUtilities::AssignRandomValues(
    array_1d<double, 3>& rValue,
    const std::string& rSeed,
    const int DomainSize,
    const double MinValue,
    const double MaxValue)
{
    AssignRandomValues(rValue[0], rSeed + "_X", DomainSize, MinValue, MaxValue);
    AssignRandomValues(rValue[1], rSeed + "_Y", DomainSize, MinValue, MaxValue);

    if (DomainSize == 3) {
        AssignRandomValues(rValue[2], rSeed + "_Z", DomainSize, MinValue, MaxValue);
    } else {
        rValue[2] = 0.0;
    }
}

ModelPart& FluidTestUtilities::CreateTestModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(PropertiesType&)>& rSetElementProperties,
    const std::function<void(PropertiesType&)>& rSetConditionProperties,
    const std::function<void(ModelPart&)>& rAddNodalSolutionStepVariablesFuncion,
    const std::function<void(NodeType&)>& rAddDofsFunction,
    const int BufferSize)
{
    auto& r_model_part = rModel.CreateModelPart(rModelPartName, BufferSize);
    rAddNodalSolutionStepVariablesFuncion(r_model_part);

    r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 1.0, 0.0);
    r_model_part.CreateNewNode(3, 2.0, 1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        rAddDofsFunction(r_node);
    }

    IndexType eq_id = 1;
    for (auto& r_node : r_model_part.Nodes()) {
        for (auto& p_dof : r_node.GetDofs()) {
            p_dof->SetEquationId(eq_id++);
        }
    }

    using nid_list = std::vector<ModelPart::IndexType>;

    auto p_elem_prop = r_model_part.CreateNewProperties(1);
    rSetElementProperties(*p_elem_prop);
    auto p_element = r_model_part.CreateNewElement(rElementName, 1, nid_list{3, 2, 1}, p_elem_prop);

    auto p_cond_prop = r_model_part.CreateNewProperties(2);
    rSetConditionProperties(*p_cond_prop);
    auto p_condition_1 = r_model_part.CreateNewCondition(rConditionName, 1, nid_list{2, 1}, p_cond_prop);
    auto p_condition_2 = r_model_part.CreateNewCondition(rConditionName, 2, nid_list{3, 2}, p_cond_prop);
    auto p_condition_3 = r_model_part.CreateNewCondition(rConditionName, 3, nid_list{1, 3}, p_cond_prop);

    p_condition_1->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({p_element}));
    p_condition_2->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({p_element}));
    p_condition_3->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({p_element}));

    return r_model_part;
}

template<class TDataType>
void FluidTestUtilities::RandomFillHistoricalVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const std::string& rSeedExtension,
    const double MinValue,
    const double MaxValue,
    const int Step)
{
    for (auto& node : rModelPart.Nodes()) {
        std::stringstream seed;
        if (rVariable.Name() == "DISTANCE")
            seed << node.Id() << "_HistoricalV_" << "WALL_DISTANCE";
        else
            seed << node.Id() << "_HistoricalV_" << rVariable.Name();
        TDataType& value = node.FastGetSolutionStepValue(rVariable, Step);
        AssignRandomValues(value, seed.str(), rModelPart.GetProcessInfo()[DOMAIN_SIZE], MinValue, MaxValue);
    }
}

template<class TContainerType, class TDataType>
void FluidTestUtilities::RandomFillNonHistoricalVariable(
    TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const std::string& rSeedExtension,
    const IndexType DomainSize,
    const double MinValue,
    const double MaxValue)
{
    for (auto& item : rContainer) {
        std::stringstream seed;
        seed << item.Id() << "_NonHistoricalV_" << rSeedExtension;
        TDataType value = rVariable.Zero();
        FluidTestUtilities::AssignRandomValues(value, seed.str(), DomainSize, MinValue, MaxValue);
        item.SetValue(rVariable, value);
    }
}

template<class TContainerType>
void FluidTestUtilities::RunEntityGetDofListTest(
    const TContainerType& rContainer,
    const ProcessInfo& rProcessInfo,
    const std::vector<const Variable<double>*>& rDofVariablesList)
{
    KRATOS_TRY

    for (const auto& r_entity : rContainer) {

        DofsVectorType dofs;
        r_entity.GetDofList(dofs, rProcessInfo);

        auto& r_geometry = r_entity.GetGeometry();

        KRATOS_CHECK_EQUAL(dofs.size(), r_geometry.PointsNumber() * rDofVariablesList.size());

        IndexType index = 0;
        for (const auto& r_node : r_geometry) {
            for (auto p_variable : rDofVariablesList) {
                KRATOS_CHECK_EQUAL(r_node.pGetDof(*p_variable), dofs[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void FluidTestUtilities::RunEntityEquationIdVectorTest(
    const TContainerType& rContainer,
    const ProcessInfo& rProcessInfo,
    const std::vector<const Variable<double>*>& rDofVariablesList)
{
    KRATOS_TRY

    for (const auto& r_entity : rContainer) {

        EquationIdVectorType equation_ids;
        r_entity.EquationIdVector(equation_ids, rProcessInfo);

        auto& r_geometry = r_entity.GetGeometry();

        KRATOS_CHECK_EQUAL(equation_ids.size(),r_geometry.PointsNumber() * rDofVariablesList.size());

        IndexType index = 0;
        for (const auto& r_node : r_geometry) {
            for (auto p_variable : rDofVariablesList) {
                KRATOS_CHECK_EQUAL(r_node.GetDof(*p_variable).EquationId(), equation_ids[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void FluidTestUtilities::RunEntityGetValuesVectorTest(
    const TContainerType& rContainer,
    const std::vector<const Variable<double>*>& rDofVariablesList)
{
    KRATOS_TRY

    for (const auto& r_entity : rContainer) {

        Vector values;
        r_entity.GetValuesVector(values);

        auto& r_geometry = r_entity.GetGeometry();

        KRATOS_CHECK_EQUAL(values.size(),r_geometry.PointsNumber() * rDofVariablesList.size());

        IndexType index = 0;
        for (const auto& r_node : r_geometry) {
            for (auto p_variable : rDofVariablesList) {
                KRATOS_CHECK_EQUAL(r_node.FastGetSolutionStepValue(*p_variable), values[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void FluidTestUtilities::RunEntityGetFirstDerivativesVectorTest(
    const TContainerType& rContainer,
    const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod)
{
    KRATOS_TRY

    for (const auto& r_entity : rContainer) {

        Vector values;
        r_entity.GetFirstDerivativesVector(values);

        IndexType index = 0;
        for (const auto& r_node : r_entity.GetGeometry()) {
            const Vector& ref_values = rValueRetrievalMethod(r_node);
            for (IndexType i = 0; i < ref_values.size(); ++i) {
                KRATOS_CHECK_EQUAL(ref_values[i], values[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void FluidTestUtilities::RunEntityGetSecondDerivativesVectorTest(
    const TContainerType& rContainer,
    const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod)
{
    KRATOS_TRY

    for (const auto& r_entity : rContainer) {

        Vector values;
        r_entity.GetSecondDerivativesVector(values);

        IndexType index = 0;
        for (const auto& r_node : r_entity.GetGeometry()) {
            const Vector& ref_values = rValueRetrievalMethod(r_node);
            for (IndexType i = 0; i < ref_values.size(); ++i) {
                KRATOS_CHECK_EQUAL(ref_values[i], values[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template <>
ModelPart::ElementsContainerType& FluidTestUtilities::GetContainer<ModelPart::ElementsContainerType>::operator()(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
ModelPart::ConditionsContainerType& FluidTestUtilities::GetContainer<ModelPart::ConditionsContainerType>::operator()(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

//// Static Operations Template Instantiations //////////////////////////

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetDofListTest<ModelPart::ConditionsContainerType>(const ModelPart::ConditionsContainerType&, const ProcessInfo&, const std::vector<const Variable<double>*>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetDofListTest<ModelPart::ElementsContainerType>(const ModelPart::ElementsContainerType&, const ProcessInfo&, const std::vector<const Variable<double>*>&);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityEquationIdVectorTest<ModelPart::ConditionsContainerType>(const ModelPart::ConditionsContainerType&, const ProcessInfo&, const std::vector<const Variable<double>*>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityEquationIdVectorTest<ModelPart::ElementsContainerType>(const ModelPart::ElementsContainerType&, const ProcessInfo&, const std::vector<const Variable<double>*>&);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetValuesVectorTest<ModelPart::ConditionsContainerType>(const ModelPart::ConditionsContainerType&, const std::vector<const Variable<double>*>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetValuesVectorTest<ModelPart::ElementsContainerType>(const ModelPart::ElementsContainerType&, const std::vector<const Variable<double>*>&);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetFirstDerivativesVectorTest<ModelPart::ConditionsContainerType>(const ModelPart::ConditionsContainerType&, const std::function<Vector(const ModelPart::NodeType&)>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetFirstDerivativesVectorTest<ModelPart::ElementsContainerType>(const ModelPart::ElementsContainerType&, const std::function<Vector(const ModelPart::NodeType&)>&);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetSecondDerivativesVectorTest<ModelPart::ConditionsContainerType>(const ModelPart::ConditionsContainerType&, const std::function<Vector(const ModelPart::NodeType&)>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RunEntityGetSecondDerivativesVectorTest<ModelPart::ElementsContainerType>(const ModelPart::ElementsContainerType&, const std::function<Vector(const ModelPart::NodeType&)>&);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillHistoricalVariable<double>(ModelPart&, const Variable<double>&, const double, const double, const int);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillHistoricalVariable<double>(ModelPart&, const Variable<double>&, const std::string&, const double, const double, const int);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillHistoricalVariable<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&, const double, const double, const int);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillHistoricalVariable<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&, const std::string&, const double, const double, const int);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::NodesContainerType, double>(ModelPart::NodesContainerType&, const Variable<double>&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::NodesContainerType, double>(ModelPart::NodesContainerType&, const Variable<double>&, const std::string&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::NodesContainerType, array_1d<double, 3>>(ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::NodesContainerType, array_1d<double, 3>>(ModelPart::NodesContainerType&, const Variable<array_1d<double, 3>>&, const std::string&, const std::size_t, const double, const double);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ConditionsContainerType, double>(ModelPart::ConditionsContainerType&, const Variable<double>&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ConditionsContainerType, double>(ModelPart::ConditionsContainerType&, const Variable<double>&, const std::string&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>(ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>(ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const std::string&, const std::size_t, const double, const double);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ElementsContainerType, double>(ModelPart::ElementsContainerType&, const Variable<double>&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ElementsContainerType, double>(ModelPart::ElementsContainerType&, const Variable<double>&, const std::string&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ElementsContainerType, array_1d<double, 3>>(ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const std::size_t, const double, const double);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidTestUtilities::RandomFillNonHistoricalVariable<ModelPart::ElementsContainerType, array_1d<double, 3>>(ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const std::string&, const std::size_t, const double, const double);

template class FluidTestUtilities::GetContainer<ModelPart::ConditionsContainerType>;
template class FluidTestUtilities::GetContainer<ModelPart::ElementsContainerType>;

} // namespace Kratos