//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Suneth Warnakulasruiya
//

// System includes
#include <algorithm>
#include <type_traits>

// Project includes
#include "containers/array_1d.h"
#include "includes/kratos_components.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/parallel_environment.h"
#include "testing/testing.h"

// Application includes
#include "custom_io/hdf5_file.h"

// Include base h
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{


namespace HDF5TestUtilities
{

template<class TDataType>
void AssignValue(TDataType& rValue)
{
    if constexpr(std::is_same_v<TDataType, int>) {
        rValue = 12345;
    } else if constexpr(std::is_same_v<TDataType, double>) {
        rValue = 1.2345;
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        rValue = array_1d<double, 3>(3, 1.2345);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 4>>) {
        rValue = array_1d<double, 4>(4, 1.2345);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 6>>) {
        rValue = array_1d<double, 6>(6, 1.2345);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 9>>) {
        rValue = array_1d<double, 9>(9, 1.2345);
    } else if constexpr(std::is_same_v<TDataType, Kratos::Vector>) {
        rValue = Kratos::Vector(2, 1.2345);
    } else if constexpr(std::is_same_v<TDataType, Kratos::Matrix>) {
        rValue = Kratos::Matrix(2, 2, 1.2345);
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
    }
}

template<class TDataType>
void CompareValues(
    const TDataType& rValue1,
    const TDataType& rValue2)
{
    if constexpr(std::is_same_v<TDataType, int>) {
        KRATOS_EXPECT_EQ(rValue1, rValue2);
    } else if constexpr(std::is_same_v<TDataType, double>) {
        KRATOS_EXPECT_EQ(rValue1, rValue2);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        KRATOS_EXPECT_VECTOR_EQ(rValue1, rValue2);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 4>>) {
        KRATOS_EXPECT_VECTOR_EQ(rValue1, rValue2);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 6>>) {
        KRATOS_EXPECT_VECTOR_EQ(rValue1, rValue2);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 9>>) {
        KRATOS_EXPECT_VECTOR_EQ(rValue1, rValue2);
    } else if constexpr(std::is_same_v<TDataType, Kratos::Vector>) {
        KRATOS_EXPECT_VECTOR_EQ(rValue1, rValue2);
    } else if constexpr(std::is_same_v<TDataType, Kratos::Matrix>) {
        KRATOS_EXPECT_MATRIX_EQ(rValue1, rValue2);
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
    }
}

template<typename TVariable>
bool AddNodalVariable(
    const std::string& rVariableName,
    ModelPart& rModelPart)
{
    if (KratosComponents<TVariable>::Has(rVariableName)) {
        const auto& r_variable = KratosComponents<TVariable>::Get(rVariableName);
        rModelPart.AddNodalSolutionStepVariable(r_variable);
        return true;
    } else {
        return false;
    }
}

template<class... TVariableType>
void AddNodalVariablesWrapper(
    const std::string& rNodalVariableName,
    ModelPart& rModelPart)
{
    KRATOS_TRY

    const bool is_added = (... || AddNodalVariable<TVariableType>(rNodalVariableName, rModelPart));

    KRATOS_ERROR_IF_NOT(is_added)
        << "Variable \"" << rNodalVariableName
        << "\" is not found in registered variables list.";

    KRATOS_CATCH("");
}

template<class TVariableType>
bool AddNodalTestData(
    const std::string& rVariableName,
    HDF5::NodesContainerType& rNodes)
{
    if (KratosComponents<TVariableType>::Has(rVariableName)) {
        const auto& r_variable = KratosComponents<TVariableType>::Get(rVariableName);
        for (auto& r_node : rNodes) {
            AssignValue(r_node.FastGetSolutionStepValue(r_variable));
        }
        return true;
    } else {
        return false;
    }
}

template<class... TVariableType>
void AddNodalTestDataWrapper(
    const std::string& rNodalVariableName,
    HDF5::NodesContainerType& rNodes)
{
    const bool is_added = (... || AddNodalTestData<TVariableType>(rNodalVariableName, rNodes));

    KRATOS_ERROR_IF_NOT(is_added)
        << "Variable \"" << rNodalVariableName
        << "\" is not found in registered variables list.";
}

template<class TComponentType>
bool AssignDataValueContainer(
    const std::string& rComponentName,
    DataValueContainer& rData,
    Flags& rFlags)
{
    if (KratosComponents<TComponentType>::Has(rComponentName)) {
        const auto& r_component = KratosComponents<TComponentType>::Get(rComponentName);
        if constexpr(std::is_same_v<TComponentType, Flags>) {
            rFlags.Set(r_component);
        } else {
            AssignValue(rData[r_component]);
        }
        return true;
    } else {
        return false;
    }
}

template<class... TComponentType>
void AssignDataValueContainerWrapper(
    const std::string& rComponentName,
    DataValueContainer& rData,
    Flags& rFlags)
{
    const bool is_added = (... || AssignDataValueContainer<TComponentType>(rComponentName, rData, rFlags));

    KRATOS_ERROR_IF_NOT(is_added)
        << "Component \"" << rComponentName
        << "\" is not found in registered variables list.";
}

template<class TComponentType>
bool CompareComponent(
    const std::string& rComponentName,
    const DataValueContainer& rData1,
    const Flags& rFlags1,
    const DataValueContainer& rData2,
    const Flags& rFlags2)
{
    if (KratosComponents<TComponentType>::Has(rComponentName)) {
        const auto& r_component = KratosComponents<TComponentType>::Get(rComponentName);
        if constexpr(std::is_same_v<TComponentType, Flags>) {
            KRATOS_EXPECT_TRUE(rFlags1.Is(r_component) == rFlags2.Is(r_component));
        } else {
            KRATOS_EXPECT_TRUE(rData1.Has(r_component) && rData2.Has(r_component));
            CompareValues(rData1[r_component], rData2[r_component]);
        }
        return true;
    } else {
        return false;
    }
}

template<class... TComponentType>
void CompareComponentWrapper(
    const std::string& rComponentName,
    const DataValueContainer& rData1,
    const Flags& rFlags1,
    const DataValueContainer& rData2,
    const Flags& rFlags2)
{
    const bool is_added = (... || CompareComponent<TComponentType>(rComponentName, rData1, rFlags1, rData2, rFlags2));

    KRATOS_ERROR_IF_NOT(is_added)
        << "Component \"" << rComponentName
        << "\" is not found in registered variables list.";
}

} // namespace HDF5TestUtilities

template<>
HDF5::File::Vector<array_1d<double,3>> TestVector(std::size_t n)
{
    HDF5::File::Vector<array_1d<double,3>> vec(n);
    for (std::size_t i = 0; i < n; ++i)
        vec(i) = array_1d<double,3>(3, i);
    return vec;
}

void TestModelPartFactory::CreateModelPart(ModelPart& rTestModelPart,
                                           std::vector<std::string> const& rElements,
                                           std::vector<std::string> const& rConditions,
                                           std::vector<std::string> const& rNodalVariables)
{
    const std::size_t num_elems_and_conds{2};
    std::size_t num_nodes{2};
    for (const auto& r_name : rElements)
    {
        const Element& r_elem = KratosComponents<Element>::Get(r_name);
        num_nodes = std::max(num_nodes, num_elems_and_conds + r_elem.GetGeometry().size());
    }
    for (const auto& r_name : rConditions)
    {
        const Condition& r_cond = KratosComponents<Condition>::Get(r_name);
        num_nodes = std::max(num_nodes, num_elems_and_conds + r_cond.GetGeometry().size());
    }
    TestModelPartFactory model_part_factory(rTestModelPart);
    model_part_factory.AddNodalVariables(rNodalVariables);
    model_part_factory.AddNodes(num_nodes);
    for (const auto& r_name : rElements)
        model_part_factory.AddElements(r_name, num_elems_and_conds);
    for (const auto& r_name : rConditions)
        model_part_factory.AddConditions(r_name, num_elems_and_conds);
    model_part_factory.AddSubModelParts();
    model_part_factory.SetBufferSize(2);
    model_part_factory.AssignNodalTestData(rNodalVariables);
}

TestModelPartFactory::TestModelPartFactory(ModelPart& rTestModelPart)
    : mrTestModelPart(rTestModelPart)
{
}

void TestModelPartFactory::AddNodalVariables(std::vector<std::string> const& rNodalVariables)
{
    for (const auto& r_variable_name : rNodalVariables) {
        HDF5TestUtilities::AddNodalVariablesWrapper<
                                                Variable<int>,
                                                Variable<double>,
                                                Variable<array_1d<double, 3>>,
                                                Variable<array_1d<double, 4>>,
                                                Variable<array_1d<double, 6>>,
                                                Variable<array_1d<double, 9>>,
                                                Variable<Kratos::Vector>,
                                                Variable<Kratos::Matrix>>(r_variable_name, mrTestModelPart);
    }
}

std::size_t TestModelPartFactory::AddNodes(std::size_t NumNodes)
{
    std::size_t start_index{0};
    if (mrTestModelPart.NumberOfNodes() > 0)
        start_index = mrTestModelPart.Nodes().back().Id();
    for (std::size_t i = start_index; i < NumNodes; ++i)
        mrTestModelPart.CreateNewNode(i + 1, 1.23, 4.56, 7.89);
    return mrTestModelPart.NumberOfNodes();
}

void TestModelPartFactory::SetBufferSize(std::size_t BufferSize)
{
    mrTestModelPart.SetBufferSize(BufferSize);
}

void TestModelPartFactory::AssignNodalTestData(std::vector<std::string> const& rNodalVariables)
{
    if (rNodalVariables.size() == 0)
        return;

    for (auto& r_name : rNodalVariables) {
        HDF5TestUtilities::AddNodalTestDataWrapper<
                                                Variable<int>,
                                                Variable<double>,
                                                Variable<array_1d<double, 3>>,
                                                Variable<array_1d<double, 4>>,
                                                Variable<array_1d<double, 6>>,
                                                Variable<array_1d<double, 9>>,
                                                Variable<Kratos::Vector>,
                                                Variable<Kratos::Matrix>>(r_name, mrTestModelPart.Nodes());
    }
}

void TestModelPartFactory::AssignNonHistoricalNodalTestData(ModelPart& rTestModelPart,
                                                            std::vector<std::string> const& rNodalVariables)
{
    if (rNodalVariables.size() == 0)
        return;

    for (auto& r_node : rTestModelPart.Nodes())
        AssignDataValueContainer(r_node.GetData(), r_node, rNodalVariables);
}

void TestModelPartFactory::AssignDataValueContainer(DataValueContainer& rData, Flags& rFlags, std::vector<std::string> const& rVariables)
{
    for (auto& r_name : rVariables) {
        HDF5TestUtilities::AssignDataValueContainerWrapper<
                                                    Flags,
                                                    Variable<int>,
                                                    Variable<double>,
                                                    Variable<array_1d<double, 3>>,
                                                    Variable<array_1d<double, 4>>,
                                                    Variable<array_1d<double, 6>>,
                                                    Variable<array_1d<double, 9>>,
                                                    Variable<Kratos::Vector>,
                                                    Variable<Kratos::Matrix>>(r_name, rData, rFlags);
    }
}

std::size_t TestModelPartFactory::AddElements(std::string const& rElement, std::size_t NumElems)
{
    const auto& r_elem = KratosComponents<Element>::Get(rElement);
    std::vector<ModelPart::IndexType> node_ids(r_elem.GetGeometry().size());
    std::size_t start_index = (mrTestModelPart.NumberOfElements() > 0) ? mrTestModelPart.Elements().back().Id() : 0;
    if (!mrTestModelPart.HasProperties(1))
    {
        mrTestModelPart.CreateNewProperties(1);
    }
    for (std::size_t i = 0; i < NumElems; ++i)
    {
        for (std::size_t j = 0; j < r_elem.GetGeometry().size(); ++j)
            node_ids.at(j) = i + j + 1;
        mrTestModelPart.CreateNewElement(rElement, start_index + i + 1, node_ids, mrTestModelPart.pGetProperties(1));
    }
    return mrTestModelPart.NumberOfElements();
}

std::size_t TestModelPartFactory::AddConditions(std::string const& rCondition, std::size_t NumConds)
{
    const auto& r_cond = KratosComponents<Condition>::Get(rCondition);
    std::vector<ModelPart::IndexType> node_ids(r_cond.GetGeometry().size());
    std::size_t start_index = (mrTestModelPart.NumberOfConditions() > 0) ? mrTestModelPart.Conditions().back().Id() : 0;
    if (!mrTestModelPart.HasProperties(1))
    {
        mrTestModelPart.CreateNewProperties(1);
    }
    for (std::size_t i = 0; i < NumConds; ++i)
    {
        for (std::size_t j = 0; j < r_cond.GetGeometry().size(); ++j)
            node_ids.at(j) = i + j + 1;
        mrTestModelPart.CreateNewCondition(rCondition, start_index + i + 1, node_ids, mrTestModelPart.pGetProperties(1));
    }
    return mrTestModelPart.NumberOfConditions();
}

void TestModelPartFactory::AddSubModelParts()
{
    AddEmptySubModelPart();
    AddElementsSubModelPart();
    AddConditionsSubModelPart();
    AddElementsAndConditionsSubModelPart();
}

void TestModelPartFactory::AddEmptySubModelPart()
{
    mrTestModelPart.CreateSubModelPart("EmptySubModelPart");
}

void TestModelPartFactory::AddElementsSubModelPart()
{
    ModelPart& sub_model_part =
        mrTestModelPart.CreateSubModelPart("ElementsSubModelPart");
    if (mrTestModelPart.NumberOfElements() > 0)
    {
        auto p_elem = *mrTestModelPart.Elements().ptr_begin();
        sub_model_part.AddElement(p_elem);
        for (auto it = p_elem->GetGeometry().ptr_begin();
             it != p_elem->GetGeometry().ptr_end(); ++it)
            sub_model_part.AddNode(*it);
    }
}

void TestModelPartFactory::AddConditionsSubModelPart()
{
    ModelPart& sub_model_part =
        mrTestModelPart.CreateSubModelPart("ConditionsSubModelPart");
    if (mrTestModelPart.NumberOfConditions() > 0)
    {
        auto p_cond = *mrTestModelPart.Conditions().ptr_begin();
        sub_model_part.AddCondition(p_cond);
        for (auto it = p_cond->GetGeometry().ptr_begin();
             it != p_cond->GetGeometry().ptr_end(); ++it)
            sub_model_part.AddNode(*it);
    }
}

void TestModelPartFactory::AddElementsAndConditionsSubModelPart()
{
    ModelPart& sub_model_part =
        mrTestModelPart.CreateSubModelPart("ElementsAndConditionsSubModelPart");
    if (mrTestModelPart.NumberOfElements() > 0 && mrTestModelPart.NumberOfConditions() > 0)
    {
        auto p_elem = *mrTestModelPart.Elements().ptr_begin();
        sub_model_part.AddElement(p_elem);
        for (auto it = p_elem->GetGeometry().ptr_begin();
             it != p_elem->GetGeometry().ptr_end(); ++it)
            sub_model_part.AddNode(*it);
        auto p_cond = *mrTestModelPart.Conditions().ptr_begin();
        sub_model_part.AddCondition(p_cond);
        for (auto it = p_cond->GetGeometry().ptr_begin();
             it != p_cond->GetGeometry().ptr_end(); ++it)
            sub_model_part.AddNode(*it);
    }
}

void CompareNodes(HDF5::NodesContainerType& rNodes1, HDF5::NodesContainerType& rNodes2)
{
    KRATOS_EXPECT_TRUE(rNodes1.size() == rNodes2.size());
    KRATOS_EXPECT_TRUE(rNodes1.size() > 0);
    for (const auto& r_node_1 : rNodes1)
    {
        HDF5::NodeType& r_node_2 = rNodes2[r_node_1.Id()];
        for (unsigned i = 0; i < 3; ++i)
            KRATOS_EXPECT_TRUE(r_node_1.Coordinates()[i] == r_node_2.Coordinates()[i]);
    }
}

void CompareElements(HDF5::ElementsContainerType& rElements1, HDF5::ElementsContainerType& rElements2)
{
    KRATOS_EXPECT_TRUE(rElements1.size() == rElements2.size());
    KRATOS_EXPECT_TRUE(rElements1.size() > 0);
    for (auto it = rElements1.begin(); it != rElements1.end(); ++it)
    {
        HDF5::ElementType& r_elem_1 = *it;
        HDF5::ElementType& r_elem_2 = rElements2[r_elem_1.Id()];
        KRATOS_EXPECT_TRUE(GeometricalObject::IsSame(r_elem_1, r_elem_2));
        for (unsigned i = 0; i < r_elem_1.GetGeometry().size(); ++i)
        {
            KRATOS_EXPECT_TRUE(r_elem_1.GetGeometry()[i].Id() == r_elem_2.GetGeometry()[i].Id());
            KRATOS_EXPECT_TRUE(r_elem_1.GetProperties().Id() == r_elem_2.GetProperties().Id());
        }
    }
}

void CompareConditions(HDF5::ConditionsContainerType& rConditions1, HDF5::ConditionsContainerType& rConditions2)
{
    KRATOS_EXPECT_TRUE(rConditions1.size() == rConditions2.size());
    KRATOS_EXPECT_TRUE(rConditions1.size() > 0);
    for (auto it = rConditions1.begin(); it != rConditions1.end(); ++it)
    {
        HDF5::ConditionType& r_cond_1 = *it;
        HDF5::ConditionType& r_cond_2 = rConditions2[r_cond_1.Id()];
        KRATOS_EXPECT_TRUE(GeometricalObject::IsSame(r_cond_1, r_cond_2));
        for (unsigned i = 0; i < r_cond_1.GetGeometry().size(); ++i)
        {
            KRATOS_EXPECT_TRUE(r_cond_1.GetGeometry()[i].Id() == r_cond_2.GetGeometry()[i].Id());
            KRATOS_EXPECT_TRUE(r_cond_1.GetProperties().Id() == r_cond_2.GetProperties().Id());
        }
    }
}

void CompareModelParts(ModelPart& rModelPart1, ModelPart& rModelPart2)
{
    KRATOS_EXPECT_TRUE(rModelPart1.IsSubModelPart() == rModelPart2.IsSubModelPart());
    if (!rModelPart1.IsSubModelPart() || rModelPart1.NumberOfNodes() > 0)
        CompareNodes(rModelPart1.Nodes(), rModelPart2.Nodes());
    if (!rModelPart1.IsSubModelPart() || rModelPart1.NumberOfElements() > 0)
        CompareElements(rModelPart1.Elements(), rModelPart2.Elements());
    if (!rModelPart1.IsSubModelPart() || rModelPart1.NumberOfConditions() > 0)
       CompareConditions(rModelPart1.Conditions(), rModelPart2.Conditions());
    KRATOS_EXPECT_TRUE(rModelPart1.NumberOfSubModelParts() == rModelPart2.NumberOfSubModelParts());
    for (ModelPart& sub_model_part_1 : rModelPart1.SubModelParts())
    {
        ModelPart& sub_model_part_2 = rModelPart2.GetSubModelPart(sub_model_part_1.Name());
        CompareModelParts(sub_model_part_1, sub_model_part_2);
    }
}

void CompareNonHistoricalNodalData(
    const std::vector<std::string>& rFlagNames,
    HDF5::NodesContainerType& rNodes1,
    HDF5::NodesContainerType& rNodes2)
{
    KRATOS_EXPECT_TRUE(rNodes1.size() == rNodes2.size());
    for (auto& r_node1 : rNodes1) {
        auto& r_node2 = rNodes2[r_node1.Id()];
        CompareDataValueContainers(rFlagNames, r_node1.GetData(), r_node1, r_node2.GetData(), r_node2);
    }
}

void CompareDataValueContainers(
    const std::vector<std::string>& rFlagNames,
    DataValueContainer const& rData1,
    Flags const& rFlags1,
    DataValueContainer const& rData2,
    Flags const& rFlags2)
{
    for (const auto& r_value1 : rData1) {
        HDF5TestUtilities::CompareComponentWrapper<
                                                Variable<int>,
                                                Variable<double>,
                                                Variable<array_1d<double, 3>>,
                                                Variable<array_1d<double, 4>>,
                                                Variable<array_1d<double, 6>>,
                                                Variable<array_1d<double, 9>>,
                                                Variable<Kratos::Vector>,
                                                Variable<Kratos::Matrix>>(r_value1.first->Name(), rData1, rFlags1, rData2, rFlags2);
    }

    for (const auto& r_flag_name : rFlagNames) {
        HDF5TestUtilities::CompareComponentWrapper<
                                                Flags>(r_flag_name, rData1, rFlags1, rData2, rFlags2);
    }
}

HDF5::File::Pointer pGetTestSerialFile()
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    return HDF5::File::Pointer(new HDF5::File(Testing::GetDefaultDataCommunicator(), file_params));
}

HDF5::File GetTestFile()
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    return HDF5::File(Testing::GetDefaultDataCommunicator(), file_params);
}

HDF5::File GetTestSerialFile()
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    return HDF5::File(Testing::GetDefaultDataCommunicator(), file_params);
}

} // namespace Testing
} // namespace Kratos.
