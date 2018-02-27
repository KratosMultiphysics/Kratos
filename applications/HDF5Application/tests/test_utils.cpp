#include "tests/test_utils.h"

#include <algorithm>

#include "containers/array_1d.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/kratos_components.h"
#include "testing/testing.h"
#include "custom_io/hdf5_file_serial.h"

namespace Kratos
{
namespace Testing
{

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

std::size_t TestModelPartFactory::AddNodes(std::size_t NumNodes)
{
    std::size_t start_index{0};
    if (mrTestModelPart.NumberOfNodes() > 0)
        start_index = mrTestModelPart.Nodes().back().Id();
    for (std::size_t i = start_index; i < NumNodes; ++i)
        mrTestModelPart.CreateNewNode(i + 1, 1.23, 4.56, 7.89);
    return mrTestModelPart.NumberOfNodes();
}

void TestModelPartFactory::AddNodalVariables(std::vector<std::string> const& rNodalVariables)
{
    for (const auto& r_variable_name : rNodalVariables)
    {
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name))
        {
            mrTestModelPart.AddNodalSolutionStepVariable(
                KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name));
        }
        else if (KratosComponents<Variable<double>>::Has(r_variable_name))
        {
            mrTestModelPart.AddNodalSolutionStepVariable(
                KratosComponents<Variable<double>>::Get(r_variable_name));
        }
        else if (KratosComponents<Variable<int>>::Has(r_variable_name))
        {
            mrTestModelPart.AddNodalSolutionStepVariable(
                KratosComponents<Variable<int>>::Get(r_variable_name));
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_variable_name << std::endl;
        }
    }
}

void TestModelPartFactory::SetBufferSize(std::size_t BufferSize)
{
    mrTestModelPart.SetBufferSize(BufferSize);
}

void TestModelPartFactory::AssignNodalTestData(std::vector<std::string> const& rNodalVariables)
{
    if (rNodalVariables.size() == 0)
        return;

    for (auto& r_node : mrTestModelPart.Nodes())
    {
        for (VariableData const& r_variable_data : *r_node.pGetVariablesList())
        {
            if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_data.Name()))
            {
                r_node.FastGetSolutionStepValue(
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(
                        r_variable_data.Name())) = array_1d<double, 3>(3, 1.23);
            }
            else if (KratosComponents<Variable<double>>::Has(r_variable_data.Name()))
            {
                r_node.FastGetSolutionStepValue(
                    KratosComponents<Variable<double>>::Get(r_variable_data.Name())) = 1.23;
            }
            else if (KratosComponents<Variable<int>>::Has(r_variable_data.Name()))
            {
                r_node.FastGetSolutionStepValue(
                    KratosComponents<Variable<int>>::Get(r_variable_data.Name())) = 123;
            }
            else
            {
                KRATOS_ERROR << "Unsupported variable type: " << r_variable_data.Name()
                             << std::endl;
            }
        }
    }
}

std::size_t TestModelPartFactory::AddElements(std::string const& rElement, std::size_t NumElems)
{
    const auto& r_elem = KratosComponents<Element>::Get(rElement);
    std::vector<ModelPart::IndexType> node_ids(r_elem.GetGeometry().size());
    std::size_t start_index = (mrTestModelPart.NumberOfElements() > 0) ? mrTestModelPart.Elements().back().Id() : 0;
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
    mrTestModelPart.CreateSubModelPart("SubModelPart0"); // Empty sub model part.
    ModelPart& sub_model_part1 =
        *mrTestModelPart.CreateSubModelPart("SubModelPart1");
    if (mrTestModelPart.NumberOfElements() > 0)
    {
        // Sub model part 1 contains elements but not conditions.
        auto p_elem = *mrTestModelPart.Elements().ptr_begin();
        sub_model_part1.AddElement(p_elem);
        for (auto it = p_elem->GetGeometry().ptr_begin();
             it != p_elem->GetGeometry().ptr_end(); ++it)
            sub_model_part1.AddNode(*it);
    }
    ModelPart& sub_model_part2 =
        *mrTestModelPart.CreateSubModelPart("SubModelPart2");
    if (mrTestModelPart.NumberOfConditions() > 0)
    {
        // Sub model part 2 contains conditions but not elements.
        auto p_cond = *mrTestModelPart.Conditions().ptr_begin();
        sub_model_part2.AddCondition(p_cond);
        for (auto it = p_cond->GetGeometry().ptr_begin();
             it != p_cond->GetGeometry().ptr_end(); ++it)
            sub_model_part2.AddNode(*it);
    }
    ModelPart& sub_model_part3 =
        *mrTestModelPart.CreateSubModelPart("SubModelPart3");
    if (mrTestModelPart.NumberOfElements() > 0 && mrTestModelPart.NumberOfConditions() > 0)
    {
        // Sub model part 3 contains elements and conditions.
        auto p_elem = *mrTestModelPart.Elements().ptr_begin();
        sub_model_part3.AddElement(p_elem);
        for (auto it = p_elem->GetGeometry().ptr_begin();
             it != p_elem->GetGeometry().ptr_end(); ++it)
            sub_model_part3.AddNode(*it);
        auto p_cond = *mrTestModelPart.Conditions().ptr_begin();
        sub_model_part3.AddCondition(p_cond);
        for (auto it = p_cond->GetGeometry().ptr_begin();
             it != p_cond->GetGeometry().ptr_end(); ++it)
            sub_model_part3.AddNode(*it);
    }
}

void CompareNodes(HDF5::NodesContainerType& rNodes1, HDF5::NodesContainerType& rNodes2)
{
    KRATOS_CHECK(rNodes1.size() == rNodes2.size());
    KRATOS_CHECK(rNodes1.size() > 0);
    for (const auto& r_node_1 : rNodes1)
    {
        HDF5::NodeType& r_node_2 = rNodes2[r_node_1.Id()];
        for (unsigned i = 0; i < 3; ++i)
            KRATOS_CHECK(r_node_1.Coordinates()[i] == r_node_2.Coordinates()[i]);
    }
}

void CompareElements(HDF5::ElementsContainerType& rElements1, HDF5::ElementsContainerType& rElements2)
{
    KRATOS_CHECK(rElements1.size() == rElements2.size());
    KRATOS_CHECK(rElements1.size() > 0);
    for (auto it = rElements1.begin(); it != rElements1.end(); ++it)
    {
        HDF5::ElementType& r_elem_1 = *it;
        HDF5::ElementType& r_elem_2 = rElements2[r_elem_1.Id()];
        KRATOS_CHECK(GeometricalObject::IsSame(r_elem_1, r_elem_2));
        for (unsigned i = 0; i < r_elem_1.GetGeometry().size(); ++i)
        {
            KRATOS_CHECK(r_elem_1.GetGeometry()[i].Id() == r_elem_2.GetGeometry()[i].Id());
            KRATOS_CHECK(r_elem_1.GetProperties().Id() == r_elem_2.GetProperties().Id());
        }
    }
}

void CompareConditions(HDF5::ConditionsContainerType& rConditions1, HDF5::ConditionsContainerType& rConditions2)
{
    KRATOS_CHECK(rConditions1.size() == rConditions2.size());
    KRATOS_CHECK(rConditions1.size() > 0);
    for (auto it = rConditions1.begin(); it != rConditions1.end(); ++it)
    {
        HDF5::ConditionType& r_cond_1 = *it;
        HDF5::ConditionType& r_cond_2 = rConditions2[r_cond_1.Id()];
        KRATOS_CHECK(GeometricalObject::IsSame(r_cond_1, r_cond_2));
        for (unsigned i = 0; i < r_cond_1.GetGeometry().size(); ++i)
        {
            KRATOS_CHECK(r_cond_1.GetGeometry()[i].Id() == r_cond_2.GetGeometry()[i].Id());
            KRATOS_CHECK(r_cond_1.GetProperties().Id() == r_cond_2.GetProperties().Id());
        }
    }
}

void CompareModelParts(ModelPart& rModelPart1, ModelPart& rModelPart2)
{
    KRATOS_CHECK(rModelPart1.IsSubModelPart() == rModelPart2.IsSubModelPart());
    if (!rModelPart1.IsSubModelPart() || rModelPart1.NumberOfNodes() > 0)
        CompareNodes(rModelPart1.Nodes(), rModelPart2.Nodes());
    if (!rModelPart1.IsSubModelPart() || rModelPart1.NumberOfElements() > 0)
        CompareElements(rModelPart1.Elements(), rModelPart2.Elements());
    if (!rModelPart1.IsSubModelPart() || rModelPart1.NumberOfConditions() > 0)
       CompareConditions(rModelPart1.Conditions(), rModelPart2.Conditions());
    KRATOS_CHECK(rModelPart1.NumberOfSubModelParts() == rModelPart2.NumberOfSubModelParts());
    for (ModelPart& sub_model_part_1 : rModelPart1.SubModelParts())
    {
        ModelPart& sub_model_part_2 = rModelPart2.GetSubModelPart(sub_model_part_1.Name());
        CompareModelParts(sub_model_part_1, sub_model_part_2);
    }
}

HDF5::File::Pointer pGetTestSerialFile()
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    return HDF5::File::Pointer(new HDF5::FileSerial(file_params));
}

HDF5::File GetTestFile()
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    return HDF5::File(file_params);
}

HDF5::FileSerial GetTestSerialFile()
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    return HDF5::FileSerial(file_params);
}

} // namespace Testing
} // namespace Kratos.
