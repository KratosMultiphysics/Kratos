#include "tests/test_utils.h"

#include <algorithm>

#include "containers/array_1d.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/kratos_components.h"
#include "testing/testing.h"
#include "custom_io/hdf5_file_serial.h"
#include "custom_utilities/registered_component_lookup.h"

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

namespace {
template <typename TVariable>
class AddNodalVariableFunctor;
}

void TestModelPartFactory::AddNodalVariables(std::vector<std::string> const& rNodalVariables)
{
    for (const auto& r_name : rNodalVariables)
        RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>>(r_name)
            .Execute<AddNodalVariableFunctor>(mrTestModelPart);
}

namespace{
template <typename TVariable>
class AddNodalVariableFunctor
{
public:
    void operator()(TVariable const& rVariable, ModelPart& rModelPart)
    {
        rModelPart.AddNodalSolutionStepVariable(rVariable);
    }
};
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

namespace {
template <typename TVariable>
class AssignNodalSolutionStepValueFunctor;
}

void TestModelPartFactory::AssignNodalTestData(std::vector<std::string> const& rNodalVariables)
{
    if (rNodalVariables.size() == 0)
        return;

    for (auto& r_name : rNodalVariables)
    {
        RegisteredComponentLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>>(r_name)
            .Execute<AssignNodalSolutionStepValueFunctor>(mrTestModelPart.Nodes());
    }
}

namespace {
void AssignValue(double&);
void AssignValue(int&);
void AssignValue(array_1d<double, 3>&);

template <typename TVariable>
class AssignNodalSolutionStepValueFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    HDF5::NodesContainerType& rNodes)
    {
        for (auto& r_node : rNodes)
            AssignValue(r_node.FastGetSolutionStepValue(rVariable));
    }
};


void AssignValue(double& d)
{
    d = 1.2345;
}

void AssignValue(int& i)
{
    i = 12345;
}

void AssignValue(array_1d<double, 3>& v)
{
    v = array_1d<double, 3>(3, 1.2345);
}
}

void TestModelPartFactory::AssignNonHistoricalNodalTestData(ModelPart& rTestModelPart,
                                                            std::vector<std::string> const& rNodalVariables)
{
    if (rNodalVariables.size() == 0)
        return;

    for (auto& r_node : rTestModelPart.Nodes())
        AssignDataValueContainer(r_node.Data(), r_node, rNodalVariables);
}

namespace {
template <typename TVariable>
class AssignDataValueContainerFunctor;
}

void TestModelPartFactory::AssignDataValueContainer(DataValueContainer& rData, Flags& rFlags, std::vector<std::string> const& rVariables)
{
    for (auto& r_name : rVariables)
        RegisteredComponentLookup<Flags, Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                                 Variable<HDF5::Vector<double>>, Variable<HDF5::Matrix<double>>>(r_name)
            .Execute<AssignDataValueContainerFunctor>(rData, rFlags);
}

namespace {
void AssignValue(HDF5::Vector<double>&);
void AssignValue(HDF5::Matrix<double>&);

template <typename TVariable>
class AssignDataValueContainerFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    DataValueContainer& rData,
                    Flags&)
    {
        AssignValue(rData[rVariable]);
    }
};

template <>
class AssignDataValueContainerFunctor<Flags>
{
public:
    void operator()(Flags const& rVariable,
                    DataValueContainer& rData,
                    Flags& rFlags)
    {
        rFlags.Set(rVariable, static_cast<bool>((counter++) % 2));
    }
private:
    int counter = 0;
};

void AssignValue(HDF5::Vector<double>& v)
{
    const std::size_t dim = 2;
    v.resize(dim);
    std::fill(v.begin(), v.end(), 1.2345);
}

void AssignValue(HDF5::Matrix<double>& m)
{
    const std::size_t dim = 2;
    m.resize(dim, dim, false);
    std::fill(m.data().begin(), m.data().end(), 1.2345);
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

void CompareNonHistoricalNodalData(HDF5::NodesContainerType& rNodes1,
                                   HDF5::NodesContainerType& rNodes2)
{
    KRATOS_CHECK(rNodes1.size() == rNodes2.size());
    for (auto& r_node1 : rNodes1)
    {
        auto& r_node2 = rNodes2[r_node1.Id()];
        CompareDataValueContainers(r_node1.Data(), r_node1, r_node2.Data(), r_node2);
    }
}

namespace {
template <typename TVariable>
class CompareVariableFunctor;
}

void CompareDataValueContainers(DataValueContainer const& rData1, Flags const& rFlags1, DataValueContainer const& rData2, Flags const& rFlags2)
{
    for (const auto& r_value1 : rData1)
        RegisteredComponentLookup<Flags, Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                                 Variable<HDF5::Vector<double>>, Variable<HDF5::Matrix<double>>>(
            r_value1.first->Name())
            .Execute<CompareVariableFunctor>(rData1, rFlags1, rData2, rFlags2);
}

namespace {
void CompareValues(int, int);
void CompareValues(double, double);
void CompareValues(HDF5::Vector<double> const&, HDF5::Vector<double> const&);
void CompareValues(HDF5::Matrix<double> const&, HDF5::Matrix<double> const&);

template <typename TVariable>
class CompareVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    DataValueContainer const& rData1,
                    Flags const&,
                    DataValueContainer const& rData2,
                    Flags const&)
    {
        KRATOS_CHECK(rData1.Has(rVariable) && rData2.Has(rVariable));
        CompareValues(rData1[rVariable], rData2[rVariable]);
    }
};

template <>
class CompareVariableFunctor<Flags>
{
public:
    void operator()(Flags const& rVariable,
                    DataValueContainer const&,
                    Flags const& rFlags1,
                    DataValueContainer const&,
                    Flags const& rFlags2)
    {
        KRATOS_CHECK(rFlags1.Is(rVariable) == rFlags2.Is(rVariable));
    }
};

void CompareValues(int i1, int i2)
{
    KRATOS_CHECK(i1 == i2);
}

void CompareValues(double d1, double d2)
{
    KRATOS_CHECK(d1 == d2);
}

void CompareValues(HDF5::Vector<double> const& v1, HDF5::Vector<double> const& v2)
{
    KRATOS_CHECK(v1.size() == v2.size());
    for (std::size_t i = 0; i < v1.size(); ++i)
        KRATOS_CHECK(v1(i) == v2(i));
}

void CompareValues(HDF5::Matrix<double> const& m1, HDF5::Matrix<double> const& m2)
{
    KRATOS_CHECK(m1.size1() == m2.size1());
    KRATOS_CHECK(m1.size2() == m2.size2());
    for (std::size_t i = 0; i < m1.size1(); ++i)
        for (std::size_t j = 0; j < m1.size2(); ++j)
            KRATOS_CHECK(m1(i,j) == m2(i,j));
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
