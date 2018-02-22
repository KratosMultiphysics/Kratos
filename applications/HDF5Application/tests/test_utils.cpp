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

void CreateTestModelPart(ModelPart& rModelPart,
                         std::string const& rElementType,
                         std::string const& rConditionType)
{
    rModelPart.SubModelParts().clear();
    rModelPart.Elements().clear();
    rModelPart.Conditions().clear();
    rModelPart.Nodes().clear();
    rModelPart.rProperties().clear();
    const unsigned num_elems_and_conds = 3;
    const Element& r_elem = KratosComponents<Element>::Get(rElementType);
    const Condition& r_cond = KratosComponents<Condition>::Get(rConditionType);
    const unsigned num_nodes =
        num_elems_and_conds +
        std::max(r_elem.GetGeometry().size(), r_cond.GetGeometry().size());
    // Create nodes.
    for (unsigned i = 0; i < num_nodes; ++i)
        rModelPart.CreateNewNode(i + 1, 1.0, 2.0, 3.0);
    // Create elements.
    for (unsigned i = 0; i < num_elems_and_conds; ++i)
    {
        std::vector<ModelPart::IndexType> elem_node_ids(r_elem.GetGeometry().size());
        for (unsigned j = 0; j < r_elem.GetGeometry().size(); ++j)
            elem_node_ids.at(j) = i + j + 1;
        rModelPart.CreateNewElement(rElementType, i + 1, elem_node_ids, rModelPart.pGetProperties(1));
        std::vector<ModelPart::IndexType> cond_node_ids(r_cond.GetGeometry().size());
        for (unsigned j = 0; j < r_cond.GetGeometry().size(); ++j)
            cond_node_ids.at(j) = i + j + 1;
        rModelPart.CreateNewCondition(rConditionType, i + 1, cond_node_ids, rModelPart.pGetProperties(1));
    }
}

void CompareNodes(HDF5::NodesContainerType& rNodes1, HDF5::NodesContainerType& rNodes2)
{
    KRATOS_CHECK(rNodes1.size() == rNodes2.size());
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
    CompareNodes(rModelPart1.Nodes(), rModelPart2.Nodes());
    CompareElements(rModelPart1.Elements(), rModelPart2.Elements());
    CompareConditions(rModelPart1.Conditions(), rModelPart2.Conditions());
    for (ModelPart& sub_model_part_1 : rModelPart1.SubModelParts())
    {
        ModelPart& sub_model_part_2 = rModelPart2.GetSubModelPart(sub_model_part_1.Name());
        CompareModelParts(sub_model_part_1, sub_model_part_2);
    }
}

HDF5::File::Pointer pGetFile()
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    return boost::make_shared<HDF5::FileSerial>(file_params);
}

} // namespace Testing
} // namespace Kratos.
