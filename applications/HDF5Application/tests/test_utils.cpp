#include "tests/test_utils.h"

#include "containers/array_1d.h"

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

void CreateTestMesh(HDF5::NodesContainerType& rNodes,
                    HDF5::PropertiesContainerType& rProperties,
                    HDF5::ElementsContainerType& rElements,
                    HDF5::Vector<int>& rElementIds,
                    HDF5::Vector<int>& rPropertiesIds,
                    HDF5::Matrix<int>& rConnectivities)
{
    const unsigned num_elems = 2;
    const unsigned num_nodes = num_elems + 3;

    rNodes.clear();
    rProperties.clear();
    rElements.clear();
    rElementIds.resize(num_elems, false);
    rPropertiesIds.resize(num_elems, false);
    rConnectivities.resize(num_elems, 3, false);

    const HDF5::ElementType& Element2D3N = KratosComponents<HDF5::ElementType>::Get("Element2D3N");
    // Create nodes.
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        auto p_node = boost::make_shared<HDF5::NodeType>(
            i + 1, 1.0, 2.0, 3.0);
        rNodes.push_back(p_node);
    }
    // Create elements.
    Element::NodesArrayType geom_nodes(3);
    for (unsigned i = 0; i < num_elems; ++i)
    {
        rElementIds[i] = i + 1;
        rPropertiesIds[i] = 1;
        for (unsigned j = 0; j < 3; ++j)
        {
            rConnectivities(i, j) = i + j + 1;
            geom_nodes(j) = rNodes(i + j + 1);
        }
        auto p_elem = Element2D3N.Create(i + 1, geom_nodes, rProperties(1));
        rElements.push_back(p_elem);
    }
}

void CreateTestMesh(HDF5::NodesContainerType& rNodes,
                    HDF5::PropertiesContainerType& rProperties,
                    HDF5::ConditionsContainerType& rConditions,
                    HDF5::Vector<int>& rConditionIds,
                    HDF5::Vector<int>& rPropertiesIds,
                    HDF5::Matrix<int>& rConnectivities)
{
    const unsigned num_conds = 2;
    const unsigned num_nodes = num_conds + 3;

    rNodes.clear();
    rProperties.clear();
    rConditions.clear();
    rConditionIds.resize(num_conds, false);
    rPropertiesIds.resize(num_conds, false);
    rConnectivities.resize(num_conds, 3, false);

    const HDF5::ConditionType& SurfaceCondition3D3N = KratosComponents<HDF5::ConditionType>::Get("SurfaceCondition3D3N");
    // Create nodes.
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        auto p_node = boost::make_shared<HDF5::NodeType>(
            i + 1, 1.0, 2.0, 3.0);
        rNodes.push_back(p_node);
    }
    // Create conditions.
    Condition::NodesArrayType geom_nodes(3);
    for (unsigned i = 0; i < num_conds; ++i)
    {
        rConditionIds[i] = i + 1;
        rPropertiesIds[i] = 1;
        for (unsigned j = 0; j < 3; ++j)
        {
            rConnectivities(i, j) = i + j + 1;
            geom_nodes(j) = rNodes(i + j + 1);
        }
        auto p_cond = SurfaceCondition3D3N.Create(i + 1, geom_nodes, rProperties(1));
        rConditions.push_back(p_cond);
    }
}

} // namespace Testing
} // namespace Kratos.
