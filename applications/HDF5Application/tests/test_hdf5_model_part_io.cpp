//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

// System includes
#include <sstream>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_components.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_io/hdf5_file_serial.h"
#include "custom_io/hdf5_model_part_io.h"
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadNodes, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    const unsigned num_nodes = 10;
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        double xyz = i;
        write_model_part.CreateNewNode(i + 1, xyz, xyz, xyz);
    }
    HDF5::File::Pointer p_file = pGetFile();
    HDF5::ModelPartIO model_part_io(p_file, "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    CompareModelParts(read_model_part, write_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements1, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    const unsigned num_elems = 5;
    const unsigned num_nodes = num_elems + 2;
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        double xyz = i;
        write_model_part.CreateNewNode(i + 1, xyz, xyz, xyz);
    }
    const HDF5::ElementType& Element2D3N = KratosComponents<HDF5::ElementType>::Get("Element2D3N");
    HDF5::ElementsContainerType& r_elems = write_model_part.Elements();
    HDF5::ElementType::NodesArrayType geom_nodes(3);
    ModelPart::PropertiesType::Pointer p_prop = write_model_part.pGetProperties(1);
    for (unsigned i = 0; i < num_elems; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
            geom_nodes(j) = write_model_part.pGetNode(i + j + 1);
        r_elems.push_back(Element2D3N.Create(i + 1, geom_nodes, p_prop));
    }
    HDF5::File::Pointer p_file = pGetFile();
    HDF5::ModelPartIO model_part_io(p_file, "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteElements(write_model_part.Elements());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadElements(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Elements());
    CompareModelParts(read_model_part, write_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements2, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    const unsigned num_tri_elems = 10;
    const unsigned num_quad_elems = 15;
    const unsigned num_nodes = (num_tri_elems + 2) + (num_quad_elems + 2);
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        double xyz = i;
        write_model_part.CreateNewNode(i + 1, xyz, xyz, xyz);
    }
    const HDF5::ElementType& Element2D3N = KratosComponents<HDF5::ElementType>::Get("Element2D3N");
    const HDF5::ElementType& Element2D4N = KratosComponents<HDF5::ElementType>::Get("Element2D4N");
    HDF5::ElementsContainerType& r_elems = write_model_part.Elements();
    HDF5::ElementType::NodesArrayType tri_nodes(3), quad_nodes(4);
    ModelPart::PropertiesType::Pointer p_prop = write_model_part.pGetProperties(1);
    for (unsigned i = 0; i < num_tri_elems; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
            tri_nodes(j) = write_model_part.pGetNode(i + j + 1);
        r_elems.push_back(Element2D3N.Create(i + 1, tri_nodes, p_prop));
    }
    for (unsigned i = 0; i < num_quad_elems; ++i)
    {
        for (unsigned j = 0; j < 4; ++j)
            quad_nodes(j) = write_model_part.pGetNode(num_tri_elems + i + j + 1);
        r_elems.push_back(Element2D4N.Create(num_tri_elems + i + 1, quad_nodes, p_prop));
    }
    HDF5::File::Pointer p_file = pGetFile();
    HDF5::ModelPartIO model_part_io(p_file, "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteElements(write_model_part.Elements());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadElements(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Elements());
    CompareModelParts(read_model_part, write_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions1, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    const unsigned num_conds = 5;
    const unsigned num_nodes = num_conds + 2;
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        double xyz = i;
        write_model_part.CreateNewNode(i + 1, xyz, xyz, xyz);
    }
    const HDF5::ConditionType& SurfaceCondition3D3N = KratosComponents<HDF5::ConditionType>::Get("SurfaceCondition3D3N");
    HDF5::ConditionsContainerType& r_conds = write_model_part.Conditions();
    HDF5::ConditionType::NodesArrayType geom_nodes(3);
    ModelPart::PropertiesType::Pointer p_prop = write_model_part.pGetProperties(1);
    for (unsigned i = 0; i < num_conds; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
            geom_nodes(j) = write_model_part.pGetNode(i + j + 1);
        r_conds.push_back(SurfaceCondition3D3N.Create(i + 1, geom_nodes, p_prop));
    }
    HDF5::File::Pointer p_file = pGetFile();
    HDF5::ModelPartIO model_part_io(p_file, "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteConditions(write_model_part.Conditions());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadConditions(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Conditions());
    CompareModelParts(read_model_part, write_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions2, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    const unsigned num_tri_conds = 10;
    const unsigned num_quad_conds = 15;
    const unsigned num_nodes = (num_tri_conds + 2) + (num_quad_conds + 2);
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        double xyz = i;
        write_model_part.CreateNewNode(i + 1, xyz, xyz, xyz);
    }
    const HDF5::ConditionType& SurfaceCondition3D3N = KratosComponents<HDF5::ConditionType>::Get("SurfaceCondition3D3N");
    const HDF5::ConditionType& SurfaceCondition3D4N = KratosComponents<HDF5::ConditionType>::Get("SurfaceCondition3D4N");
    HDF5::ConditionsContainerType& r_conds = write_model_part.Conditions();
    HDF5::ConditionType::NodesArrayType tri_nodes(3), quad_nodes(4);
    ModelPart::PropertiesType::Pointer p_prop = write_model_part.pGetProperties(1);
    for (unsigned i = 0; i < num_tri_conds; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
            tri_nodes(j) = write_model_part.pGetNode(i + j + 1);
        r_conds.push_back(SurfaceCondition3D3N.Create(i + 1, tri_nodes, p_prop));
    }
    for (unsigned i = 0; i < num_quad_conds; ++i)
    {
        for (unsigned j = 0; j < 4; ++j)
            quad_nodes(j) = write_model_part.pGetNode(num_tri_conds + i + j + 1);
        r_conds.push_back(SurfaceCondition3D4N.Create(num_tri_conds + i + 1, quad_nodes, p_prop));
    }
    HDF5::File::Pointer p_file = pGetFile();
    HDF5::ModelPartIO model_part_io(p_file, "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteConditions(write_model_part.Conditions());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadConditions(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Conditions());
    CompareModelParts(read_model_part, write_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_Properties, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    ModelPart read_model_part("test_read");
    HDF5::PropertiesContainerType& r_write_properties = write_model_part.rProperties();
    HDF5::PropertiesContainerType& r_read_properties = read_model_part.rProperties();
    HDF5::File::Pointer p_file = pGetFile();
    HDF5::ModelPartIO model_part_io(p_file, "/Step");

    r_write_properties[1][DOMAIN_SIZE] = 2;
    r_write_properties[3][TIME] = 1.2345;
    r_write_properties[3][STRAIN] = HDF5::Vector<double>(3, 1.234567);
    r_write_properties[4][LOCAL_AXES_MATRIX] = HDF5::Matrix<double>(3, 3, 1.23);
    model_part_io.WriteProperties(r_write_properties);
    model_part_io.ReadProperties(r_read_properties);
    KRATOS_CHECK(read_model_part.NumberOfProperties() == write_model_part.NumberOfProperties());
    KRATOS_CHECK(r_read_properties[1][DOMAIN_SIZE] == r_write_properties[1][DOMAIN_SIZE]);
    KRATOS_CHECK(r_read_properties[3][TIME] == r_write_properties[3][TIME]);
    KRATOS_CHECK(r_read_properties[3][STRAIN].size() == r_write_properties[3][STRAIN].size());
    for (unsigned i = 0; i < r_read_properties[3][STRAIN].size(); ++i)
        KRATOS_CHECK(r_read_properties[3][STRAIN](i) == r_write_properties[3][STRAIN](i));
    KRATOS_CHECK(r_read_properties[4][LOCAL_AXES_MATRIX].size1() == r_write_properties[4][LOCAL_AXES_MATRIX].size1());
    KRATOS_CHECK(r_read_properties[4][LOCAL_AXES_MATRIX].size2() == r_write_properties[4][LOCAL_AXES_MATRIX].size2());
    for (unsigned i = 0; i < r_read_properties[4][LOCAL_AXES_MATRIX].size1(); ++i)
        for (unsigned j = 0; j < r_read_properties[4][LOCAL_AXES_MATRIX].size2(); ++j)
            KRATOS_CHECK(r_read_properties[4][LOCAL_AXES_MATRIX](i, j) == r_write_properties[4][LOCAL_AXES_MATRIX](i, j));
}

} // namespace Testing
} // namespace Kratos.
