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

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"

// Application includes
#include "tests/test_utils.h"
#include "custom_io/hdf5_model_part_io.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadNodes, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    TestModelPartFactory::CreateModelPart(write_model_part);
    KRATOS_CHECK(write_model_part.NumberOfNodes() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    CompareNodes(read_model_part.Nodes(), write_model_part.Nodes());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements1, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    TestModelPartFactory::CreateModelPart(write_model_part, {{"Element2D3N"}});
    KRATOS_CHECK(write_model_part.NumberOfElements() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteElements(write_model_part.Elements());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadElements(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Elements());
    CompareElements(read_model_part.Elements(), write_model_part.Elements());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements2, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    TestModelPartFactory::CreateModelPart(write_model_part,
                                          {{"Element2D3N"}, {"Element2D4N"}});
    KRATOS_CHECK(write_model_part.NumberOfElements() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteElements(write_model_part.Elements());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadElements(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Elements());
    CompareElements(read_model_part.Elements(), write_model_part.Elements());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements3, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    TestModelPartFactory::CreateModelPart(write_model_part, {}, {});
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteElements(write_model_part.Elements());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadElements(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Elements());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions1, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    TestModelPartFactory::CreateModelPart(write_model_part, {}, {{"SurfaceCondition3D3N"}});
    KRATOS_CHECK(write_model_part.NumberOfConditions() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteConditions(write_model_part.Conditions());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadConditions(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Conditions());
    CompareConditions(read_model_part.Conditions(), write_model_part.Conditions());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions2, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    TestModelPartFactory::CreateModelPart(
        write_model_part, {},
        {{"SurfaceCondition3D3N"}, {"SurfaceCondition3D4N"}});
    KRATOS_CHECK(write_model_part.NumberOfConditions() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteConditions(write_model_part.Conditions());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadConditions(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Conditions());
    CompareConditions(read_model_part.Conditions(), write_model_part.Conditions());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions3, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    TestModelPartFactory::CreateModelPart(write_model_part, {}, {});
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(write_model_part.Nodes());
    model_part_io.WriteConditions(write_model_part.Conditions());
    ModelPart read_model_part("test_read");
    model_part_io.ReadNodes(read_model_part.Nodes());
    model_part_io.ReadConditions(read_model_part.Nodes(), read_model_part.rProperties(), read_model_part.Conditions());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_Properties, KratosHDF5TestSuite)
{
    ModelPart write_model_part("test_write");
    ModelPart read_model_part("test_read");
    HDF5::PropertiesContainerType& r_write_properties = write_model_part.rProperties();
    HDF5::PropertiesContainerType& r_read_properties = read_model_part.rProperties();
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");

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

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadModelPart1, KratosHDF5TestSuite)
{
    ModelPart model_part_out;
    TestModelPartFactory::CreateModelPart(model_part_out, {{"Element2D3N"}},
                                          {{"SurfaceCondition3D3N"}});
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/ModelData");
    model_part_io.WriteModelPart(model_part_out);
    ModelPart model_part_in;
    model_part_io.ReadModelPart(model_part_in);
    CompareModelParts(model_part_in, model_part_out);
}

} // namespace Testing
} // namespace Kratos.
