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
#include "containers/model.h"

// Application includes
#include "tests/test_utils.h"
#include "custom_io/hdf5_model_part_io.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadNodes, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part);
    KRATOS_CHECK(r_write_model_part.NumberOfNodes() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    CompareNodes(r_read_model_part.Nodes(), r_write_model_part.Nodes());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements1, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part, {{"Element2D3N"}});
    KRATOS_CHECK(r_write_model_part.NumberOfElements() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    model_part_io.WriteElements(r_write_model_part.Elements());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    model_part_io.ReadElements(r_read_model_part.Nodes(), r_read_model_part.rProperties(), r_read_model_part.Elements());
    CompareElements(r_read_model_part.Elements(), r_write_model_part.Elements());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements2, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part,
                                          {{"Element2D3N"}, {"Element2D4N"}});
    KRATOS_CHECK(r_write_model_part.NumberOfElements() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    model_part_io.WriteElements(r_write_model_part.Elements());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    model_part_io.ReadElements(r_read_model_part.Nodes(), r_read_model_part.rProperties(), r_read_model_part.Elements());
    CompareElements(r_read_model_part.Elements(), r_write_model_part.Elements());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadElements3, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part, {}, {});
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    model_part_io.WriteElements(r_write_model_part.Elements());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    model_part_io.ReadElements(r_read_model_part.Nodes(), r_read_model_part.rProperties(), r_read_model_part.Elements());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions1, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part, {}, {{"SurfaceCondition3D3N"}});
    KRATOS_CHECK(r_write_model_part.NumberOfConditions() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    model_part_io.WriteConditions(r_write_model_part.Conditions());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    model_part_io.ReadConditions(r_read_model_part.Nodes(), r_read_model_part.rProperties(), r_read_model_part.Conditions());
    CompareConditions(r_read_model_part.Conditions(), r_write_model_part.Conditions());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions2, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(
        r_write_model_part, {},
        {{"SurfaceCondition3D3N"}, {"SurfaceCondition3D4N"}});
    KRATOS_CHECK(r_write_model_part.NumberOfConditions() > 0);
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    model_part_io.WriteConditions(r_write_model_part.Conditions());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    model_part_io.ReadConditions(r_read_model_part.Nodes(), r_read_model_part.rProperties(), r_read_model_part.Conditions());
    CompareConditions(r_read_model_part.Conditions(), r_write_model_part.Conditions());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadConditions3, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part, {}, {});
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    model_part_io.WriteConditions(r_write_model_part.Conditions());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    model_part_io.ReadConditions(r_read_model_part.Nodes(), r_read_model_part.rProperties(), r_read_model_part.Conditions());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_Properties, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    HDF5::PropertiesContainerType& r_write_properties = r_write_model_part.rProperties();
    TestModelPartFactory::AssignDataValueContainer(r_write_properties[1].Data(),
                                                   {{"DOMAIN_SIZE"}});
    TestModelPartFactory::AssignDataValueContainer(r_write_properties[3].Data(),
                                                   {{"TIME"}, {"STRAIN"}});
    TestModelPartFactory::AssignDataValueContainer(r_write_properties[4].Data(),
                                                   {{"LOCAL_AXES_MATRIX"}});
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/Step");
    model_part_io.WriteProperties(r_write_properties);
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    HDF5::PropertiesContainerType& r_read_properties = r_read_model_part.rProperties();
    model_part_io.ReadProperties(r_read_properties);
    KRATOS_CHECK(r_read_model_part.NumberOfProperties() == r_write_model_part.NumberOfProperties());
    CompareDataValueContainers(r_read_properties[1].Data(),
                               r_write_properties[1].Data());
    CompareDataValueContainers(r_read_properties[3].Data(),
                               r_write_properties[3].Data());
    CompareDataValueContainers(r_read_properties[4].Data(),
                               r_write_properties[4].Data());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_ModelPartIO_ReadModelPart1, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_model_part_out = this_model.CreateModelPart("test_out");
    TestModelPartFactory::CreateModelPart(r_model_part_out, {{"Element2D3N"}},
                                          {{"SurfaceCondition3D3N"}});
    HDF5::ModelPartIO model_part_io(pGetTestSerialFile(), "/ModelData");
    model_part_io.WriteModelPart(r_model_part_out);
    ModelPart& r_model_part_in = this_model.CreateModelPart("test_in");
    model_part_io.ReadModelPart(r_model_part_in);
    CompareModelParts(r_model_part_in, r_model_part_out);
}

} // namespace Testing
} // namespace Kratos.
