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
#include "includes/kratos_parameters.h"
#include "includes/communicator.h"

// Application includes
#include "tests/test_utils.h"
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_container_component_io.h"
#include "custom_utilities/container_io_utils.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_ReadNodalResults2, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");

    std::vector<std::string> variables_list {
        "DISPLACEMENT",
        "PRESSURE",
        "REFINEMENT_LEVEL"
    };

    // "shuffle" the list of variables to check whether it's handled
    // without deadlocks.
    std::rotate(
        variables_list.begin(),
        variables_list.begin() + (r_read_model_part.GetCommunicator().GetDataCommunicator().Rank() % variables_list.size()),
        variables_list.end()
    );

    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    auto p_test_file = Kratos::make_shared<HDF5::File>(r_read_model_part.GetCommunicator().GetDataCommunicator(), file_params);

    r_read_model_part.AddNodalSolutionStepVariable(DISPLACEMENT); // array_1d
    r_read_model_part.AddNodalSolutionStepVariable(PRESSURE); // double
    r_read_model_part.AddNodalSolutionStepVariable(REFINEMENT_LEVEL); // int
    r_write_model_part.AddNodalSolutionStepVariable(DISPLACEMENT); // array_1d
    r_write_model_part.AddNodalSolutionStepVariable(PRESSURE); // double
    r_write_model_part.AddNodalSolutionStepVariable(REFINEMENT_LEVEL); // int
    for (unsigned i = 0; i < 15; ++i)
    {
        r_read_model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
        r_write_model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
    }
    r_read_model_part.SetBufferSize(2);
    r_write_model_part.SetBufferSize(2);

    for (auto& r_node : r_write_model_part.Nodes())
    {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>(3, r_node.Id() + 0.2345);
        r_node.FastGetSolutionStepValue(PRESSURE) = r_node.Id() + 0.2345;
        r_node.FastGetSolutionStepValue(REFINEMENT_LEVEL) = r_node.Id() + 4;
    }

    Parameters io_params(R"({
        "prefix": "/Step",
        "list_of_variables": []
    })");
    io_params["list_of_variables"].SetStringArray(variables_list);

    HDF5::ContainerComponentIO<ModelPart::NodesContainerType, HDF5::Internals::HistoricalIO, Variable<int>, Variable<double>, Variable<array_1d<double, 3>>, Variable<array_1d<double, 4>>, Variable<array_1d<double, 6>>, Variable<array_1d<double, 9>>, Variable<Kratos::Vector>, Variable<Kratos::Matrix>> data_io(io_params, p_test_file);
    data_io.Write(r_write_model_part.Nodes(), HDF5::Internals::HistoricalIO(0), Parameters("""{}"""));
    data_io.Read(r_read_model_part.Nodes(), HDF5::Internals::HistoricalIO(0), r_read_model_part.GetCommunicator());

    for (unsigned i = 0; i < 15; ++i)
    {
        HDF5::NodeType& r_read_node = r_read_model_part.Nodes()[i + 1];
        HDF5::NodeType& r_write_node = r_write_model_part.Nodes()[i + 1];
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_X) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_X));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Y) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Y));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Z) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Z));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(PRESSURE) ==
                     r_write_node.FastGetSolutionStepValue(PRESSURE));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(REFINEMENT_LEVEL) ==
                     r_write_node.FastGetSolutionStepValue(REFINEMENT_LEVEL));
    }
}

KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_ReadNodalResults, KratosHDF5TestSuite)
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");

    Model this_model;
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");

    auto p_test_file = Kratos::make_shared<HDF5::File>(r_read_model_part.GetCommunicator().GetDataCommunicator(), file_params);

    r_read_model_part.AddNodalSolutionStepVariable(DISPLACEMENT); // array_1d
    r_read_model_part.AddNodalSolutionStepVariable(PRESSURE); // double
    r_read_model_part.AddNodalSolutionStepVariable(REFINEMENT_LEVEL); // int
    r_write_model_part.AddNodalSolutionStepVariable(DISPLACEMENT); // array_1d
    r_write_model_part.AddNodalSolutionStepVariable(PRESSURE); // double
    r_write_model_part.AddNodalSolutionStepVariable(REFINEMENT_LEVEL); // int
    for (unsigned i = 0; i < 15; ++i)
    {
        r_read_model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
        r_write_model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
    }
    r_read_model_part.SetBufferSize(2);
    r_write_model_part.SetBufferSize(2);

    for (auto& r_node : r_write_model_part.Nodes())
    {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>(3, r_node.Id() + 0.2345);
        r_node.FastGetSolutionStepValue(PRESSURE) = r_node.Id() + 0.2345;
        r_node.FastGetSolutionStepValue(REFINEMENT_LEVEL) = r_node.Id() + 4;
    }

    Parameters io_params(R"(
        {
            "prefix": "/Step",
            "list_of_variables": ["DISPLACEMENT", "PRESSURE", "REFINEMENT_LEVEL"]
        })");
    HDF5::ContainerComponentIO<ModelPart::NodesContainerType, HDF5::Internals::HistoricalIO, Variable<int>, Variable<double>, Variable<array_1d<double, 3>>, Variable<array_1d<double, 4>>, Variable<array_1d<double, 6>>, Variable<array_1d<double, 9>>, Variable<Kratos::Vector>, Variable<Kratos::Matrix>> data_io(io_params, p_test_file);
    data_io.Write(r_write_model_part.Nodes(), HDF5::Internals::HistoricalIO(0), Parameters("""{}"""));
    data_io.Read(r_read_model_part.Nodes(), HDF5::Internals::HistoricalIO(0), r_read_model_part.GetCommunicator());
    for (unsigned i = 0; i < 15; ++i)
    {
        HDF5::NodeType& r_read_node = r_read_model_part.Nodes()[i + 1];
        HDF5::NodeType& r_write_node = r_write_model_part.Nodes()[i + 1];
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_X) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_X));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Y) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Y));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Z) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Z));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(PRESSURE) ==
                     r_write_node.FastGetSolutionStepValue(PRESSURE));
        KRATOS_EXPECT_TRUE(r_read_node.FastGetSolutionStepValue(REFINEMENT_LEVEL) ==
                     r_write_node.FastGetSolutionStepValue(REFINEMENT_LEVEL));
    }
}

} // namespace Testing
} // namespace Kratos.
