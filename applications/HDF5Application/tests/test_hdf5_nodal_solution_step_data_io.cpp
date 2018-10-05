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
#include "includes/kratos_parameters.h"
#include "includes/communicator.h"

// Application includes
#include "tests/test_utils.h"
#include "custom_io/hdf5_file_serial.h"
#include "custom_io/hdf5_nodal_solution_step_data_io.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_ReadNodalResults2, KratosHDF5TestSuite)
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    auto p_test_file = Kratos::make_shared<HDF5::FileSerial>(file_params);

    Model this_model;
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
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
    HDF5::NodalSolutionStepDataIO data_io(io_params, p_test_file);
    data_io.WriteNodalResults(r_write_model_part.Nodes());
    data_io.ReadNodalResults(r_read_model_part.Nodes(), r_read_model_part.GetCommunicator());
    for (unsigned i = 0; i < 15; ++i)
    {
        HDF5::NodeType& r_read_node = r_read_model_part.Nodes()[i + 1];
        HDF5::NodeType& r_write_node = r_write_model_part.Nodes()[i + 1];
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_X) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_X));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Y) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Y));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Z) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Z));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(PRESSURE) ==
                     r_write_node.FastGetSolutionStepValue(PRESSURE));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(REFINEMENT_LEVEL) ==
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
    auto p_test_file = Kratos::make_shared<HDF5::FileSerial>(file_params);

    Model this_model;
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
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
    HDF5::NodalSolutionStepDataIO data_io(io_params, p_test_file);
    data_io.WriteNodalResults(r_write_model_part.Nodes());
    data_io.ReadNodalResults(r_read_model_part.Nodes(), r_read_model_part.GetCommunicator());
    for (unsigned i = 0; i < 15; ++i)
    {
        HDF5::NodeType& r_read_node = r_read_model_part.Nodes()[i + 1];
        HDF5::NodeType& r_write_node = r_write_model_part.Nodes()[i + 1];
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_X) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_X));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Y) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Y));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(DISPLACEMENT_Z) ==
                     r_write_node.FastGetSolutionStepValue(DISPLACEMENT_Z));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(PRESSURE) ==
                     r_write_node.FastGetSolutionStepValue(PRESSURE));
        KRATOS_CHECK(r_read_node.FastGetSolutionStepValue(REFINEMENT_LEVEL) ==
                     r_write_node.FastGetSolutionStepValue(REFINEMENT_LEVEL));
    }
}

} // namespace Testing
} // namespace Kratos.
