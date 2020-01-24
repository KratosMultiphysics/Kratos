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
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "testing/testing.h"

// Application includes
#include "custom_io/hdf5_model_part_io.h"
#include "custom_io/hdf5_nodal_flag_value_io.h"
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(HDF5NodalFlagValueIO_WriteNodalFlags1, KratosHDF5TestSuite)
{
    Parameters settings(R"({
        "prefix": "/Results",
        "list_of_variables": ["SLIP",
                              "ACTIVE",
                              "STRUCTURE"]
        })");
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part);
    TestModelPartFactory::AssignNonHistoricalNodalTestData(
        r_write_model_part, {{"SLIP"}, {"ACTIVE"}, {"STRUCTURE"}});
    auto p_file = pGetTestSerialFile();
    HDF5::ModelPartIO model_part_io(p_file, "/ModelData");
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    HDF5::NodalFlagValueIO nodal_value_io(settings, p_file);
    nodal_value_io.WriteNodalFlags(r_write_model_part.Nodes());
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    nodal_value_io.ReadNodalFlags(r_read_model_part.Nodes(),
                                  r_read_model_part.GetCommunicator());
    CompareNonHistoricalNodalData(r_read_model_part.Nodes(), r_write_model_part.Nodes());
}

} // namespace Testing
} // namespace Kratos.
