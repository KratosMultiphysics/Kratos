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
#include "custom_io/hdf5_container_component_io.h"
#include "custom_utilities/container_io_utils.h"
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(HDF5NodalFlagValueIO_WriteNodalFlags1, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    TestModelPartFactory::CreateModelPart(r_write_model_part);

    std::vector<std::string> variables_list {
        "SLIP",
        "ACTIVE",
        "STRUCTURE"
    };

    // "shuffle" the list of variables to check whether it's handled
    // without deadlocks.
    std::rotate(
        variables_list.begin(),
        variables_list.begin() + (r_read_model_part.GetCommunicator().GetDataCommunicator().Rank() % variables_list.size()),
        variables_list.end()
    );

    Parameters settings(R"({
        "prefix": "/Results",
        "list_of_variables": []
    })");
    settings["list_of_variables"].SetStringArray(variables_list);

    TestModelPartFactory::AssignNonHistoricalNodalTestData(
        r_write_model_part, variables_list);
    auto p_file = pGetTestSerialFile();
    HDF5::ModelPartIO model_part_io(Parameters(R"({"prefix":"/ModelData"})"), p_file);
    model_part_io.WriteNodes(r_write_model_part.Nodes());
    HDF5::ContainerComponentIO<ModelPart::NodesContainerType, HDF5::Internals::FlagIO, Flags> nodal_value_io(settings, p_file);
    nodal_value_io.Write(r_write_model_part.Nodes(), HDF5::Internals::FlagIO{}, Parameters("""{}"""));
    model_part_io.ReadNodes(r_read_model_part.Nodes());
    nodal_value_io.Read(r_read_model_part.Nodes(), HDF5::Internals::FlagIO{}, r_read_model_part.GetCommunicator());
    CompareNonHistoricalNodalData({}, r_read_model_part.Nodes(), r_write_model_part.Nodes());
}

} // namespace Testing
} // namespace Kratos.
