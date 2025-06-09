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
#include "custom_io/hdf5_points_data.h"
#include "custom_utilities/container_io_utils.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_PointsData1, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_test_model_part = this_model.CreateModelPart("TestModelPart");
    TestModelPartFactory::CreateModelPart(r_test_model_part);
    KRATOS_EXPECT_TRUE(r_test_model_part.NumberOfNodes() > 0);

    HDF5::Internals::PointsData<HDF5::Internals::NodesIO> data("/Nodes", pGetTestSerialFile());
    data.Write(r_test_model_part.Nodes(), HDF5::Internals::NodesIO{}, Parameters(R"({})"));

    HDF5::NodesContainerType new_nodes;
    data.Read(new_nodes, HDF5::Internals::NodesIO{});
    CompareNodes(new_nodes, r_test_model_part.Nodes());
}

} // namespace Testing
} // namespace Kratos.
