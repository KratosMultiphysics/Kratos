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

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_PointsData1, KratosHDF5TestSuite)
{
    Model this_model;
    ModelPart& r_test_model_part = this_model.CreateModelPart("TestModelPart");
    TestModelPartFactory::CreateModelPart(r_test_model_part);
    KRATOS_CHECK(r_test_model_part.NumberOfNodes() > 0);
    HDF5::Internals::PointsData data;
    data.SetData(r_test_model_part.Nodes());
    auto test_file = GetTestSerialFile();
    HDF5::WriteInfo info;
    data.WriteData(test_file, "/Nodes", info);
    data.Clear();
    KRATOS_CHECK(data.size() == 0);
    data.ReadData(test_file, "/Nodes", info.StartIndex, info.BlockSize);
    HDF5::NodesContainerType new_nodes;
    data.CreateNodes(new_nodes);
    CompareNodes(new_nodes, r_test_model_part.Nodes());
}

} // namespace Testing
} // namespace Kratos.
