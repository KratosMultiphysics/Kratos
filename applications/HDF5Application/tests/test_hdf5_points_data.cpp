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

// Application includes
#include "custom_io/hdf5_points_data.h"
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_PointsData1, KratosHDF5TestSuite)
{
    ModelPart test_model_part;
    TestModelPartFactory::CreateModelPart(test_model_part);
    KRATOS_CHECK(test_model_part.NumberOfNodes() > 0);
    HDF5::Internals::PointsData data;
    data.SetData(test_model_part.Nodes());
    HDF5::File::Pointer p_file = pGetFile();
    HDF5::WriteInfo info;
    data.WriteData(*p_file, "/Nodes", info);
    data.Clear();
    KRATOS_CHECK(data.size() == 0);
    data.ReadData(*p_file, "/Nodes", info.StartIndex, info.BlockSize);
    HDF5::NodesContainerType new_nodes;
    data.CreateNodes(new_nodes);
    CompareNodes(new_nodes, test_model_part.Nodes());
}

} // namespace Testing
} // namespace Kratos.
