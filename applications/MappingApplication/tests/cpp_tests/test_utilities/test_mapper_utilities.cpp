//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos {
namespace Testing {

void CreateNodesForMapping(ModelPart& rModelPart, const int NumNodes)
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    const int size = rModelPart.GetCommunicator().TotalProcesses();

    const int start_id = NumNodes * rank + 1;

    // creating nodes with random coordinates
    for (int i=0; i< NumNodes; ++i)
        rModelPart.CreateNewNode(i+start_id, i*0.1*rank*size+0.134,
                                             i*0.2+rank*3.48*size,
                                             i*0.3*rank*6.13*size);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_AssignInterfaceEquationIds, KratosMappingApplicationSerialTestSuite)
{
    const int num_nodes = 11;
    ModelPart model_part("ForTest");

    CreateNodesForMapping(model_part, num_nodes);

    MapperUtilities::AssignInterfaceEquationIds(model_part.GetCommunicator());

    int idx = 0;

    for (const auto& r_node : model_part/*.GetCommunicator().LocalMesh()*/.Nodes())
    {
        KRATOS_CHECK_EQUAL(idx, r_node.GetValue(INTERFACE_EQUATION_ID));
        idx += 1;
    }
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_ComputeBoundingBox, KratosMappingApplicationSerialTestSuite)
{
    ModelPart model_part("ForTest");
    model_part.CreateNewNode(1, 0.2, 5.3, -8.3);
    model_part.CreateNewNode(2, 8.2, 25.3, 16.4);
    model_part.CreateNewNode(3, -9.2, -17.13, 1.5);
    model_part.CreateNewNode(4, 12.6, 5.3, -8.3);

    const auto bbox = MapperUtilities::ComputeLocalBoundingBox(model_part);

    // std::cout << MapperUtilities::BoundingBoxStringStream(bbox) << std::endl;

    KRATOS_CHECK_EQUAL(bbox.size(), 6);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[0], 12.6);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[1], -9.2);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[2], 25.3);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[3], -17.13);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[4], 16.4);
    KRATOS_CHECK_DOUBLE_EQUAL(bbox[5], -8.3);
}

KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_ComputeBoundingBoxWithTol, KratosMappingApplicationSerialTestSuite)
{
    std::vector<double> bboxes_wrong_size(5);
    std::vector<double> bboxes_with_tol;

    KRATOS_CHECK_EXCEPTION_IS_THROWN(MapperUtilities::ComputeBoundingBoxesWithTolerance(bboxes_wrong_size, 1.235, bboxes_with_tol),
        "Error: Bounding Boxes size has to be a multiple of 6!");

    // Cretae a vector containing the fake bboxes
    const std::size_t num_entries = 24;
    std::vector<double> bboxes(num_entries);

    const double factor = 1.2589;
    const double offset = 8.4;

    for (std::size_t i=0; i<num_entries; ++i)
        bboxes[i] = i*factor - offset;

    const double tolerance = 5.478;

    MapperUtilities::ComputeBoundingBoxesWithTolerance(bboxes,
                                                       tolerance,
                                                       bboxes_with_tol);

    for (std::size_t i=0; i<num_entries; i+=2)
        KRATOS_CHECK_DOUBLE_EQUAL(bboxes_with_tol[i], (i*factor-offset + tolerance));

    for (std::size_t i=1; i<num_entries; i+=2)
        KRATOS_CHECK_DOUBLE_EQUAL(bboxes_with_tol[i], (i*factor-offset - tolerance));
}




}  // namespace Testing
}  // namespace Kratos