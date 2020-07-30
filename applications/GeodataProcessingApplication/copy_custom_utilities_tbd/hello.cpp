//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//                   Simon Wenczowski
//

#include "includes/variables.h"
#include "hello.h"

// Project includes
#include "utilities/intersection_utilities.h"

#include <iostream>

namespace Kratos {

void Hello::Greet()
{
    KRATOS_INFO("Hello") << "A first greeting from the GeodataProcessingApplication: Hello World!" << std::endl;
}

// test function
void Hello::test_nicola()
{
    KRATOS_INFO("Nicola") << "test_nicola" << std::endl;
}

// function to check if a point is internal
/**
 * VolumeModelPart     (Elements3D4N)
 * GeometryModelPart   (Elements3D3N)
*/
ModelPart& Hello::CheckIfInternal( ModelPart& VolumeModelPart, ModelPart& GeometryModelPart ){

    /*
    // TEST THAT WORKS

    std::vector<std::size_t> index_node;
    for (int i = 0; i < static_cast<int>(VolumeModelPart.NumberOfNodes()); ++i){
        auto p_node = VolumeModelPart.NodesBegin() + i;
        index_node.push_back( p_node->Id() );
    }
    KRATOS_INFO("node") << index_node << std::endl;
    // return index_node;
    */

    std::vector<std::size_t> coord_z;
    for (auto& node : VolumeModelPart.Nodes()){
        coord_z.push_back(node.Coordinates()[2]);
    }
    const double max_z = *std::max_element(coord_z.begin(), coord_z.end());    // max value of z
    KRATOS_INFO("max_z") << max_z << std::endl;

    int num_nodes = VolumeModelPart.NumberOfNodes();
    KRATOS_INFO("num_nodes") << num_nodes << std::endl;

    /*
    // we set all node with +inf DISTANCE
    // UGLY... I can call an existing function
    for (auto& node : VolumeModelPart.Nodes()){
        node.GetSolutionStepValue(DISTANCE) = std::numeric_limits<double>::max();
    }
    */

    // just for test. DELETE it
    int num_neg_dist = 0;

    // loop in the node of Volume
    for (auto& node : VolumeModelPart.Nodes()){
        const array_1d<double,3> p_coords = node.Coordinates();
        
        double x = p_coords[0];
        double y = p_coords[1];
        double z = p_coords[2] + max_z;

        array_1d<double,3> t_coords;
        t_coords[0] = x; t_coords[1] = y; t_coords[2] = z;

        // Call the line-triangle intersection util
        array_1d<double,3> pnt_inter = ZeroVector(3);
        const double epsilon = 1e-12;

        // loop in the elements of Geometry
        for (auto& elem : GeometryModelPart.Elements()){

            // auto &r_geom = elem->GetGeometry();     // Geometry object
            auto &r_geom = elem.GetGeometry();     // Geometry object

            /**
             * is_intersected:
             *     -1 (the triangle is degenerate)
             *     0 (disjoint - no intersection)
             *     1 (intersect in a unique point)
             *     2 (are in the same plane)
            */
            const int is_intersected = IntersectionUtilities::ComputeTriangleLineIntersection(
                r_geom,
                p_coords,
                t_coords,
                pnt_inter,
                epsilon);
            
            if (is_intersected == 1)
            {
                num_neg_dist++;
                // we set the negative DISTANCE
                node.GetSolutionStepValue(DISTANCE) = -1;
                
                // KRATOS_INFO("pnt_inter") << pnt_inter << std::endl;
                // KRATOS_INFO("DISTANCE") << node.GetSolutionStepValue(DISTANCE) << std::endl;
            }
            
            // KRATOS_INFO("is_intersected") << is_intersected << std::endl;
        }

    }
    KRATOS_INFO("num_neg_dist") << num_neg_dist << std::endl;
    // return intersection;

    // std::cin.get();      // input
}

}
