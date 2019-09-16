//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Simon Wenczowski
//
//

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "containers/model.h"

// Application includes
#include "extrusion_height_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void ExtrusionHeightUtilities::SetExtrusionHeight( const double radius, const double max_height, const double free_board ){

        int delete_counter = 0;
        // brute-force iteration over all nodes in model part
        for ( int i_node = 0; i_node < static_cast<int>( mrModelPart.NumberOfNodes() ); ++i_node )
        {
            auto node_center = mrModelPart.NodesBegin() + i_node;
            array_1d<double,3>& coord_center = node_center->Coordinates();
            double min_height = 1.0e10;

            #pragma omp parallel for reduction(min:min_height)
            for ( int j_node = 0; j_node < static_cast<int>( mrModelPart.NumberOfNodes() ); ++j_node ){

                // nothhing must be modified about the surrounding nodes
                const auto node_other = mrModelPart.NodesBegin() + j_node;
                const array_1d<double,3>& coord_other = node_other->Coordinates();
                const double dist = std::sqrt( (coord_center[0]-coord_other[0]) * (coord_center[0]-coord_other[0])
                                             + (coord_center[1]-coord_other[1]) * (coord_center[1]-coord_other[1]));

                if( dist < radius && coord_other[2] < min_height ){

                    min_height = coord_other[2];
                }
            }
            // computation of the desired height of the topper of the domain
            const double extrusion_height = min_height + max_height;

            if( coord_center[2] > extrusion_height - free_board ){
                // the node itself is located at a higher level than the desired height (*)
                node_center->Set( TO_ERASE, true );
                delete_counter++;

            } else {
                // the node is assigned an extrusion height based on the lowest neighbor
                node_center->Set( TO_ERASE, false );
                node_center->SetValue( EXTRUSION_HEIGHT, extrusion_height );
            }
        }
        // remove the nodes that were above the level (*)
        KRATOS_INFO_IF("ExtrusionHeightUtilities", delete_counter > 0) << "In total " << delete_counter << " points are deleted."  << std::endl;
        mrModelPart.RemoveNodes( TO_ERASE );
    }


    void ExtrusionHeightUtilities::SmoothExtrusionHeight( const double radius, const int iterations, const double free_board ){

        for (int iter = 0; iter < iterations; ++iter){

            // brute-force iteration over all nodes in model part
            for ( int i_node = 0; i_node < static_cast<int>( mrModelPart.NumberOfNodes() ); ++i_node )
            {
                auto node_center = mrModelPart.NodesBegin() + i_node;
                array_1d<double,3>& r_coord_center = node_center->Coordinates();
                double& r_ext_h_center = node_center->GetValue( EXTRUSION_HEIGHT );

                int neighbor_counter = 0;
                double total_height = 0.0;

                #pragma omp parallel for reduction(+:total_height,neighbor_counter)
                for ( int j_node = 0; j_node < static_cast<int>( mrModelPart.NumberOfNodes() ); ++j_node ){

                    // nothhing must be modified about the surrounding nodes
                    const auto node_other = mrModelPart.NodesBegin() + j_node;
                    const array_1d<double,3>& r_coord_other = node_other->Coordinates();
                    const double dist = std::sqrt( (r_coord_center[0]-r_coord_other[0]) * (r_coord_center[0]-r_coord_other[0])
                                                 + (r_coord_center[1]-r_coord_other[1]) * (r_coord_center[1]-r_coord_other[1]));

                    if( dist < radius ){
                        neighbor_counter++;
                        total_height += node_other->GetValue( EXTRUSION_HEIGHT );
                    }
                }

                if ( neighbor_counter > 0 ){
                    const double average_neighbors = total_height / static_cast<double>(neighbor_counter);
                    if ( average_neighbors > ( r_coord_center[2] + free_board ) ){
                        // smoothing the filed with its neighbor values
                        r_ext_h_center = 0.5 * r_ext_h_center + 0.5 * average_neighbors;
                    }
                }

            }
        }
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ExtrusionHeightUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
