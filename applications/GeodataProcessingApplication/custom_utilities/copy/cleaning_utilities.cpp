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
#include "containers/model.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "cleaning_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void CleaningUtilities::CleanIsolatedNodes(){

        const int initial_num = mrModelPart.Nodes().size();
        auto& r_nodes_array = mrModelPart.Nodes();
        const auto& r_elem_array = mrModelPart.Elements();

        // marking all nodes as "superfluous"
        #pragma omp parallel for
        for( int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node ){
            auto node = r_nodes_array.begin() + i_node;
            node->Set(TO_ERASE, true);
        }

        // saving the nodes that belong to an element
        #pragma omp parallel for
        for( int i_elem = 0; i_elem < static_cast<int>(r_elem_array.size()); ++i_elem ){
            const auto elem = r_elem_array.begin() + i_elem;
            auto& r_geom = elem->GetGeometry();

            for (unsigned int i = 0; i < r_geom.size(); ++i){
                r_geom[i].Set(TO_ERASE, false);
            }
        }

        mrModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        const int final_num = mrModelPart.Nodes().size();
        KRATOS_INFO("CleaningUtilities") << "In total " << (initial_num - final_num) <<" superfluous nodes were cleared" << std::endl;

    }


    ModelPart& CleaningUtilities::HardCopyBeforeSurfaceDiscretization( ModelPart& OriginalModelPart, ModelPart& NewModelPart ){

        auto& r_sub_model_part = NewModelPart.GetSubModelPart("AuxSubModelPart");
        Properties::Pointer p_prop = NewModelPart.pGetProperties(0);
        // Properties::Pointer p_prop = OriginalModelPart.pGetProperties(0);    // [NG]

        // copying every node to the auxiliary model part
        // #pragma omp parallel for
        for( int i_node = 0; i_node < static_cast<int>( OriginalModelPart.NumberOfNodes() ); ++i_node ){

            auto p_node = OriginalModelPart.NodesBegin() + i_node;
            auto node = NewModelPart.CreateNewNode( p_node->Id(),
                                                    p_node->Coordinates()[0],
                                                    p_node->Coordinates()[1],
                                                    p_node->Coordinates()[2] );
            // #pragma omp critical
            r_sub_model_part.AddNode( node );
        }

        // copying all the elements into the auxiliary part
        // #pragma omp parallel for
        for( int i_elem = 0; i_elem < static_cast<int>( OriginalModelPart.NumberOfElements() ); ++i_elem ){

            auto p_elem= OriginalModelPart.ElementsBegin() + i_elem;
            auto& r_geom = p_elem->GetGeometry();

            std::vector<ModelPart::IndexType> elem_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id() };
            auto elem = NewModelPart.CreateNewElement(  "Element3D4N",
                                                        p_elem->Id(),
                                                        elem_nodes,
                                                        p_prop );
            // #pragma omp critical
            r_sub_model_part.AddElement( elem );
        }


        /* [NG]
        // copying all the conditions into the auxiliary part
        // #pragma omp parallel for
        for( int i_cond = 0; i_cond < static_cast<int>( OriginalModelPart.NumberOfConditions() ); ++i_cond ){

            auto p_cond = OriginalModelPart.ConditionsBegin() + i_cond;
            auto& r_geom = p_cond->GetGeometry();

            std::vector<ModelPart::IndexType> cond_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id() };
            auto cond = NewModelPart.CreateNewCondition(  "SurfaceCondition3D3N",
                                                        p_cond->Id(),
                                                        cond_nodes,
                                                        p_prop );
            // // #pragma omp critical
            // r_sub_model_part.AddCondition( cond );
        }
        [NG] */


        // transferring the distance field
        #pragma omp parallel for
        for( int i_node = 0; i_node < static_cast<int>( OriginalModelPart.NumberOfNodes() ); ++i_node ){

            auto p_node_recv = NewModelPart.NodesBegin() + i_node;
            auto p_node_send = OriginalModelPart.NodesBegin() + i_node;

            p_node_recv->FastGetSolutionStepValue( DISTANCE ) = p_node_send->FastGetSolutionStepValue( DISTANCE );
        }

        return NewModelPart;
    }


    ModelPart& CleaningUtilities::HardCopyAfterSurfaceDiscretization( ModelPart& OriginalModelPart, ModelPart& NewModelPart ){

        auto& r_sub_model_part_bound = NewModelPart.GetSubModelPart("Complete_Boundary");
        auto& r_sub_model_part_fluid = NewModelPart.GetSubModelPart("Parts_Fluid");
        Properties::Pointer p_prop = NewModelPart.pGetProperties(0);

        std::vector<std::size_t> index_node;
        std::vector<std::size_t> index_element;
        std::vector<std::size_t> index_condition;

        // copying every node to the auxiliary model part
        // #pragma omp parallel for
        for( int i_node = 0; i_node < static_cast<int>( OriginalModelPart.NumberOfNodes() ); ++i_node ){

            auto p_node = OriginalModelPart.NodesBegin() + i_node;
            auto node = NewModelPart.CreateNewNode( p_node->Id(),
                                                    p_node->Coordinates()[0],
                                                    p_node->Coordinates()[1],
                                                    p_node->Coordinates()[2] );
            // #pragma omp critical
            index_node.push_back( p_node->Id() );
        }
        r_sub_model_part_fluid.AddNodes( index_node );


        // copying all the elements into the auxiliary part
        // #pragma omp parallel for
        for( int i_elem = 0; i_elem < static_cast<int>( OriginalModelPart.NumberOfElements() ); ++i_elem ){

            auto p_elem= OriginalModelPart.ElementsBegin() + i_elem;
            auto& r_geom = p_elem->GetGeometry();

            std::vector<ModelPart::IndexType> elem_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id() };
            auto elem = NewModelPart.CreateNewElement(   "Element3D4N",
                                                            p_elem->Id(),
                                                            elem_nodes,
                                                            p_prop );
            // #pragma omp critical
            index_element.push_back( p_elem->Id() );
        }
        r_sub_model_part_fluid.AddElements( index_element );

        // clearing the nodes vector
        index_node.clear();

        // copying all conditions
        // #pragma omp parallel for
        for( int i_cond = 0; i_cond < static_cast<int>( OriginalModelPart.NumberOfConditions() ); ++i_cond ){
            auto p_cond = OriginalModelPart.ConditionsBegin() + i_cond;
            auto& r_geom = p_cond->GetGeometry();

            std::vector<ModelPart::IndexType> cond_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id() };
            auto cond = NewModelPart.CreateNewCondition(  "Condition3D3N",
                                                            p_cond->Id(),
                                                            cond_nodes,
                                                            p_prop );

            // #pragma omp critical
            index_condition.push_back( p_cond->Id() );
            // #pragma omp critical
            index_node.insert(index_node.end(), cond_nodes.begin(), cond_nodes.end());
        }

        r_sub_model_part_bound.AddNodes( index_node );
        r_sub_model_part_bound.AddConditions( index_condition );

        return NewModelPart;
    }


    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const CleaningUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
