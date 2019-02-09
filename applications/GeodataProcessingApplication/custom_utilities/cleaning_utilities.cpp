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


    ModelPart& CleaningUtilities::ReBuildModelPart(){

        Model& owner_model = mrModelPart.GetModel();
        auto& aux_model_part = owner_model.CreateModelPart("AuxModelPartForCopy");
        Properties::Pointer p_prop = aux_model_part.pGetProperties(0);

        // copying everything to the auxiliary model part
        for( int i_node = 0; i_node < static_cast<int>( mrModelPart.NumberOfNodes() ); ++i_node ){
            auto p_node = mrModelPart.NodesBegin() + i_node;

            auto node = aux_model_part.CreateNewNode(   p_node->Id(),
                                                        p_node->Coordinates()[0],
                                                        p_node->Coordinates()[1],
                                                        p_node->Coordinates()[2] );
        }

        for( int i_cond = 0; i_cond < static_cast<int>( mrModelPart.NumberOfConditions() ); ++i_cond ){
            auto p_cond = mrModelPart.ConditionsBegin() + i_cond;
            auto& r_geom = p_cond->GetGeometry();

            std::vector<ModelPart::IndexType> cond_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id() };
            auto cond = aux_model_part.CreateNewCondition(  "Condition3D3N",
                                                            p_cond->Id(),
                                                            cond_nodes,
                                                            p_prop );
        }

        for( int i_elem = 0; i_elem < static_cast<int>( mrModelPart.NumberOfElements() ); ++i_elem ){
            auto p_elem= mrModelPart.ElementsBegin() + i_elem;
            auto& r_geom = p_elem->GetGeometry();

            std::vector<ModelPart::IndexType> elem_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id() };
            auto elem = aux_model_part.CreateNewElement(   "Element3D4N",
                                                            p_elem->Id(),
                                                            elem_nodes,
                                                            p_prop );
        }

        // delete and re-create
        std::string model_part_name = mrModelPart.Name();
        owner_model.DeleteModelPart( model_part_name );

        ModelPart& new_model_part = owner_model.CreateModelPart( model_part_name );

        // copying everything to the newly generated model part
        for( int i_node = 0; i_node < static_cast<int>( aux_model_part.NumberOfNodes() ); ++i_node ){
            auto p_node = aux_model_part.NodesBegin() + i_node;

            auto node = new_model_part.CreateNewNode(   p_node->Id(),
                                                        p_node->Coordinates()[0],
                                                        p_node->Coordinates()[1],
                                                        p_node->Coordinates()[2] );
        }

        for( int i_cond = 0; i_cond < static_cast<int>( aux_model_part.NumberOfConditions() ); ++i_cond ){
            auto p_cond = aux_model_part.ConditionsBegin() + i_cond;
            auto& r_geom = p_cond->GetGeometry();

            std::vector<ModelPart::IndexType> cond_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id() };
            auto cond = new_model_part.CreateNewCondition(  "Condition3D3N",
                                                            p_cond->Id(),
                                                            cond_nodes,
                                                            p_prop );
        }

        for( int i_elem = 0; i_elem < static_cast<int>( aux_model_part.NumberOfElements() ); ++i_elem ){
            auto p_elem= aux_model_part.ElementsBegin() + i_elem;
            auto& r_geom = p_elem->GetGeometry();

            std::vector<ModelPart::IndexType> elem_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id() };
            auto elem = new_model_part.CreateNewElement(   "Element3D4N",
                                                            p_elem->Id(),
                                                            elem_nodes,
                                                            p_prop );
        }

        owner_model.DeleteModelPart("AuxModelPartForCopy");
        return new_model_part;
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
