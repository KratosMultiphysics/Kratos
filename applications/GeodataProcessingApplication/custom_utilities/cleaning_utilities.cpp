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
//                   Nicola Germano
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


    // // [NG] function to clear conditions
    // void CleaningUtilities::CleanConditions(){
    //     const int initial_cond = mrModelPart.Conditions().size();
    //     auto& r_conditions_array = mrModelPart.Conditions();
    //     const auto& r_nodes_array = mrModelPart.Nodes();        // list_node
    //     std::set<int> nodes_set(begin(r_nodes_array), end(r_nodes_array));    // we conver the array into a set

    //     for (int i_cond = 0; i_cond < static_cast<int>(r_conditions_array.size()); ++i_cond){
    //         const auto cond = r_conditions_array.begin() + i_cond;
    //         auto& r_geom = cond->GetGeometry();  // nodes

    //         for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
    //             // bool exist = std::find(std::begin(r_nodes_array), std::end(r_nodes_array), r_geom[i_node]) != std::end(r_nodes_array);
    //             bool exist = std::find(std::begin(nodes_set), std::end(nodes_set), r_geom[i_node]) != std::end(nodes_set);
    //             if (!exist){
    //                 // we delete the condition (and break from loop) if at least one node it is not in main model part
    //                 cond->Set(TO_ERASE, true);
    //                 break;
    //             }
    //         }
    //     }

    //     mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
    //     const int final_cond = mrModelPart.Conditions().size();
    //     KRATOS_INFO("\nCleaningUtilities") << "In total " << (initial_cond - final_cond) <<" superfluous conditions were cleared" << std::endl;
    // }


    // void CleaningUtilities::CleanConditions(){
    //     const int initial_cond = mrModelPart.Conditions().size();
    //     // auto& r_conditions_array = mrModelPart.Conditions();

    //     // nodes in ModelPart
    //     const auto& r_nodes_array = mrModelPart.Nodes();        // list_node
    //     std::set<int> nodes_set(begin(r_nodes_array), end(r_nodes_array));    // we conver the array into a set

    //     KRATOS_INFO("\nCleaningUtilities") << "nodes_set.size(): " << nodes_set.size() << std::endl;

    //     // we get Conditions from BottomModelPart
    //     auto& r_bottom_model_part = mrModelPart.GetSubModelPart("BottomModelPart");
    //     // const int initial_bottom_cond = r_bottom_model_part.Conditions().size();
    //     auto& r_bottom_conditions_array = r_bottom_model_part.Conditions();

    //     KRATOS_INFO("\nCleaningUtilities") << "r_bottom_conditions_array.size(): " << r_bottom_conditions_array.size() << std::endl;

    //     for (int i_cond = 0; i_cond < static_cast<int>(r_bottom_conditions_array.size()); ++i_cond) {
    //         const auto cond = r_bottom_conditions_array.begin() + i_cond;
    //         auto& r_geom = cond->GetGeometry();  // nodes

    //         KRATOS_INFO("\tCondition") << cond->Id() << std::endl;

    //         for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
    //             KRATOS_INFO("\t\tNode") << r_geom[i_node].Id() << std::endl;
    //             // bool exist = std::find(std::begin(r_nodes_array), std::end(r_nodes_array), r_geom[i_node]) != std::end(r_nodes_array);
    //             // bool exist = std::find(std::begin(nodes_set), std::end(nodes_set), r_geom[i_node]) != std::end(nodes_set);
    //             // if (!exist){
    //             //     // we delete the condition (and break from loop) if at least one node it is not in main model part
    //             //     KRATOS_INFO("\t CleanConditions") << "The Condition " << cond->Id() <<" will be delete!" << std::endl;
    //             //     cond->Set(TO_ERASE, true);
    //             //     break;
    //             // }

    //             // auto existing_node_it = r_nodes_array.find(r_geom[i_node].Id());
    //             if (nodes_set.find(r_geom[i_node].Id()) != nodes_set.end()) {
    //             // if (existing_node_it == r_nodes_array.end()) {
    //                 // we delete the condition (and break from loop) if at least one node it is not in main model part
    //                 KRATOS_INFO("\t CleanConditions") << "The Condition " << cond->Id() <<" will be delete!" << std::endl;
    //                 cond->Set(TO_ERASE, true);
    //                 break;
    //             }

    //             // /*** TRY WITH AT FUNCTION ***/
    //             // // TODO: solve the segmentation fault that appears when these lines are performed
    //             // auto existing_node_it = r_nodes_array.find(r_geom[i_node].Id());
    //             // if (existing_node_it == mrModelPart.NodesEnd()) {
    //             //     // we delete the condition (and break from loop) if at least one node it is not in main model part
    //             //     KRATOS_INFO("\t* CleanConditions") << "The Condition " << cond->Id() <<" will be delete!" << std::endl;
    //             //     cond->Set(TO_ERASE, true);
    //             //     break;
    //             // }
    //         }
    //     }

    //     mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
    //     const int final_cond = mrModelPart.Conditions().size();
    //     KRATOS_INFO("\nCleaningUtilities") << "In total " << (initial_cond - final_cond) <<" superfluous conditions were cleared" << std::endl;
    // }

    void CleaningUtilities::CleanConditions(){
        const int initial_cond = mrModelPart.Conditions().size();

        // nodes in ModelPart
        const auto& r_nodes_array = mrModelPart.Nodes();        // list_node

        // we get Nodes and Conditions from BottomModelPart
        auto& r_bottom_model_part = mrModelPart.GetSubModelPart("BottomModelPart");
        auto& r_bottom_nodes_array = r_bottom_model_part.Nodes();
        auto& r_bottom_conditions_array = r_bottom_model_part.Conditions();

        KRATOS_INFO("\t* CleanConditions") << "size r_bottom_nodes_array " << r_bottom_nodes_array.size() << std::endl;

        // all Nodes in BottomModelPart will set as "TO_ERASE true"
        for (int i_node = 0; i_node < static_cast<int>(r_bottom_nodes_array.size()); ++i_node) {
            auto node = r_nodes_array.begin() + i_node;
            node->Set(TO_ERASE, true);
        }

        // all Node in mrModelPart will set as "TO_ERASE false"
        // in this case only superfluous nodes remain as "TO_ERASE true"
        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            auto node = r_nodes_array.begin() + i_node;
            node->Set(TO_ERASE, false);
        }

        for (int i_cond = 0; i_cond < static_cast<int>(r_bottom_conditions_array.size()); ++i_cond) {
            const auto cond = r_bottom_conditions_array.begin() + i_cond;
            auto& r_geom = cond->GetGeometry();  // Nodes

            for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
                if (r_geom[i_node].Is(TO_ERASE)) {
                    // if at least one Node in current Condition is as TO_ERASE true, the current Condition will be delete
                    KRATOS_INFO("\t* CleanConditions") << "The Condition " << cond->Id() <<" will be delete!" << std::endl;
                    cond->Set(TO_ERASE, true);
                    break;
                }
            }
        }

        mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
        mrModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        const int final_cond = mrModelPart.Conditions().size();
        KRATOS_INFO("\nCleaningUtilities") << "In total " << (initial_cond - final_cond) <<" superfluous conditions were cleared" << std::endl;
    }


    // [NG] function to clear conditions in the angles
    void CleaningUtilities::CleanConditionsAngles() {
		const int initial_cond = mrModelPart.Conditions().size();

		auto& r_skin_model_part = mrModelPart.GetSubModelPart("SKIN_ISOSURFACE");
		auto& r_skin_nodes_array = r_skin_model_part.Nodes();
		// std::set<int> skin_nodes_set(begin(r_skin_nodes_array), end(r_skin_nodes_array));    // we conver the array into a set

		auto& r_bottom_model_part = mrModelPart.GetSubModelPart("BottomModelPart");
		auto& r_bottom_conditions_array = r_bottom_model_part.Conditions();

		for (int i_cond = 0; i_cond < static_cast<int>(r_bottom_conditions_array.size()); ++i_cond) {
			const auto cond = r_bottom_conditions_array.begin() + i_cond;
			auto& r_geom = cond->GetGeometry();  // nodes

			int count = 0;
			for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
                // bool exist_in_skin = std::find(std::begin(r_skin_nodes_array), std::end(r_skin_nodes_array), r_geom[i_node]) != std::end(r_skin_nodes_array);
				// // bool exist_in_skin = std::find(std::begin(skin_nodes_set), std::end(skin_nodes_set), r_geom[i_node]) != std::end(skin_nodes_set);
				// if (exist_in_skin) {
				// 	count++;
				// }
                auto existing_node_in_skin_it = r_skin_nodes_array.find(r_geom[i_node].Id());
                if (existing_node_in_skin_it != r_skin_model_part.NodesEnd()) {
                    count++;
                }
			}
			if (count == 3) {
				KRATOS_INFO("\t CleanConditionsAngles") << "The Condition " << cond->Id() <<" will be delete!" << std::endl;
				cond->Set(TO_ERASE, true);
			}
		}
		mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
		const int final_cond = mrModelPart.Conditions().size();
        KRATOS_INFO("\nCleanConditionsAngles") << "In total " << (initial_cond - final_cond) <<" superfluous conditions were cleared" << std::endl;
	}


    // [NG] fill sub model part with new conditions from the skin of building
    // CHECK IF IT NEED TO BE MOVED TO ANOTHER PART
    void CleaningUtilities::FillBottom(){

        std::vector<std::size_t> id_condition;
        std::vector<std::size_t> id_condition_new;
        std::vector<std::size_t> id_node_new;

        for (ModelPart::SubModelPartIterator submodelpart_it = mrModelPart.SubModelPartsBegin(); submodelpart_it != mrModelPart.SubModelPartsEnd(); ++submodelpart_it){
            ModelPart& submodelpart = *submodelpart_it;

            // we fill all id conditions into id_condition
            for (int i_cond = 0; i_cond < static_cast<int>(submodelpart.NumberOfConditions()); ++i_cond){
                auto p_cond = submodelpart.ConditionsBegin() + i_cond;
                id_condition.push_back(p_cond->Id());
            }
        }

        // loop in the main model part and we take only the conditions that are not in id_condition
        auto& r_conditions_array = mrModelPart.Conditions();
        for (int i_cond = 0; i_cond < static_cast<int>(r_conditions_array.size()); ++i_cond){
            const auto cond = r_conditions_array.begin() + i_cond;
            bool exist = std::find(std::begin(id_condition), std::end(id_condition), cond->Id()) != std::end(id_condition);
            if (!exist){
                id_condition_new.push_back(cond->Id());
                auto& r_geom = cond->GetGeometry();  // nodes
                for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
                    id_node_new.push_back(r_geom[i_node].Id());
                }
            }
        }

        // we delete the duplicate nodes
        std::sort(id_node_new.begin(), id_node_new.end());
        id_node_new.erase(unique(id_node_new.begin(), id_node_new.end()), id_node_new.end());

        // we add nodes and conditions into bottom sub model part
        auto& r_sub_model_part_bottom = mrModelPart.GetSubModelPart("bottom");
        r_sub_model_part_bottom.AddNodes(id_node_new);
        r_sub_model_part_bottom.AddConditions(id_condition_new);
        KRATOS_INFO("CleaningUtilities") << "bottom sub model part filled!" << std::endl;
    }


    ModelPart& CleaningUtilities::HardCopyBeforeSurfaceDiscretization( ModelPart& OriginalModelPart, ModelPart& NewModelPart ){

        auto& r_sub_model_part = NewModelPart.GetSubModelPart("AuxSubModelPart");
        // Properties::Pointer p_prop = NewModelPart.pGetProperties(0);
        Properties::Pointer p_prop = OriginalModelPart.pGetProperties(0);    // [NG]

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

            auto p_elem = OriginalModelPart.ElementsBegin() + i_elem;
            auto& r_geom = p_elem->GetGeometry();

            std::vector<ModelPart::IndexType> elem_nodes{ r_geom[0].Id() , r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id() };
            auto elem = NewModelPart.CreateNewElement(  "Element3D4N",
                                                        p_elem->Id(),
                                                        elem_nodes,
                                                        p_prop );
            // #pragma omp critical
            r_sub_model_part.AddElement( elem );
        }

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
            auto cond = NewModelPart.CreateNewCondition(  "SurfaceCondition3D3N" /*"Condition3D3N"*/,
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
