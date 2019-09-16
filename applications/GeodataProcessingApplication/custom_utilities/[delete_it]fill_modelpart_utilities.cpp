//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "fill_modelpart_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    // Function to fill "Parts_Fluid" sub model part from main_model_part
    void FillModelpartUtilities::FillFluid(){

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
        KRATOS_INFO("FillModelpartUtilities") << "bottom sub model part filled!" << std::endl;
    }


    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const FillModelpartUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
