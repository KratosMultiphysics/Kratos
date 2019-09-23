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
#include "fill_cfd_modelpart_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    ModelPart& FillCfdModelpartUtilities::FillPartsFluid(ModelPart& CfdModelPart, ModelPart& OriginModelPart, const std::string& ElementModelName) {
        
        auto& r_cfd_nodes_array = CfdModelPart.Nodes();
        std::set<int> cfd_nodes_set(begin(r_cfd_nodes_array), end(r_cfd_nodes_array));    // we conver the array into a set
        auto& r_cfd_elems_array = CfdModelPart.Elements();

        // CHECK IF THE SUB MODEL PART ALREADY EXISTS; OTHERWISE IT MUST BE CREATED
        if (CfdModelPart.HasModelPart("Parts_Fluid"))
            auto& r_fluid_model_part = CfdModelPart.GetSubModelPart("Parts_Fluid");
        else
            auto& r_fluid_model_part = CfdModelPart.CreateSubModelPart("Parts_Fluid");

        auto& r_elem_model_part = OriginModelPart.GetSubModelPart(ElementModelName);
        auto& r_elements_array = r_elem_model_part.Elements();

        Properties::Pointer p_prop = OriginModelPart.pGetProperties(0);
        std::vector<ModelPart::IndexType> elem_ids;

        for (int i_elem = 0; i_elem < static_cast<int>(r_elements_array.size()); ++i_elem) {
            const auto elem = r_elements_array.begin() + i_elem;
            auto& r_geom = elem->GetGeometry();  // nodes

            for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
                // bool node_exist = std::find(std::begin(r_cfd_nodes_array), std::end(r_cfd_nodes_array), r_geom[i_node]) != std::end(r_cfd_nodes_array);
                bool node_exist = std::find(std::begin(cfd_nodes_set), std::end(cfd_nodes_set), r_geom[i_node]) != std::end(cfd_nodes_set);

                if (!node_exist) {
                    auto new_node = CfdModelPart.CreateNewNode( r_geom[i_node].Id(),
                                                                r_geom[i_node].Coordinates()[0],
                                                                r_geom[i_node].Coordinates()[1],
                                                                r_geom[i_node].Coordinates()[2]);
                    r_fluid_model_part.AddNode(new_node);
                }
            }
            
            bool elem_exist = std::find(std::begin(r_cfd_elems_array), std::end(r_cfd_elems_array), *elem) != std::end(r_cfd_elems_array);
            if (!elem_exist) {
                auto new_elem = CfdModelPart.CreateNewElement(  "Element3D4N",
                                                                elem->Id(),
                                                                {r_geom[0].Id(), r_geom[1].Id(), r_geom[2].Id(), r_geom[3].Id()},
                                                                p_prop);
            }

            // r_fluid_model_part.AddElements(elem.Id());
            elem_ids.push_back(elem->Id());
        }

        r_fluid_model_part.AddElements(elem_ids);

        return CfdModelPart;
    }


    void FillCfdModelpartUtilities::FillInlet(){}


    void FillCfdModelpartUtilities::FillOutlet() {}


    void FillCfdModelpartUtilities::FillSlip() {}


    void FillCfdModelpartUtilities::FillNoslip() {}



    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const FillCfdModelpartUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
