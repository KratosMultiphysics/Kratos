// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/skin_detection_process.h"

namespace Kratos
{
template<SizeType TDim>
SkinDetectionProcess<TDim>::SkinDetectionProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mThisParameters(ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part" : "SkinModelPart",
        "name_auxiliar_condition"  : "Condition"
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::Execute()
{
    KRATOS_TRY;

    const SizeType number_of_elements = mrModelPart.Elements().size();

    // We clean the current database of neighbours
    for(IndexType i = 0; i < number_of_elements; ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;

        const SizeType reserve_size = ComputeReserveSize(it_elem);

        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
        r_neighbour_elements.reserve(reserve_size);
        r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end() );
    }

    /* NEIGHBOUR ELEMENTS */
    // Create the HashMapVectorIntElementPointerType
    HashMapVectorIntElementPointerType face_map;

    for(IndexType i = 0; i < number_of_elements; ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;

        GeometryType& geom = it_elem->GetGeometry();

        const SizeType reserve_size = ComputeReserveSize(it_elem);

        // Insert a pointer to the element identified by the hash value ids if it doesn't exist
        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);

        for (IndexType i_face = 0; i_face < reserve_size; ++i_face) {

            GeometryType& sub_geom = GetSubGeometry(geom, i_face);

            /* FACES/EDGES */
            const SizeType number_nodes = sub_geom.PointsNumber();
            VectorIndexType vector_ids(number_nodes);

            /* FACE/EDGE */
            for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                vector_ids[i_node] = sub_geom[i_node].Id();
            }

            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(vector_ids.begin(), vector_ids.end());
            // Check if the elements already exist in the HashMapVectorIntElementPointerType
            HashMapVectorIntElementPointerIteratorType it_check = face_map.find(vector_ids);

            if(it_check != face_map.end() ) {
                // If it exists the node is added as a neighbour, reciprocally
                r_neighbour_elements.push_back(it_check->second);
                ElementPointerVector& neighbour_elements = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
                neighbour_elements.push_back(*it_elem.base());
            } else {
                // If it doesn't exist it is added to the database
                face_map.insert( HashMapVectorIntElementPointerType::value_type(vector_ids, *it_elem.base()) );
            }
        }
    }

    // We create the auxiliar ModelPart
    ModelPart::Pointer p_auxiliar_model_part = mrModelPart.CreateSubModelPart(mThisParameters["name_auxiliar_model_part"].GetString());

    // The auxiliar name of the condition
    const std::string& name_condition = mThisParameters["name_auxiliar_model_part"].GetString();

    // The number of conditions
    IndexType condition_id = mrModelPart.Conditions().size();

    // Compute number the neighbour elements
//     #pragma omp parallel for // This can not be OMP because we are creating new conditions. We can create a vector in OMP and then fill the model part...
    for(int i = 0; i < static_cast<int>(number_of_elements); ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;

        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);

        GeometryType& geom = it_elem->GetGeometry();

        const SizeType reserve_size = ComputeReserveSize(it_elem);

        // If the real size is smaller than the
        if (r_neighbour_elements.size() < reserve_size) {
            for (IndexType i_face = 0; i_face < reserve_size; ++i_face) {

                GeometryType& sub_geom = GetSubGeometry(geom, i_face);

                /* FACES/EDGES */
                const SizeType number_nodes = sub_geom.PointsNumber();
                VectorIndexType vector_ids(number_nodes);

                /* FACE/EDGE */
                for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                    vector_ids[i_node] = sub_geom[i_node].Id();
                }

                /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
                std::sort(vector_ids.begin(), vector_ids.end());
                // Check if the elements already exist in the HashMapVectorIntElementPointerType
                HashMapVectorIntElementPointerIteratorType it_check = face_map.find(vector_ids);

                if(it_check == face_map.end() ) { // If it doesn't exist we add the skin to the auxiliar ModelPart
//                     #pragma omp atomic
                    condition_id += 1;

                    const std::string complete_name = name_condition + std::to_string(TDim) + "D" + std::to_string(number_nodes) + "N"; // If the condition doesn't follow this structure...sorry, we then need to modify this...
                    auto p_cond = p_auxiliar_model_part->CreateNewCondition(complete_name, condition_id, sub_geom, it_elem->pGetProperties());
                    p_cond->Set(INTERFACE, true);
                }
            }
        }
    }

    // Now we set the falg on the nodes. The list of nodes of the auxiliar model part
    auto& nodes_array = p_auxiliar_model_part->Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->Set(INTERFACE, true);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::ExecuteFinalize()
{
    ClearNeighbours();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::ClearNeighbours()
{
    for(IndexType i = 0; i < mrModelPart.Elements().size(); ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;
        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
        r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType SkinDetectionProcess<2>::ComputeReserveSize(ElementsIteratorType itElem)
{
    const auto& geometry = itElem->GetGeometry();
    return geometry.EdgesNumber();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType SkinDetectionProcess<3>::ComputeReserveSize(ElementsIteratorType itElem)
{
    const auto& geometry = itElem->GetGeometry();
    return geometry.FacesNumber();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
GeometryType& SkinDetectionProcess<2>::GetSubGeometry(
    GeometryType& BaseGeometry,
    const IndexType IndexSubGeometry
    )
{
    return BaseGeometry.Edges()[IndexSubGeometry];
}

/***********************************************************************************/
/***********************************************************************************/

template<>
GeometryType& SkinDetectionProcess<3>::GetSubGeometry(
    GeometryType& BaseGeometry,
    const IndexType IndexSubGeometry
    )
{
    return BaseGeometry.Faces()[IndexSubGeometry];
}

/***********************************************************************************/
/***********************************************************************************/

template class SkinDetectionProcess<2>;
template class SkinDetectionProcess<3>;
// class SkinDetectionProcess

} // namespace Kratos
