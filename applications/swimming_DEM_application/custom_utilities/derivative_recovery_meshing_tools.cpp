//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2016-03-08 08:56:42 $
//
//

#include "derivative_recovery_meshing_tools.h"

namespace Kratos
{
template <std::size_t TDim>
void DerivativeRecoveryMeshingTools<TDim>::FillUpEdgesModelPartFromSimplicesModelPart(ModelPart& r_edges_model_part, ModelPart& r_tetra_model_part, std::string element_type)
{
    std::set< std::set<int> > set_of_all_edges; // actually, this is a set of pairs of Nodes Ids

    for (int i = 0; i < (int)r_tetra_model_part.NumberOfElements(); ++i){
        ElementIteratorType it_tetra = r_tetra_model_part.ElementsBegin() + i;
        Geometry<Node<3> >& geom = it_tetra->GetGeometry();

        for (unsigned int i_first_node = 0; i_first_node < TDim; ++i_first_node){
            int first_ID = geom[i_first_node].Id();

            for (unsigned int i_second_node = i_first_node + 1; i_second_node < TDim + 1; ++i_second_node){
                int second_ID = geom[i_second_node].Id();
                std::set<int> pair;
                pair.insert(first_ID);
                pair.insert(second_ID);
                set_of_all_edges.insert(pair);
            }
        }
    }

    Properties::Pointer p_properties = r_tetra_model_part.pGetProperties(0);
    std::vector<ModelPart::IndexType> pair;
    pair.resize(2);
	ModelPart::IndexType elem_id = 0;
	ModelPart::IndexType default_index = 0;

    for (auto i_edge : set_of_all_edges) {
        int i_node = 0;

        for (std::set<int>::iterator it_node = i_edge.begin(); it_node != i_edge.end(); ++it_node){
            pair[i_node] = *it_node;
            ++i_node;
        }

        r_edges_model_part.CreateNewElement(element_type, elem_id, pair, p_properties, default_index);
        ++elem_id;
    }
}

// Explicit instantiations
template class DerivativeRecoveryMeshingTools<2>;
template class DerivativeRecoveryMeshingTools<3>;

}  // namespace Kratos.
