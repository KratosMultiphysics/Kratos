//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "fluid_mesh_utilities.h"

namespace Kratos
{

bool FluidMeshUtilities::AllElementsAreSimplex(const ModelPart& rModelPart)
{
    // Check that the provided model part is not empty
    const unsigned int gl_n_elems = rModelPart.GetCommunicator().GlobalNumberOfElements();
    KRATOS_ERROR_IF(gl_n_elems == 0) << "There are no elements in model part '" << rModelPart.FullName() << "'." << std::endl;

    // Initialize the counter in case current partition is empty
    unsigned int n_simplex = 0;

    // Compute the total number of simplex elements in the mesh
    const auto& r_comm = rModelPart.GetCommunicator();
    const auto& r_elems = r_comm.LocalMesh().Elements();
    n_simplex = block_for_each<SumReduction<IndexType>>(r_elems, [](const auto& rElement){
        const auto& r_geometry = rElement.GetGeometry();
        return r_geometry.LocalSpaceDimension() + 1 == r_geometry.PointsNumber() ? 1 : 0;
    });
    n_simplex = r_comm.GetDataCommunicator().SumAll(n_simplex);

    // Check if the total number of simplex matches the total number of elements
    return n_simplex == gl_n_elems;
}

void FluidMeshUtilities::AssignNeighbourElementsToConditions(
    ModelPart& rModelPart,
    const bool CheckRepeatedConditions)
{
    // Create a map of the conditions with their connectivity as key
    // Note that the map data type is an array of conditions (not a single one)
    // This is intentionally done to track repeated conditions and for the eventual case in which we allow for repeated conditions
    ConditionsConnectivityMapType faces_map;
    for (auto it_cond = rModelPart.ConditionsBegin(); it_cond != rModelPart.ConditionsEnd(); ++it_cond) {
        // Reset the condition VISITED flag
        it_cond->Set(VISITED, false);

        // Creat the condition array of ids
        auto& r_geometry = it_cond->GetGeometry();
        const SizeType n_nodes = r_geometry.PointsNumber();
        DenseVector<int> cond_ids(n_nodes);
        for(IndexType i = 0; i < n_nodes; ++i) {
            r_geometry[i].Set(BOUNDARY,true); //FIXME: Check the purpose of this
            cond_ids[i] = r_geometry[i].Id();
        }

        // Sort the array of ids as this will be the key of the map
        std::sort(cond_ids.begin(), cond_ids.end());

        // Insert a pointer to the condition identified by the hash value ids
        auto aux_it = faces_map.find(cond_ids);
        if(aux_it == faces_map.end()) {
            faces_map.insert(std::make_pair(cond_ids, std::vector<Condition::Pointer>({*it_cond.base()})));
        } else {
            KRATOS_ERROR_IF(CheckRepeatedConditions)
                << "Condition " << it_cond->Id() << " shares the connectivity with condition " << aux_it->second[0]->Id() << ". Repeated conditions are not allowed in fluid meshes." << std::endl;
        }
    }

    // Now loop for all the elements and for each face of the element check if it is in the conditions map
    for (auto it_elem = rModelPart.ElementsBegin(); it_elem != rModelPart.ElementsEnd(); ++it_elem) {
        // Get the boundary entities (edges or faces)
        auto& r_geometry = it_elem->GetGeometry();
        const auto& r_boundary_entities = r_geometry.GenerateBoundariesEntities();

        // Loop the boundary entities to check their connectivities in the faces map
        DenseVector<int> r_face_ids(r_boundary_entities[0].PointsNumber());
        for (const auto& r_face : r_boundary_entities) {
            // Set the sorted array of ids
            IndexType i_node = 0;
            for (const auto& r_node : r_face) {
                r_face_ids[i_node++] = r_node.Id();
            }
            std::sort(r_face_ids.begin(), r_face_ids.end());

            // Search for current face in the conditions map
            auto it_face = faces_map.find(r_face_ids);
            if (it_face != faces_map.end()) {
                // Mark the found condition(s) as visited for the check below
                auto& r_list_conditions = it_face->second;
                for (const auto& rp_cond : r_list_conditions) {
                    rp_cond->Set(VISITED, true);
                }

                // Assign current element as condition parent (first neighbour entry)
                GlobalPointersVector<Element> vector_of_neighbours;
                vector_of_neighbours.resize(1);
                vector_of_neighbours(0) = Element::WeakPointer(*it_elem.base());
                for (const auto& rp_cond : r_list_conditions) {
                    rp_cond->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
                }
            }
        }
    }

    // Check that all of the conditions belong to at least an element
    // This is particularly useful in MPI to ensure that conditions lie in the same partition that the parents do
    if (rModelPart.GetCommunicator().LocalMesh().NumberOfConditions() != 0) {
        block_for_each(rModelPart.Conditions(), [](Condition& rCondition){
            KRATOS_ERROR_IF_NOT(rCondition.Is(VISITED)) << "Found a condition without any corresponding element. Id of condition " << rCondition.Id() << "." << std::endl;
        });
    }
}

} // namespace Kratos
