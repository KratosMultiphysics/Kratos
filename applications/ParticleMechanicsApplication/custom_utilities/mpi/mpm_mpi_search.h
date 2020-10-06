//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#ifndef KRATOS_MPM_MPI_SEARCH_IMPORT_H
#define KRATOS_MPM_MPI_SEARCH_IMPORT_H

// Project inlcudes
#include "includes/model_part.h"

namespace Kratos {

typedef ModelPart::ElementsContainerType ElementsArrayType;
typedef ModelPart::ConditionsContainerType ConditionsArrayType;

namespace MPM_MPI_SEARCH {

/**
 * @brief Search element connectivity for each particle in all mpi-partitions.
 * @details A search is performed to know in which grid element the material point falls.
 * If one or more material points fall in the grid element, the grid element is
 * set to be active and its connectivity is associated to the material point
 * element.
 * STEPS:
 * 1) All the elements are set to be INACTIVE
 * 2) For each particle:
 *      - Search element connectivity on current partition
 *      - If connecitiviy is found, update parent geometry of particle and set corresponding grid element active.
 *      - If connectivity could not be found:
 *              1. Send particle to all other partitions and remove from current model_part.
 *              2. Search particle on all other ranks (BinBasedSearch).
 *              3. Gather search results from all ranks for each particle.
 *              4. Add particle to model_part in which particle has been found (Thereby, lower ranks
 *                 have priority. Example: If rank 4 and 7 find the same particle, it will
 *                 only be added to the model_part in rank 4)
 *              5. Update parent geometry of particle and set corresponding grid element active.
 * 3) Sync nodal displacements at interface. (This is done to avoid the following problem: Imagine, one arbitrary particle
 *    leaves its current grid element and enters a new one. One node, associated to the new grid element, is an interface node.
 *    Furthermore, in the previous time step, this interface node was only active in an adjacent partition (not in the current
 *    one). Therefore, the displacement of this node is zero in the current partition and leads to wrong results.)
 **/
void SearchElementMPI(ModelPart& rBackgroundGridModelPart,
                      ModelPart& rMPMModelPart,
                      const std::size_t MaxNumberOfResults,
                      const double Tolerance);

/**
 * @brief Search particle connectivity in other partitions.
 * @details This function is called when a particle could not be found inside the current partition.
 * STEPS:
 * 1. Send particle to all other partitions and remove from current model_part.
 * 2. Search particle on all other ranks (BinBasedSearch).
 * 3. Gather search results from all ranks for each particle.
 * 4. Add particle to model_part in which particle has been found (Thereby, lower ranks have priority. Example: If rank 4 and 7
 *    find the same particle, it will only be added to the model_part in rank 4)
 * 5. Update parent geometry of particle and set corresponding grid element active.
 * @param rMissingElements contains missing particle elements
 * @param rMissingConditions contains missing particle conditions
 * */
void SearchElementsInOtherPartitions(ModelPart& rMPMModelPart, ModelPart& rBackgroundGridModelPart,
                                     std::vector<Element::Pointer>& rMissingElements,
                                     std::vector<Condition::Pointer>& rMissingConditions,
                                     const std::size_t MaxNumberOfResults, const double Tolerance);

/**
 * @brief Search element using a bin based approach.
 * @param rMissingElements Particle elements to be searched. If particle is found, it is removed from the vector.
 * @param rMissingConditions Particle condition to be searched. If particle is found, it is removed from the vector.
 * @param element_search_results Vector containing search results: [0, 1, 0] -> [not found, found, not found]
 * @param condition_search_results Vector containing search results: [0, 1, 0] -> [not found, found, not found]
 **/
void BinBasedSearchElementsAndConditionsMPI(ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<int>& element_search_results,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        std::vector<int>& condition_search_results,
        const std::size_t MaxNumberOfResults, const double Tolerance);

} // namespace MPM_MPI_SEARCH
} // namespace Kratos

#endif // KRATOS_MPM_MPI_SEARCH_IMPORT_H