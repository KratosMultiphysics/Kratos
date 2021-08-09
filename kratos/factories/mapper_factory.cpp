//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "mapper_factory.h"
#include "mappers/mapper_define.h"


namespace Kratos {

typedef typename MapperDefinitions::SparseSpaceType SparseSpaceType;
typedef typename MapperDefinitions::DenseSpaceType  DenseSpaceType;

template<>
std::unordered_map<std::string, typename Mapper<SparseSpaceType, DenseSpaceType>::Pointer>& MapperFactory<SparseSpaceType,
    DenseSpaceType>::GetRegisteredMappersList()
{
    static std::unordered_map<std::string, typename Mapper<SparseSpaceType, DenseSpaceType>::Pointer> registered_mappers;

    return registered_mappers;
}

    /* // CommRank is used as input bcs the MyPID function of the non-MPI MapperCommunicator is used
    // since this function is called before the MapperMPICommunicato is initialized
    void CheckInterfaceModelParts(const int CommRank)
    {
        const int num_nodes_origin = MapperUtilities::ComputeNumberOfNodes(mrModelPartOrigin);
        const int num_conditions_origin = MapperUtilities::ComputeNumberOfConditions(mrModelPartOrigin);
        const int num_elements_origin = MapperUtilities::ComputeNumberOfElements(mrModelPartOrigin);

        const int num_nodes_destination = MapperUtilities::ComputeNumberOfNodes(mrModelPartDestination);
        const int num_conditions_destination = MapperUtilities::ComputeNumberOfConditions(mrModelPartDestination);
        const int num_elements_destination = MapperUtilities::ComputeNumberOfElements(mrModelPartDestination);

        // Check if the ModelPart contains entities
        KRATOS_ERROR_IF(num_nodes_origin + num_conditions_origin + num_elements_origin < 1)
            << "Neither Nodes nor Conditions nor Elements found "
            << "in the Origin ModelPart" << std::endl;

        KRATOS_ERROR_IF(num_nodes_destination + num_conditions_destination + num_elements_destination < 1)
            << "Neither Nodes nor Conditions nor Elements found "
            << "in the Destination ModelPart" << std::endl;

        // Check if the inpt ModelParts contain both Elements and Conditions
        // This is NOT possible, bcs the InterfaceObjects are constructed
        // with whatever exists in the Modelpart (see the InterfaceObjectManagerBase,
        // function "InitializeInterfaceGeometryObjectManager")
        KRATOS_ERROR_IF(num_conditions_origin > 0 && num_elements_origin > 0)
            << "Origin ModelPart contains both Conditions and Elements "
            << "which is not permitted" << std::endl;

        KRATOS_ERROR_IF(num_conditions_destination > 0 && num_elements_destination > 0)
            << "Destination ModelPart contains both Conditions and Elements "
            << "which is not permitted" << std::endl;

        if (mEchoLevel >= 2) {
            std::vector<double> model_part_origin_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartOrigin);
            std::vector<double> model_part_destination_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartDestination);

            bool bbox_overlapping = MapperUtilities::ComputeBoundingBoxIntersection(
                                                        model_part_origin_bbox,
                                                        model_part_destination_bbox);
            if(CommRank == 0)
            {
                if (!bbox_overlapping) {
                    std::cout << "MAPPER WARNING, the bounding boxes of the "
                              << "Modelparts do not overlap! "
                              << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
                                                                              model_part_destination_bbox)
                              << std::endl;
                } else if (mEchoLevel >= 3)
                {
                    std::cout << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
                                                                              model_part_destination_bbox)
                              << std::endl;
                }
            }
        }
    }
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MapperFactory< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.

