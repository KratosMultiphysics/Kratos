//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "custom_utilities/mapper_typedefs.h"

#include "mapper_factory.h"


namespace Kratos
{

template<> KRATOS_API(MAPPING_APPLICATION) std::unordered_map<std::string, typename Mapper<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>();

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<> KRATOS_API(MAPPING_APPLICATION) std::unordered_map<std::string, typename Mapper<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>();
#endif

template<>
std::unordered_map<std::string, typename Mapper<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>()
{
    static std::unordered_map<std::string, typename Mapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>::Pointer> registered_mappers;

    return registered_mappers;
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
std::unordered_map<std::string, typename Mapper<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>()
{
    static std::unordered_map<std::string, typename Mapper<MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType>::Pointer> registered_mappers;

    return registered_mappers;
}
#endif

ModelPart& MapperFactory::ReadInterfaceModelPart(ModelPart& rModelPart,
                                                    Parameters InterfaceParameters,
                                                    const std::string& InterfaceSide)
{
    int echo_level = 0;
    // read the echo_level temporarily, bcs the mJsonParameters have not yet been validated and defaults assigned
    if (InterfaceParameters.Has("echo_level"))
    {
        echo_level = std::max(echo_level, InterfaceParameters["echo_level"].GetInt());
    }

    int comm_rank = rModelPart.GetCommunicator().MyPID();

    std::string key_sub_model_part = "interface_submodel_part_";
    key_sub_model_part.append(InterfaceSide);

    if (InterfaceParameters.Has(key_sub_model_part))
    {
        const std::string name_interface_submodel_part = InterfaceParameters[key_sub_model_part].GetString();

        if (echo_level >= 3 && comm_rank == 0)
        {
            std::cout << "Mapper: SubModelPart used for " << InterfaceSide << "-ModelPart" << std::endl;
        }

        return rModelPart.GetSubModelPart(name_interface_submodel_part);
    }
    else
    {
        if (echo_level >= 3 && comm_rank == 0)
        {
            std::cout << "Mapper: Main ModelPart used for " << InterfaceSide << "-ModelPart" << std::endl;
        }

        return rModelPart;
    }
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

}  // namespace Kratos.

