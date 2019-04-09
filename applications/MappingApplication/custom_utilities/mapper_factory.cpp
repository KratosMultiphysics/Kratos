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
#include "spaces/ublas_space.h"

#include "mapper_factory.h"
// #include "custom_utilities/mapper_utilities.h"


namespace Kratos
{
//     namespace MapperTypedefs
//     {
//         typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

//         typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//         typedef typename Mapper<SparseSpaceType, DenseSpaceType>::Pointer MapperPointer;
//         typedef std::unordered_map<std::string, MapperPointer> MapperPointerMap;

// #ifdef KRATOS_USING_MPI // mpi-parallel compilation
//         typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> MPISparseSpaceType;
//         typedef typename Mapper<MPISparseSpaceType, DenseSpaceType>::Pointer MPIMapperPointer;
//         typedef std::unordered_map<std::string, MPIMapperPointer> MPIMapperPointerMap;
// #endif
//     }

//     // Template instantiations
//     template<MapperTypedefs::SparseSpaceType, MapperTypedefs::DenseSpaceType>
//     typename Mapper<MapperTypedefs::SparseSpaceType, MapperTypedefs::DenseSpaceType>::Pointer MapperFactory::CreateMapper(
//         ModelPart& rModelPartOrigin,
//         ModelPart& rModelPartDestination,
//         Parameters MapperSettings)
//     {

//     }

//     typename Mapper<MapperTypedefs::MPISparseSpaceType, MapperTypedefs::DenseSpaceType>::Pointer MapperFactory::CreateMapper(
//         ModelPart& rModelPartOrigin,
//         ModelPart& rModelPartDestination,
//         Parameters MapperSettings)
//     {

//     }

//     void MapperFactory::Register(const std::string& rMapperName,
//                                  typename Mapper<MapperTypedefs::SparseSpaceType,
//                                     MapperTypedefs::DenseSpaceType>::Pointer pMapperPrototype)
//     {
//         GetRegisteredMappersList().insert(
//             make_pair(rMapperName, pMapperPrototype)); // TODO std::?
//     }

//     void MapperFactory::Register(const std::string& rMapperName,
//                                  typename Mapper<MapperTypedefs::MPISparseSpaceType,
//                                     MapperTypedefs::DenseSpaceType>::Pointer pMapperPrototype)
//     {
//         GetRegisteredMappersList().insert(
//             make_pair(rMapperName, pMapperPrototype)); // TODO std::?
//     }

//     std::unordered_map<std::string, typename Mapper<MapperTypedefs::SparseSpaceType,
//         MapperTypedefs::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList()
//     {
//         static std::unordered_map<std::string, typename Mapper<MapperTypedefs::SparseSpaceType,
//         MapperTypedefs::DenseSpaceType>::Pointer> registered_mappers;

//         return registered_mappers;
//     }

//     std::unordered_map<std::string, typename Mapper<MapperTypedefs::MPISparseSpaceType,
//         MapperTypedefs::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList()
//     {
//         static std::unordered_map<std::string, typename Mapper<MapperTypedefs::MPISparseSpaceType,
//         MapperTypedefs::DenseSpaceType>::Pointer> registered_mappers;

//         return registered_mappers;
//     }

    // template<class TSparseSpace, class TDenseSpace>
    // typename Mapper<TSparseSpace, TDenseSpace>::Pointer MapperFactory::CreateMapper(
    //     ModelPart& rModelPartOrigin,
    //     ModelPart& rModelPartDestination,
    //     Parameters MapperSettings)
    // {
    //     ModelPart& r_interface_model_part_origin = ReadInterfaceModelPart(rModelPartOrigin, MapperSettings, "origin");
    //     ModelPart& r_interface_model_part_destination = ReadInterfaceModelPart(rModelPartDestination, MapperSettings, "destination");

    //     const std::string mapper_name = MapperSettings["mapper_type"].GetString();

    //     const auto& mapper_list = GetRegisteredMappersList<TSparseSpace, TDenseSpace>();

    //     if (mapper_list.find(mapper_name) != mapper_list.end())
    //     {
    //         const bool is_mpi_execution = GetIsMPIExecution();

    //         // Removing Parameters that are not needed by the Mapper
    //         MapperSettings.RemoveValue("mapper_type");
    //         MapperSettings.RemoveValue("interface_submodel_part_origin");
    //         MapperSettings.RemoveValue("interface_submodel_part_destination");

    //         return mapper_list.at(mapper_name)->Clone(r_interface_model_part_origin,
    //                                                   r_interface_model_part_destination,
    //                                                   MapperSettings,
    //                                                   is_mpi_execution);
    //     }
    //     else
    //     {
    //         std::stringstream err_msg;
    //         err_msg << "The requested Mapper \"" << mapper_name <<"\" is not not available!\n"
    //                 << "The following Mappers are available:" << std::endl;

    //         for (auto const& registered_mapper : mapper_list)
    //             err_msg << "\t" << registered_mapper.first << "\n";

    //         KRATOS_ERROR << err_msg.str() << std::endl;
    //     }
    // }

    // template<class TSparseSpace, class TDenseSpace>
    // void MapperFactory::Register(const std::string& rMapperName,
    //                              typename Mapper<TSparseSpace, TDenseSpace>::Pointer pMapperPrototype)
    // {
    //     GetRegisteredMappersList<TSparseSpace, TDenseSpace>().insert(
    //         make_pair(rMapperName, pMapperPrototype)); // TODO std::?
    // }

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
    // template<class TSparseSpace, class TDenseSpace>
    // std::unordered_map<std::string, typename Mapper<TSparseSpace,
    //     TDenseSpace>::Pointer>& MapperFactory::GetRegisteredMappersList()
    // {
    //     static std::unordered_map<std::string, typename Mapper<TSparseSpace, TDenseSpace>::Pointer> registered_mappers;

    //     return registered_mappers;
    // }

    // TODO is doing this ok or do I pollute the namespace?

//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

//     template<>
//     std::unordered_map<std::string, typename Mapper<SparseSpaceType,
//         DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<SparseSpaceType, DenseSpaceType>()
//     {
//         static std::unordered_map<std::string, typename Mapper<SparseSpaceType, DenseSpaceType>::Pointer> registered_mappers;

//         return registered_mappers;
//     }

// #ifdef KRATOS_USING_MPI // mpi-parallel compilation
//     typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> MPISparseSpaceType;

//     template<>
//     std::unordered_map<std::string, typename Mapper<MPISparseSpaceType,
//         DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<MPISparseSpaceType, DenseSpaceType>()
//     {
//         static std::unordered_map<std::string, typename Mapper<MPISparseSpaceType, DenseSpaceType>::Pointer> registered_mappers;

//         return registered_mappers;
//     }

// #endif


}  // namespace Kratos.

