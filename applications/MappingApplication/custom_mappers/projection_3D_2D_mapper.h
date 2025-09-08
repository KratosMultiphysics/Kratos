//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"
#include "custom_mappers/barycentric_mapper.h"
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

// Geometric definitions
using NodeType = Node;
using GeometryType = Geometry<NodeType>;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @brief This method retrieves the first geometry from a model part
 * @param rModelPart The input model part
 * @return The first geometry of the model part
 */
GeometryType::Pointer GetGeometryFromModelPart(const ModelPart& rModelPart)
{
    return (rModelPart.NumberOfElements() > 0 ? rModelPart.ElementsBegin()->pGetGeometry() : rModelPart.NumberOfConditions() > 0 ? rModelPart.ConditionsBegin()->pGetGeometry() : nullptr);
}

/**
 * @brief This method determines the partition where the model part has at least one entity
 * @param rModelPart The input model part where determine the partition with at least one entity
 * @return A partition with at least one partition
 */
int DeterminePartitionWithEntities(const ModelPart& rModelPart)
{
    // Get partition with entities
    auto p_geometry = GetGeometryFromModelPart(rModelPart);
    const int partition_entity = (p_geometry != nullptr) ? rModelPart.GetCommunicator().GetDataCommunicator().Rank() : -1;
    return rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(partition_entity);
}

/**
 * @brief This method determines maximum local space dimension of a ModelPart
 * @param rModelPart The model part to be considered
 * @return The maximum local space dimension
 */
unsigned int DetermineModelPartMaximumLocalDimension(ModelPart& rModelPart)
{
    // Getting local space dimension
    auto p_geometry = GetGeometryFromModelPart(rModelPart);
    const unsigned int local_space_dimension = (p_geometry == nullptr) ? 0 : p_geometry->LocalSpaceDimension();
    return rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_space_dimension);
}

/**
 * @brief This method determines the 2D model part
 * @param rFirstModelPart The first ModelPart
 * @param rSecondModelPart The second ModelPart
 * @return The 2D model part
 */
ModelPart& Determine2DModelPart(
    ModelPart& rFirstModelPart,
    ModelPart& rSecondModelPart
    )
{
    // Getting the maximum local space dimension
    const unsigned int max_local_space_dimension_1 = DetermineModelPartMaximumLocalDimension(rFirstModelPart);
    const unsigned int max_local_space_dimension_2 = DetermineModelPartMaximumLocalDimension(rSecondModelPart);
    KRATOS_ERROR_IF(max_local_space_dimension_1 == 3 && max_local_space_dimension_2 == 3) << "Both model parts are 3D" << std::endl;
    KRATOS_ERROR_IF(max_local_space_dimension_1 == 1 || max_local_space_dimension_2 == 1) << "One model part is 1D, not compatible" << std::endl;
    KRATOS_ERROR_IF(max_local_space_dimension_1 == 0 || max_local_space_dimension_2 == 0) << "Impossible to determine local space dimension in at least one model part" << std::endl;
    if (max_local_space_dimension_1 == 2) {
        return rFirstModelPart;
    } else if (max_local_space_dimension_2 == 2) {
        return rSecondModelPart;
    } else { // Corner case a priori impossible
        KRATOS_ERROR << "Impossible to detect 2D model part" << std::endl;
    }
    return rFirstModelPart;
}

/**
 * @brief This method determines the 3D model part
 * @details The counter part of the previous method
 * @param rFirstModelPart The first ModelPart
 * @param rSecondModelPart The second ModelPart
 * @return The 3D model part
 */
ModelPart& Determine3DModelPart(
    ModelPart& rFirstModelPart,
    ModelPart& rSecondModelPart
    )
{
    // Getting the maximum local space dimension
    const unsigned int max_local_space_dimension_1 = DetermineModelPartMaximumLocalDimension(rFirstModelPart);
    const unsigned int max_local_space_dimension_2 = DetermineModelPartMaximumLocalDimension(rSecondModelPart);
    KRATOS_ERROR_IF(max_local_space_dimension_1 == 3 && max_local_space_dimension_2 == 3) << "Both model parts are 3D" << std::endl;
    KRATOS_ERROR_IF(max_local_space_dimension_1 == 1 || max_local_space_dimension_2 == 1) << "One model part is 1D, not compatible" << std::endl;
    KRATOS_ERROR_IF(max_local_space_dimension_1 == 0 || max_local_space_dimension_2 == 0) << "Impossible to determine local space dimension in at least one model part" << std::endl;
    if (max_local_space_dimension_1 == 3) {
        return rFirstModelPart;
    } else if (max_local_space_dimension_2 == 3) {
        return rSecondModelPart;
    } else { // Corner case a priori impossible
        KRATOS_ERROR << "Impossible to detect 3D model part" << std::endl;
    }
    return rFirstModelPart;
}

///@}
///@name Kratos Classes
///@{

/**
 * @ingroup MapingApplication
 * @class Projection3D2DMapper
 * @brief This mapper simplifies the mapping between two model parts thanks to the projection over a reference plane
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class Projection3D2DMapper
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Projection3D2DMapper
    KRATOS_CLASS_POINTER_DEFINITION(Projection3D2DMapper);

    /// BaseType definitions
    using BaseType = InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>;
    using BaseMapperUniquePointerType = Kratos::unique_ptr<BaseType>;
    using TMappingMatrixType = typename BaseType::TMappingMatrixType;
    using MapperUniquePointerType = typename BaseType::MapperUniquePointerType;

    /// Interface definitions
    using InterfaceCommunicatorType = typename TMapperBackend::InterfaceCommunicatorType;
    using MapperInterfaceInfoUniquePointerType = typename InterfaceCommunicatorType::MapperInterfaceInfoUniquePointerType;

    /// Other mappers definition
    using NearestNeighborMapperType = NearestNeighborMapper<TSparseSpace, TDenseSpace, TMapperBackend>;
    using NearestElementMapperType = NearestElementMapper<TSparseSpace, TDenseSpace, TMapperBackend>;
    using BarycentricMapperType = BarycentricMapper<TSparseSpace, TDenseSpace, TMapperBackend>;

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief Entity type mesh considered
     */
    enum class EntityTypeMesh
    {
        NONE,
        CONDITIONS,
        ELEMENTS
    };

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    Projection3D2DMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination
        ) : BaseType(rModelPartOrigin, rModelPartDestination),
            // NOTE: TO AVOID ERROR AT REGISTER
            mr2DModelPart(rModelPartOrigin),
            mr3DModelPart(rModelPartDestination)
    {
    }

    // Constructor with settings
    Projection3D2DMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters
        ) : BaseType(rModelPartOrigin, rModelPartDestination, JsonParameters),
            mr2DModelPart(Determine2DModelPart(rModelPartOrigin, rModelPartDestination)),
            mr3DModelPart(Determine3DModelPart(rModelPartOrigin, rModelPartDestination))
    {
        KRATOS_TRY;

        // Validate input
        this->ValidateInput();

        // Copying parameters
        mCopiedParameters = JsonParameters.Clone();

        // Checking if the 2D modelpart is the origin model part
        CheckOriginIs2D();
        if (mOriginIs2D) {
            KRATOS_INFO("Projection3D2DMapper") << "The 2D model part is the origin model part" << std::endl;
        } else {
            KRATOS_INFO("Projection3D2DMapper") << "The 3D model part is the origin model part" << std::endl;
        }

        // Type of metamapper considered
        mMetaMapperType = mCopiedParameters["base_mapper"].GetString();

        // Pre mapper creation operations when model part of origin is 2D
        if (mOriginIs2D) {
            // Type of mesh entity considered
            GetEntityMeshType();

            // Getting the normal and reference plane
            GetNormalAndReferencePlane();

            // Move mesh
            MoveModelParts();
        }

        // Cleaning the parameters
        mCopiedParameters.RemoveValue("base_mapper");

        // Initializing the base mapper
        CreateBaseMapper();

        // Post mapper creation operations when model part of origin is 2D
        if (mOriginIs2D) {
            // Unmove mesh
            UnMoveModelParts();
        }

        // Calling initialize
        this->Initialize();

        // Now we copy the mapping matrix
        BaseType::mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(mpBaseMapper->GetMappingMatrix());

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~Projection3D2DMapper() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Updates the mapping-system after the geometry/mesh has changed
     * @details After changes in the topology (e.g. remeshing or sliding interfaces)
     * the relations for the mapping have to be recomputed. This means that
     * the search has to be conducted again and the mapping-system has to be
     * rebuilt, hence this is expensive
     * @param MappingOptions flags used to specify how the update has to be done
     * @param SearchRadius search radius used for the search
     */
    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius
        ) override
    {
        KRATOS_TRY;

        // Pre mapper creation operations when model part of origin is 2D
        if (mOriginIs2D) {
            // Move mesh
            MoveModelParts();
        }

        // Initializing the base mapper
        CreateBaseMapper();

        // Update interface base mapper
        mpBaseMapper->UpdateInterface(MappingOptions, SearchRadius);

        // Post mapper creation operations when model part of origin is 2D
        if (mOriginIs2D) {
            // Unmove mesh
            UnMoveModelParts();
        }

        // Calling initialize
        BaseType::UpdateInterface(MappingOptions, SearchRadius);

        // Now we copy the mapping matrix
        BaseType::mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(mpBaseMapper->GetMappingMatrix());

        KRATOS_CATCH("");
    }

    /**
    * @brief Cloning the Mapper
    * returns a clone of the current Mapper
    * pure virtual, has to be implemented in every derived mapper,
    * used in the creation of the Mappers
    * @see MapperFactory
    */
    MapperUniquePointerType Clone(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters
        ) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<Projection3D2DMapper<TSparseSpace, TDenseSpace, TMapperBackend>>(rModelPartOrigin, rModelPartDestination, JsonParameters);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the 2D model part
     * @return The 2D model part
     */
    ModelPart& Get2DModelPart()
    {
        return mr2DModelPart;
    }

    /**
     * @brief This method returns the 3D model part
     * @return The 3D model part
     */
    ModelPart& Get3DModelPart()
    {
        return mr3DModelPart;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Projection3D2DMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Projection3D2DMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mr2DModelPart;                   /// The 2D model part
    ModelPart& mr3DModelPart;                   /// The 3D model part
    BaseMapperUniquePointerType mpBaseMapper;   /// Pointer to the base mapper
    array_1d<double, 3> mNormalPlane;           /// The normal defining the plane to project
    Point mPointPlane;                          /// The coordinates of the plane to project
    Parameters mCopiedParameters;               /// The copied parameters. We copy the parameters to avoid conflicts in inverse mapping
    std::string mMetaMapperType;                /// The meta mapper type
    EntityTypeMesh mEntityTypeMesh;             /// The type of mesh considered
    bool mOriginIs2D;                           /// If the 2D modelpart is the origin model part

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Checking if the 2D modelpart is the origin model part
     */
    void CheckOriginIs2D()
    {
        KRATOS_TRY;

        const auto* p_origin_model_part = &(this->GetOriginModelPart());
        const auto* p_2d_model_part = &(this->Get2DModelPart());
        mOriginIs2D = p_origin_model_part == p_2d_model_part ? true : false;

        KRATOS_CATCH("");
    }

    /**
     * @brief Getting the mesh entity type considered
     */
    void GetEntityMeshType()
    {
        KRATOS_TRY;

        // A priori in the constructor we have already checked that the 2D is not empty
        const auto& r_2d_model_part = this->Get2DModelPart();
        if (r_2d_model_part.NumberOfConditions() > 0) {
            mEntityTypeMesh = EntityTypeMesh::CONDITIONS;
        } else if (r_2d_model_part.NumberOfElements() > 0) {
            mEntityTypeMesh = EntityTypeMesh::ELEMENTS;
        } else {
            mEntityTypeMesh = EntityTypeMesh::NONE; // Enum is zero, will be used to check if not elements/conditions in MPI
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Getting normal and reference plane
     */
    void GetNormalAndReferencePlane()
    {
        KRATOS_TRY;

        // We retrieve the values of interest
        const auto& r_2d_model_part = this->Get2DModelPart();
        const bool is_distributed = r_2d_model_part.IsDistributed();
        auto p_geometry = GetGeometryFromModelPart(r_2d_model_part);

        // MPI data
        const auto& r_communicator = r_2d_model_part.GetCommunicator();
        const auto& r_data_communicator = r_communicator.GetDataCommunicator();
        const int mpi_rank = r_data_communicator.Rank();
        const int mpi_size = r_data_communicator.Size();

        // Getting the partition with entities index
        const int partition_entity = DeterminePartitionWithEntities(r_2d_model_part);

        // Define the send tag
        const int tag_send_normal = 1;
        const int tag_send_point = 2;

        // Getting from parameters if not elements or conditions
        if (partition_entity != mpi_rank) {
            // Now transfer the normal plane and the point plane between partitions
            if (is_distributed) {
                // The partitions that receive
                r_data_communicator.Recv(mNormalPlane, partition_entity, tag_send_normal);
                r_data_communicator.Recv(mPointPlane.Coordinates(), partition_entity, tag_send_point);
            }
        } else {
            GeometryType::CoordinatesArrayType aux_coords;
            noalias(mPointPlane.Coordinates()) = p_geometry->Center();
            p_geometry->PointLocalCoordinates(aux_coords, mPointPlane);
            const bool is_pure_2d_geometry = p_geometry->LocalSpaceDimension() == p_geometry->WorkingSpaceDimension() ? true : false;
            if (is_pure_2d_geometry) {
                mNormalPlane[0] = 0.0;
                mNormalPlane[1] = 0.0;
                mNormalPlane[2] = 1.0;
            } else {
                noalias(mNormalPlane) = p_geometry->UnitNormal(aux_coords);
            }

            // Doing a check that all normals are aligned
            if (!is_pure_2d_geometry) {
                std::size_t check_normal;
                const double numerical_limit = std::numeric_limits<double>::epsilon() * 1.0e4;
                struct normal_check {
                    normal_check(array_1d<double, 3>& rNormal) : reference_normal(rNormal) {};
                    array_1d<double, 3> reference_normal;
                    GeometryType::CoordinatesArrayType aux_coords;
                };
                if (mEntityTypeMesh == EntityTypeMesh::CONDITIONS) {
                    check_normal = block_for_each<SumReduction<std::size_t>>(r_2d_model_part.Conditions(), normal_check(mNormalPlane), [&numerical_limit](auto& r_cond, normal_check& nc) {
                        auto& r_geom = r_cond.GetGeometry();
                        r_geom.PointLocalCoordinates(nc.aux_coords, r_geom.Center());
                        const auto normal = r_geom.UnitNormal(nc.aux_coords);
                        return (norm_2(normal - nc.reference_normal) > numerical_limit);
                    });
                } else {
                    check_normal = block_for_each<SumReduction<std::size_t>>(r_2d_model_part.Elements(), normal_check(mNormalPlane), [&numerical_limit](auto& r_elem, normal_check& nc) {
                        auto& r_geom = r_elem.GetGeometry();
                        r_geom.PointLocalCoordinates(nc.aux_coords, r_geom.Center());
                        const auto normal = r_geom.UnitNormal(nc.aux_coords);
                        return (norm_2(normal - nc.reference_normal) > numerical_limit);
                    });
                }
                KRATOS_ERROR_IF_NOT(check_normal == 0) << "The 2D reference model part has not consistent normals. Please check that is properly aligned" << std::endl;
            }

            // The partition that sends
            if (is_distributed) {
                const auto& r_point_coordinates = mPointPlane.Coordinates();
                for (int i_rank = 0; i_rank < mpi_size; ++i_rank) {
                    if (i_rank != partition_entity) {
                        r_data_communicator.Send(mNormalPlane, i_rank, tag_send_normal);
                        r_data_communicator.Send(r_point_coordinates, i_rank, tag_send_point);
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Create the base mapper
     */
    void CreateBaseMapper()
    {
        KRATOS_TRY;

        // Model parts
        auto& r_origin_model_part = this->GetOriginModelPart();
        auto& r_destination_model_part = this->GetDestinationModelPart();

        // Creating the base mapper
        if (mMetaMapperType == "nearest_neighbor") {
            if (mCopiedParameters.Has("interpolation_type")) mCopiedParameters.RemoveValue("interpolation_type");
            if (mCopiedParameters.Has("local_coord_tolerance")) mCopiedParameters.RemoveValue("local_coord_tolerance");
            mpBaseMapper = Kratos::make_unique<NearestNeighborMapperType>(r_origin_model_part, r_destination_model_part, mCopiedParameters);
        } else if (mMetaMapperType == "nearest_element") {
            if (mCopiedParameters.Has("interpolation_type")) mCopiedParameters.RemoveValue("interpolation_type");
            mpBaseMapper = Kratos::make_unique<NearestElementMapperType>(r_origin_model_part, r_destination_model_part, mCopiedParameters);
        } else if (mMetaMapperType == "barycentric") {
            mpBaseMapper = Kratos::make_unique<BarycentricMapperType>(r_origin_model_part, r_destination_model_part, mCopiedParameters);
        } else {
            KRATOS_ERROR << "Mapper " << mCopiedParameters["base_mapper"].GetString() << " is not available as base mapper for projection" << std::endl;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Moves the model parts
     * @details The 3D model part is projected to the 2D plane
     */
    void MoveModelParts()
    {
        KRATOS_TRY;

        // The 3D model part
        auto& r_3d_model_part = this->Get3DModelPart();

        // Save current configuration
        MapperUtilities::SaveCurrentConfiguration(r_3d_model_part);

        // Definition of projection variables
        struct ProjectionVariables
        {
            ProjectionVariables(array_1d<double, 3>& rNormal, Point& rPoint) : reference_normal(rNormal), reference_point(rPoint) {};
            array_1d<double, 3> reference_normal;
            Point reference_point;
            double distance;
            array_1d<double, 3> projected_point_coordinates;
        };

        // Iterate over the existing nodes
        block_for_each(r_3d_model_part.Nodes(), ProjectionVariables(mNormalPlane, mPointPlane), [&](auto& r_node, ProjectionVariables& p) {
            noalias(p.projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(p.reference_point, r_node, p.reference_normal, p.distance).Coordinates();
            noalias(r_node.Coordinates()) = p.projected_point_coordinates;
        });

        KRATOS_CATCH("");
    }

    /**
     * @brief Unmoves the model parts
     */
    void UnMoveModelParts()
    {
        KRATOS_TRY;

        // The 3D model part
        auto& r_3d_model_part = this->Get3DModelPart();

        // Restore configuration
        MapperUtilities::RestoreCurrentConfiguration(r_3d_model_part);

        KRATOS_CATCH("");
    }

    // Functions for customizing the behavior of this Mapper
    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems
        ) override
    {
        // Calling base mapper method. But not sure if this must be changed
        AccessorInterpolativeMapperBase<TMapperBackend>::CreateMapperLocalSystems(*mpBaseMapper, rModelPartCommunicator, rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return AccessorInterpolativeMapperBase<TMapperBackend>::GetMapperInterfaceInfo(*mpBaseMapper);
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"                    : {},
            "echo_level"                         : 0,
            "interpolation_type"                 : "unspecified",
            "local_coord_tolerance"              : 0.25,
            "use_initial_configuration"          : false,
            "print_pairing_status_to_file"       : false,
            "pairing_status_file_path"           : "",
            "base_mapper"                        : "nearest_neighbor"
        })");
    }

    ///@}

}; // Class Projection3D2DMapper

///@} addtogroup block
}  // namespace Kratos.
