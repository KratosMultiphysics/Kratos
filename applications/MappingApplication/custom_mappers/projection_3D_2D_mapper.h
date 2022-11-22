//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
#include "custom_mappers/coupling_geometry_mapper.h"
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/tessellation_utilities/delaunator_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

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
    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
    typedef Kratos::unique_ptr<BaseType> BaseMapperUniquePointerType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;

    /// Interface definitions
    typedef typename TMapperBackend::InterfaceCommunicatorType InterfaceCommunicatorType;
    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    /// Other mappers definition
    typedef NearestNeighborMapper<TSparseSpace, TDenseSpace, TMapperBackend> NearestNeighborMapperType;
    typedef NearestElementMapper<TSparseSpace, TDenseSpace, TMapperBackend>   NearestElementMapperType;
    typedef BarycentricMapper<TSparseSpace, TDenseSpace, TMapperBackend>         BarycentricMapperType;
    // typedef CouplingGeometryMapper<TSparseSpace, TDenseSpace> CouplingGeometryMapperType; // Not compatible with base class InterpolativeMapperBase

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    Projection3D2DMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination
        ) : BaseType(rModelPartOrigin, rModelPartDestination) 
    {

    }

    // Constructor with settings
    Projection3D2DMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters
        ) : BaseType(rModelPartOrigin, rModelPartDestination, JsonParameters)
    {
        KRATOS_TRY;

        // Validate input
        this->ValidateInput();

        // We generate retrieve the values of interest
        mNormalPlane = JsonParameters["normal_plane"].GetVector();
        mPointPlane.Coordinates() = JsonParameters["reference_plane_coordinates"].GetVector();

        // Create the base mapper
        const std::string& r_mapper_name = JsonParameters["base_mapper"].GetString();

        // TODO: Before creating the mappers, we must generate the projected modelparts
        /* Origin model part */
        auto& r_origin_model = rModelPartOrigin.GetModel();
        auto& r_projected_origin_modelpart = r_origin_model.CreateModelPart("projected_origin_modelpart");

        // Iterate over the existing nodes
        double distance;
        array_1d<double, 3> projected_point_coordinates;
        for (auto& r_node : rModelPartOrigin.Nodes()) {
            noalias(projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(mPointPlane, r_node, mNormalPlane, distance).Coordinates();
            r_projected_origin_modelpart.CreateNewNode(r_node.Id(), projected_point_coordinates[0], projected_point_coordinates[1], projected_point_coordinates[2]); // TODO: This is assuming the plane is always XY, to fix after this works
        }

        // In case of nearest_element we generate "geometries" to be able to interpolate
        if (r_mapper_name == "nearest_element" || r_mapper_name == "barycentric") { // || r_mapper_name == "coupling_geometry") {
            DelaunatorUtilities::CreateTriangleMeshFromNodes(r_projected_origin_modelpart);
        }

        /* Destination model part */
        auto& r_destination_model = rModelPartDestination.GetModel();
        auto& r_projected_destination_modelpart = r_destination_model.CreateModelPart("projected_destination_modelpart");
        for (auto& r_node : rModelPartDestination.Nodes()) {
            noalias(projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(mPointPlane, r_node, mNormalPlane, distance).Coordinates();
            r_projected_destination_modelpart.CreateNewNode(r_node.Id(), projected_point_coordinates[0], projected_point_coordinates[1], projected_point_coordinates[2]); // TODO: This is assuming the plane is always XY, to fix after this works
        }

        // In case of nearest_element or barycentric or coupling_geometry we generate "geometries" to be able to interpolate
        if (r_mapper_name == "nearest_element" || r_mapper_name == "barycentric") { // || r_mapper_name == "coupling_geometry") {
            DelaunatorUtilities::CreateTriangleMeshFromNodes(r_projected_destination_modelpart);
        }

        // Initializing the base mapper
        if (r_mapper_name == "nearest_neighbor") {
            mpBaseMapper = Kratos::make_unique<NearestNeighborMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        } else if (r_mapper_name == "nearest_element") {
            mpBaseMapper = Kratos::make_unique<NearestElementMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        } else if (r_mapper_name == "barycentric") {
            mpBaseMapper = Kratos::make_unique<BarycentricMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        // } else if (r_mapper_name == "coupling_geometry") {
        //    mpBaseMapper = Kratos::make_unique<CouplingGeometryMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        } else {
            KRATOS_ERROR << "ERROR:: Mapper " << r_mapper_name << " is not available as base mapper for projection" << std::endl;
        }

        // Now we copy the mapping matrix
        this->GetMappingMatrix() = mpBaseMapper->GetMappingMatrix();

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

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    int AreMeshesConforming() const override
    {
        KRATOS_WARNING_ONCE("Mapper") << "Developer-warning: \"AreMeshesConforming\" is deprecated and will be removed in the future" << std::endl;
        return BaseType::mMeshesAreConforming;
    }

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

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    BaseMapperUniquePointerType mpBaseMapper; /// Pointer to the base mapper
    array_1d<double, 3> mNormalPlane;         /// The normal defining the plane to project
    Point mPointPlane;                        /// The coordinates of the plane to project

    ///@}
    ///@name Private Operations
    ///@{

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
            "search_settings"              : {},
            "echo_level"                   : 0,
            "interpolation_type"           : "unspecified",
            "local_coord_tolerance"        : 0.25,
            "use_initial_configuration"    : false,
            "print_pairing_status_to_file" : false,
            "pairing_status_file_path"     : "",
            "dual_mortar"                   : false,
            "precompute_mapping_matrix"     : false,
            "modeler_name"                  : "UNSPECIFIED",
            "modeler_parameters"            : {},
            "consistency_scaling"           : true,
            "row_sum_tolerance"             : 1e-12,
            "destination_is_slave"          : true,
            "linear_solver_settings"        : {},
            "base_mapper"                  : "nearest_neighbor",
            "normal_plane"                 : [0.0,0.0,1.0],
            "reference_plane_coordinates"  : [0.0,0.0,0.0]
        })");
    }

    ///@}

}; // Class Projection3D2DMapper

///@} addtogroup block
}  // namespace Kratos.