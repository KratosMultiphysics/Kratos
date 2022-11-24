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
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/tessellation_utilities/delaunator_utilities.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/parallel_utilities.h"

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

namespace ModelPartUtility
{
/**
 * @brief This function returns the corresponding model part considering the parameters
 * @param rOriginModelPart The model part to be considered
 * @param ThisParameters The configuration parameters
*/
ModelPart& GetOriginModelPart(
    ModelPart& rOriginModelPart,
    Parameters ThisParameters = Parameters(R"({})")
    )
{
    // We get the model part name
    const std::string origin_2d_sub_model_part_name = ThisParameters.Has("origin_2d_sub_model_part_name") ? ThisParameters["origin_2d_sub_model_part_name"].GetString() : "";

    // We check if the submodelpart exists
    if (origin_2d_sub_model_part_name != "") {
        // We check if the submodelpart is the modelpart
        if (rOriginModelPart.Name() == origin_2d_sub_model_part_name) {
            return rOriginModelPart;
        } else { // We retrieve the orginal modelpart if not
            KRATOS_ERROR_IF_NOT(rOriginModelPart.HasSubModelPart(origin_2d_sub_model_part_name)) << "Submodelpart " << rOriginModelPart.FullName() << "." << origin_2d_sub_model_part_name << " does not exist" << std::endl;
            return rOriginModelPart.GetSubModelPart(origin_2d_sub_model_part_name);
        }
    } else { // We return the origin model part if not defined
        return rOriginModelPart;
    }
}
} // namespace Utility

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
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;

    /// Interface definitions
    typedef typename TMapperBackend::InterfaceCommunicatorType InterfaceCommunicatorType;
    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    /// Other mappers definition
    typedef NearestNeighborMapper<TSparseSpace, TDenseSpace, TMapperBackend> NearestNeighborMapperType;
    typedef NearestElementMapper<TSparseSpace, TDenseSpace, TMapperBackend>   NearestElementMapperType;
    typedef BarycentricMapper<TSparseSpace, TDenseSpace, TMapperBackend>         BarycentricMapperType;

    /// Geometric definitions
    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

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
        ) : BaseType(ModelPartUtility::GetOriginModelPart(rModelPartOrigin, JsonParameters), rModelPartDestination, JsonParameters)
    {
        KRATOS_TRY;

        // Validate input
        this->ValidateInput();

        // We copy the parameters to avoid conflicts in inverse mapping
        Parameters copied_parameters = JsonParameters.Clone();

        // Create the base mapper
        const std::string mapper_name = copied_parameters["base_mapper"].GetString();

        // Origin model
        auto& r_origin_model = rModelPartOrigin.GetModel();

        // Destination model part
        auto& r_destination_model = rModelPartDestination.GetModel();

        // 2D model parts name if any
        const std::string origin_2d_sub_model_part_name = copied_parameters["origin_2d_sub_model_part_name"].GetString();
        const bool is_2d_origin = (origin_2d_sub_model_part_name != "") ? true : false;

        // We retrieve the values of interest
        if (is_2d_origin) {
            noalias(mNormalPlane) = copied_parameters["normal_plane"].GetVector();
            noalias(mPointPlane.Coordinates()) = copied_parameters["reference_plane_coordinates"].GetVector();
        } else {
            const auto& r_origin_sub_model_part = this->GetOriginModelPart();
            GeometryType::Pointer p_geometry = nullptr;
            int entity = 0;
            if (r_origin_sub_model_part.NumberOfElements() > 0) {
                const auto first_element = r_origin_sub_model_part.ElementsBegin();
                p_geometry = first_element->pGetGeometry();
            } else if (r_origin_sub_model_part.NumberOfConditions() > 0) {
                const auto first_condition = r_origin_sub_model_part.ConditionsBegin();
                p_geometry = first_condition->pGetGeometry();
                entity = 1;
            }
            // Getting from parameters if not elements or conditions
            if (p_geometry == nullptr) {
                KRATOS_ERROR_IF_NOT(mapper_name == "nearest_neighbor") << "The mapper \"nearest_element\"  or \"barycentric\" cannot be used without elements or conditions" << std::endl;
                noalias(mNormalPlane) = copied_parameters["normal_plane"].GetVector();
                noalias(mPointPlane.Coordinates()) = copied_parameters["reference_plane_coordinates"].GetVector();
            } else {
                GeometryType::CoordinatesArrayType aux_coords;
                noalias(mPointPlane.Coordinates()) = p_geometry->Center();
                p_geometry->PointLocalCoordinates(aux_coords, mPointPlane);
                noalias(mNormalPlane) = p_geometry->UnitNormal(aux_coords);

                // Doing a check that all normals are aligned
                std::size_t check_normal;
                const double numerical_limit = std::numeric_limits<double>::epsilon() * 1.0e4;
                struct normal_check {
                    normal_check(array_1d<double, 3>& rNormal) : reference_normal(rNormal) {};
                    array_1d<double, 3> reference_normal;
                    GeometryType::CoordinatesArrayType aux_coords;
                };
                if (entity == 0) {
                    check_normal = block_for_each<SumReduction<std::size_t>>(r_origin_sub_model_part.Elements(), normal_check(mNormalPlane), [&numerical_limit](auto& r_elem, normal_check& nc) {
                        auto& r_geom = r_elem.GetGeometry();
                        r_geom.PointLocalCoordinates(nc.aux_coords, r_geom.Center());
                        const auto normal = r_geom.UnitNormal(nc.aux_coords);
                        if (norm_2(normal - nc.reference_normal) > numerical_limit) { return 1; } else { return 0; }
                    });
                } else {
                    check_normal = block_for_each<SumReduction<std::size_t>>(r_origin_sub_model_part.Conditions(), normal_check(mNormalPlane), [&numerical_limit](auto& r_cond, normal_check& nc) {
                        auto& r_geom = r_cond.GetGeometry();
                        r_geom.PointLocalCoordinates(nc.aux_coords, r_geom.Center());
                        const auto normal = r_geom.UnitNormal(nc.aux_coords);
                        if (norm_2(normal - nc.reference_normal) > numerical_limit) { return 1; } else { return 0; }
                    });
                }
                KRATOS_ERROR_IF_NOT(check_normal == 0) << "The 2D reference model part has not consistent normals. Please check that is properly aligned" << std::endl;
            }
        }

        // Cleaning the parameters
        copied_parameters.RemoveValue("normal_plane");
        copied_parameters.RemoveValue("reference_plane_coordinates");
        copied_parameters.RemoveValue("base_mapper");
        copied_parameters.RemoveValue("origin_2d_sub_model_part_name");
        copied_parameters.RemoveValue("destination_2d_sub_model_part_name");

        // Projected origin model part. If the origin model part is not 2D we project it
        if (is_2d_origin) {
            auto& r_projected_origin_modelpart = r_origin_model.CreateModelPart("projected_origin_modelpart");

            // Iterate over the existing nodes
            double distance;
            array_1d<double, 3> projected_point_coordinates;
            for (auto& r_node : rModelPartOrigin.Nodes()) {
                noalias(projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(mPointPlane, r_node, mNormalPlane, distance).Coordinates();
                r_projected_origin_modelpart.CreateNewNode(r_node.Id(), projected_point_coordinates[0], projected_point_coordinates[1], projected_point_coordinates[2]);
            }

            // In case of nearest_element we generate "geometries" to be able to interpolate
            if (mapper_name == "nearest_element" || mapper_name == "barycentric") {
                DelaunatorUtilities::CreateTriangleMeshFromNodes(r_projected_origin_modelpart);
            }
        }

        // Projected destination model part
        {
            auto& r_projected_destination_modelpart = r_destination_model.CreateModelPart("projected_destination_modelpart");

            // Iterate over the existing nodes
            double distance;
            array_1d<double, 3> projected_point_coordinates;
            for (auto& r_node : rModelPartDestination.Nodes()) {
                noalias(projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(mPointPlane, r_node, mNormalPlane, distance).Coordinates();
                r_projected_destination_modelpart.CreateNewNode(r_node.Id(), projected_point_coordinates[0], projected_point_coordinates[1], projected_point_coordinates[2]);
            }
        }

        // Initializing the base mapper
        {
            auto& r_projected_origin_modelpart = r_origin_model.HasModelPart("projected_origin_modelpart") ? r_origin_model.GetModelPart("projected_origin_modelpart") : this->GetOriginModelPart();
            auto& r_projected_destination_modelpart = r_destination_model.GetModelPart("projected_destination_modelpart");
            if (mapper_name == "nearest_neighbor") {
                copied_parameters.RemoveValue("interpolation_type");
                copied_parameters.RemoveValue("local_coord_tolerance");
                mpBaseMapper = Kratos::make_unique<NearestNeighborMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, copied_parameters);
            } else if (mapper_name == "nearest_element") {
                copied_parameters.RemoveValue("interpolation_type");
                mpBaseMapper = Kratos::make_unique<NearestElementMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, copied_parameters);
            } else if (mapper_name == "barycentric") {
                mpBaseMapper = Kratos::make_unique<BarycentricMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, copied_parameters);
            } else {
                KRATOS_ERROR << "ERROR:: Mapper " << mapper_name << " is not available as base mapper for projection" << std::endl;
            }
        }

        // Calling initialize
        this->Initialize();

        // Now we copy the mapping matrix
        BaseType::mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(mpBaseMapper->GetMappingMatrix());

        // Remove the created model parts to not crash in InvertedMap
        if (r_origin_model.HasModelPart("projected_origin_modelpart")) {
            r_origin_model.DeleteModelPart("projected_origin_modelpart");
        }
        r_destination_model.DeleteModelPart("projected_destination_modelpart");

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

    /**
     * @brief This function creates the inverted mapping parameters if they are required to be differemt from the forward mapping parameters
     * @details This function has to be implemented in the derived classes in case the inverted mapping parameters are required to be different from the forward mapping parameters
     * @return The inverted mapping parameters
     */
    Parameters GetInvertedMappingParameters(Parameters ForwardMappingParameters) override
    {
        // Copy the parameters
        Parameters inverted_parameters = ForwardMappingParameters.Clone();

        // Invserse 2D submodelparts mapping parameters
        const std::string origin_2d_sub_model_part_name = inverted_parameters["origin_2d_sub_model_part_name"].GetString();
        const std::string destination_2d_sub_model_part_name = inverted_parameters["destination_2d_sub_model_part_name"].GetString();
        inverted_parameters["origin_2d_sub_model_part_name"].SetString(destination_2d_sub_model_part_name);
        inverted_parameters["destination_2d_sub_model_part_name"].SetString(origin_2d_sub_model_part_name);
        
        return inverted_parameters;
    }

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
            "search_settings"                    : {},
            "echo_level"                         : 0,
            "interpolation_type"                 : "unspecified",
            "origin_2d_sub_model_part_name"      : "",
            "destination_2d_sub_model_part_name" : "",
            "local_coord_tolerance"              : 0.25,
            "use_initial_configuration"          : false,
            "print_pairing_status_to_file"       : false,
            "pairing_status_file_path"           : "",
            "base_mapper"                        : "nearest_neighbor",
            "normal_plane"                       : [0.0,0.0,1.0],
            "reference_plane_coordinates"        : [0.0,0.0,0.0]
        })");
    }

    ///@}

}; // Class Projection3D2DMapper

///@} addtogroup block
}  // namespace Kratos.