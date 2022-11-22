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
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
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
        const std::string mapper_name = JsonParameters["base_mapper"].GetString();

        // Cleaning the parameters
        JsonParameters.RemoveValue("normal_plane");
        JsonParameters.RemoveValue("reference_plane_coordinates");
        JsonParameters.RemoveValue("base_mapper");

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
        if (mapper_name == "nearest_element" || mapper_name == "barycentric") { // || mapper_name == "coupling_geometry") {
            DelaunatorUtilities::CreateTriangleMeshFromNodes(r_projected_origin_modelpart);
        }

        /* Destination model part */
        auto& r_destination_model = rModelPartDestination.GetModel();
        auto& r_projected_destination_modelpart = r_destination_model.CreateModelPart("projected_destination_modelpart");
        for (auto& r_node : rModelPartDestination.Nodes()) {
            noalias(projected_point_coordinates) = GeometricalProjectionUtilities::FastProject(mPointPlane, r_node, mNormalPlane, distance).Coordinates();
            r_projected_destination_modelpart.CreateNewNode(r_node.Id(), projected_point_coordinates[0], projected_point_coordinates[1], projected_point_coordinates[2]); // TODO: This is assuming the plane is always XY, to fix after this works
        }

        // In destination we only consider nodes, so no Delaunay is required
        // // In case of nearest_element or barycentric or coupling_geometry we generate "geometries" to be able to interpolate
        // if (mapper_name == "nearest_element" || mapper_name == "barycentric") { // || mapper_name == "coupling_geometry") {
        //     DelaunatorUtilities::CreateTriangleMeshFromNodes(r_projected_destination_modelpart);
        // }

        // Initializing the base mapper
        if (mapper_name == "nearest_neighbor") {
            JsonParameters.RemoveValue("interpolation_type");
            JsonParameters.RemoveValue("local_coord_tolerance");
            mpBaseMapper = Kratos::make_unique<NearestNeighborMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        } else if (mapper_name == "nearest_element") {
            JsonParameters.RemoveValue("interpolation_type");
            mpBaseMapper = Kratos::make_unique<NearestElementMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        } else if (mapper_name == "barycentric") {
            mpBaseMapper = Kratos::make_unique<BarycentricMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        // } else if (mapper_name == "coupling_geometry") {
        //    mpBaseMapper = Kratos::make_unique<CouplingGeometryMapperType>(r_projected_origin_modelpart, r_projected_destination_modelpart, JsonParameters);
        } else {
            KRATOS_ERROR << "ERROR:: Mapper " << mapper_name << " is not available as base mapper for projection" << std::endl;
        }

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
    * After changes in the topology (e.g. remeshing or sliding interfaces)
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
        mpBaseMapper->UpdateInterface(MappingOptions, SearchRadius);
    }

    /**
    * @brief Mapping from Origin to Destination, Scalar Variable
    * Data is exchanged on the Interface, from the Origin-Modelpart
    * to the Destination-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see InverseMap
    */
    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions
        ) override
    {
        mpBaseMapper->Map(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /**
    * @brief Mapping from Origin to Destination, Vector Variable
    * Same as Map, but maps an array3-variable
    * @see Map
    */
    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions
        ) override
    {
        mpBaseMapper->Map(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /**
    * @brief Mapping from Destination to Origin, Scalar Variable
    * Data is exchanged on the Interface, from the Destination-Modelpart
    * to the Origin-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * It does the opposite of Map
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see Map
    */
    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions
        ) override
    {
        mpBaseMapper->InverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /**
    * @brief Mapping from Destination to Origin, Vector Variable
    * Same as InveseMap, but maps an array3-variable
    * @see InverseMap
    */
    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions
        ) override
    {
        mpBaseMapper->InverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
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
            "base_mapper"                  : "nearest_neighbor",
            "normal_plane"                 : [0.0,0.0,1.0],
            "reference_plane_coordinates"  : [0.0,0.0,0.0]
        })");
    }

    ///@}

}; // Class Projection3D2DMapper

///@} addtogroup block
}  // namespace Kratos.