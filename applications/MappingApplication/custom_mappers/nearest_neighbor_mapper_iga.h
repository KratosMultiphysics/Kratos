//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti


#pragma once

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(MAPPING_APPLICATION) NearestNeighborInterfaceInfoIGA : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    NearestNeighborInterfaceInfoIGA()
        : mShapeFunctionValues(),
        mNearestNeighborId(),
        mNumberOfNearestNeighbors(0),
        mNearestNeighborDistance(std::numeric_limits<double>::max())
    {}

    explicit NearestNeighborInterfaceInfoIGA(const CoordinatesArrayType& rCoordinates,
                                             const IndexType SourceLocalSystemIndex,
                                             const IndexType SourceRank)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank),  
        mShapeFunctionValues(),
        mNearestNeighborId(),
        mNumberOfNearestNeighbors(0),
        mNearestNeighborDistance(std::numeric_limits<double>::max()) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfoIGA>();
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfoIGA>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Geometry_Center;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    // Shape function values for the local system evaluated in the coordinates of the gauss point
    void GetValue(std::vector<double>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mShapeFunctionValues;
    }
    
    // Nearest neighbor ids (control point ids) for the local system
    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNearestNeighborId;
    }

    void GetValue(double& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNearestNeighborDistance;
    }

    void GetValue(IndexType& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNumberOfNearestNeighbors;
    }

private:
    std::vector<double> mShapeFunctionValues;
    std::vector<int> mNearestNeighborId;
    IndexType mNumberOfNearestNeighbors;
    double mNearestNeighborDistance;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("ShapeFunctionValues", mShapeFunctionValues);
        rSerializer.save("NearestNeighborId", mNearestNeighborId);
        rSerializer.save("NumberofNearestNeighbors", mNumberOfNearestNeighbors);
        rSerializer.save("NearestNeighborDistance", mNearestNeighborDistance);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("ShapeFunctionValues", mShapeFunctionValues);
        rSerializer.load("NearestNeighborId", mNearestNeighborId);
        rSerializer.load("NumberofNearestNeighbors", mNumberOfNearestNeighbors);
        rSerializer.load("NearestNeighborDistance", mNearestNeighborDistance);
    }
};

class KRATOS_API(MAPPING_APPLICATION) NearestNeighborLocalSystemIGA : public MapperLocalSystem
{
public:

    explicit NearestNeighborLocalSystemIGA(NodePointerType pNode) : mpNode(pNode) {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;
        return mpNode->Coordinates();
    }

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<NearestNeighborLocalSystemIGA>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

    void SetPairingStatusForPrinting() override;

private:
    NodePointerType mpNode;

};

/// Nearest Neighbor Mapper for IGA-FEM simulations 
/** This class implements the Nearest Neighbor Mapping technique for IGA-FEM partitioned simulations.
* Each node on the destination side (FEM side) gets assigned is's closest integration point (neighbor) on the other side of the interface (IGA side).
* Once the nearest integration point is found, the shape function values of the integration point are used to interpolate the values from the IGA side to the FEM side.
* For information about the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class KRATOS_API(MAPPING_APPLICATION) NearestNeighborMapperIGA
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapperIGA
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapperIGA);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    NearestNeighborMapperIGA(ModelPart& rModelPartOrigin,
                             ModelPart& rModelPartDestination)
                             : BaseType(rModelPartOrigin, rModelPartDestination) {}

    NearestNeighborMapperIGA(ModelPart& rModelPartOrigin,
                             ModelPart& rModelPartDestination,
                             Parameters JsonParameters)
                             : BaseType(rModelPartOrigin,
                                       rModelPartDestination,
                                       JsonParameters)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(!JsonParameters.Has("is_origin_iga") || !JsonParameters["is_origin_iga"].GetBool())
            << "NearestNeighborMapperIGA expects the origin model part to be IGA.\n"
            << "Please set \"is_origin_iga\": true in the mapper settings." << std::endl;

        if (rModelPartOrigin.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) {
            KRATOS_ERROR_IF(rModelPartOrigin.GetCommunicator().GlobalNumberOfConditions() == 0)
                << "No conditions exist in ModelPart \"" << rModelPartOrigin.FullName() << "\"" << std::endl;
        }

        if (rModelPartDestination.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) {
            KRATOS_ERROR_IF(rModelPartDestination.GetCommunicator().GlobalNumberOfNodes() == 0)
                << "No nodes exist in ModelPart \"" << rModelPartDestination.FullName() << "\"" << std::endl;
        }

        this->ValidateInput();
        this->Initialize();

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<NearestNeighborMapperIGA<TSparseSpace, TDenseSpace, TMapperBackend>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "NearestNeighborMapperIGA";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestNeighborMapperIGA";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:

    ///@name Private Operations
    ///@{

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) override
    {
        MapperUtilities::CreateMapperLocalSystemsFromNodes(
            NearestNeighborLocalSystemIGA(nullptr),
            rModelPartCommunicator,
            rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return Kratos::make_unique<NearestNeighborInterfaceInfoIGA>();
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"              : {},
            "is_origin_iga"                : true,
            "use_initial_configuration"    : false,
            "echo_level"                   : 0,
            "print_pairing_status_to_file" : false,
            "pairing_status_file_path"     : ""
        })");
    }

    ///@}

}; // Class NearestNeighborMapperIGA

///@} addtogroup block
}  // namespace Kratos.