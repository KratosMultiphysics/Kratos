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

#pragma once

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(MAPPING_APPLICATION) BeamMapperInterfaceInfo : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    BeamMapperInterfaceInfo() {}

    explicit BeamMapperInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                 const IndexType SourceLocalSystemIndex,
                                 const IndexType SourceRank)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<BeamMapperInterfaceInfo>();
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<BeamMapperInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Node_Coords;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mBeamMapperId;
    }

    void GetValue(double& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mBeamMapperDistance;
    }

private:

    std::vector<int> mBeamMapperId = {};
    double mBeamMapperDistance = std::numeric_limits<double>::max();

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("BeamMapperId", mBeamMapperId);
        rSerializer.save("BeamMapperDistance", mBeamMapperDistance);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("BeamMapperId", mBeamMapperId);
        rSerializer.load("BeamMapperDistance", mBeamMapperDistance);
    }
};

class KRATOS_API(MAPPING_APPLICATION) BeamMapperLocalSystem : public MapperLocalSystem
{
public:

    explicit BeamMapperLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

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
        return Kratos::make_unique<BeamMapperLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

    void SetPairingStatusForPrinting() override;

private:
    NodePointerType mpNode;

};

/// Nearest Neighbor Mapper
/** This class implements the Nearest Neighbor Mapping technique.
* Each node on the destination side gets assigned is's closest neighbor on the other side of the interface.
* In the mapping phase every node gets assigned the value of it's neighbor
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class KRATOS_API(MAPPING_APPLICATION) BeamMapper
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of BeamMapper
    KRATOS_CLASS_POINTER_DEFINITION(BeamMapper);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    BeamMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination)
                          : BaseType(rModelPartOrigin, rModelPartDestination) {}

    BeamMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters)
                          : BaseType(rModelPartOrigin,
                                     rModelPartDestination,
                                     JsonParameters)
    {
        KRATOS_TRY;

        auto check_has_nodes = [](const ModelPart& rModelPart){
            if (rModelPart.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()) {
                KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfNodes() == 0) << "No nodes exist in ModelPart \"" << rModelPart.FullName() << "\"" << std::endl;
            }
        };
        check_has_nodes(rModelPartOrigin);
        check_has_nodes(rModelPartDestination);

        this->ValidateInput();
        this->Initialize();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~BeamMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<BeamMapper<TSparseSpace, TDenseSpace, TMapperBackend>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);

        KRATOS_CATCH("");
    }

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
        return "BeamMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BeamMapper";
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
            BeamMapperLocalSystem(nullptr),
            rModelPartCommunicator,
            rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return Kratos::make_unique<BeamMapperInterfaceInfo>();
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"              : {},
            "use_initial_configuration"    : false,
            "echo_level"                   : 0,
            "print_pairing_status_to_file" : false,
            "pairing_status_file_path"     : ""
        })");
    }

    ///@}

}; // Class BeamMapper

///@} addtogroup block
}  // namespace Kratos.