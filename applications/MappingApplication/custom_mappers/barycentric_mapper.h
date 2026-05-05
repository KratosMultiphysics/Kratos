//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"
#include "custom_utilities/closest_points.h"
#include "custom_utilities/projection_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

enum class BarycentricInterpolationType {
    LINE,
    TRIANGLE,
    TETRAHEDRA
};

class KRATOS_API(MAPPING_APPLICATION) BarycentricInterfaceInfo : public MapperInterfaceInfo
{
public:
    explicit BarycentricInterfaceInfo(const BarycentricInterpolationType InterpolationType);

    explicit BarycentricInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                      const IndexType SourceLocalSystemIndex,
                                      const IndexType SourceRank,
                                      const BarycentricInterpolationType InterpolationType);

    MapperInterfaceInfo::Pointer Create() const override;

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override;

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override;

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    BarycentricInterpolationType GetInterpolationType() const { return mInterpolationType; }

    const ClosestPointsContainer& GetClosestPoints() const { return mClosestPoints; }

    std::size_t GetNumSearchResults() const { return mNumSearchResults; }

private:
    BarycentricInterpolationType mInterpolationType;
    ClosestPointsContainer mClosestPoints;
    std::size_t mNumSearchResults = 0;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};

class KRATOS_API(MAPPING_APPLICATION) BarycentricLocalSystem : public MapperLocalSystem
{
public:

    explicit BarycentricLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

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
        return Kratos::make_unique<BarycentricLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

    void SetPairingStatusForPrinting() override;

    bool IsDoneSearching() const override;

private:
    NodePointerType mpNode;
    mutable ProjectionUtilities::PairingIndex mPairingIndex = ProjectionUtilities::PairingIndex::Unspecified;

};

/// Barycentric Mapper
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class KRATOS_API(MAPPING_APPLICATION) BarycentricMapper
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of BarycentricMapper
    KRATOS_CLASS_POINTER_DEFINITION(BarycentricMapper);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    BarycentricMapper(ModelPart& rModelPartOrigin,
                      ModelPart& rModelPartDestination)
                      : BaseType(rModelPartOrigin, rModelPartDestination) {}

    BarycentricMapper(ModelPart& rModelPartOrigin,
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

        const std::string interpolation_type = JsonParameters["interpolation_type"].GetString();
        if (interpolation_type == "line") {
            mInterpolationType = BarycentricInterpolationType::LINE;
        } else if (interpolation_type == "triangle") {
            mInterpolationType = BarycentricInterpolationType::TRIANGLE;
        } else if (interpolation_type == "tetrahedra") {
            mInterpolationType = BarycentricInterpolationType::TETRAHEDRA;
        } else {
            KRATOS_ERROR << "BarycentricMapper: No \"interpolation_type\" was specified, please select \"line\", \"triangle\" or \"tetrahedra\"" << std::endl;
        }

        this->Initialize();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~BarycentricMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<BarycentricMapper<TSparseSpace, TDenseSpace, TMapperBackend>>(
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
        return "BarycentricMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BarycentricMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:
    ///@name Member Variables
    ///@{

    BarycentricInterpolationType mInterpolationType;

    ///@}

    ///@name Private Operations
    ///@{

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) override
    {
        MapperUtilities::CreateMapperLocalSystemsFromNodes(
            BarycentricLocalSystem(nullptr),
            rModelPartCommunicator,
            rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return Kratos::make_unique<BarycentricInterfaceInfo>(mInterpolationType);
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"              : {},
            "interpolation_type"           : "unspecified",
            "local_coord_tolerance"        : 0.25,
            "use_initial_configuration"    : false,
            "echo_level"                   : 0,
            "print_pairing_status_to_file" : false,
            "pairing_status_file_path"     : ""
        })");
    }

    ///@}

}; // Class BarycentricMapper

///@} addtogroup block
}  // namespace Kratos.