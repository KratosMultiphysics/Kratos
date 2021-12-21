//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#if !defined(KRATOS_BARYCENTRIC_MAPPER_H_INCLUDED )
#define  KRATOS_BARYCENTRIC_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"


namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(MAPPING_APPLICATION) BarycentricInterfaceInfo : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    explicit BarycentricInterfaceInfo(const std::size_t NumInterpolationNodes)
    {
        Initialize(NumInterpolationNodes);
    }

    explicit BarycentricInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                      const IndexType SourceLocalSystemIndex,
                                      const IndexType SourceRank,
                                      const std::size_t NumInterpolationNodes)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank)
    {
        Initialize(NumInterpolationNodes);
    }

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<BarycentricInterfaceInfo>(mNodeIds.size());
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<BarycentricInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank,
            mNodeIds.size());
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Node_Coords;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNodeIds;
    }

    void GetValue(std::vector<double>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNeighborCoordinates;
    }

private:

    std::vector<int> mNodeIds;
    std::vector<double> mNeighborCoordinates;

    void Initialize(const std::size_t NumInterpolationNodes)
    {
        mNodeIds.resize(NumInterpolationNodes);
        mNeighborCoordinates.resize(3*NumInterpolationNodes);
        std::fill(mNodeIds.begin(), mNodeIds.end(), -1);
        std::fill(mNeighborCoordinates.begin(), mNeighborCoordinates.end(), std::numeric_limits<double>::max());
    }

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("NodeIds", mNodeIds);
        rSerializer.save("NeighborCoords", mNeighborCoordinates);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("NodeIds", mNodeIds);
        rSerializer.load("NeighborCoords", mNeighborCoordinates);
    }
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
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        return mpNode->Coordinates();
    }

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<BarycentricLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

    void SetPairingStatusForPrinting() override;

private:
    NodePointerType mpNode;

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

        this->ValidateInput();

        const std::string interpolation_type = JsonParameters["interpolation_type"].GetString();
        if (interpolation_type == "line") {
            mNumInterpolationNodes = 2;
        } else if (interpolation_type == "triangle") {
            mNumInterpolationNodes = 3;
        } else if (interpolation_type == "tetrahedra") {
            mNumInterpolationNodes = 4;
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

    std::size_t mNumInterpolationNodes;

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
        return Kratos::make_unique<BarycentricInterfaceInfo>(mNumInterpolationNodes);
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

#endif // KRATOS_BARYCENTRIC_MAPPER_H_INCLUDED  defined