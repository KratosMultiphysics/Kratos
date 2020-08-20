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

#if !defined(KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED )
#define  KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"


namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(MAPPING_APPLICATION) NearestNeighborInterfaceInfo : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    NearestNeighborInterfaceInfo() {}

    explicit NearestNeighborInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                 const IndexType SourceLocalSystemIndex,
                                 const IndexType SourceRank)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfo>();
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Node_Coords;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject,
                             const double NeighborDistance) override;

    void GetValue(int& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNearestNeighborId;
    }

    void GetValue(double& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNearestNeighborDistance;
    }

private:

    int mNearestNeighborId = -1;
    double mNearestNeighborDistance = std::numeric_limits<double>::max();

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("NearestNeighborId", mNearestNeighborId);
        rSerializer.save("NearestNeighborDistance", mNearestNeighborDistance);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("NearestNeighborId", mNearestNeighborId);
        rSerializer.load("NearestNeighborDistance", mNearestNeighborDistance);
    }
};

class KRATOS_API(MAPPING_APPLICATION) NearestNeighborLocalSystem : public MapperLocalSystem
{
public:

    explicit NearestNeighborLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        return mpNode->Coordinates();
    }

    /// Turn back information as a string.
    std::string PairingInfo(const int EchoLevel) const override;

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
template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) NearestNeighborMapper : public InterpolativeMapperBase<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapper);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination)
                          : BaseType(rModelPartOrigin, rModelPartDestination) {}

    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters)
                          : BaseType(rModelPartOrigin,
                                     rModelPartDestination,
                                     JsonParameters)
    {
        this->ValidateInput();
        this->Initialize();
    }

    /// Destructor.
    ~NearestNeighborMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<NearestNeighborMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "NearestNeighborMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestNeighborMapper";
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
        MapperUtilities::CreateMapperLocalSystemsFromNodes<NearestNeighborLocalSystem>(
            rModelPartCommunicator,
            rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return Kratos::make_unique<NearestNeighborInterfaceInfo>();
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_radius"                : -1.0,
            "search_iterations"            : 3,
            "echo_level"                   : 0,
            "print_pairing_status_to_file" : true,
            "pairing_status_file_path"     : ""
        })");
    }

    ///@}

}; // Class NearestNeighborMapper

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED  defined