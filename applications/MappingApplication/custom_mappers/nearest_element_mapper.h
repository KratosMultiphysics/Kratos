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

#if !defined(KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED )
#define  KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interpolative_mapper_base.h"


namespace Kratos
{
///@name Kratos Classes
///@{

class NearestElementInterfaceInfo : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    NearestElementInterfaceInfo() {}

    explicit NearestElementInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                const IndexType SourceLocalSystemIndex,
                                const IndexType SourceRank)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>();
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Geometry_Center;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject,
                             const double NeighborDistance) override;

    void ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject,
                                             const double NeighborDistance) override;

    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNodeIds;
    }

    void GetValue(std::vector<double>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mShapeFunctionValues;
    }

    void GetValue(double& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mClosestProjectionDistance;
    }

private:

    std::vector<int> mNodeIds;
    std::vector<double> mShapeFunctionValues;
    double mClosestProjectionDistance = std::numeric_limits<double>::max();

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("NodeIds", mNodeIds);
        rSerializer.save("SFValues", mShapeFunctionValues);
        rSerializer.save("ClosestProjectionDistance", mClosestProjectionDistance);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("NodeIds", mNodeIds);
        rSerializer.load("SFValues", mShapeFunctionValues);
        rSerializer.load("ClosestProjectionDistance", mClosestProjectionDistance);
    }

};

class NearestElementLocalSystem : public MapperLocalSystem
{
public:

    explicit NearestElementLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        return mpNode->Coordinates();
    }

    /// Turn back information as a string.
    std::string PairingInfo(const int EchoLevel, const int CommRank) const override;

private:
    NodePointerType mpNode;

};

/// Interpolative Mapper
/** This class implements the Nearest Element Mapping technique.
* Each node on the destination side gets assigned is's closest condition or element (distance to center)
* on the other side of the interface.
* In the mapping phase every node gets assigned the interpolated value of the condition/element.
* The interpolation is done with the shape funcitons
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
template<class TSparseSpace, class TDenseSpace>
class NearestElementMapper : public InterpolativeMapperBase<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestElementMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace> BaseType;
    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : BaseType(rModelPartOrigin,rModelPartDestination) {}

    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters)
                         : BaseType(rModelPartOrigin,
                                    rModelPartDestination,
                                    JsonParameters)
    {
        this->Initialize();
    }

    /// Destructor.
    ~NearestElementMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<NearestElementMapper<TSparseSpace, TDenseSpace>>(
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
        return "NearestElementMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestElementMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    ///@name Private Operations
    ///@{

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) override
    {
        MapperUtilities::CreateMapperLocalSystemsFromNodes<NearestElementLocalSystem>(
            rModelPartCommunicator,
            rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return Kratos::make_unique<NearestElementInterfaceInfo>();
    }

    ///@}

}; // Class NearestElementMapper

}  // namespace Kratos.

#endif // KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED  defined
