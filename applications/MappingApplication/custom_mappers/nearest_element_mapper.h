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
#include "custom_utilities/projection_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// @brief Options for configuring the behavior of the Mapper
struct NearestElementOptions
{
    bool UseApproximation = true; // Whether or not an approximation should be computed, if the full projection fails
    double LocalCoordTol = 0.25; // tolerance for the local coordinates, until which the projection is considered
};

class KRATOS_API(MAPPING_APPLICATION) NearestElementInterfaceInfo : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    explicit NearestElementInterfaceInfo(const NearestElementOptions Options={}) : mOptions(Options) {}

    explicit NearestElementInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                const IndexType SourceLocalSystemIndex,
                                const IndexType SourceRank,
                                         const NearestElementOptions Options={})
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank), mOptions(Options) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>(mOptions);
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank,
            mOptions);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Geometry_Center;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    void ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject) override;

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

    void GetValue(int& rValue,
                  const InfoType ValueType) const override
    {
        rValue = (int)mPairingIndex;
    }

    std::size_t GetNumSearchResults() const { return mNumSearchResults; }

private:

    std::vector<int> mNodeIds;
    std::vector<double> mShapeFunctionValues;
    double mClosestProjectionDistance = std::numeric_limits<double>::max();
    ProjectionUtilities::PairingIndex mPairingIndex = ProjectionUtilities::PairingIndex::Unspecified;
    NearestElementOptions mOptions; // this is not needed after searching, hence no need to serialize it
    std::size_t mNumSearchResults = 0;

    void SaveSearchResult(const InterfaceObject& rInterfaceObject,
                          const bool ComputeApproximation);

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("NodeIds", mNodeIds);
        rSerializer.save("SFValues", mShapeFunctionValues);
        rSerializer.save("ClosestProjectionDistance", mClosestProjectionDistance);
        rSerializer.save("PairingIndex", (int)mPairingIndex);
        rSerializer.save("NumSearchResults", mNumSearchResults);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("NodeIds", mNodeIds);
        rSerializer.load("SFValues", mShapeFunctionValues);
        rSerializer.load("ClosestProjectionDistance", mClosestProjectionDistance);
        int temp;
        rSerializer.load("PairingIndex", temp);
        mPairingIndex = (ProjectionUtilities::PairingIndex)temp;
        rSerializer.load("NumSearchResults", mNumSearchResults);
    }

};

class KRATOS_API(MAPPING_APPLICATION) NearestElementLocalSystem : public MapperLocalSystem
{
public:

    explicit NearestElementLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

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
        return Kratos::make_unique<NearestElementLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

    void SetPairingStatusForPrinting() override;

    bool IsDoneSearching() const override;

private:
    NodePointerType mpNode;
    mutable ProjectionUtilities::PairingIndex mPairingIndex = ProjectionUtilities::PairingIndex::Unspecified;

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
template<class TSparseSpace, class TDenseSpace, class TMapperBackend>
class KRATOS_API(MAPPING_APPLICATION) NearestElementMapper
    : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestElementMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, TMapperBackend> BaseType;
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
        KRATOS_TRY;

        this->ValidateInput();

        const bool use_approximation = JsonParameters["use_approximation"].GetBool();
        const double local_coord_tol = JsonParameters["local_coord_tolerance"].GetDouble();
        KRATOS_ERROR_IF(local_coord_tol < 0.0) << "The local-coord-tolerance cannot be negative" << std::endl;

        mOptions = NearestElementOptions{use_approximation, local_coord_tol};

        this->Initialize();

        KRATOS_CATCH("");
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
        KRATOS_TRY;

        return Kratos::make_unique<NearestElementMapper<TSparseSpace, TDenseSpace, TMapperBackend>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);

        KRATOS_CATCH("");
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
        BaseType::PrintData(rOStream);
    }

private:
    ///@name Member Variables
    ///@{

    NearestElementOptions mOptions;;

    ///@}

    ///@name Private Operations
    ///@{

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) override
    {
        MapperUtilities::CreateMapperLocalSystemsFromNodes(
            NearestElementLocalSystem(nullptr),
            rModelPartCommunicator,
            rLocalSystems);
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return Kratos::make_unique<NearestElementInterfaceInfo>(mOptions);
    }

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters( R"({
            "search_settings"              : {},
            "local_coord_tolerance"        : 0.25,
            "use_approximation"            : true,
            "use_initial_configuration"    : false,
            "echo_level"                   : 0,
            "print_pairing_status_to_file" : false,
            "pairing_status_file_path"     : ""
        })");
    }

    ///@}

}; // Class NearestElementMapper

}  // namespace Kratos.
