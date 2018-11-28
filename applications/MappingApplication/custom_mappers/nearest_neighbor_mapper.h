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
#include "mapper.h"
#include "custom_utilities/mapping_matrix_builder.h"
#include "custom_searching/interface_communicator.h"
#include "custom_utilities/interface_preprocessor.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"


namespace Kratos
{

///@name Kratos Classes
///@{

class NearestNeighborInterfaceInfo : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    NearestNeighborInterfaceInfo() {}

    NearestNeighborInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                 const IndexType SourceLocalSystemIndex,
                                 const IndexType SourceRank)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfo>();
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank=0) const override
    {
        return Kratos::make_shared<NearestNeighborInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank);
    }

    void ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                             const double NeighborDistance) override;

    void GetValue(int& rValue,
                  const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
    {
        rValue = mNearestNeighborId;
    }

    void GetValue(double& rValue,
                  const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
    {
        rValue = mNearestNeighborDistance;
    }

private:

    int mNearestNeighborId;
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

class NearestNeighborLocalSystem : public MapperLocalSystem
{
public:

    NearestNeighborLocalSystem() {}

    NearestNeighborLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<NearestNeighborLocalSystem>(pNode);
    }

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    bool UseNodesAsBasis() const override { return true; }

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        return mpNode->Coordinates();
    }

    /// Turn back information as a string.
    std::string PairingInfo(const int EchoLevel, const int CommRank) const override;

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
class NearestNeighborMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapper);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef Kratos::unique_ptr<InterfacePreprocessor> InterfacePreprocessorPointerType;
    typedef Kratos::unique_ptr<InterfaceCommunicator> InterfaceCommunicatorPointerType;
    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    typedef MappingMatrixBuilder<TSparseSpace, TDenseSpace> MappingMatrixBuilderType;
    typedef Kratos::unique_ptr<MappingMatrixBuilderType> MappingMatrixBuilderPointerType;
    typedef typename MappingMatrixBuilderType::MapperLocalSystemPointer MapperLocalSystemPointer;
    typedef typename MappingMatrixBuilderType::MapperLocalSystemPointerVector MapperLocalSystemPointerVector;
    typedef typename MappingMatrixBuilderType::MapperLocalSystemPointerVectorPointer MapperLocalSystemPointerVectorPointer;
    typedef typename MappingMatrixBuilderType::InterfaceVectorContainerType InterfaceVectorContainerType;
    typedef typename MappingMatrixBuilderType::InterfaceVectorContainerPointerType InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination)
                          : mrModelPartOrigin(rModelPartOrigin),
                            mrModelPartDestination(rModelPartDestination) {}

    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters)
                          : mrModelPartOrigin(rModelPartOrigin),
                            mrModelPartDestination(rModelPartDestination),
                            mMapperSettings(JsonParameters)
    {
        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartOrigin);
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartDestination);

        mpMapperLocalSystems = Kratos::make_shared<MapperLocalSystemPointerVector>();

        mpInterfacePreprocessor = Kratos::make_unique<InterfacePreprocessor>(mrModelPartDestination,
                                                                             mpMapperLocalSystems);
        InitializeMappingMatrixBuilder();

        ValidateInput(mMapperSettings);
        InitializeInterfaceCommunicator();

        InitializeInterface();
    }

    /// Destructor.
    ~NearestNeighborMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        // Set the Flags according to the type of remeshing
        if (MappingOptions.Is(MapperFlags::REMESHED)) {
            InitializeInterface(MappingOptions);
        }
        else {
            BuildMappingMatrix(MappingOptions);
        }

        if (mpInverseMapper) {
            mpInverseMapper->UpdateInterface(MappingOptions,
                                             SearchRadius);
        }
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            GetInverseMapper()->Map(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            GetInverseMapper()->Map(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

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
    ///@name Access
    ///@{

    TMappingMatrixType& GetMappingMatrix() override
    {
        return mpMappingMatrixBuilder->GetMappingMatrix();
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
    }

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    MappingMatrixBuilderPointerType mpMappingMatrixBuilder;
    InterfacePreprocessorPointerType mpInterfacePreprocessor;
    InterfaceCommunicatorPointerType mpIntefaceCommunicator;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;
    MapperLocalSystemPointerVectorPointer mpMapperLocalSystems;

    ///@}
    ///@name Private Operations
    ///@{

    void ValidateInput(Parameters AllMapperSettings);

    void ValidateParameters(Parameters AllMapperSettings)
    {
        Parameters default_settings = Parameters( R"({
            "search_radius"            : -1.0,
            "search_iterations"        : 3,
            "echo_level"               : 0
        })");

        AllMapperSettings.ValidateAndAssignDefaults(default_settings);
    }

    void InitializeInterfaceCommunicator();

    void InitializeMappingMatrixBuilder();

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void BuildMappingMatrix(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartOrigin.GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartDestination.GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<double>& rOriginVariable,
                              const Variable<double>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
                     const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<array_1d<double, 3>>& rOriginVariable,
                              const Variable<array_1d<double, 3>>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    ///@}
    ///@name Private  Access
    ///@{

    MapperUniquePointerType& GetInverseMapper()
    {
        if (!mpInverseMapper) {
            InitializeInverseMapper();
        }
        return mpInverseMapper;
    }

    void InitializeInverseMapper()
    {
        mpInverseMapper = Clone(mrModelPartDestination,
                                mrModelPartOrigin,
                                mMapperSettings);
    }

    ///@}

}; // Class NearestNeighborMapper

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED  defined