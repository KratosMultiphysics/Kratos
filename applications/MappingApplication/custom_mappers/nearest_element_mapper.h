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

class NearestElementInterfaceInfo : public MapperInterfaceInfo
{
public:

    /// Default constructor.
    NearestElementInterfaceInfo() {}

    NearestElementInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                const IndexType SourceLocalSystemIndex,
                                const IndexType SourceRank)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>();
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank=0) const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank);
    }

    void ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                             const double NeighborDistance) override;

    void ProcessSearchResultForApproximation(const InterfaceObject::Pointer& rpInterfaceObject,
                                             const double NeighborDistance) override;

    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
    {
        rValue = mNodeIds;
    }

    void GetValue(std::vector<double>& rValue,
                  const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
    {
        rValue = mShapeFunctionValues;
    }

    void GetValue(double& rValue,
                  const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
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

    NearestElementLocalSystem() {}

    NearestElementLocalSystem(NodePointerType pNode) : mpNode(pNode) {}

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<NearestElementLocalSystem>(pNode);
    }

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    bool UseNodesAsBasis() const override { return true; }

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
class NearestElementMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestElementMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

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
    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination) {}

    NearestElementMapper(ModelPart& rModelPartOrigin,
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
    ~NearestElementMapper() override = default;

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
        return Kratos::make_unique<NearestElementMapper<TSparseSpace, TDenseSpace>>(
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

}; // Class NearestElementMapper

}  // namespace Kratos.

#endif // KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED  defined
