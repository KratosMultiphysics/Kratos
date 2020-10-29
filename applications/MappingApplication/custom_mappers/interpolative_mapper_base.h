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

#if !defined(KRATOS_INTERPOLATIVE_MAPPER_BASE_H_INCLUDED )
#define  KRATOS_INTERPOLATIVE_MAPPER_BASE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"
#include "custom_searching/interface_communicator.h"
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"


namespace Kratos
{
///@name Kratos Classes
///@{

template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) InterpolativeMapperBase : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterpolativeMapperBase
    KRATOS_CLASS_POINTER_DEFINITION(InterpolativeMapperBase);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef Kratos::unique_ptr<InterfaceCommunicator> InterfaceCommunicatorPointerType;
    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef Kratos::unique_ptr<TMappingMatrixType> TMappingMatrixUniquePointerType;

    typedef Variable<double> ComponentVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    InterpolativeMapperBase(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination) {}

    InterpolativeMapperBase(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters)
                         : mrModelPartOrigin(rModelPartOrigin),
                           mrModelPartDestination(rModelPartDestination),
                           mMapperSettings(JsonParameters)
    {
        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartOrigin);
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartDestination);
    }

    /// Destructor.
    ~InterpolativeMapperBase() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        KRATOS_WARNING_IF("Mapper", mMapperSettings["use_initial_configuration"].GetBool()) << "Updating the interface while using the initial configuration for mapping!" << std::endl;

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
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
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
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
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

    ///@}
    ///@name Access
    ///@{

    TMappingMatrixType& GetMappingMatrix() override
    {
        return *mpMappingMatrix;
    }

    ///@}
    ///@name Inquiry
    ///@{

    int AreMeshesConforming() const override
    {
        return mpIntefaceCommunicator->AreMeshesConforming();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "InterpolativeMapperBase";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterpolativeMapperBase";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

protected:

   /**
    * @brief Initializing the Mapper
    * This has to be called in the constructor of the
    * derived classes, since it involves calls to
    * pure virtual functions
    */
    void Initialize()
    {
        InitializeInterfaceCommunicator();
        InitializeInterface();
    }

    void ValidateInput();

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    TMappingMatrixUniquePointerType mpMappingMatrix;

    MapperLocalSystemPointerVector mMapperLocalSystems;

    InterfaceCommunicatorPointerType mpIntefaceCommunicator;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;

    ///@}
    ///@name Private Operations
    ///@{



    void InitializeInterfaceCommunicator();

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

    void PrintPairingInfo(const int EchoLevel);

    // functions for customizing the behavior of this Mapper
    virtual void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) = 0;

    virtual MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const = 0;

    virtual Parameters GetMapperDefaultSettings() const = 0;

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
        mpInverseMapper = this->Clone(mrModelPartDestination,
                                      mrModelPartOrigin,
                                      mMapperSettings);
    }

}; // Class InterpolativeMapperBase

}  // namespace Kratos.

#endif // KRATOS_INTERPOLATIVE_MAPPER_BASE_H_INCLUDED  defined
