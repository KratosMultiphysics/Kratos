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

#if !defined(KRATOS_MAPPER_H_INCLUDED )
#define  KRATOS_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/matrix_based_mapping_operation_utility.h"
#include "custom_utilities/interface_preprocessor.h"
#include "custom_utilities/mapper_local_system.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_searching/interface_search_structure_base.h"
#include "mapping_application_variables.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Base Class for all Mappers
/** This is the base class for every mapper.
* It contains the three pure virtual functions that have to be implemented by every mapper:
* - Map: Basic function that maps a field from one ModelPart to another Modelpart
*        Mapping Direction: Origin => Destionation
* - InverseMap: This function does the opposite of the "Map" function
*               Mapping Direction: Destination => Origin
* - UpdateInterface: Called when the interface is changed. It recomputes the neighbors and
*   other information related to the relations btw entities (node, elements,...) on the interfaces
* It is also responsible for initializing the MapperCommunicator or the MapperMPICommuniator
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/

template<class TSparseSpace, class TDenseSpace>
class Mapper
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of Mapper
    KRATOS_CLASS_POINTER_DEFINITION(Mapper);

    typedef MappingOperationUtility<TSparseSpace, TDenseSpace> MappingOperationUtilityType;
    typedef typename Kratos::unique_ptr<MappingOperationUtilityType> MappingOperationUtilityPointerType;
    using InterfacePreprocessorPointerType = Kratos::unique_ptr<InterfacePreprocessor>;
    using SearchStructureBaseType = InterfaceSearchStructureBase;
    using SearchStructurePointerType = Kratos::unique_ptr<SearchStructureBaseType>;
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    using MapperUniquePointerType = Kratos::unique_ptr<Mapper>;

    using MapperInterfaceInfoUniquePointerType = typename SearchStructureBaseType::MapperInterfaceInfoUniquePointerType;

    using MapperLocalSystemPointer = typename MappingOperationUtilityType::MapperLocalSystemPointer;
    using MapperLocalSystemPointerVector = typename MappingOperationUtilityType::MapperLocalSystemPointerVector;
    using MapperLocalSystemPointerVectorPointer = typename MappingOperationUtilityType::MapperLocalSystemPointerVectorPointer;

    using TSystemMatrixType = typename MappingOperationUtilityType::TSystemMatrixType;
    using TSystemVectorType = typename MappingOperationUtilityType::TSystemVectorType;

    using TSystemMatrixUniquePointerType = typename MappingOperationUtilityType::TSystemMatrixUniquePointerType;
    using TSystemVectorUniquePointerType = typename MappingOperationUtilityType::TSystemVectorUniquePointerType;


    ///@}
    ///@name Life Cycle
    ///@{

    Mapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination) :
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mGeneralMapperSettings(Parameters(R"({})")) {}

    /// Destructor.
    virtual ~Mapper()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius)
    {
        std::cout << "\n\n";
        if (MappingOptions.IsDefined(MapperFlags::REMESHED))
        {
            KRATOS_INFO("\tFlagsCheck") << "IsDefined" << std::endl;
            if (MappingOptions.Is(MapperFlags::REMESHED))
            {
                KRATOS_INFO("\t\tFlagsCheck") << "IsDefined AND Is" << std::endl;
            }
            else
            {
                KRATOS_INFO("\t\tFlagsCheck") << "IsDefined AND NOT Is" << std::endl;
            }
        }
        else
        {
            KRATOS_INFO("\tFlagsCheck") << "NOT IsDefined" << std::endl;
                        KRATOS_INFO("\tFlagsCheck") << "IsDefined" << std::endl;
            if (MappingOptions.Is(MapperFlags::REMESHED))
            {
                KRATOS_INFO("\t\tFlagsCheck") << "NOT IsDefined AND Is" << std::endl;
            }
            else
            {
                KRATOS_INFO("\t\tFlagsCheck") << "NOT IsDefined AND NOT Is" << std::endl;
            }
        }

        std::cout << "\n";

        if (MappingOptions.Is(MapperFlags::REMESHED))
        {
            KRATOS_INFO("\tFlagsCheck") << "Is" << std::endl;
        }
        else
        {
            KRATOS_INFO("\tFlagsCheck") << "NOT Is" << std::endl;
        }

        std::cout << "\n\n";

        UpdateInterfaceInternal(MappingOptions, SearchRadius);
        if (mInverseMapperIsInitialized)
            mpInverseMapper->UpdateInterface(MappingOptions, SearchRadius);
    }

    /* This function maps from Origin to Destination */
    virtual void Map(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions)
    {
        TMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /* This function maps from Origin to Destination */
    virtual void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
                     const Variable< array_1d<double, 3> >& rDestinationVariable,
                     Kratos::Flags MappingOptions)
    {
        TMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /* This function maps from Destination to Origin */
    virtual void InverseMap(const Variable<double>& rOriginVariable,
                            const Variable<double>& rDestinationVariable,
                            Kratos::Flags MappingOptions)
    {
        TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /* This function maps from Destination to Origin */
    virtual void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                            const Variable< array_1d<double, 3> >& rDestinationVariable,
                            Kratos::Flags MappingOptions)
    {
        TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    virtual MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) = 0;

    // Developer function only being used if this class needs to be accessed from outside!
    // can be overridden in case it is needed (e.g. for Mortar where Mdd != I !)
    virtual double GetMappingMatrixEntry(const IndexType RowIndex, const IndexType ColumnIndex)
    {
        // return mpMappingOperationUtility->GetMappingMatrixEntry(RowIndex, ColumnIndex);
        KRATOS_ERROR << "Not implemented" << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Mapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Mapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;

    Parameters mGeneralMapperSettings = Parameters(R"({})");
    Parameters mGeneralMapperSettingsDefaults = Parameters( R"(
    {
        "search_radius"     : 0.0,
        "search_iterations" : 3,
        "echo_level"        : 0
    }  )" );

    MappingOperationUtilityPointerType mpMappingOperationUtility;
    InterfacePreprocessorPointerType mpInterfacePreprocessor;
    SearchStructurePointerType mpSearchStructure;
    MapperLocalSystemPointerVectorPointer mpMapperLocalSystems;

    // The mapping matrix and the corresponding vectors
    TSystemMatrixUniquePointerType mpMdo;
    TSystemVectorUniquePointerType mpQo;
    TSystemVectorUniquePointerType mpQd;

    // global, aka of the entire submodel-parts
    // int mNumConditionsOrigin;
    // int mNumConditionsDestination;

    // int mNumNodesOrigin;
    // int mNumNodesDestination;

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    // Constructor, can only be called by derived classes (actual mappers)
    Mapper(ModelPart& rModelPartOrigin,
           ModelPart& rModelPartDestination,
           Parameters MapperSettings);


    /**
     * This function can be overridden by derived Mappers to do sth different
     * */
    virtual void Initialize();

    virtual void InitializeSearchStructure();

    virtual void InitializeMappingOperationUtility();

    virtual void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    virtual void BuildMappingMatrix(Kratos::Flags MappingOptions = Kratos::Flags());

    virtual void UpdateInterfaceInternal(Kratos::Flags MappingOptions, double SearchRadius);


    template< typename TDataType >
    void TMap(const Variable<TDataType>& rOriginVariable,
              const Variable<TDataType>& rDestinationVariable,
              Kratos::Flags MappingOptions,
              const bool UseTranspose = false)
    {
        CheckForConservative(MappingOptions);

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE))
        {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE); // TODO test this!!!
            const bool use_transpose = true;
            TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions, use_transpose);
        }
        else
        {
            KRATOS_DEBUG_ERROR_IF_NOT(mpMappingOperationUtility)<< "mpMappingOperationUtility "
                << "is a nullptr" << std::endl;

            mpMappingOperationUtility->ExecuteMapping(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                rOriginVariable, rDestinationVariable,
                MappingOptions, UseTranspose);
        }
    }

    template< typename TDataType >
    void TInverseMap(const Variable<TDataType>& rOriginVariable,
                     const Variable<TDataType>& rDestinationVariable,
                     Kratos::Flags MappingOptions,
                     const bool UseTranspose = false)
    {
        CheckForConservative(MappingOptions);

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE))
        {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE); // TODO test this!!!
            const bool use_transpose = true;
            TMap(rOriginVariable, rDestinationVariable, MappingOptions, use_transpose);
        }
        else
            GetInverseMapper()->TMap(rDestinationVariable, rOriginVariable, MappingOptions, UseTranspose);
    }

    MapperUniquePointerType& GetInverseMapper()
    {
        InitializeInverseMapper(); // Checks if it was initialized
        return mpInverseMapper;
    }

    // template< typename T>
    // void TestFunction(T someParam);

    // TSystemMatrixType& GetMdo()
    // {
    //     TSystemMatrixType& rMdo = *mpMdo;

    //     return rMdo;
    // }

    // TSystemVectorType& GetQo()
    // {
    //     TSystemVectorType& rQo = *mpQo;

    //     return rQo;
    // }

    // TSystemVectorType& GetQd()
    // {
    //     TSystemVectorType& rQd = *mpQd;

    //     return rQd;
    // }

    /**
     * This function can be overridden by derived Mappers if they need some
     * Mapper-specific settings. After initializing/validating the mapper-specific
     * settings this function should set the mGeneralMapperParameters
     * */
    virtual void ValidateMapperSpecificSettings(Parameters AllMapperSettings)
    {
        mGeneralMapperSettings = AllMapperSettings;
    }

    virtual MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const = 0;
    virtual MapperLocalSystemPointer GetMapperLocalSystem() const = 0;

    virtual InterfaceObject::ConstructionType GetInterfaceObjectConstructionTypeOrigin() const = 0;


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    bool mInverseMapperIsInitialized = false;

    MapperUniquePointerType mpInverseMapper;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeInverseMapper()
    {
        if (!mInverseMapperIsInitialized)
        {
            mpInverseMapper = Clone(mrModelPartDestination,
                                    mrModelPartOrigin,
                                    mGeneralMapperSettings); // TODO how to handle this ...? => some parameters wil be validated in the derived clases (Mappers)

            mInverseMapperIsInitialized = true;
        }
    }

    void ValidateParameters(Parameters AllMapperSettings)
    {
        ValidateMapperSpecificSettings(AllMapperSettings);
        mGeneralMapperSettings.RecursivelyValidateAndAssignDefaults(mGeneralMapperSettingsDefaults);
    }

    // From outside the user might specify CONSERVATIVE
    // This is translated to USE_TRANSPOSE for internal use
    // Note that if the user would specify USE_TRANSPOSE it would have the same effect
    void CheckForConservative(Kratos::Flags& rMappingOptions)
    {
        if (rMappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            rMappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!!
            rMappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!
        }
    }

    void AssignInterfaceEquationIds();
    void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // Mapper& operator=(Mapper const& rOther) {}

    /// Copy constructor.
    //Mapper(Mapper const& rOther);

    ///@}

}; // Class Mapper

}  // namespace Kratos.

#endif // KRATOS_MAPPER_H_INCLUDED  defined