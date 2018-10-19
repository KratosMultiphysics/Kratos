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
#include "custom_searching/interface_communicator.h"
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
/** This is the base class for every mapper. This is the equivalent to a Kratos-SolvingStrategy.
 * By default it solves the mapping-system Qd = Mdo * Qo, but it can be extended, using the same
 * idea as the Solving-Strategies in Python
 * It holds the basic helper-objects needed to do the mapping, e.g. searching, building the mapping-
 * system and preparing the Interfaces
 * Furthermore it manages the InverseMapper, which is the same mapper but with reversed input of
 * ModelParts to map in the opposite direction
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
    typedef Kratos::unique_ptr<InterfacePreprocessor> InterfacePreprocessorPointerType;
    typedef Kratos::unique_ptr<InterfaceCommunicator> InterfaceCommunicatorPointerType;
    typedef std::size_t IndexType;

    typedef Kratos::unique_ptr<Mapper> MapperUniquePointerType;

    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    typedef typename MappingOperationUtilityType::MapperLocalSystemPointer MapperLocalSystemPointer;
    typedef typename MappingOperationUtilityType::MapperLocalSystemPointerVector MapperLocalSystemPointerVector;
    typedef typename MappingOperationUtilityType::MapperLocalSystemPointerVectorPointer MapperLocalSystemPointerVectorPointer;

    typedef typename MappingOperationUtilityType::TSystemMatrixType TSystemMatrixType;
    typedef typename MappingOperationUtilityType::TSystemVectorType TSystemVectorType;

    typedef typename MappingOperationUtilityType::TSystemMatrixUniquePointerType TSystemMatrixUniquePointerType;
    typedef typename MappingOperationUtilityType::TSystemVectorUniquePointerType TSystemVectorUniquePointerType;

    typedef typename MappingOperationUtilityType::ComponentVariableType ComponentVariableType;


    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
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

    /**
    * @brief Updates the mapping-system after the geometry/mesh has changed
    * After changes in the topology (e.g. remeshing or sliding interfaces)
    * the relations for the mapping have to be recomputed. This means that
    * the search has to be conducted again and the mapping-system has to be
    * rebuilt, hence this is expensive
    * @param MappingOptions flags used to specify how the update has to be done
    * @param SearchRadius search radius used for the search
    * @author Philipp Bucher
    */
    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius)
    {
        // std::cout << "\n\n";
        // if (MappingOptions.IsDefined(MapperFlags::REMESHED))
        // {
        //     KRATOS_INFO("\tFlagsCheck") << "IsDefined" << std::endl;
        //     if (MappingOptions.Is(MapperFlags::REMESHED))
        //     {
        //         KRATOS_INFO("\t\tFlagsCheck") << "IsDefined AND Is" << std::endl;
        //     }
        //     else
        //     {
        //         KRATOS_INFO("\t\tFlagsCheck") << "IsDefined AND NOT Is" << std::endl;
        //     }
        // }
        // else
        // {
        //     KRATOS_INFO("\tFlagsCheck") << "NOT IsDefined" << std::endl;
        //                 KRATOS_INFO("\tFlagsCheck") << "IsDefined" << std::endl;
        //     if (MappingOptions.Is(MapperFlags::REMESHED))
        //     {
        //         KRATOS_INFO("\t\tFlagsCheck") << "NOT IsDefined AND Is" << std::endl;
        //     }
        //     else
        //     {
        //         KRATOS_INFO("\t\tFlagsCheck") << "NOT IsDefined AND NOT Is" << std::endl;
        //     }
        // }

        // std::cout << "\n";

        // if (MappingOptions.Is(MapperFlags::REMESHED))
        // {
        //     KRATOS_INFO("\tFlagsCheck") << "Is" << std::endl;
        // }
        // else
        // {
        //     KRATOS_INFO("\tFlagsCheck") << "NOT Is" << std::endl;
        // }

        // std::cout << "\n\n";

        UpdateInterfaceInternal(MappingOptions, SearchRadius);
        if (mpInverseMapper) {
            mpInverseMapper->UpdateInterface(MappingOptions, SearchRadius);
        }
    }

    /**
    * @brief Mapping from Origin to Destination, Scalar Variable
    * Data is exchanged on the Interface, from the Origin-Modelpart
    * to the Destination-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see InverseMap
    * @author Philipp Bucher
    */
    virtual void Map(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions)
    {
        TMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /**
    * @brief Mapping from Origin to Destination, Vector Variable
    * Same as Map, but maps an array3-variable
    * @see Map
    * @author Philipp Bucher
    */
    virtual void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
                     const Variable< array_1d<double, 3> >& rDestinationVariable,
                     Kratos::Flags MappingOptions)
    {
        TMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /**
    * @brief Mapping from Destination to Origin, Scalar Variable
    * Data is exchanged on the Interface, from the Destination-Modelpart
    * to the Origin-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * It does the opposite of Map
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see Map
    * @author Philipp Bucher
    */
    virtual void InverseMap(const Variable<double>& rOriginVariable,
                            const Variable<double>& rDestinationVariable,
                            Kratos::Flags MappingOptions)
    {
        TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /**
    * @brief Mapping from Destination to Origin, Vector Variable
    * Same as InveseMap, but maps an array3-variable
    * @see InverseMap
    * @author Philipp Bucher
    */
    virtual void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                            const Variable< array_1d<double, 3> >& rDestinationVariable,
                            Kratos::Flags MappingOptions)
    {
        TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
    }

    /**
    * @brief Cloning the Mapper
    * returns a clone of the current Mapper
    * pure virtual, has to be implemented in every derived mapper,
    * used in the creation of the Mappers
    * @see MapperFactory
    * @author Philipp Bucher
    */
    virtual MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                          ModelPart& rModelPartDestination,
                                          Parameters JsonParameters) = 0;

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

    Parameters mGeneralMapperSettings;

    // utilities needed by the mapper
    MappingOperationUtilityPointerType mpMappingOperationUtility;
    InterfacePreprocessorPointerType mpInterfacePreprocessor;
    InterfaceCommunicatorPointerType mpIntefaceCommunicator;
    MapperLocalSystemPointerVectorPointer mpMapperLocalSystems;

    // The mapping matrix and the corresponding vectors
    // size of mpQo: num-nodes origin
    // size of mpQd: num-nodes destination
    // size of mpMdo: num-nodes destination x num-nodes origin
    TSystemMatrixUniquePointerType mpMdo;
    TSystemVectorUniquePointerType mpQo;
    TSystemVectorUniquePointerType mpQd;

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    // Constructor, can only be called by derived classes (actual mappers)
    // Mappers should only be constructed through the MapperFactory
    Mapper(ModelPart& rModelPartOrigin,
           ModelPart& rModelPartDestination,
           Parameters MapperSettings);

    /**
     * This function can be overridden by derived Mappers to do sth different
     * */
    virtual void Initialize();

    virtual void InitializeInterfaceCommunicator();

    virtual void InitializeMappingOperationUtility();

    virtual void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    virtual void BuildMappingMatrix(Kratos::Flags MappingOptions = Kratos::Flags());

    virtual void UpdateInterfaceInternal(Kratos::Flags MappingOptions, double SearchRadius);

    MapperUniquePointerType& GetInverseMapper()
    {
        InitializeInverseMapper(); // Checks if it was initialized
        return mpInverseMapper;
    }

        /* This function maps from Destination to Origin */
    virtual void MapInternal(const Variable<double>& rOriginVariable,
                             const Variable<double>& rDestinationVariable,
                             Kratos::Flags MappingOptions,
                             const bool UseTranspose)
    {
        mpMappingOperationUtility->InitializeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                rOriginVariable, rDestinationVariable,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->ExecuteMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                rOriginVariable, rDestinationVariable,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->FinalizeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                rOriginVariable, rDestinationVariable,
                MappingOptions, UseTranspose);
    }

    /* This function maps from Destination to Origin */
    virtual void MapInternal(const Variable< array_1d<double, 3> >& rOriginVariable,
                             const Variable< array_1d<double, 3> >& rDestinationVariable,
                             Kratos::Flags MappingOptions,
                             const bool UseTranspose)
    {
        const auto& var_x_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_X");
        const auto& var_y_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Y");
        const auto& var_z_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Z");

        const auto& var_x_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_X");
        const auto& var_y_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Y");
        const auto& var_z_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Z");

        // X-Component
        mpMappingOperationUtility->InitializeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_x_origin, var_x_destination,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->ExecuteMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_x_origin, var_x_destination,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->FinalizeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_x_origin, var_x_destination,
                MappingOptions, UseTranspose);

        // Y-Component
        mpMappingOperationUtility->InitializeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_y_origin, var_y_destination,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->ExecuteMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_y_origin, var_y_destination,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->FinalizeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_y_origin, var_y_destination,
                MappingOptions, UseTranspose);

        // Z-Component
        mpMappingOperationUtility->InitializeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_z_origin, var_z_destination,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->ExecuteMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_z_origin, var_z_destination,
                MappingOptions, UseTranspose);
        mpMappingOperationUtility->FinalizeMappingStep(
                *mpMdo, *mpQo, *mpQd,
                mrModelPartOrigin, mrModelPartDestination,
                var_z_origin, var_z_destination,
                MappingOptions, UseTranspose);
    }

    /**
     * This function can be overridden by derived Mappers if they need some
     * Mapper-specific settings. After initializing/validating the mapper-specific
     * settings this function should set the mGeneralMapperParameters
     * */
    virtual void ValidateMapperSpecificSettings(Parameters AllMapperSettings)
    {
        // mGeneralMapperSettings = AllMapperSettings;
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

    MapperUniquePointerType mpInverseMapper = nullptr;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeInverseMapper()
    {
        if (!mpInverseMapper) {
            mpInverseMapper = Clone(mrModelPartDestination,
                                    mrModelPartOrigin,
                                    mGeneralMapperSettings); // TODO how to handle this ...? => some parameters wil be validated in the derived clases (Mappers)
        }
    }

    void ValidateInput(Parameters AllMapperSettings);

    void ValidateParameters(Parameters AllMapperSettings)
    {
        ValidateMapperSpecificSettings(AllMapperSettings);

        Parameters default_settings = Parameters( R"({
            "search_radius"            : -1.0,
            "search_iterations"        : 3,
            "echo_level"               : 0,
            "mapper_specific_settings" : {}
        })");

        mGeneralMapperSettings.ValidateAndAssignDefaults(default_settings);
    }

    template< typename TDataType >
    void TMap(const Variable<TDataType>& rOriginVariable,
              const Variable<TDataType>& rDestinationVariable,
              Kratos::Flags MappingOptions,
              const bool UseTranspose = false)
    {
        CheckForConservative(MappingOptions);

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE); // TODO test this!!!
            const bool use_transpose = true;
            TInverseMap(rOriginVariable, rDestinationVariable, MappingOptions, use_transpose);
        }
        else {
            KRATOS_DEBUG_ERROR_IF_NOT(mpMappingOperationUtility)<< "mpMappingOperationUtility "
                << "is a nullptr" << std::endl;

            MapInternal(rOriginVariable,
                        rDestinationVariable,
                        MappingOptions,
                        UseTranspose);
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

    void PrintPairingInfo();

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

/// output stream function
template<class TSparseSpace, class TDenseSpace>
inline std::ostream & operator << (
    std::ostream& rOStream,
    const Mapper<TSparseSpace, TDenseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_MAPPER_H_INCLUDED  defined