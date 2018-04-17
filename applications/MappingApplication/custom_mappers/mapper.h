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
#include "custom_utilities/interface_communicator.h"
#include "custom_utilities/matrix_based_mapping_operation_utility.h"
#include "custom_utilities/interface_preprocessor.h"
#include "includes/kratos_parameters.h"
// #include "custom_utilities/mapper_utilities.h"
#include "custom_utilities/mapper_flags.h"

// For MPI-parallel Mapper
#ifdef KRATOS_USING_MPI
#include "mpi.h" // TODO needed here?
#include "custom_utilities/interface_communicator_mpi.h"
#endif


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

    using InterfaceCommunicatorPointerType = InterfaceCommunicator::Pointer;
    typedef MappingOperationUtility<TSparseSpace, TDenseSpace> MappingOperationUtilityType;
    typedef typename MappingOperationUtilityType::Pointer MappingOperationUtilityPointerType;
    using InterfacePreprocessorPointerType = InterfacePreprocessor::Pointer;
    using ModelPartPointerType = ModelPart::Pointer;
    using SizeType = std::size_t;
    using IndexType = std::size_t;

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

    virtual void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius);

    /* This function maps from Origin to Destination */
    virtual void Map(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    /* This function maps from Origin to Destination */
    virtual void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
                     const Variable< array_1d<double, 3> >& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    /* This function maps from Destination to Origin */
    virtual void InverseMap(const Variable<double>& rOriginVariable,
                            const Variable<double>& rDestinationVariable,
                            Kratos::Flags MappingOptions);

    /* This function maps from Destination to Origin */
    virtual void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                            const Variable< array_1d<double, 3> >& rDestinationVariable,
                            Kratos::Flags MappingOptions);

    virtual Mapper<TSparseSpace, TDenseSpace>::Pointer Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) = 0;

    // Developer function only being used if this class needs to be accessed from outside!
    InterfaceCommunicatorPointerType pGetInterfaceCommunicator()
    {
        return mpInterfaceCommunicator;
    }

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

    bool mInverseMapperIsInitialized = false;

    InterfaceCommunicatorPointerType mpInterfaceCommunicator;
    MappingOperationUtilityPointerType mpMappingOperationUtility;
    InterfacePreprocessorPointerType mpInterfacePreprocessor;
    ModelPartPointerType mpInterfaceModelPart;

    Mapper::Pointer mpInverseMapper;

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

    // template< class TVarType>
    // void InitializeMappingStep(const TVarType& rVarOrigin,
    //                            const TVarType& rVarDestination,
    //                            const Kratos::Flags& MappingOptions)
    // {
    //     if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) // TODO replace this with USE_TRANSPOSE!
    //     {
    //         mpMappingMatrixBuilder->UpdateSystemVector(mrModelPartDestination, *mpQd, rVarDestination);
    //     }
    //     else
    //     {
    //         mpMappingMatrixBuilder->UpdateSystemVector(mrModelPartOrigin, *mpQo, rVarOrigin);
    //     }
    // }

    // virtual void ExecuteMappingStep(const Kratos::Flags& MappingOptions) // Override this class in Mortar
    // {
    //     if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE))
    //     {
    //         const bool transpose_flag = true; // constexpr?
    //         mpMappingMatrixBuilder->Multiply(*mpMdo, *mpQd, *mpQo, transpose_flag);
    //     }
    //     else
    //     {
    //         mpMappingMatrixBuilder->Multiply(*mpMdo, *mpQo, *mpQd);
    //     }
    // }

    // template< class TVarType>
    // void FinalizeMappingStep(const TVarType& rVarOrigin,
    //                          const TVarType& rVarDestination,
    //                          const Kratos::Flags& MappingOptions)
    // {
    //     double factor = 1.0f;
    //     ProcessMappingOptions(MappingOptions, factor);

    //     if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE))
    //     {

    //         mpMappingMatrixBuilder->Update(mrModelPartOrigin, *mpQo, rVarOrigin, MappingOptions, factor);
    //     }
    //     else
    //     {
    //         mpMappingMatrixBuilder->Update(mrModelPartDestination, *mpQd, rVarDestination, MappingOptions, factor);
    //     }
    // }

    template< typename TDataType >
    void TMap(const Variable<TDataType>& rOriginVariable,
              const Variable<TDataType>& rDestinationVariable,
              Kratos::Flags MappingOptions)
    {
        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            MappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!!
            MappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!

            InverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else
        {
            mpMappingOperationUtility->ExecuteMapping(rOriginVariable,
                                                      rDestinationVariable,
                                                      MappingOptions);
        }
    }

    template< typename TDataType >
    void TInverseMap(const Variable<TDataType>& rOriginVariable,
                     const Variable<TDataType>& rDestinationVariable,
                     Kratos::Flags MappingOptions)
    {
        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            MappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!!
            MappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!

            Map(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else
        {
            // Construct the inverse mapper if it hasn't been done before
            if (!mInverseMapperIsInitialized) InitializeInverseMapper();

            mpInverseMapper->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    /**
     * This function can be overridden by derived Mappers if they need some
     * Mapper-specific settings. After initializing/validating the mapper-specific
     * settings this function should set the mGeneralMapperParameters
     * */
    virtual void ValidateMapperSpecificSettings(Parameters AllMapperSettings)
    {
        mGeneralMapperSettings = AllMapperSettings;
    }

    virtual InterfaceCommunicatorPointerType CreateInterfaceCommunicator(
        ModelPart& rModelPartOrigin, ModelPartPointerType pInterfaceModelPart) const
    {
        return Kratos::make_shared<InterfaceCommunicator>(rModelPartOrigin, pInterfaceModelPart);
    }

    virtual MappingOperationUtilityPointerType CreateMappingOperationUtility(
        ModelPartPointerType pInterfaceModelPart) const
    {   // here we could return the MatrixFree variant in the future
        return Kratos::make_shared<MatrixBasedMappingOperationUtility<TSparseSpace, TDenseSpace>>(pInterfaceModelPart);
    }

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    virtual InterfaceCommunicatorPointerType CreateMPIInterfaceCommunicator(
        ModelPart& rModelPartOrigin, ModelPartPointerType pInterfaceModelPart) const
    {
        return Kratos::make_shared<InterfaceCommunicatorMPI>(rModelPartOrigin, pInterfaceModelPart);
    }
#endif

    void InitializeInverseMapper()
    {
        mpInverseMapper = Clone(mrModelPartDestination, // TODO needs "this->" ?
                                mrModelPartOrigin,
                                mGeneralMapperSettings); // TODO how to handle this ...?

        mInverseMapperIsInitialized = true;
    }

    /**
    This function return information abt how the interface has to be constructed.
    This information is specific to every mapper!
    */
    virtual Parameters GetInterfaceParameters() = 0;



    // void ComputeNumberOfNodesAndConditions()
    // {
    //     // Compute the quantities of the local model_parts
    //     mNumConditionsOrigin = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfConditions();
    //     mNumConditionsDestination = mrModelPartDestination.GetCommunicator().LocalMesh().NumberOfConditions();

    //     mNumNodesOrigin = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
    //     mNumNodesDestination = mrModelPartDestination.GetCommunicator().LocalMesh().NumberOfNodes();

    //     // Compute the quantities of the global model_parts
    //     mrModelPartOrigin.GetCommunicator().SumAll(mNumConditionsOrigin);
    //     mrModelPartDestination.GetCommunicator().SumAll(mNumConditionsDestination);

    //     mrModelPartOrigin.GetCommunicator().SumAll(mNumNodesOrigin);
    //     mrModelPartDestination.GetCommunicator().SumAll(mNumNodesDestination);
    // }

    // void ProcessMappingOptions(const Kratos::Flags& rMappingOptions,
    //                            double& Factor)
    // {
    //     if (rMappingOptions.Is(MapperFlags::SWAP_SIGN))
    //     {
    //         Factor *= (-1);
    //     }
    // }

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ValidateParameters(Parameters AllMapperSettings)
    {
        ValidateMapperSpecificSettings(AllMapperSettings);
        mGeneralMapperSettings.RecursivelyValidateAndAssignDefaults(mGeneralMapperSettingsDefaults);
    }

    void GenerateInterfaceModelPart()
    {
        mpInterfacePreprocessor->GenerateInterfaceModelPart(GetInterfaceParameters());
    }

    void InitializeInterfaceCommunicator();

    void InitializeMappingOperationUtility();


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