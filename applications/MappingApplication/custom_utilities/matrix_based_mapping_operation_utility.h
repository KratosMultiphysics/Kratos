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

#if !defined(KRATOS_MATRIX_BASED_MAPPING_OPERATION_UTILITY_H )
#define  KRATOS_MATRIX_BASED_MAPPING_OPERATION_UTILITY_H

// System includes

// External includes

// Project includes
#include "mapping_operation_utility.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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

/// Short class definition.
/** Detail class definition.
*/
template<class TSparseSpace, class TDenseSpace>
class MatrixBasedMappingOperationUtility
    : public MappingOperationUtility<TSparseSpace, TDenseSpace>
{
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixBasedMappingOperationUtility
    KRATOS_CLASS_POINTER_DEFINITION(MatrixBasedMappingOperationUtility);

    using BaseType = MappingOperationUtility<TSparseSpace, TDenseSpace>;

    using SizeType = typename BaseType::SizeType;
    using IndexType = typename BaseType::IndexType;

    using MapperLocalSystemPointerType = typename BaseType::MapperLocalSystemPointer;
    using MapperLocalSystemPointerVector = typename BaseType::MapperLocalSystemPointerVector;
    using MapperLocalSystemPointerVectorPointer = typename BaseType::MapperLocalSystemPointerVectorPointer;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::TSystemMatrixUniquePointerType TSystemMatrixUniquePointerType;
    typedef typename BaseType::TSystemVectorUniquePointerType TSystemVectorUniquePointerType;

    using DoubleVariableType = typename BaseType::DoubleVariableType;
    using ComponentVariableType = typename BaseType::ComponentVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixBasedMappingOperationUtility(Parameters Settings);

    /// Destructor.
    virtual ~MatrixBasedMappingOperationUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ResizeAndInitializeVectors(
        TSystemMatrixUniquePointerType& rpMdo,
        TSystemVectorUniquePointerType& rpQo,
        TSystemVectorUniquePointerType& rpQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        MapperLocalSystemPointerVector& rMapperLocalSystems) const override;

    // The "Build" function
    void BuildMappingMatrix(const MapperLocalSystemPointerVector& rMapperLocalSystems,
                            TSystemMatrixType& rMdo) const override;

    // The "Solve" function
    void InitializeMappingStep(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const DoubleVariableType& rOriginVariable,
        const DoubleVariableType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const override;

    // The "Solve" function
    void InitializeMappingStep(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const ComponentVariableType& rOriginVariable,
        const ComponentVariableType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const override;

    // The "Solve" function
    void ExecuteMappingStep(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const DoubleVariableType& rOriginVariable,
        const DoubleVariableType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const override;

    // The "Solve" function
    void ExecuteMappingStep(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const ComponentVariableType& rOriginVariable,
        const ComponentVariableType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const override;

    // The "Solve" function
    void FinalizeMappingStep(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const DoubleVariableType& rOriginVariable,
        const DoubleVariableType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const override;

    // The "Solve" function
    void FinalizeMappingStep(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const ComponentVariableType& rOriginVariable,
        const ComponentVariableType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const override;

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








            //     VectorComponentType var_component_x_origin =
//         KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_X"));
//     VectorComponentType var_component_y_origin =
//         KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_Y"));
//     VectorComponentType var_component_z_origin =
//         KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_Z"));

//     VectorComponentType var_component_x_destination =
//         KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_X"));
//     VectorComponentType var_component_y_destination =
//         KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_Y"));
//     VectorComponentType var_component_z_destination =
//         KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_Z"));

//     // X-Component
//     InitializeMappingStep<VectorComponentType>(var_component_x_origin,
//                                                var_component_x_destination,
//                                                MappingOptions);

//     ExecuteMappingStep(MappingOptions);

//     FinalizeMappingStep<VectorComponentType>(var_component_x_origin,
//                                              var_component_x_destination,
//                                              MappingOptions);

//     // Y-Component
//     InitializeMappingStep<VectorComponentType>(var_component_y_origin,
//                                                var_component_y_destination,
//                                                MappingOptions);

//     ExecuteMappingStep(MappingOptions);

//     FinalizeMappingStep<VectorComponentType>(var_component_y_origin,
//                                              var_component_y_destination,
//                                              MappingOptions);

//     // Z-Component
//     InitializeMappingStep<VectorComponentType>(var_component_z_origin,
//                                                var_component_z_destination,
//                                                MappingOptions);

//     ExecuteMappingStep(MappingOptions);

//     FinalizeMappingStep<VectorComponentType>(var_component_z_origin,
//                                              var_component_z_destination,
//                                              MappingOptions);


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
    std::string Info() const override
    {
        return "MatrixBasedMappingOperationUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


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
    // MatrixBasedMappingOperationUtility& operator=(MatrixBasedMappingOperationUtility const& rOther) {}

    /// Copy constructor.
    // MatrixBasedMappingOperationUtility(MatrixBasedMappingOperationUtility const& rOther) {}


    ///@}

}; // Class MatrixBasedMappingOperationUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MATRIX_BASED_MAPPING_OPERATION_UTILITY_H  defined
