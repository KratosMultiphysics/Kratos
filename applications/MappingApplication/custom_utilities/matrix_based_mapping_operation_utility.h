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
#include "mapper_utilities.h"


namespace Kratos
{
///@addtogroup MappingApplication
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

/// Object for constructing the mapping-matrices of the mapping-system
/** This class assembles the Mapping Matrix
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

    typedef MappingOperationUtility<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::MapperLocalSystemPointerVector MapperLocalSystemPointerVector;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::TSystemMatrixUniquePointerType TSystemMatrixUniquePointerType;
    typedef typename BaseType::TSystemVectorUniquePointerType TSystemVectorUniquePointerType;

    typedef typename BaseType::DoubleVariableType DoubleVariableType;
    typedef typename BaseType::ComponentVariableType ComponentVariableType;

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

    void BuildMappingSystem(
        TSystemMatrixUniquePointerType& rpMdo,
        TSystemVectorUniquePointerType& rpQo,
        TSystemVectorUniquePointerType& rpQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        MapperLocalSystemPointerVector& rMapperLocalSystems) const override;

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
        const bool UseTranspose) const override
    {
        ExecuteMapping(rMdo, rQo, rQd, UseTranspose);
    }

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
        const bool UseTranspose) const override
    {
        ExecuteMapping(rMdo, rQo, rQd, UseTranspose);
    }

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

    template< class TVectorType, class TVarType >
    void TInitializeMappingStep(
        TVectorType& rQo,
        TVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const TVarType& rOriginVariable,
        const TVarType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const
    {
        if (UseTranspose) {
            MapperUtilities::UpdateSystemVectorFromModelPart(
                rQd,
                rModelPartDestination,
                rDestinationVariable,
                MappingOptions);
        }
        else {
            MapperUtilities::UpdateSystemVectorFromModelPart(
                rQo,
                rModelPartOrigin,
                rOriginVariable,
                MappingOptions);
        }
    }

    template< class TVectorType, class TVarType >
    void TFinalizeMappingStep(
        TVectorType& rQo,
        TVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const TVarType& rOriginVariable,
        const TVarType& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const
    {
        if (UseTranspose) {
            MapperUtilities::UpdateModelPartFromSystemVector(
                rQo,
                rModelPartOrigin,
                rOriginVariable,
                MappingOptions);
        }
        else {
            MapperUtilities::UpdateModelPartFromSystemVector(
                rQd,
                rModelPartDestination,
                rDestinationVariable,
                MappingOptions);
        }
    }

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

    void ExecuteMapping(TSystemMatrixType& rMdo,
                        TSystemVectorType& rQo,
                        TSystemVectorType& rQd,
                        const bool UseTranspose) const
    {
        if (UseTranspose) {
            TSparseSpace::TransposeMult(rMdo, rQd, rQo); // rQo = rMdo^T * rQo
        }
        else {
            TSparseSpace::Mult(rMdo, rQo, rQd); // rQd = rMdo * rQo
        }
    }

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
