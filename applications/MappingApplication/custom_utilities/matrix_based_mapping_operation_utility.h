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
    using MapperLocalSystemPointerType = typename BaseType::MapperLocalSystemPointer;
    using MapperLocalSystemPointerVector = typename BaseType::MapperLocalSystemPointerVector;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixBasedMappingOperationUtility();

    /// Destructor.
    virtual ~MatrixBasedMappingOperationUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ResizeAndInitializeVectors(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination) override { }

    // The "Build" function
    void BuildMappingMatrix(MapperLocalSystemPointerVector& rMapperLocalSystems,
                                    TSystemMatrixType& rMdo) override { }


    void UpdateInterface() override { }

    // The "Solve" function
    void ExecuteMapping(const Variable<double>& rOriginVariable,
                                const Variable<double>& rDestinationVariable,
                                Kratos::Flags MappingOptions) override { }

    // The "Solve" function
    void ExecuteMapping(const Variable<array_1d<double, 3>>& rOriginVariable,
                                const Variable<array_1d<double, 3>>& rDestinationVariable,
                                Kratos::Flags MappingOptions) override { }


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
