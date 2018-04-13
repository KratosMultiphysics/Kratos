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

#if !defined(KRATOS_MATRIX_BASED_MAPPING_OPERATION_UTILITY_MPI_H)
#define  KRATOS_MATRIX_BASED_MAPPING_OPERATION_UTILITY_MPI_H

// System includes

// External includes
// #include "mpi.h"// included in "trilinos_space"
#include "Epetra_FEVector.h"
// #include "Epetra_FECrsMatrix.h" // included in "trilinos_space"

// Project includes
#include "includes/define.h"
#include "trilinos_space.h"
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
class MatrixBasedMappingOperationUtilityMPI : public MappingOperationUtility
{
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixBasedMappingOperationUtilityMPI
    KRATOS_CLASS_POINTER_DEFINITION(MatrixBasedMappingOperationUtilityMPI);

    using BaseType = MappingOperationUtility;
    using ModelPartPointerType = BaseType::ModelPartPointerType;

    using SparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixBasedMappingOperationUtilityMPI(ModelPartPointerType pInterfaceModelPart);

    /// Destructor.
    virtual ~MatrixBasedMappingOperationUtilityMPI() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
        return "MatrixBasedMappingOperationUtilityMPI";
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
    // MatrixBasedMappingOperationUtilityMPI& operator=(MatrixBasedMappingOperationUtilityMPI const& rOther) {}

    /// Copy constructor.
    // MatrixBasedMappingOperationUtilityMPI(MatrixBasedMappingOperationUtilityMPI const& rOther) {}

    ///@}

    }; // Class MatrixBasedMappingOperationUtilityMPI

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MATRIX_BASED_MAPPING_OPERATION_UTILITY_MPI_H  defined
