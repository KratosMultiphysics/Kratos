//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//

#if !defined(KRATOS_REORDERER_H_INCLUDED )
#define  KRATOS_REORDERER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

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

/// Base class for all reorderer objects in Kratos used in linear solvers.
/** Reorderer is the base class for all reorderer objects in Kratos used in linear solvers.
    two template parameter:
    - TSparseSpaceType which specify type
      of the unknowns, coefficients, sparse matrix, vector of
unknowns, right hand side vector and their respective operators.
    - TDenseMatrixType which specify type of the
      matrices used as temporary matrices or multi solve unknowns and
right hand sides and their operators.
*/
template<class TSparseSpaceType, class TDenseSpaceType>
class Reorderer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Reorderer
    KRATOS_CLASS_POINTER_DEFINITION(Reorderer);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef unsigned int IndexType;

    typedef unsigned int SizeType;

    typedef std::vector<IndexType> IndexVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Reorderer() {}

    /// Copy constructor.
    Reorderer(Reorderer const& rOther) {}

    /// Destructor.
    virtual ~Reorderer() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Reorderer& operator=(Reorderer const& rOther)
    {
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        CalculateIndexPermutation(rA);
    }

    virtual void Reorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    virtual void InverseReorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    virtual IndexVectorType& CalculateIndexPermutation(SparseMatrixType& rA, IndexType InitialIndex = IndexType())
    {
        SizeType size = TSparseSpaceType::Size1(rA);

        if(mIndexPermutation.size() != size)
            mIndexPermutation.resize(size);

        for(SizeType i = 0 ; i < size ; i++)
            mIndexPermutation[i] = i;

        return mIndexPermutation;
    }

    ///@}
    ///@name Access
    ///@{

    virtual IndexVectorType& GetIndexPermutation()
    {
        return mIndexPermutation;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Reorderer";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Reorderer";
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

    IndexVectorType mIndexPermutation;


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


    ///@}

}; // Class Reorderer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Reorderer<TSparseSpaceType, TDenseSpaceType>& rThis);

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Reorderer<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REORDERER_H_INCLUDED  defined 


