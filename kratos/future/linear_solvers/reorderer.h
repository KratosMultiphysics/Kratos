//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos::Future
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
    - TMatrixType which specify type
      of the unknowns, coefficients, sparse matrix, vector of
unknowns, right hand side vector and their respective operators.
    - TDenseMatrixType which specify type of the
      matrices used as temporary matrices or multi solve unknowns and
right hand sides and their operators.
*/
template<class TMatrixType, class TVectorType>
class Reorderer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Reorderer
    KRATOS_CLASS_POINTER_DEFINITION(Reorderer);

    using VectorType = TVectorType;

    using SparseMatrixType = TMatrixType;

    using DataType = typename SparseMatrixType::DataType;

    using SizeType = typename SparseMatrixType::SizeType;

    using IndexType = typename SparseMatrixType::IndexType;

    using DenseMatrixType = DenseMatrix<DataType>;

    using IndexVectorType = std::vector<IndexType>;

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

    virtual void Initialize(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
        CalculateIndexPermutation(rA);
    }

    virtual void Reorder(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
    }

    virtual void InverseReorder(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB)
    {
    }

    virtual IndexVectorType& CalculateIndexPermutation(
        SparseMatrixType& rA,
        IndexType InitialIndex = IndexType())
    {
        SizeType size = rA.size1();

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
template<class TMatrixType, class TVectorType>
inline std::istream& operator >> (
    std::istream& rIStream,
    Reorderer<TMatrixType, TVectorType>& rThis)
{
}

/// output stream function
template<class TMatrixType, class TVectorType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const Reorderer<TMatrixType, TVectorType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos::Future.
