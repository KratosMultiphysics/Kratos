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
//                   Riccardo Rossi
//





#if !defined(KRATOS_MATRIX_VECTOR_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_MATRIX_VECTOR_OPERATOR_PYTHON_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
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

/// Short class definition.
/** Detail class definition.
*/
template<class TMatrixType, class TVectorType, class TResultType = TVectorType>
class MatrixVectorOperatorPython : public def_visitor<MatrixVectorOperatorPython<TMatrixType, TVectorType, TResultType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixVectorOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(MatrixVectorOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixVectorOperatorPython() {}

    /// Copy constructor.
    MatrixVectorOperatorPython(const MatrixVectorOperatorPython& rOther) {}


    /// Destructor.
    virtual ~MatrixVectorOperatorPython() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template <class TClassType>
    void visit(TClassType& ThisClass) const
    {
        ThisClass
        .def("__mul__", &mul)
        .def("__rmul__", &rmul)
        ;
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


    ///@}
    ///@name Friends
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

    static
    TResultType
    mul(TMatrixType& ThisMatrix, TVectorType ThisVector)
    {
        return prod(ThisMatrix, ThisVector);
    }

    static
    TResultType
    rmul(TMatrixType& ThisMatrix, TVectorType ThisVector)
    {
        return prod(ThisVector, ThisMatrix);
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
    MatrixVectorOperatorPython& operator=(const MatrixVectorOperatorPython& rOther);



    ///@}

}; // Class MatrixVectorOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_MATRIX_VECTOR_OPERATOR_PYTHON_H_INCLUDED  defined 


