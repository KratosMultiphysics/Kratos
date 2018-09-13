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



#if !defined(KRATOS_MATRIX_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_MATRIX_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED



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
template<class TMatrixType, class TScalarType>
class MatrixScalarAssignmentOperatorPython : public def_visitor<MatrixScalarAssignmentOperatorPython<TMatrixType, TScalarType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixScalarAssignmentOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(MatrixScalarAssignmentOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixScalarAssignmentOperatorPython() {}

    /// Copy constructor.
    MatrixScalarAssignmentOperatorPython(const MatrixScalarAssignmentOperatorPython& rOther);

    /// Destructor.
    virtual ~MatrixScalarAssignmentOperatorPython() {}


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
        .def("__iadd__", &iadd, return_internal_reference<1>())
        .def("__isub__", &isub, return_internal_reference<1>())
        .def("__imul__", &imul, return_internal_reference<1>())
        .def("__idiv__", &idiv, return_internal_reference<1>())
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
    TMatrixType&
    iadd(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix += scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
    }

    static
    TMatrixType&
    isub(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix += scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
    }

    static
    TMatrixType&
    imul(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix *= ThisScalar;
    }

    static
    TMatrixType&
    idiv(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix /= ThisScalar;
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
    MatrixScalarAssignmentOperatorPython& operator=(const MatrixScalarAssignmentOperatorPython& rOther);


    ///@}

}; // Class MatrixScalarAssignmentOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_MATRIX_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED  defined 


