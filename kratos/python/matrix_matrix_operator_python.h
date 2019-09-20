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



#if !defined(KRATOS_MATRIX_MATRIX_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_MATRIX_MATRIX_OPERATOR_PYTHON_H_INCLUDED



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
template<class TMatrixType, class TOtherMatrixType, class TResultType>
class MatrixMatrixOperatorPython : public def_visitor<MatrixMatrixOperatorPython<TMatrixType, TOtherMatrixType, TResultType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixMatrixOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(MatrixMatrixOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixMatrixOperatorPython() {}


    /// Copy constructor.
    MatrixMatrixOperatorPython(const MatrixMatrixOperatorPython& rOther) {}

    /// Destructor.
    virtual ~MatrixMatrixOperatorPython() {}


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
        .def("__add__", &add)
        .def("__sub__", &sub)
        .def("__mul__", &mul)
        .def("__radd__", &radd)
        .def("__rsub__", &rsub)
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
    add(TMatrixType& ThisContainer, TOtherMatrixType const& OtherContainer)
    {
        return ThisContainer + OtherContainer;
    }

    static
    TResultType
    sub(TMatrixType& ThisContainer, TOtherMatrixType const& OtherContainer)
    {
        return ThisContainer - OtherContainer;
    }

    static
    TResultType
    mul(TMatrixType& ThisContainer, TOtherMatrixType const& OtherContainer)
    {
        return prod(ThisContainer, OtherContainer);
    }

    static
    TResultType
    radd(TMatrixType& ThisContainer, TOtherMatrixType const& OtherContainer)
    {
        return OtherContainer + ThisContainer;
    }

    static
    TResultType
    rsub(TMatrixType& ThisContainer, TOtherMatrixType const& OtherContainer)
    {
        return OtherContainer - ThisContainer;
    }

    static
    TResultType
    rmul(TMatrixType& ThisContainer, TOtherMatrixType const& OtherContainer)
    {
        return prod(OtherContainer, ThisContainer);
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
    MatrixMatrixOperatorPython& operator=(const MatrixMatrixOperatorPython& rOther);


    ///@}

}; // Class MatrixMatrixOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_MATRIX_MATRIX_OPERATOR_PYTHON_H_INCLUDED  defined 


