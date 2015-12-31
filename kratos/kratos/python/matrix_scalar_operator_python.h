// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#if !defined(KRATOS_MATRIX_SCALAR_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_MATRIX_SCALAR_OPERATOR_PYTHON_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


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
template<class TMatrixType, class TScalarType, class TAddResultType = TMatrixType, class TMultResultType = TMatrixType>
class MatrixScalarOperatorPython : public def_visitor<MatrixScalarOperatorPython<TMatrixType, TScalarType, TAddResultType, TMultResultType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixScalarOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(MatrixScalarOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MatrixScalarOperatorPython() {}

    /// Copy constructor.
    MatrixScalarOperatorPython(const MatrixScalarOperatorPython& rOther) {}

    /// Destructor.
    virtual ~MatrixScalarOperatorPython() {}


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
        .def("__div__", &div)
        .def("__radd__", &radd)
        .def("__rsub__", &rsub)
        .def("__rmul__", &rmul)
        //.def("__rdiv__", &rdiv)
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
    TAddResultType
    add(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix + scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
    }

    static
    TAddResultType
    sub(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix - scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
    }

    static
    TMultResultType
    mul(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix * ThisScalar;
    }

    static
    TMultResultType
    div(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisMatrix / ThisScalar;
    }

    static
    TAddResultType
    radd(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar) + ThisMatrix;
    }

    static
    TAddResultType
    rsub(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar) - ThisMatrix;
    }

    static
    TMultResultType
    rmul(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    {
        return ThisScalar * ThisMatrix;
    }

    //static
    //	TMultResultType
    //	rdiv(TMatrixType& ThisMatrix, TScalarType ThisScalar)
    //{
    //	return scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar) / ThisMatrix;
    //}



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
    MatrixScalarOperatorPython& operator=(const MatrixScalarOperatorPython& rOther);


    ///@}

}; // Class MatrixScalarOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_MATRIX_SCALAR_OPERATOR_PYTHON_H_INCLUDED  defined 


