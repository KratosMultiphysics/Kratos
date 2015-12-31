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



#if !defined(KRATOS_MATRIX_MATRIX_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_MATRIX_MATRIX_OPERATOR_PYTHON_H_INCLUDED



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


