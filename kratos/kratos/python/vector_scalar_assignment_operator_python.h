/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_VECTOR_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_VECTOR_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED



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
template<class TContainerType, class TScalarType>
class VectorScalarAssignmentOperatorPython : public def_visitor<VectorScalarAssignmentOperatorPython<TContainerType, TScalarType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VectorScalarAssignmentOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(VectorScalarAssignmentOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VectorScalarAssignmentOperatorPython() {}

    /// Copy constructor.
    VectorScalarAssignmentOperatorPython(const VectorScalarAssignmentOperatorPython& rOther);

    /// Destructor.
    virtual ~VectorScalarAssignmentOperatorPython() {}


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
    TContainerType&
    iadd(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            ThisContainer[i] += ThisScalar;
        return ThisContainer;
    }

    static
    TContainerType&
    isub(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            ThisContainer[i] -= ThisScalar;
        return ThisContainer;
    }

    static
    TContainerType&
    imul(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            ThisContainer[i] *= ThisScalar;
        return ThisContainer;
    }

    static
    TContainerType&
    idiv(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            ThisContainer[i] /= ThisScalar;
        return ThisContainer;
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
    VectorScalarAssignmentOperatorPython& operator=(const VectorScalarAssignmentOperatorPython& rOther);


    ///@}

}; // Class VectorScalarAssignmentOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_VECTOR_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED  defined 


