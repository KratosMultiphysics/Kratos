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




#if !defined(KRATOS_VECTOR_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_VECTOR_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED



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


