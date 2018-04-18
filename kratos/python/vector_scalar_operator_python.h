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




#if !defined(KRATOS_VECTOR_SCALAR_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_VECTOR_SCALAR_OPERATOR_PYTHON_H_INCLUDED



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
template<class TContainerType, class TScalarType, class TResultType>
class VectorScalarOperatorPython : public def_visitor<VectorScalarOperatorPython<TContainerType, TScalarType, TResultType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VectorScalarOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(VectorScalarOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VectorScalarOperatorPython() {}

    /// Copy constructor.
    VectorScalarOperatorPython(const VectorScalarOperatorPython& rOther);

    /// Destructor.
    virtual ~VectorScalarOperatorPython() {}


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
        .def("__rdiv__", &rdiv)
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
    add(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        TResultType result(ThisContainer);
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            result[i] += ThisScalar;
        return result;
    }

    static
    TResultType
    sub(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        TResultType result(ThisContainer);
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            result[i] -= ThisScalar;
        return result;
    }

    static
    TResultType
    mul(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        return ThisContainer * ThisScalar;
    }

    static
    TResultType
    div(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        return ThisContainer / ThisScalar;
    }

    static
    TResultType
    radd(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        TResultType result(ThisContainer);
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            result[i] += ThisScalar;
        return result;
    }

    static
    TResultType
    rsub(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        TResultType result(ThisContainer);
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            result[i] = ThisScalar - result[i];
        return result;
    }

    static
    TResultType
    rmul(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        return ThisScalar * ThisContainer;
    }

    static
    TResultType
    rdiv(TContainerType& ThisContainer, TScalarType ThisScalar)
    {
        TResultType result(ThisContainer);
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            result[i] = ThisScalar / result[i];
        return result;
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
    VectorScalarOperatorPython& operator=(const VectorScalarOperatorPython& rOther);


    ///@}

}; // Class VectorScalarOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_VECTOR_SCALAR_OPERATOR_PYTHON_H_INCLUDED  defined 


