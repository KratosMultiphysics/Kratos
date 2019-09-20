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


#if !defined(KRATOS_VECTOR_VECTOR_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_VECTOR_VECTOR_OPERATOR_PYTHON_H_INCLUDED



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
template<class TContainerType, class TOtherContainerType, class TResultType>
class VectorVectorOperatorPython : public def_visitor<VectorVectorOperatorPython<TContainerType, TOtherContainerType, TResultType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VectorVectorOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(VectorVectorOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VectorVectorOperatorPython() {}

    /// Copy constructor.
    VectorVectorOperatorPython(const VectorVectorOperatorPython& rOther);

    /// Destructor.
    virtual ~VectorVectorOperatorPython() {}


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
    add(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        return ThisContainer + OtherContainer;
    }

    static
    TResultType
    sub(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        return ThisContainer - OtherContainer;
    }

    static
    typename TResultType::value_type
    mul(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        return inner_prod(ThisContainer, OtherContainer);
    }

    static
    TResultType
    radd(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        return OtherContainer + ThisContainer;
    }

    static
    TResultType
    rsub(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        return OtherContainer - ThisContainer;
    }

    static
    typename TResultType::value_type
    rmul(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        return inner_prod(OtherContainer, ThisContainer);
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
    VectorVectorOperatorPython& operator=(const VectorVectorOperatorPython& rOther);


    ///@}

}; // Class VectorVectorOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_VECTOR_VECTOR_OPERATOR_PYTHON_H_INCLUDED  defined 


